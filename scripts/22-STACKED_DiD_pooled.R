# =============================================================================
# CAZ DiD ANALYSIS — CLEAN PRIMARY PIPELINE
# =============================================================================
#
# Aim:
#   Estimate the effect of CAZ schemes on road traffic injuries using a
#   road-link × quarter panel.
#
# Main model:
#   Stacked event-study Difference-in-Differences estimated using PPML.
#
# Why:
#   - Injury outcome is a rare count with many zeros.
#   - CAZ timing differs across schemes.
#   - Effects may vary over time.
#   - Effects may differ by scheme.
#   - Treated roads had different COVID lockdown/recovery patterns.
#
# Main outputs:
#   1. Overall PPML IRR
#   2. Overall PPML dynamic event-study plot
#   3. Scheme-specific PPML event-study plots
#   4. C&S robustness event-study plot
#   5. Final summary table
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================
# Load only packages needed for the cleaned analysis.

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(lubridate)
library(fixest)
library(did)
library(patchwork)


# =============================================================================
# 1. SETTINGS
# =============================================================================
# Define the outcome, event-study window, and output folder.

outcome_var <- "total_inj_adj_All"

K <- 8L   # quarters before CAZ implementation
L <- 8L   # quarters after CAZ implementation

outdir <- here("output", "pooled", "All_clean")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
# These functions keep the rest of the code short and readable.

# Extract event-study coefficients from fixest models.
# This avoids memory-heavy broom::tidy(conf.int = TRUE).
extract_fixest_event_study <- function(model) {
  
  ct <- coeftable(model)
  
  tibble(
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"]
  ) %>%
    filter(str_detect(term, "event_time_f::")) %>%
    mutate(
      event_time = str_extract(term, "(?<=event_time_f::)-?\\d+") %>% as.numeric(),
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se,
      irr = exp(estimate),
      irr_lo = exp(ci_lo),
      irr_hi = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo = 100 * (irr_lo - 1),
      pct_hi = 100 * (irr_hi - 1)
    ) %>%
    arrange(event_time)
}


# Plot event-study coefficients on the log-IRR or ATT scale.
plot_event_study <- function(df, title, subtitle, ylab, colour = "#E74C3C") {
  
  ggplot(df, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = colour) +
    geom_line(linewidth = 0.8, colour = colour) +
    geom_point(size = 1.8, colour = colour) +
    scale_x_continuous(breaks = seq(-K, L, by = 2)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Quarters relative to CAZ implementation",
      y = ylab
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


# Run C&S as a robustness check.
run_cs <- function(data, xformla = NULL) {
  
  att <- att_gt(
    yname         = "outcome_raw",
    tname         = "qtr_int",
    idname        = "uid_int",
    gname         = "g",
    data          = data,
    control_group = "notyettreated",
    xformla       = xformla,
    bstrap        = TRUE,
    anticipation  = 0,
    panel         = TRUE
  )
  
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  dyn <- aggte(att, type = "dynamic", na.rm = TRUE)
  
  list(att = att, agg = agg, dyn = dyn)
}


# Extract C&S dynamic effects into a plotting table.
extract_cs_event_study <- function(dyn, trim = c(-12, 12)) {
  
  tibble(
    event_time = dyn$egt,
    estimate   = dyn$att.egt,
    se         = dyn$se.egt
  ) %>%
    mutate(
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se
    ) %>%
    filter(event_time >= trim[1], event_time <= trim[2])
}


# =============================================================================
# 3. LOAD ROAD-LINK PANEL
# =============================================================================
# Load matched road-link × quarter data and convert quarter variables.

road_panel <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

cat("Loaded rows:", nrow(road_panel), "\n")
cat("Road links:", n_distinct(road_panel$identifier), "\n")
cat("Quarters:", n_distinct(road_panel$quarter_year), "\n")


# =============================================================================
# 4. ADJUST CAZ START QUARTERS
# =============================================================================
# If a CAZ started in the second half of a quarter, assign treatment to the
# next quarter. This avoids treating a mostly pre-policy quarter as post-policy.

road_caz_props <- readRDS(here("data", "processed", "roads_caz_props.rds"))

scheme_start <- road_caz_props %>%
  distinct(scheme, startDt, caz_start_q) %>%
  filter(!is.na(startDt)) %>%
  mutate(
    start_date = dmy(startDt),
    raw_q = as.yearqtr(start_date),
    q_start = as.Date(raw_q),
    q_end = as.Date(raw_q + 0.25) - 1,
    q_mid = q_start + as.integer(difftime(q_end, q_start, units = "days")) / 2,
    caz_start_q_adj = if_else(start_date > q_mid, raw_q + 0.25, raw_q)
  ) %>%
  select(scheme, caz_start_q_adj)

road_panel <- road_panel %>%
  left_join(scheme_start, by = "scheme") %>%
  mutate(caz_start_q = coalesce(caz_start_q_adj, caz_start_q)) %>%
  select(-caz_start_q_adj)

scheme_timing <- road_panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q)

print(scheme_timing)

rm(road_caz_props, scheme_start)


# =============================================================================
# 5. LOAD OA COVARIATES
# =============================================================================
# These are time-invariant covariates used only for C&S robustness adjustment.

matched_covars <- readRDS(
  here("data", "processed", "OA_matched_full_pooled.rds")
) %>%
  mutate(
    log1p_pop_density = log1p(pmax(pop_density, 0)),
    log1p_road_density_m_km2 = log1p(pmax(road_density_m_km2, 0))
  ) %>%
  select(OA, log1p_pop_density, IMD, log1p_road_density_m_km2) %>%
  distinct(OA, .keep_all = TRUE)


# =============================================================================
# 6. BUILD MODEL PANEL
# =============================================================================
# Create the main analysis panel:
#   - one road-link × scheme × quarter row
#   - treatment timing
#   - event-time variables
#   - COVID treated-road interactions

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

model_panel <- road_panel %>%
  select(
    panel_id, identifier, OA, scheme, quarter_year,
    caz_start_q, treat_group, all_of(outcome_var)
  ) %>%
  rename(outcome_raw = all_of(outcome_var)) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  left_join(matched_covars, by = "OA") %>%
  mutate(
    uid = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    ),
    
    post_flag = as.integer(treat_group == 1 & quarter_year >= ref_start),
    
    covid_lockdown = as.integer(
      quarter_year >= as.yearqtr("2020 Q1") &
        quarter_year <= as.yearqtr("2021 Q1")
    ),
    
    covid_recovery = as.integer(
      quarter_year >= as.yearqtr("2021 Q2") &
        quarter_year <= as.yearqtr("2021 Q4")
    ),
    
    covid_lockdown_treated = covid_lockdown * treat_group,
    covid_recovery_treated = covid_recovery * treat_group,
    
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  ) %>%
  filter(!is.na(outcome_raw))

rm(road_panel)

schemes_all <- sort(unique(model_panel$scheme))


# =============================================================================
# 7. BASIC CHECKS
# =============================================================================
# Confirm sample size, sparsity, treatment/control counts, and duplicated
# controls from matching multiplicity.

summary_by_group <- model_panel %>%
  group_by(group) %>%
  summarise(
    units = n_distinct(uid_int),
    observations = n(),
    mean_injury = mean(outcome_raw),
    pct_zero = 100 * mean(outcome_raw == 0),
    .groups = "drop"
  )

print(summary_by_group)

summary_by_scheme <- model_panel %>%
  distinct(scheme, uid_int, treat_group) %>%
  group_by(scheme) %>%
  summarise(
    treated_units = sum(treat_group == 1),
    control_units = sum(treat_group == 0),
    .groups = "drop"
  )

print(summary_by_scheme)

duplicate_controls <- model_panel %>%
  count(identifier, scheme, quarter_year, treat_group, name = "n_rows") %>%
  filter(n_rows > 1)

duplicate_summary <- duplicate_controls %>%
  group_by(treat_group) %>%
  summarise(
    duplicated_cells = n(),
    mean_reuse = mean(n_rows),
    max_reuse = max(n_rows),
    .groups = "drop"
  )

print(duplicate_summary)


# =============================================================================
# 8. RAW TREND PLOT
# =============================================================================
# Shows treated/control raw injury trends and highlights COVID periods.

trend_df <- model_panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop")

p_trend <- ggplot(trend_df, aes(x = as.Date(quarter_year), y = mean_injury, colour = group)) +
  annotate(
    "rect",
    xmin = as.Date(as.yearqtr("2020 Q2")),
    xmax = as.Date(as.yearqtr("2021 Q1") + 0.25),
    ymin = -Inf, ymax = Inf,
    alpha = 0.08,
    fill = "grey50"
  ) +
  annotate(
    "rect",
    xmin = as.Date(as.yearqtr("2021 Q2")),
    xmax = as.Date(as.yearqtr("2022 Q1") + 0.25),
    ymin = -Inf, ymax = Inf,
    alpha = 0.06,
    fill = "grey70"
  ) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Mean injuries over time",
    subtitle = "Shaded periods: COVID lockdown and recovery",
    x = NULL,
    y = "Mean injuries per road-link-quarter",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_trend)
ggsave(file.path(outdir, "01_raw_trends.png"), p_trend, width = 10, height = 6, dpi = 300)


# =============================================================================
# 9. BUILD STACKED EVENT-STUDY DATA
# =============================================================================
# For each scheme:
#   - keep treated roads and matched controls for that scheme
#   - restrict to K quarters before and L quarters after implementation
#   - define event time relative to implementation quarter

stacked <- map_dfr(schemes_all, function(sc) {
  
  sc_start <- scheme_timing %>%
    filter(scheme == sc) %>%
    pull(caz_start_q)
  
  if (length(sc_start) == 0 || is.na(sc_start)) return(NULL)
  
  sc_start_int <- as.integer(round((as.numeric(sc_start) - min_qtr) * 4)) + 1L
  
  model_panel %>%
    filter(
      scheme == sc,
      qtr_int >= sc_start_int - K,
      qtr_int <= sc_start_int + L
    ) %>%
    mutate(
      stack_scheme = sc,
      event_time = qtr_int - sc_start_int,
      event_time_f = relevel(factor(event_time), ref = "-1"),
      uid_stack = paste0(uid_int, "_", sc),
      post_stack = as.integer(treat_group == 1 & event_time >= 0)
    )
})

cat("Stacked rows:", nrow(stacked), "\n")
cat("Stacked units:", n_distinct(stacked$uid_stack), "\n")


# =============================================================================
# 10. PRIMARY MODEL SELECTION: WITHOUT VS WITH COVID INTERACTIONS
# =============================================================================
# We first estimate the stacked PPML model without COVID interaction terms.
# Then we estimate the same model with treated × COVID lockdown/recovery terms.
# The comparison checks whether the CAZ estimate is sensitive to differential
# COVID disruption among treated roads.

# Model A: no COVID interactions
m_ppml_no_covid <- feglm(
  outcome_raw ~ post_stack |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

# Model B: with COVID interactions
m_ppml_covid <- feglm(
  outcome_raw ~ post_stack +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

# Compare models in regression table
etable(
  m_ppml_no_covid,
  m_ppml_covid,
  headers = c("No COVID interactions", "With COVID interactions"),
  dict = c(
    "post_stack" = "CAZ post-treatment",
    "covid_lockdown_treated" = "Treated × COVID lockdown",
    "covid_recovery_treated" = "Treated × COVID recovery"
  )
)

# Extract comparable IRRs
compare_covid_models <- tibble(
  model = c("No COVID interactions", "With COVID interactions"),
  estimate_log_irr = c(
    coef(m_ppml_no_covid)["post_stack"],
    coef(m_ppml_covid)["post_stack"]
  ),
  se = c(
    se(m_ppml_no_covid)["post_stack"],
    se(m_ppml_covid)["post_stack"]
  )
) %>%
  mutate(
    ci_lo = estimate_log_irr - 1.96 * se,
    ci_hi = estimate_log_irr + 1.96 * se,
    irr = exp(estimate_log_irr),
    irr_lo = exp(ci_lo),
    irr_hi = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo = 100 * (irr_lo - 1),
    pct_hi = 100 * (irr_hi - 1)
  )

print(compare_covid_models)

write_csv(
  compare_covid_models,
  file.path(outdir, "covid_model_comparison.csv")
)

# Plot comparison
p_covid_compare <- ggplot(
  compare_covid_models,
  aes(x = pct_change, y = fct_reorder(model, pct_change))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Sensitivity to COVID adjustment",
    subtitle = "Stacked PPML estimates with and without treated-road COVID interactions",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_covid_compare)

ggsave(
  file.path(outdir, "covid_model_comparison.png"),
  p_covid_compare,
  width = 8,
  height = 4.5,
  dpi = 300
)

# Use the COVID-adjusted model as the primary model for the rest of the analysis.
# This is preferred because the trend plot showed treated roads had a different
# COVID lockdown/recovery pattern from controls.
m_ppml <- m_ppml_covid

overall_ppml <- tibble(
  model = "Stacked PPML, COVID-adjusted",
  estimate_log_irr = coef(m_ppml)["post_stack"],
  se = se(m_ppml)["post_stack"]
) %>%
  mutate(
    ci_lo = estimate_log_irr - 1.96 * se,
    ci_hi = estimate_log_irr + 1.96 * se,
    irr = exp(estimate_log_irr),
    irr_lo = exp(ci_lo),
    irr_hi = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo = 100 * (irr_lo - 1),
    pct_hi = 100 * (irr_hi - 1)
  )

print(overall_ppml)

# =============================================================================
# 11. PRIMARY MODEL: OVERALL EVENT STUDY
# =============================================================================
# Estimates dynamic effects by event time, averaged across all schemes.
# Negative event times test anticipation/pre-trends.
# Positive event times show how the effect evolves after implementation.

m_ppml_es <- feglm(
  outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

ppml_es <- extract_fixest_event_study(m_ppml_es)

p_ppml_es <- plot_event_study(
  ppml_es,
  title = "Stacked PPML event study",
  subtitle = "Log incidence rate ratios; reference = quarter -1",
  ylab = "Log incidence rate ratio"
)

print(p_ppml_es)
ggsave(file.path(outdir, "02_ppml_event_study_overall.png"), p_ppml_es, width = 10, height = 7, dpi = 300)


# =============================================================================
# 12. SCHEME-SPECIFIC EVENT STUDIES
# =============================================================================
# Runs the same PPML event-study separately for each CAZ scheme.
# This is memory-efficient and directly answers whether effects differ by scheme.

run_scheme_ppml <- function(sc) {
  
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
        covid_lockdown_treated + covid_recovery_treated |
        uid_stack + qtr_int,
      data = d,
      family = "poisson",
      cluster = ~OA,
      lean = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  extract_fixest_event_study(fit) %>%
    mutate(scheme = sc)
}

scheme_es <- map_dfr(schemes_all, run_scheme_ppml)

p_scheme_es <- ggplot(scheme_es, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~scheme, scales = "free_y", ncol = 2) +
  labs(
    title = "Scheme-specific PPML event studies",
    subtitle = "Each panel estimates dynamic effects for one CAZ scheme",
    x = "Quarters relative to CAZ implementation",
    y = "Log incidence rate ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

print(p_scheme_es)
ggsave(file.path(outdir, "03_ppml_event_study_by_scheme.png"), p_scheme_es, width = 12, height = 10, dpi = 300)


# =============================================================================
# 13. SCHEME-SPECIFIC SUMMARY EFFECTS
# =============================================================================
# Summarises each scheme's average post-treatment effect from the scheme-specific
# event-study coefficients.

scheme_summary <- scheme_es %>%
  filter(event_time >= 0) %>%
  group_by(scheme) %>%
  summarise(
    mean_log_irr = mean(estimate, na.rm = TRUE),
    se_log_irr = sd(estimate, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ci_lo = mean_log_irr - 1.96 * se_log_irr,
    ci_hi = mean_log_irr + 1.96 * se_log_irr,
    irr = exp(mean_log_irr),
    irr_lo = exp(ci_lo),
    irr_hi = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo = 100 * (irr_lo - 1),
    pct_hi = 100 * (irr_hi - 1)
  )

print(scheme_summary)

p_scheme_summary <- ggplot(scheme_summary, aes(x = pct_change, y = fct_reorder(scheme, pct_change))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Average post-CAZ effect by scheme",
    subtitle = "Based on average post-treatment PPML event-study coefficients",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_scheme_summary)
ggsave(file.path(outdir, "04_scheme_summary_effects.png"), p_scheme_summary, width = 9, height = 6, dpi = 300)




# =============================================================================
# 14. C&S ROBUSTNESS CHECK
# =============================================================================
# C&S is used as a robustness check on the additive injury-count scale.
# Zero-only road-link series are removed because they provide no identifying
# information for this estimator.

cs_data <- model_panel %>%
  group_by(uid_int) %>%
  filter(sum(outcome_raw, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  mutate(uid_int = as.integer(factor(uid_int)))

# =============================================================================
# FORMAL SCHEME-SPECIFIC STACKED PPML EFFECTS
# =============================================================================
# This estimates one average post-treatment effect for each CAZ scheme.
# Coefficients are log incidence rate ratios.
# exp(coef) gives the IRR; (exp(coef)-1)*100 gives % change in injuries.

stacked <- stacked %>%
  mutate(
    stack_scheme = factor(stack_scheme),
    scheme_post = if_else(treat_group == 1 & post_stack == 1,
                          as.character(stack_scheme),
                          "control"),
    scheme_post = factor(scheme_post, levels = c("control", schemes_all))
  )

m_ppml_scheme <- feglm(
  outcome_raw ~ i(scheme_post, ref = "control") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  lean    = TRUE
)

etable(m_ppml_scheme)

cs_fit <- run_cs(
  cs_data,
  xformla = ~ log1p_pop_density + IMD + log1p_road_density_m_km2
)

summary(cs_fit$agg)
summary(cs_fit$att)

cs_es <- extract_cs_event_study(cs_fit$dyn)

p_cs_es <- plot_event_study(
  cs_es,
  title = "Callaway & Sant'Anna robustness event study",
  subtitle = "ATT on additive injury-count scale",
  ylab = "ATT, injury count",
  colour = "#2E6FAB"
)

print(p_cs_es)
ggsave(file.path(outdir, "05_cs_event_study.png"), p_cs_es, width = 10, height = 7, dpi = 300)


# =============================================================================
# SCHEME-SPECIFIC STACKED PPML EFFECTS
# =============================================================================
# This estimates one average post-treatment effect for each CAZ scheme.
# Coefficients are log incidence rate ratios.
# exp(coef) gives the IRR; (exp(coef)-1)*100 gives % change in injuries.

stacked <- stacked %>%
  mutate(
    stack_scheme = factor(stack_scheme),
    scheme_post = if_else(treat_group == 1 & post_stack == 1,
                          as.character(stack_scheme),
                          "control"),
    scheme_post = factor(scheme_post, levels = c("control", schemes_all))
  )

m_ppml_scheme <- feglm(
  outcome_raw ~ i(scheme_post, ref = "control") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  lean    = TRUE
)

etable(m_ppml_scheme)



# =============================================================================
# EXTRACT SCHEME-SPECIFIC IRRs
# =============================================================================

ct <- coeftable(m_ppml_scheme)

scheme_ppml_table <- tibble(
  term = rownames(ct),
  estimate_log_irr = ct[, "Estimate"],
  se = ct[, "Std. Error"],
  p_value = ct[, "Pr(>|z|)"]
) %>%
  filter(str_detect(term, "scheme_post::")) %>%
  mutate(
    scheme = str_remove(term, "scheme_post::"),
    ci_lo = estimate_log_irr - 1.96 * se,
    ci_hi = estimate_log_irr + 1.96 * se,
    irr = exp(estimate_log_irr),
    irr_lo = exp(ci_lo),
    irr_hi = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo = 100 * (irr_lo - 1),
    pct_hi = 100 * (irr_hi - 1),
    result = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE ~ ""
    )
  ) %>%
  select(
    scheme,
    estimate_log_irr,
    se,
    irr,
    irr_lo,
    irr_hi,
    pct_change,
    pct_lo,
    pct_hi,
    p_value,
    result
  ) %>%
  arrange(pct_change)

print(scheme_ppml_table)


#Sheffield changed most:
# the above table averaged the post-event-study coefficients, which can be distorted by noisy event-time estimates 
#and does not account properly for covariance across event-time coefficients. 
#This directly estimates the average post-period contrast, somore stable.

# =============================================================================
# 15. FINAL COMPARISON TABLE
# =============================================================================
# Presents the main PPML estimate and C&S robustness estimate side by side.

cs_summary <- tibble(
  model = "Callaway & Sant'Anna",
  estimate = cs_fit$agg$overall.att,
  se = cs_fit$agg$overall.se,
  ci_lo = estimate - 1.96 * se,
  ci_hi = estimate + 1.96 * se,
  scale = "Additive injury-count scale"
)

ppml_summary <- overall_ppml %>%
  transmute(
    model,
    estimate = estimate_log_irr,
    se,
    ci_lo,
    ci_hi,
    scale = "Log-IRR; exponentiate for IRR"
  )

final_results <- bind_rows(ppml_summary, cs_summary)

print(final_results)

write_csv(final_results, file.path(outdir, "final_results_overall.csv"))
write_csv(scheme_summary, file.path(outdir, "final_results_by_scheme.csv"))
write_csv(ppml_es, file.path(outdir, "ppml_event_study_overall.csv"))
write_csv(scheme_es, file.path(outdir, "ppml_event_study_by_scheme.csv"))
write_csv(cs_es, file.path(outdir, "cs_event_study.csv"))


# =============================================================================
# 16. FINAL COMBINED FIGURE
# =============================================================================
# Combines primary PPML and C&S robustness event-study plots.

p_final <- p_ppml_es / p_cs_es +
  plot_annotation(
    title = "CAZ effects on road traffic injuries",
    subtitle = "Primary stacked PPML model and C&S robustness check",
    caption = "PPML uses log incidence rate ratios; C&S uses additive injury-count ATT."
  )

print(p_final)
ggsave(file.path(outdir, "06_final_combined_event_studies.png"), p_final, width = 10, height = 12, dpi = 300)


# =============================================================================
# 17. SAVE CORE DATA AND RESULTS
# =============================================================================
# Saves only what is needed for reproducibility and reporting.

arrow::write_parquet(
  model_panel,
  here("data", "processed", "final_model_data_clean_All.parquet")
)

saveRDS(
  list(
    overall_ppml = overall_ppml,
    ppml_event_study = ppml_es,
    scheme_event_study = scheme_es,
    scheme_summary = scheme_summary,
    cs_event_study = cs_es,
    cs_summary = cs_summary,
    final_results = final_results
  ),
  here("data", "processed", "caz_did_results_clean_All.rds")
)

cat("\nAnalysis complete.\n")
cat("Outputs saved to:", outdir, "\n")