# =============================================================================
# CAZ DiD ANALYSIS — MAIN PPML ANALYSIS SCRIPT
# =============================================================================
#
# Main framework:
#   Matched stacked Difference-in-Differences estimated using PPML.
#
# Model structure:
#   1. Pooled average treatment effect across all schemes
#   2. Pooled dynamic event-study effects and parallel-trends assessment
#   3. Scheme-specific average treatment effects
#   4. Scheme-specific dynamic event-study effects
#
# Notes:
#   - Outcome is a sparse count, so models are estimated with Poisson
#     pseudo-maximum likelihood (PPML).
#   - Standard errors are clustered at OA level.
#   - COVID lockdown/recovery interactions are included in the main models.
#   - Robustness, sensitivity, C&S checks, sparsity diagnostics, and alternative
#     clustering checks should be kept in a separate script.
# =============================================================================
library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(lubridate)
library(fixest)
library(patchwork)

outcome_var <- "total_inj_adj_All"

K <- 8L   # quarters before CAZ implementation
L <- 8L   # quarters after CAZ implementation

outdir <- here("output", "pooled", "All_clean")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
# Extract event-study coefficients from fixest objects without using broom,
# then convert log-rate coefficients to IRRs and percentage changes.

extract_fixest_event_study <- function(model) {
  ct <- coeftable(model)
  
  tibble(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"]
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

# =============================================================================
# data
# =============================================================================
# Load the matched road-link × quarter panel and apply a majority-quarter rule:
# if a CAZ starts in the second half of a quarter, treatment starts next quarter.

road_panel <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q = as.yearqtr(caz_start_q)
  )

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

cat("\nScheme timing after majority-quarter adjustment:\n")
print(scheme_timing)

rm(road_caz_props, scheme_start)

# =============================================================================
# OA COVARIATES
# =============================================================================
# Time-invariant OA covariates are retained in the model panel so they are
# available for robustness analyses, especially C&S checks in the separate script.

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
#  BUILD MODEL PANEL
# =============================================================================
# Construct the road-link × scheme × quarter panel for modelling.
# COVID periods are coded from visual inspection:
#   Lockdown/disruption: 2020 Q1–2021 Q1
#   Recovery:            2021 Q2–2021 Q3

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

    # -----------------------------------------------------------------
    # COVID period: single three-level factor
    #   pre_covid  = before Q2 2020             (reference)
    #   lockdown   = Q2 2020 – Q1 2021          (first full lockdown quarter
    #                                            through end of third lockdown)
    #   recovery   = Q2 2021 – Q3 2021          (restrictions mostly lifted
    #                                            by August 2021)
    # Anything from Q4 2021 onward = pre_covid  (back to baseline)
    # -----------------------------------------------------------------

model_panel <- road_panel %>%
  select(
    panel_id, identifier, OA, scheme, quarter_year,
    caz_start_q, treat_group, all_of(outcome_var)
  ) %>%
  rename(outcome_raw = all_of(outcome_var)) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  left_join(matched_covars, by = "OA") %>%
  mutate(
    uid     = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    ),
    
    post_flag = as.integer(treat_group == 1 & quarter_year >= ref_start),
    
    covid_period = case_when(
      quarter_year >= as.yearqtr("2020 Q2") &
        quarter_year <= as.yearqtr("2021 Q1") ~ "lockdown",
      
      quarter_year >= as.yearqtr("2021 Q2") &
        quarter_year <= as.yearqtr("2021 Q3") ~ "recovery",
      
      TRUE ~ "pre_covid"
    ),
    
    covid_period = factor(
      covid_period,
      levels = c("pre_covid", "lockdown", "recovery")
    ),
    
    covid_lockdown_treated = as.integer(covid_period == "lockdown") * treat_group,
    covid_recovery_treated = as.integer(covid_period == "recovery") * treat_group,
    
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  ) %>%
  filter(!is.na(outcome_raw))
rm(road_panel)

schemes_all <- sort(unique(model_panel$scheme))

# =============================================================================
# BASIC SAMPLE CHECKS
# =============================================================================
# Confirm sample size, sparsity, number of OAs/road links, treatment/control
# counts, and duplicate control rows arising from matching with replacement.

summary_by_group <- model_panel %>%
  group_by(group) %>%
  summarise(
    units = n_distinct(uid_int),
    observations = n(),
    mean_injury = mean(outcome_raw),
    pct_zero = 100 * mean(outcome_raw == 0),
    .groups = "drop"
  )

cat("\nSummary by group:\n")
print(summary_by_group)

cat("\nNumber of OAs:", n_distinct(model_panel$OA), "\n")
cat("Number of road links:", n_distinct(model_panel$identifier), "\n")

summary_by_scheme <- model_panel %>%
  distinct(scheme, uid_int, treat_group) %>%
  group_by(scheme) %>%
  summarise(
    treated_units = sum(treat_group == 1),
    control_units = sum(treat_group == 0),
    .groups = "drop"
  )

cat("\nSummary by scheme:\n")
print(summary_by_scheme)

duplicate_summary <- model_panel %>%
  count(identifier, scheme, quarter_year, treat_group, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  group_by(treat_group) %>%
  summarise(
    duplicated_cells = n(),
    mean_reuse = mean(n_rows),
    max_reuse = max(n_rows),
    .groups = "drop"
  )

cat("\nDuplicate rows from matched-control reuse:\n")
print(duplicate_summary)

# =============================================================================
# RAW TREND PLOT
# =============================================================================
# Plot raw treated/control trends and shade the COVID periods used in models.

trend_df <- model_panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop")

p_trend <- ggplot(trend_df, aes(x = as.Date(quarter_year), y = mean_injury, colour = group)) +
  annotate(
    "rect",
    xmin = as.Date(as.yearqtr("2020 Q2")),
    xmax = as.Date(as.yearqtr("2021 Q3") + 0.25),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08,
    fill = "grey70"
  ) +
  annotate(
    "rect",
    xmin = as.Date(as.yearqtr("2021 Q2")),
    xmax = as.Date(as.yearqtr("2021 Q4") + 0.25),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.06,
    fill = "grey70"
  ) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Mean injuries over time",
    subtitle = "Shaded periods: COVID lockdown/disruption and recovery",
    x = NULL,
    y = "Mean injuries per road-link-quarter",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_trend)

ggsave(
  file.path(outdir, "00_raw_trends.png"),
  p_trend,
  width = 10,
  height = 6,
  dpi = 300
)

# =============================================================================
# BUILD STACKED DATA
# =============================================================================
# Each CAZ scheme becomes a separate matched DiD stack. The same control road
# can appear in more than one stack, but receives a stack-specific unit ID.

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
}) %>%
  mutate(stack_scheme = factor(stack_scheme))

cat("\nStacked rows:", nrow(stacked), "\n")
cat("Stacked units:", n_distinct(stacked$uid_stack), "\n")

###    stackes fixed effects 
stacked <- stacked %>%
  mutate(stack_time_fe = interaction(stack_scheme, qtr_int, drop = TRUE))


# =============================================================================
# MODEL 1 — POOLED AVERAGE TREATMENT EFFECT
# =============================================================================
# This model estimates one average post-treatment CAZ effect across all schemes.
# A no-COVID version is fitted only to show sensitivity to COVID adjustment.

m1_no_covid <- feglm(
  outcome_raw ~ post_stack |
    uid_stack + + stack_time_fe,
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

m1_covid <- feglm(
  outcome_raw ~ post_stack +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + + stack_time_fe,
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

cat("\nModel 1: pooled average effect, with and without COVID interactions\n")
etable(
  m1_no_covid,
  m1_covid,
  headers = c("No COVID interactions", "COVID-adjusted"),
  dict = c(
    "post_stack" = "CAZ post-treatment",
    "covid_lockdown_treated" = "Treated × COVID lockdown",
    "covid_recovery_treated" = "Treated × COVID recovery"
  )
)

model1_results <- tibble(
  model = c("No COVID interactions", "COVID-adjusted"),
  estimate_log_irr = c(
    coef(m1_no_covid)["post_stack"],
    coef(m1_covid)["post_stack"]
  ),
  se = c(
    se(m1_no_covid)["post_stack"],
    se(m1_covid)["post_stack"]
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

print(model1_results)

p_model1 <- ggplot(
  model1_results,
  aes(x = pct_change, y = fct_reorder(model, pct_change))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Model 1: pooled average CAZ effect",
    subtitle = "Stacked PPML estimates with and without treated-road COVID interactions",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_model1)

ggsave(
  file.path(outdir, "01_model1_pooled_average_effect.png"),
  p_model1,
  width = 8,
  height = 4.5,
  dpi = 300
)

# Main headline estimate from Model 1.
m1_main <- m1_covid
model1_main_result <- model1_results %>% filter(model == "COVID-adjusted")

# =============================================================================
#  MODEL 2 — POOLED DYNAMIC EVENT STUDY
# =============================================================================
# This is the primary dynamic model. It estimates pooled lead and lag effects
# relative to quarter -1 and is used for parallel-trends assessment.

m2_event <- feglm(
  outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + + stack_time_fe,
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

model2_results <- extract_fixest_event_study(m2_event)

p_model2 <- plot_event_study(
  model2_results,
  title = "Model 2: pooled stacked PPML event study",
  subtitle = "Primary dynamic model; reference period = quarter -1",
  ylab = "Log incidence rate ratio"
)

print(p_model2)

ggsave(
  file.path(outdir, "02_model2_pooled_event_study.png"),
  p_model2,
  width = 10,
  height = 7,
  dpi = 300
)

# Parallel-trends diagnostics from the primary PPML event-study.
# The closest pre-treatment window is the most relevant diagnostic.
pt_8_2 <- wald(
  m2_event,
  keep = "event_time_f::-(8|7|6|5|4|3|2):treat_group"
)

pt_6_2 <- wald(
  m2_event,
  keep = "event_time_f::-(6|5|4|3|2):treat_group"
)

pt_4_2 <- wald(
  m2_event,
  keep = "event_time_f::-(4|3|2):treat_group"
)

cat("\nParallel-trends tests from Model 2:\n")
print(pt_8_2)
print(pt_6_2)
print(pt_4_2)

# =============================================================================
# MODEL 3 — SCHEME-SPECIFIC AVERAGE TREATMENT EFFECTS
# =============================================================================
# This model estimates one formal average post-treatment effect for each scheme.

stacked <- stacked %>%
  mutate(
    scheme_post = if_else(
      treat_group == 1 & post_stack == 1,
      as.character(stack_scheme),
      "control"
    ),
    scheme_post = factor(scheme_post, levels = c("control", schemes_all))
  )

m3_scheme <- feglm(
  outcome_raw ~ i(scheme_post, ref = "control") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + stack_time_fe,
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  lean = TRUE
)

cat("\nModel 3: scheme-specific average effects\n")
etable(m3_scheme)

ct_m3 <- coeftable(m3_scheme)

model3_results <- tibble(
  term = rownames(ct_m3),
  estimate_log_irr = ct_m3[, "Estimate"],
  se = ct_m3[, "Std. Error"],
  p_value = ct_m3[, "Pr(>|z|)"]
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
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.10 ~ ".",
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
    sig
  ) %>%
  arrange(pct_change)

print(model3_results)

p_model3 <- ggplot(
  model3_results,
  aes(x = pct_change, y = fct_reorder(scheme, pct_change))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Model 3: scheme-specific average CAZ effects",
    subtitle = "Formal scheme-specific stacked PPML estimates",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_model3)

ggsave(
  file.path(outdir, "03_model3_scheme_average_effects.png"),
  p_model3,
  width = 9,
  height = 6,
  dpi = 300
)

# =============================================================================
# MODEL 4 — SCHEME-SPECIFIC DYNAMIC EVENT STUDIES
# =============================================================================
# This runs the event-study model separately for each scheme. These estimates
# show whether the timing and shape of effects differ across CAZ schemes.

run_scheme_event_ppml <- function(sc) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
        covid_lockdown_treated + covid_recovery_treated |
        uid_stack + stack_time_fe,
      data = d,
      family = "poisson",
      cluster = ~OA,
      lean = TRUE
    ),
    error = function(e) {
      cat("Scheme event study failed for:", sc, "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  extract_fixest_event_study(fit) %>%
    mutate(scheme = sc)
}

model4_results <- map_dfr(schemes_all, run_scheme_event_ppml)

p_model4 <- ggplot(model4_results, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~scheme, scales = "free_y", ncol = 2) +
  labs(
    title = "Model 4: scheme-specific PPML event studies",
    subtitle = "Dynamic effects estimated separately for each CAZ scheme",
    x = "Quarters relative to CAZ implementation",
    y = "Log incidence rate ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_model4)

ggsave(
  file.path(outdir, "04_model4_scheme_event_studies.png"),
  p_model4,
  width = 12,
  height = 10,
  dpi = 300
)

# =============================================================================
#  OUTPUTS
# =============================================================================
write_csv(model1_results, file.path(outdir, "model1_pooled_average_effect.csv"))
write_csv(model2_results, file.path(outdir, "model2_pooled_event_study.csv"))
write_csv(model3_results, file.path(outdir, "model3_scheme_average_effects.csv"))
write_csv(model4_results, file.path(outdir, "model4_scheme_event_studies.csv"))

arrow::write_parquet(
  model_panel,
  here("data", "processed", "final_model_panel_main.parquet")
)

arrow::write_parquet(
  stacked,
  here("data", "processed", "final_stacked_data_main.parquet")
)

saveRDS(
  list(
    model1_results = model1_results,
    model1_main_result = model1_main_result,
    model2_results = model2_results,
    model3_results = model3_results,
    model4_results = model4_results,
    pt_8_2 = pt_8_2,
    pt_6_2 = pt_6_2,
    pt_4_2 = pt_4_2
  ),
  here("data", "processed", "caz_ppml_main_results.rds")
)

cat("\nMAIN PPML ANALYSIS COMPLETE.\n")
cat("Outputs saved to:", outdir, "\n")
