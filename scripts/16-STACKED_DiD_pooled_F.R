# =============================================================================
# CAZ DiD ANALYSIS - REORGANISED MAIN PPML SCRIPT
# =============================================================================
#
# Framework:
#   Matched stacked Difference-in-Differences estimated using PPML.
#
# Main estimand:
#   Equal-weighted average effect across CAZ schemes. Within each scheme,
#   treated road-link units sum to 1 and matched controls sum to 1.
#
# Main model setup:
#   - Outcome: total_inj_adj_All
#   - Model family: PPML
#   - Unit fixed effects: road-link x scheme stack
#   - Time adjustment: scheme-specific linear quarter trend
#   - SEs: clustered by OA
#   - COVID adjustment: flexible scheme-by-treatment lockdown/recovery terms
#   - Main event-study reference: year before implementation, event times -4:-1
#
# Script organisation:
#   SECTION 1: Main analysis
#     1. Pooled average effect
#     2. Pooled dynamic event study, year-reference model
#     3. Scheme-specific average effects
#     4. Scheme-specific dynamic event studies
#
#   SECTION 2: Supporting diagnostics and sensitivity checks
#     A. Zero-exclusion diagnostics
#     B. Raw trends and sample checks
#     C. Alternative COVID/reference/weight specifications
#     D. Parallel-trends diagnostics
#     E. Event-time composition diagnostics
#     F. Bradford diagnostic investigation
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(lubridate)
library(fixest)
library(patchwork)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
count  <- dplyr::count

outcome_var <- "total_inj_adj_All"
outdir <- here("output", "pooled", "All_clean")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

add_irr_columns <- function(df, estimate_col = "estimate", se_col = "se") {
  df %>%
    mutate(
      ci_lo      = .data[[estimate_col]] - 1.96 * .data[[se_col]],
      ci_hi      = .data[[estimate_col]] + 1.96 * .data[[se_col]],
      irr        = exp(.data[[estimate_col]]),
      irr_lo     = exp(ci_lo),
      irr_hi     = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo     = 100 * (irr_lo - 1),
      pct_hi     = 100 * (irr_hi - 1)
    )
}

extract_event_study <- function(model, var_prefix) {
  ct <- coeftable(model)
  
  tibble(
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"]
  ) %>%
    filter(str_detect(term, paste0("^", var_prefix, "::"))) %>%
    mutate(
      event_time = str_match(
        term,
        paste0("^", var_prefix, "::(-?\\d+):treat_group$")
      )[, 2] %>% as.numeric()
    ) %>%
    filter(!is.na(event_time)) %>%
    add_irr_columns() %>%
    arrange(event_time)
}

extract_scheme_effects <- function(model, var_prefix = "scheme_post") {
  ct <- coeftable(model)
  tibble(
    term             = rownames(ct),
    estimate_log_irr = ct[, "Estimate"],
    se               = ct[, "Std. Error"],
    p_value          = ct[, "Pr(>|z|)"]
  ) %>%
    filter(str_detect(term, paste0("^", var_prefix, "::"))) %>%
    mutate(scheme = str_remove(term, paste0("^", var_prefix, "::"))) %>%
    add_irr_columns("estimate_log_irr", "se") %>%
    mutate(sig = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*", p_value < 0.10 ~ ".", TRUE ~ ""
    )) %>%
    select(scheme, estimate_log_irr, se, irr, irr_lo, irr_hi,
           pct_change, pct_lo, pct_hi, p_value, sig) %>%
    arrange(pct_change)
}


plot_event_study <- function(df, title, subtitle, ylab, colour = "#E74C3C") {
  ggplot(df, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = colour) +
    geom_line(linewidth = 0.8, colour = colour) +
    geom_point(size = 1.8, colour = colour) +
    scale_x_continuous(breaks = pretty(df$event_time, n = 10)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Quarters relative to CAZ implementation",
      y = ylab
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


###### 
make_clean_reference <- function(stacked_data) {
  # If the conventional year-reference period (-4:-1) overlaps COVID/recovery,
  # use the common pre-COVID year 2019 Q2-2020 Q1 as a clean reference instead.
  stacked_data %>%
    mutate(
      is_covid = covid_period %in% c("lockdown", "recovery"),
      pre_covid_ref = quarter_year >= as.yearqtr("2019 Q2") &
        quarter_year <= as.yearqtr("2020 Q1")
    ) %>%
    group_by(stack_scheme) %>%
    mutate(
      normal_ref = event_time >= -4 & event_time <= -1,
      normal_ref_has_covid = any(normal_ref & is_covid, na.rm = TRUE),
      clean_ref_year = case_when(
        !normal_ref_has_covid & normal_ref ~ TRUE,
        normal_ref_has_covid & pre_covid_ref & event_time < 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      event_time_ref_clean = if_else(
        clean_ref_year,
        "ref_year",
        as.character(event_time)
      )
    ) %>%
    ungroup() %>%
    mutate(
      event_time_ref_clean = relevel(factor(event_time_ref_clean), ref = "ref_year")
    )
}

# =============================================================================
# DATA PREPARATION
# =============================================================================

road_panel <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

road_caz_props <- readRDS(here("data", "processed", "roads_caz_props.rds"))

scheme_start <- road_caz_props %>%
  distinct(scheme, startDt, caz_start_q) %>%
  filter(!is.na(startDt)) %>%
  mutate(
    start_date      = dmy(startDt),
    raw_q           = as.yearqtr(start_date),
    q_start         = as.Date(raw_q),
    q_end           = as.Date(raw_q + 0.25) - 1,
    q_mid           = q_start + as.integer(difftime(q_end, q_start, units = "days")) / 2,
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

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

model_panel <- road_panel %>%
  select(
    panel_id, identifier, OA, scheme, quarter_year,
    caz_start_q, treat_group, all_of(outcome_var)
  ) %>%
  rename(outcome_raw = all_of(outcome_var)) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  mutate(
    uid     = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    post_flag = as.integer(treat_group == 1 & quarter_year >= ref_start),
    covid_period = factor(
      case_when(
        quarter_year >= as.yearqtr("2020 Q2") &
          quarter_year <= as.yearqtr("2021 Q1") ~ "lockdown",
        quarter_year >= as.yearqtr("2021 Q2") &
          quarter_year <= as.yearqtr("2021 Q4") ~ "recovery",
        TRUE ~ "non_pandemic"
      ),
      levels = c("non_pandemic", "lockdown", "recovery")
    ),
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  ) %>%
  filter(!is.na(outcome_raw)) %>%
  group_by(uid) %>%
  mutate(unit_total_injury_all_periods = sum(outcome_raw, na.rm = TRUE)) %>%
  ungroup()

cat("\nAll-zero road-link units before exclusion:\n")
model_panel %>%
  distinct(uid, group, unit_total_injury_all_periods) %>%
  group_by(group) %>%
  summarise(
    units                 = n(),
    zero_injury_units     = sum(unit_total_injury_all_periods == 0),
    pct_zero_injury_units = 100 * mean(unit_total_injury_all_periods == 0),
    .groups = "drop"
  ) %>%
  print()

model_panel_for_zero_diag <- model_panel

model_panel <- model_panel %>%
  filter(unit_total_injury_all_periods > 0) %>%
  select(-unit_total_injury_all_periods)

cat("\nRows after excluding all-zero road-link units:", nrow(model_panel), "\n")
cat("Units after exclusion:", n_distinct(model_panel$uid), "\n")

rm(road_panel)

schemes_all <- sort(unique(model_panel$scheme))

stacked <- map_dfr(schemes_all, function(sc) {
  sc_start <- scheme_timing %>%
    filter(scheme == sc) %>%
    pull(caz_start_q)
  
  if (length(sc_start) == 0 || is.na(sc_start)) return(NULL)
  
  sc_start_int <- as.integer(round((as.numeric(sc_start) - min_qtr) * 4)) + 1L
  
  model_panel %>%
    filter(scheme == sc) %>%
    mutate(
      stack_scheme = sc,
      event_time   = qtr_int - sc_start_int,
      event_time_f = relevel(factor(event_time), ref = "-1"),
      uid_stack    = paste0(uid_int, "_", sc),
      post_stack   = as.integer(treat_group == 1 & event_time >= 0)
    )
}) %>%
  mutate(
    stack_scheme = factor(stack_scheme),
    event_time_ref = if_else(
      event_time >= -4 & event_time <= -1,
      "ref_year",
      as.character(event_time)
    ),
    event_time_ref = relevel(factor(event_time_ref), ref = "ref_year"),
    treat_scheme = interaction(treat_group, stack_scheme, drop = TRUE)
  )

cat("\nStacked rows:", nrow(stacked), "\n")
cat("Stacked units:", n_distinct(stacked$uid_stack), "\n")

cat("\nConventional reference quarter (event_time = -1) by scheme:\n")
stacked %>%
  filter(event_time == -1) %>%
  distinct(stack_scheme, quarter_year, covid_period) %>%
  arrange(stack_scheme) %>%
  print()

analysis_weight_counts <- stacked %>%
  distinct(stack_scheme, treat_group, uid_stack) %>%
  count(stack_scheme, treat_group, name = "n_units")

stacked <- stacked %>%
  left_join(analysis_weight_counts, by = c("stack_scheme", "treat_group")) %>%
  mutate(analysis_weight = 1 / n_units) %>%
  select(-n_units)

link_weight_counts <- analysis_weight_counts %>%
  pivot_wider(
    names_from = treat_group,
    values_from = n_units,
    names_prefix = "grp_"
  ) %>%
  mutate(control_link_weight = grp_1 / grp_0) %>%
  select(stack_scheme, control_link_weight)

stacked <- stacked %>%
  left_join(link_weight_counts, by = "stack_scheme") %>%
  mutate(link_weight = if_else(treat_group == 1, 1, control_link_weight)) %>%
  select(-control_link_weight)

cat("\nAnalysis weight verification at event_time = -1:\n")
stacked %>%
  filter(event_time == -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    n_rows = n(),
    sum_analysis_weight = sum(analysis_weight),
    sum_link_weight = sum(link_weight),
    .groups = "drop"
  ) %>%
  print(n = Inf)

stacked <- stacked %>%
  mutate(
    scheme_post = if_else(
      treat_group == 1 & post_stack == 1,
      as.character(stack_scheme),
      "control"
    ),
    scheme_post = factor(scheme_post, levels = c("control", schemes_all))
  )



stacked <- make_clean_reference(stacked)

clean_reference_table <- stacked %>%
  filter(clean_ref_year) %>%
  distinct(stack_scheme, quarter_year, event_time, covid_period) %>%
  arrange(stack_scheme, quarter_year)



# =============================================================================
# SECTION 1 - MAIN ANALYSIS
# =============================================================================

# -----------------------------------------------------------------------------
# Main Model 1: pooled average effect
# -----------------------------------------------------------------------------

m1_pooled_average <- feglm(
  outcome_raw ~
    post_stack +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model1_results <- tibble(
  model = "Main pooled average: equal-scheme, flexible scheme COVID",
  estimate_log_irr = coef(m1_pooled_average)["post_stack"],
  se = se(m1_pooled_average)["post_stack"]
) %>%
  add_irr_columns("estimate_log_irr", "se")

cat("\nMain Model 1: pooled average effect\n")
etable(m1_pooled_average, dict = c("post_stack" = "CAZ post-treatment"))
print(model1_results)

p_model1 <- ggplot(model1_results, aes(x = pct_change, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Main Model 1: pooled average CAZ effect",
    subtitle = "Equal-scheme weighted stacked PPML with flexible scheme-by-treatment COVID adjustment",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "01_main_pooled_average_effect.png"),
       p_model1, width = 9, height = 4.5, dpi = 300)



# Best-available Model 1: flexible scheme x treatment COVID control (already have)
# + clean pre-COVID reference for the "pre" comparison
# + explicit reporting alongside the scheme-decomposed Model 3, not as a substitute for it


# -----------------------------------------------------------------------------
# Main Model 2: pooled dynamic event study, year-reference model
# -----------------------------------------------------------------------------

m2_pooled_event_yearref <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model2_results <- extract_event_study(m2_pooled_event_yearref, "event_time_ref")

cat("\nMain Model 2: pooled year-reference event study\n")
etable(m2_pooled_event_yearref)

p_model2 <- plot_event_study(
  model2_results %>% filter(event_time >= -28, event_time <= 10),
  title = "Main Model 2: pooled stacked PPML event study",
  subtitle = "Equal-scheme weighted; reference = year before implementation (event times -4 to -1)",
  ylab = "Log incidence rate ratio"
)

ggsave(file.path(outdir, "02_main_pooled_yearref_event_study.png"),
       p_model2, width = 10, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Main Model 3: scheme-specific average effects
# -----------------------------------------------------------------------------

m3_scheme_average <- feglm(
  outcome_raw ~
    i(scheme_post, ref = "control") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model3_results <- extract_scheme_effects(m3_scheme_average)

cat("\nMain Model 3: scheme-specific average effects\n")
etable(m3_scheme_average)
print(model3_results)

p_model3 <- ggplot(
  model3_results,
  aes(x = pct_change, y = fct_reorder(scheme, pct_change))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(
    title = "Main Model 3: scheme-specific average CAZ effects",
    subtitle = "Stacked PPML; equal-scheme weights; flexible scheme-by-treatment COVID adjustment",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "03_main_scheme_average_effects.png"),
       p_model3, width = 9, height = 6, dpi = 300)


# ==============================================================================
# MODEL 3 (REVISED): scheme-specific average effects
# Fixes vs. original: (1) comparison anchored to clean pre-COVID reference period
# only, not an unweighted pool of the entire non-stationary pre-period + controls;
# (2) common fixed post-treatment window (0-5 quarters) so all 7 schemes are
# compared on equal footing (per post_event_scheme_counts, all 7 are present
# through event_time 5).
# ==============================================================================

COMMON_POST_MAX <- 5

stacked_m3 <- stacked %>%
  filter(clean_ref_year | (event_time >= 0 & event_time <= COMMON_POST_MAX)) %>%
  mutate(
    post_common        = as.integer(treat_group == 1 &
                                      event_time >= 0 & event_time <= COMMON_POST_MAX),
    scheme_post_common  = if_else(post_common == 1, as.character(stack_scheme), "control"),
    scheme_post_common  = factor(scheme_post_common, levels = c("control", schemes_all))
  )

m3_scheme_average_clean <- feglm(
  outcome_raw ~
    i(scheme_post_common, ref = "control") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked_m3,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model3_clean_results <- extract_scheme_effects(m3_scheme_average_clean, "scheme_post_common")
print(model3_clean_results, n = Inf)


cat("\nModel 3 (revised): scheme-specific effects, common 0-", COMMON_POST_MAX,
    "Q post window, clean pre-COVID reference\n", sep = "")
etable(m3_scheme_average_clean)
print(model3_clean_results, n = Inf)

# -----------------------------------------------------------------------------
# Direct before/after comparison, to document the correction 
# -----------------------------------------------------------------------------
model3_comparison <- bind_rows(
  model3_results %>% mutate(spec = "Original: unbalanced window, uncleaned reference"),
  model3_clean_results %>% mutate(spec = paste0("Revised: common 0-", COMMON_POST_MAX,
                                                "Q window, clean reference"))
) %>%
  select(spec, scheme, pct_change, pct_lo, pct_hi, p_value, sig) %>%
  arrange(scheme, spec)

print(model3_comparison, n = Inf)

p_model3_comparison <- ggplot(
  model3_comparison,
  aes(x = pct_change, y = fct_reorder(scheme, pct_change), colour = spec)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Scheme-specific average effects: original vs. corrected specification",
    subtitle = "Corrected: common post-treatment window + clean pre-COVID reference period",
    x = "% change in injuries", y = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "03b_model3_original_vs_clean_comparison.png"),
       p_model3_comparison, width = 10, height = 6, dpi = 300)

write_csv(model3_clean_results, file.path(outdir, "model3_revised_clean_common_window.csv"))
write_csv(model3_comparison, file.path(outdir, "model3_original_vs_revised_comparison.csv"))



#####    extend same correction to model 1
m1_pooled_average_clean <- feglm(
  outcome_raw ~
    post_common +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked_m3,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model1_clean_results <- tibble(
  model             = paste0("Pooled average (0-", COMMON_POST_MAX,
                             "Q, clean reference, equal-scheme weight)"),
  estimate_log_irr  = coef(m1_pooled_average_clean)["post_common"],
  se                = se(m1_pooled_average_clean)["post_common"]
) %>%
  add_irr_columns("estimate_log_irr", "se")

print(model1_clean_results)




######      Keep all the intervening quarters in the estimation (so the trend/FE structure is identified off the full series, same as your original models), but collapse them into a neutral "other" category that doesn't enter the ref-vs-post contrast:

stacked_m3c <- stacked %>%
  mutate(
    period_bucket = case_when(
      clean_ref_year                                    ~ "ref_year",
      treat_group == 1 & event_time >= 0 & event_time <= COMMON_POST_MAX ~ "post_common",
      TRUE                                               ~ "other"
    ),
    scheme_post_bucket = if_else(period_bucket == "post_common",
                                 as.character(stack_scheme), "ref_year"),
    scheme_post_bucket = factor(scheme_post_bucket, levels = c("ref_year", schemes_all)),
    # give "other" rows their own dummy so they're absorbed, not folded into "ref"
    other_flag = as.integer(period_bucket == "other")
  )

m3_scheme_average_bucketed <- feglm(
  outcome_raw ~
    i(scheme_post_bucket, ref = "ref_year") +
    other_flag:treat_scheme +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked_m3c,        #
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model3_bucketed_results <- extract_scheme_effects(m3_scheme_average_bucketed, "scheme_post_bucket")
print(model3_bucketed_results, n = Inf)

model3_full_comparison <- bind_rows(
  model3_results %>% mutate(spec = "1_Original"),
  model3_clean_results %>% mutate(spec = "2_Clean ref, rows dropped"),
  model3_bucketed_results %>% mutate(spec = "3_Clean ref, full data retained")
) %>%
  select(spec, scheme, pct_change, pct_lo, pct_hi, sig) %>%
  arrange(scheme, spec)
print(model3_full_comparison, n = Inf)

model3_bucketed_results <- extract_scheme_effects(m3_scheme_average_bucketed, "scheme_post_bucket")
print(model3_bucketed_results, n = Inf)

# Direct 3-way comparison
model3_full_comparison <- bind_rows(
  model3_results %>% mutate(spec = "1_Original"),
  model3_clean_results %>% mutate(spec = "2_Clean ref, rows dropped"),
  model3_bucketed_results %>% mutate(spec = "3_Clean ref, full data retained")
) %>%
  select(spec, scheme, pct_change, pct_lo, pct_hi, sig) %>%
  arrange(scheme, spec)
print(model3_full_comparison, n = Inf)




# Does the bucketed spec's fitted reference level differ meaningfully from
# m2_clean_reference's implicit reference, and from the individual dummies' reference?
# Check: what's the estimated log-rate for Bradford control & treated in EACH spec's
# reference category, on the linear predictor scale, to see where the divergence enters.

# 1. Simple sanity check: does the *manually computed* raw injury ratio in the
#    two post-window definitions match?
bradford_check <- stacked %>%
  filter(stack_scheme == "Bradford") %>%
  mutate(
    is_ref  = clean_ref_year,
    is_post = treat_group == 1 & event_time >= 0 & event_time <= 5
  )

bradford_check %>%
  filter(is_ref) %>%
  group_by(treat_group) %>%
  summarise(mean_inj = mean(outcome_raw), n = n(), .groups = "drop")

bradford_check %>%
  filter(event_time >= 0, event_time <= 5) %>%
  group_by(treat_group) %>%
  summarise(mean_inj = mean(outcome_raw), n = n(), .groups = "drop")





# -----------------------------------------------------------------------------
# Main Model 4: scheme-specific dynamic event studies
# -----------------------------------------------------------------------------
# These are scheme-level diagnostics. Within individual schemes, COVID terms may
# be dropped by fixest because event time, calendar time, and COVID periods can
# be collinear. Interpret these plots as scheme-specific dynamic checks, not as
# the headline pooled causal estimand.

run_scheme_event_ppml <- function(sc) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~
        i(event_time_ref, treat_group, ref = "ref_year") +
        i(covid_period, treat_group, ref = "non_pandemic") |
        uid_stack +
        qtr_int,
      data    = d,
      family  = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean    = TRUE
    ),
    error = function(e) {
      cat("Scheme event study failed for:", sc, "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  extract_event_study(fit, "event_time_ref") %>%
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
    title = "Main Model 4: scheme-specific PPML event studies",
    subtitle = "Reference = year before implementation where estimable",
    x = "Quarters relative to CAZ implementation",
    y = "Log incidence rate ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(outdir, "04_main_scheme_event_studies.png"),
       p_model4, width = 12, height = 10, dpi = 300)



### Bradford 
# Extend Bradford's dynamic event study to the full post-treatment window,
# using the clean-reference specification, to see whether late post-treatment
# quarters show a re-emerging gap that could explain Model 3's aggregate
model4_results %>%
  filter(scheme == "Bradford") %>%
  arrange(event_time) %>%
  print(n = Inf)

# =============================================================================
# SECTION 2 - SUPPORTING DIAGNOSTICS AND SENSITIVITY CHECKS
# =============================================================================

# -----------------------------------------------------------------------------
# A. Zero-exclusion diagnostics
# -----------------------------------------------------------------------------

zero_exclusion_diag_dir <- file.path(outdir, "zero_exclusion_diagnostics")
dir.create(zero_exclusion_diag_dir, showWarnings = FALSE, recursive = TRUE)

model_panel_zero_diag <- model_panel_for_zero_diag %>%
  mutate(zero_dropped = unit_total_injury_all_periods == 0)

zero_exclusion_unit <- model_panel_zero_diag %>%
  group_by(uid, identifier, OA, scheme, treat_group, group, zero_dropped) %>%
  summarise(
    link_stack_observations = n(),
    total_inj_all_periods = sum(outcome_raw, na.rm = TRUE),
    mean_inj_per_quarter = mean(outcome_raw, na.rm = TRUE),
    n_quarters = n_distinct(quarter_year),
    .groups = "drop"
  )

zero_obs_summary <- model_panel_zero_diag %>%
  summarise(
    link_stack_observations_total = n(),
    link_stack_observations_dropped = sum(zero_dropped),
    link_stack_observations_retained = sum(!zero_dropped),
    pct_observations_dropped = round(100 * mean(zero_dropped), 1),
    pct_observations_retained = round(100 * mean(!zero_dropped), 1)
  )

zero_unit_summary <- zero_exclusion_unit %>%
  summarise(
    link_stack_units_total = n(),
    link_stack_units_dropped = sum(zero_dropped),
    link_stack_units_retained = sum(!zero_dropped),
    pct_link_stack_units_dropped = round(100 * mean(zero_dropped), 1)
  )

zero_by_group <- zero_exclusion_unit %>%
  group_by(group, zero_dropped) %>%
  summarise(
    link_stack_units = n(),
    unique_road_links = n_distinct(identifier),
    OAs = n_distinct(OA),
    mean_quarters_observed = mean(n_quarters),
    .groups = "drop"
  ) %>%
  group_by(group) %>%
  mutate(pct_within_group = round(100 * link_stack_units / sum(link_stack_units), 1)) %>%
  ungroup()

zero_by_scheme <- zero_exclusion_unit %>%
  group_by(scheme, group, zero_dropped) %>%
  summarise(
    link_stack_units = n(),
    unique_road_links = n_distinct(identifier),
    OAs = n_distinct(OA),
    .groups = "drop"
  ) %>%
  group_by(scheme, group) %>%
  mutate(pct_within_scheme_group = round(100 * link_stack_units / sum(link_stack_units), 1)) %>%
  ungroup() %>%
  arrange(scheme, group, zero_dropped)

write_csv(zero_obs_summary,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_observation_summary.csv"))
write_csv(zero_unit_summary,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_link_stack_unit_summary.csv"))
write_csv(zero_by_group,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_by_group.csv"))
write_csv(zero_by_scheme,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_by_scheme.csv"))

# -----------------------------------------------------------------------------
# B. Raw trends and sample checks
# -----------------------------------------------------------------------------

sample_summary <- model_panel %>%
  group_by(group) %>%
  summarise(
    units = n_distinct(uid),
    observations = n(),
    total_RTI = sum(outcome_raw, na.rm = TRUE),
    mean_RTI_per_road_quarter = mean(outcome_raw, na.rm = TRUE),
    pct_zero = 100 * mean(outcome_raw == 0, na.rm = TRUE),
    .groups = "drop"
  )

scheme_sample_summary <- model_panel %>%
  distinct(scheme, uid_int, treat_group) %>%
  group_by(scheme) %>%
  summarise(
    treated_units = sum(treat_group == 1),
    control_units = sum(treat_group == 0),
    .groups = "drop"
  )

trend_df <- model_panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop")

p_trend <- ggplot(
  trend_df,
  aes(x = as.Date(quarter_year), y = mean_injury, colour = group)
) +
  annotate(
    "rect",
    xmin = as.Date(as.yearqtr("2020 Q2")),
    xmax = as.Date(as.yearqtr("2021 Q4") + 0.25),
    ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey70"
  ) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Mean injuries over time",
    subtitle = "Shaded period: COVID lockdown/recovery",
    x = NULL,
    y = "Mean injuries per road-link-quarter",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "00_raw_trends.png"),
       p_trend, width = 10, height = 6, dpi = 300)

write_csv(sample_summary, file.path(outdir, "sample_summary.csv"))
write_csv(scheme_sample_summary, file.path(outdir, "scheme_sample_summary.csv"))

# -----------------------------------------------------------------------------
# C. Alternative COVID/reference/weight specifications
# -----------------------------------------------------------------------------

m1_no_covid <- feglm(
  outcome_raw ~ post_stack |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

m1_pooled_covid <- feglm(
  outcome_raw ~
    post_stack +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

m1_link_weighted <- feglm(
  outcome_raw ~
    post_stack +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~link_weight,
  lean = TRUE
)

m2_quarter_ref <- feglm(
  outcome_raw ~
    i(event_time_f, treat_group, ref = "-1") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

model1_sensitivity_results <- tibble(
  model = c(
    "No COVID adjustment",
    "Pooled COVID adjustment",
    "Flexible scheme COVID adjustment",
    "Flexible scheme COVID, link-weighted"
  ),
  model_type = c("Support", "Support", "Main", "Support"),
  estimate_log_irr = c(
    coef(m1_no_covid)["post_stack"],
    coef(m1_pooled_covid)["post_stack"],
    coef(m1_pooled_average)["post_stack"],
    coef(m1_link_weighted)["post_stack"]
  ),
  se = c(
    se(m1_no_covid)["post_stack"],
    se(m1_pooled_covid)["post_stack"],
    se(m1_pooled_average)["post_stack"],
    se(m1_link_weighted)["post_stack"]
  )
) %>%
  add_irr_columns("estimate_log_irr", "se")

model2_quarter_ref_results <- extract_event_study(m2_quarter_ref, "event_time_f")



m2_clean_reference <- feglm(
  outcome_raw ~
    i(event_time_ref_clean, treat_group, ref = "ref_year") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

model2_cleanref_results <- extract_event_study(
  m2_clean_reference,
  "event_time_ref_clean"
)

# Joint F/Wald test: are the 6 individual post-treatment coefficients
# (event_time 0-5) jointly different from zero, even if none is individually significant?
wald(m2_clean_reference, keep = "event_time_ref_clean::[0-5]:treat_group")

# The conventional -4:-1 reference window is not a neutral default — it's a specific empirical claim ("these four quarters represent each scheme's untreated steady state"), and you've now falsified that claim for at least one scheme with a hard number: Bradford's treated/control ratio in that window is 1.30, versus 0.86 in the stable pre-COVID period — a 51% relative swing. Using a contaminated reference period doesn't add noise symmetrically; it introduces a specific, signed bias into every post-treatment coefficient for that scheme, because every event-time estimate is defined relative to that reference level. A method that anchors the reference in a period you've independently verified as stable (2019 Q2–2020 Q1, pre-COVID) is a strictly better-justified choice once you know this, not just an alternative robustness cut.
# the clean reference isn't at a constant distance from treatment across schemes the way -4:-1 is (Bradford's clean reference sits at event_time −14 to −11, Bristol's at −4 to −1, since Bristol implemented soon after the pre-COVID window while Bradford implemented years later). That's a real tradeoff — you're trading "constant relative timing, contaminated level" for "clean level, inconsistent relative timing." Report both, lead with clean reference as primary, and note this tradeoff explicitly rather than pretending it's a strict improvement with no downside.


post_event_scheme_counts <- stacked %>%
  filter(event_time >= 0) %>%
  distinct(stack_scheme, event_time) %>%
  count(event_time, name = "n_schemes") %>%
  arrange(event_time)

write_csv(model1_sensitivity_results,
          file.path(outdir, "support_model1_specification_sensitivity.csv"))
write_csv(model2_quarter_ref_results,
          file.path(outdir, "support_model2_quarter_ref_event_study.csv"))
write_csv(clean_reference_table,
          file.path(outdir, "support_clean_reference_periods.csv"))
write_csv(model2_cleanref_results,
          file.path(outdir, "support_model2_cleanref_event_study.csv"))
write_csv(post_event_scheme_counts,
          file.path(outdir, "support_post_event_scheme_counts.csv"))

# -----------------------------------------------------------------------------
# D. Parallel-trends diagnostics
# -----------------------------------------------------------------------------

wald_pretrend_yearref <- function(model, max_k) {
  pattern <- paste0(
    "event_time_ref::-(",
    paste(max_k:5, collapse = "|"),
    "):treat_group"
  )
  wald(model, keep = pattern)
}

pt_8_5  <- wald_pretrend_yearref(m2_pooled_event_yearref, 8)
pt_12_5 <- wald_pretrend_yearref(m2_pooled_event_yearref, 12)
pt_16_5 <- wald_pretrend_yearref(m2_pooled_event_yearref, 16)
pt_20_5 <- wald_pretrend_yearref(m2_pooled_event_yearref, 20)
pt_24_5 <- wald_pretrend_yearref(m2_pooled_event_yearref, 24)

wald_summary <- tibble(
  window_label = c(
    "Quarters -8 to -5",
    "Quarters -12 to -5",
    "Quarters -16 to -5",
    "Quarters -20 to -5",
    "Quarters -24 to -5"
  ),
  max_k = c(8, 12, 16, 20, 24),
  df1 = c(pt_8_5$df1, pt_12_5$df1, pt_16_5$df1, pt_20_5$df1, pt_24_5$df1),
  stat = c(pt_8_5$stat, pt_12_5$stat, pt_16_5$stat, pt_20_5$stat, pt_24_5$stat),
  p_value = c(pt_8_5$p, pt_12_5$p, pt_16_5$p, pt_20_5$p, pt_24_5$p)
) %>%
  mutate(conclusion = if_else(p_value < 0.05, "Reject H0", "Fail to reject H0"))

stacked_pre <- stacked %>%
  filter(event_time < 0) %>%
  mutate(event_time_c = event_time + 1)

m_pretrend_heterog_covid <- feglm(
  outcome_raw ~
    treat_group:event_time_c:stack_scheme +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked_pre,
  family = "poisson",
  weights = ~analysis_weight,
  cluster = ~OA,
  lean = TRUE
)

pretrend_heterog_covid_wald <- wald(
  m_pretrend_heterog_covid,
  keep = "treat_group:event_time"
)

write_csv(wald_summary, file.path(outdir, "support_parallel_trends_wald_tests.csv"))

# -----------------------------------------------------------------------------
# E. Event-time composition diagnostics
# -----------------------------------------------------------------------------

scheme_composition <- stacked_pre %>%
  count(event_time, stack_scheme, name = "n_obs") %>%
  group_by(event_time) %>%
  mutate(pct_of_bin = 100 * n_obs / sum(n_obs)) %>%
  ungroup() %>%
  arrange(event_time, desc(pct_of_bin))

composition_wide <- scheme_composition %>%
  select(event_time, stack_scheme, pct_of_bin) %>%
  pivot_wider(
    names_from = stack_scheme,
    values_from = pct_of_bin,
    values_fill = 0
  ) %>%
  arrange(event_time)

scheme_event_range <- stacked_pre %>%
  group_by(stack_scheme) %>%
  summarise(
    min_event_time = min(event_time),
    max_event_time = max(event_time),
    n_quarters_pre = n_distinct(event_time),
    .groups = "drop"
  ) %>%
  arrange(min_event_time)

write_csv(composition_wide,
          file.path(outdir, "support_pretrend_scheme_composition_by_event_time.csv"))
write_csv(scheme_event_range,
          file.path(outdir, "support_pretrend_scheme_event_range.csv"))

p_composition <- ggplot(
  scheme_composition,
  aes(x = event_time, y = pct_of_bin, fill = stack_scheme)
) +
  geom_area(position = "stack") +
  labs(
    title = "Scheme composition of pre-treatment event-time bins",
    x = "Quarters relative to CAZ implementation",
    y = "% of observations in event-time bin",
    fill = "Scheme"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "support_pretrend_scheme_composition.png"),
       p_composition, width = 10, height = 6, dpi = 300)




# -----------------------------------------------------------------------------
#. Which scheme(s) drive the pooled pretrend rejection?
# -----------------------------------------------------------------------------

pretrend_decomp <- model4_results %>%
  filter(event_time %in% -8:-5) %>%
  select(scheme, event_time, estimate, se) %>%
  left_join(
    scheme_composition %>%
      filter(event_time %in% -8:-5) %>%
      rename(scheme = stack_scheme),
    by = c("scheme", "event_time")
  ) %>%
  mutate(
    z = estimate / se,
    weighted_contribution = estimate * (pct_of_bin / 100)
  ) %>%
  arrange(event_time, desc(abs(weighted_contribution)))

cat("\nScheme-level contributions to pooled pretrend coefficients, event_time -8:-5\n")
print(pretrend_decomp, n = Inf)

# Quick visual: composition-weighted scheme contribution by event time
ggplot(pretrend_decomp, aes(x = factor(event_time), y = weighted_contribution, fill = scheme)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Composition-weighted scheme contributions to pooled pretrend coefficients",
       x = "Event time", y = "Composition-weighted log-IRR contribution") +
  theme_minimal(base_size = 12)

# -----------------------------------------------------------------------------
# F. Bradford diagnostic investigation
# -----------------------------------------------------------------------------

bradford_stack <- stacked %>%
  filter(stack_scheme == "Bradford")

bradford_period_summary <- bradford_stack %>%
  mutate(period = if_else(event_time >= 0, "post", "pre")) %>%
  group_by(treat_group, period) %>%
  summarise(
    obs = n(),
    mean_injury = mean(outcome_raw),
    pct_zero = 100 * mean(outcome_raw == 0),
    total_injury = sum(outcome_raw),
    .groups = "drop"
  )

bradford_indexed <- bradford_stack %>%
  mutate(
    group = if_else(treat_group == 1, "Treated", "Control"),
    period = if_else(event_time >= 0, "post", "pre")
  ) %>%
  group_by(group, period) %>%
  summarise(mean_inj = mean(outcome_raw), .groups = "drop") %>%
  group_by(group) %>%
  mutate(index = mean_inj / mean_inj[period == "pre"] * 100) %>%
  ungroup()

same_scheme_overlap <- stacked %>%
  group_by(stack_scheme, OA) %>%
  summarise(
    appears_treated = any(treat_group == 1),
    appears_control = any(treat_group == 0),
    .groups = "drop"
  ) %>%
  filter(appears_treated & appears_control)

control_post_pre_index <- stacked %>%
  filter(treat_group == 0) %>%
  mutate(period = if_else(event_time >= 0, "post", "pre")) %>%
  group_by(stack_scheme, period) %>%
  summarise(
    injury_per_qtr = sum(outcome_raw) / n_distinct(quarter_year),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period, values_from = injury_per_qtr) %>%
  mutate(index = 100 * post / pre) %>%
  arrange(index)

reference_year_ratio <- stacked %>%
  filter(event_time >= -4, event_time <= -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    injury_per_qtr = sum(outcome_raw) / n_distinct(quarter_year),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = treat_group,
    values_from = injury_per_qtr,
    names_prefix = "group_"
  ) %>%
  mutate(control_to_treated_ratio = group_0 / group_1) %>%
  arrange(desc(control_to_treated_ratio))


# -----------------------------------------------------------------------------
#  Bradford controls : where do controls actually come from?
# -----------------------------------------------------------------------------
# whether any OA appears as BOTH treated and control within Bradford
#    (would indicate treated/control roads sitting in the same small area)
bradford_oa_overlap <- bradford_stack %>%
  distinct(OA, treat_group) %>%
  count(OA, treat_group) %>%
  pivot_wider(names_from = treat_group, values_from = n, values_fill = 0,
              names_prefix = "treat_")

bradford_oa_overlap %>%
  mutate(status = case_when(
    treat_1 > 0 & treat_0 > 0 ~ "BOTH treated and control in same OA",
    treat_1 > 0 ~ "treated only",
    TRUE ~ "control only"
  )) %>%
  count(status)
control_reuse <- stacked %>%
  filter(treat_group == 0) %>%
  distinct(identifier, stack_scheme) %>%
  count(identifier, name = "n_schemes_as_control")

cat("\nHow many control links are reused across multiple scheme stacks:\n")
print(count(control_reuse, n_schemes_as_control))

bradford_shared <- stacked %>%
  filter(stack_scheme == "Bradford", treat_group == 0,
         identifier %in% control_reuse$identifier[control_reuse$n_schemes_as_control > 1]) %>%
  distinct(identifier) %>% nrow()
cat("\nBradford control links also used as controls elsewhere:", bradford_shared, "\n")

# 6. Leverage check: is the -5/-6 spike coming from a handful of high-injury control links?
bradford_control_leverage <- bradford_stack %>%
  filter(treat_group == 0, event_time < 0) %>%
  group_by(identifier) %>%
  summarise(total_pre_injuries = sum(outcome_raw), .groups = "drop") %>%
  arrange(desc(total_pre_injuries)) %>%
  mutate(cum_share = cumsum(total_pre_injuries) / sum(total_pre_injuries))

cat("\nTop 10 Bradford control links by pre-period injuries:\n")
print(head(bradford_control_leverage, 10))
cat("Share of all Bradford control pre-period injuries from top 10 links: ",
    round(100 * bradford_control_leverage$cum_share[10], 1), "%\n", sep = "")

# 7. Robustness: drop the top-leverage control links and re-run Bradford's event study
top_leverage_ids <- head(bradford_control_leverage$identifier, 10)

bradford_trimmed <- bradford_stack %>%
  filter(!(treat_group == 0 & identifier %in% top_leverage_ids)) %>%
  droplevels()

m_bradford_trimmed <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + qtr_int,
  data = bradford_trimmed, family = "poisson",
  cluster = ~OA, weights = ~analysis_weight, lean = TRUE
)

cat("\nBradford event study excluding top-10-leverage control links:\n")
extract_event_study(m_bradford_trimmed, "event_time_ref") %>%
  filter(event_time >= -8, event_time <= 8) %>%
  print(n = Inf)

# 8. Finer COVID control: continuous ramp instead of 3-level factor,
#    to test whether the -6/-7 spike is a residual differential-recovery artifact
bradford_stack_ramp <- bradford_stack %>%
  mutate(
    qtrs_since_lockdown_start = pmax(0, qtr_int -
                                       (as.integer(round((as.numeric(as.yearqtr("2020 Q2")) - min_qtr) * 4)) + 1L)),
    recovery_ramp = if_else(covid_period != "non_pandemic", qtrs_since_lockdown_start, 0)
  )

m_bradford_ramp <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    treat_group:recovery_ramp |
    uid_stack + qtr_int,
  data = bradford_stack_ramp, family = "poisson",
  cluster = ~OA, weights = ~analysis_weight, lean = TRUE
)

cat("\nBradford event study with continuous COVID-recovery ramp control:\n")
extract_event_study(m_bradford_ramp, "event_time_ref") %>%
  filter(event_time >= -8, event_time <= 8) %>%
  print(n = Inf)





# Rebuild Bradford's scheme-specific event study using the CLEAN reference
# (to match the Model 3 spec you're trying to reconcile), then run the joint
# post-treatment Wald test within Bradford alone.

bradford_stack_clean <- stacked %>%
  filter(stack_scheme == "Bradford") %>%
  droplevels()

m_bradford_event_clean <- feglm(
  outcome_raw ~
    i(event_time_ref_clean, treat_group, ref = "ref_year") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + qtr_int,
  data    = bradford_stack_clean,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

cat("\nBradford event study, clean pre-COVID reference:\n")
extract_event_study(m_bradford_event_clean, "event_time_ref_clean") %>%
  filter(event_time >= 0, event_time <= 8) %>%
  print(n = Inf)

# The joint test that actually matches Model 3's claim
wald(m_bradford_event_clean, keep = "event_time_ref_clean::[0-5]:treat_group")



write_csv(bradford_period_summary,
          file.path(outdir, "support_bradford_period_summary.csv"))
# -----------------------------------------------------------------------------
# 9. Pinpoint exactly which calendar quarter drives the -5/-6 divergence
#    (fully saturated quarter x treatment, no event-time binning, no linear trend)
# -----------------------------------------------------------------------------

bradford_saturated <- bradford_stack %>%
  filter(quarter_year >= as.yearqtr("2020 Q1"), quarter_year <= as.yearqtr("2022 Q4")) %>%
  mutate(qtr_f = factor(quarter_year))

m_bradford_saturated <- feglm(
  outcome_raw ~ i(qtr_f, treat_group, ref = as.character(as.yearqtr("2022 Q3"))) |
    uid_stack,
  data = bradford_saturated, family = "poisson",
  cluster = ~OA, weights = ~analysis_weight, lean = TRUE
)

cat("\nBradford: fully saturated calendar-quarter x treatment (no trend, no COVID dummy)\n")
etable(m_bradford_saturated)
# Look specifically at 2021 Q2 and 2021 Q3 (= event_time -6, -5) --
# if the coefficient is concentrated in exactly one quarter rather than ramping
# smoothly, that points to a discrete event/anticipation shock rather than a
# recovery-trajectory confound.

# -----------------------------------------------------------------------------
# 10. Leverage check restricted specifically to event_time -5 and -6
#     (the earlier check pooled ALL pre-period injuries, diluting a concentrated spike)
# -----------------------------------------------------------------------------

bradford_spike_leverage <- bradford_stack %>%
  filter(event_time %in% c(-6, -5)) %>%
  group_by(identifier, treat_group) %>%
  summarise(spike_injuries = sum(outcome_raw), n_obs = n(), .groups = "drop") %>%
  arrange(desc(spike_injuries))

cat("\nTop control AND treated links driving the -5/-6 window specifically:\n")
print(head(bradford_spike_leverage, 15))

# Compare: is the -5/-6 rise concentrated in a few links, or diffuse across many?
bradford_spike_leverage %>%
  filter(treat_group == 0) %>%
  mutate(cum_share = cumsum(spike_injuries) / sum(spike_injuries)) %>%
  slice_head(n = 10)

# -----------------------------------------------------------------------------
# 11. Anticipation check: was Bradford's CAZ start date announced/delayed
#     well before caz_start_q? Test whether treated roads diverge from an
#     EARLIER hypothetical "announcement" cutoff, not just the go-live date.
#     (fill in the actual announcement date once confirmed)
# -----------------------------------------------------------------------------

bradford_announce_test <- bradford_stack %>%
  mutate(
    # placeholder - replace with actual publicly announced/confirmed start date
    announce_q = as.yearqtr("2021 Q3"),
    post_announce = as.integer(treat_group == 1 & quarter_year >= announce_q)
  )

m_bradford_announce <- feglm(
  outcome_raw ~ post_announce + post_stack |
    uid_stack + qtr_int,
  data = bradford_announce_test, family = "poisson",
  cluster = ~OA, weights = ~analysis_weight, lean = TRUE
)
etable(m_bradford_announce)
# If post_announce absorbs a chunk of what was previously loading onto the
# -5/-6 pretrend coefficients, that's evidence of genuine anticipation
# rather than a COVID artifact, and argues for redefining "treatment start"
# using the announcement date rather than go-live for Bradford (and possibly
# other delayed schemes).

# -----------------------------------------------------------------------------
# 12. Same diagnostic run scheme-by-scheme, since the composition chart shows
#     Newcastle and Birmingham also drive large -6/-5 contributions
# -----------------------------------------------------------------------------

run_saturated_check <- function(sc, lo_q = "2020 Q1", hi_q = "2022 Q4", ref_q = NULL) {
  d <- stacked %>%
    filter(
      stack_scheme == sc,
      quarter_year >= as.yearqtr(lo_q),
      quarter_year <= as.yearqtr(hi_q)
    ) %>%
    mutate(qtr_chr = as.character(quarter_year))
  
  ref <- if (is.null(ref_q)) {
    sort(unique(d$qtr_chr))[1]
  } else {
    as.character(as.yearqtr(ref_q))
  }
  
  if (!ref %in% d$qtr_chr) {
    stop(
      "Reference quarter ", ref, " is not present for scheme ", sc,
      ". Available quarters: ", paste(sort(unique(d$qtr_chr)), collapse = ", ")
    )
  }
  
  d <- d %>%
    mutate(
      qtr_f_ref = if_else(qtr_chr == ref, "REF_QTR", qtr_chr),
      qtr_f_ref = factor(
        qtr_f_ref,
        levels = c("REF_QTR", sort(setdiff(unique(qtr_chr), ref)))
      )
    )
  
  feglm(
    outcome_raw ~ i(qtr_f_ref, treat_group, ref = "REF_QTR") | uid_stack,
    data = d,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = TRUE
  )
}

for (sc in c("Bradford", "Newcastle", "Birmingham")) {
  cat("\n===", sc, "===\n")
  print(etable(run_saturated_check(sc)))
}


# -----------------------------------------------------------------------------
# Decompose: is the -6/-5 divergence driven by treated levels rising,
# or control levels falling?
# -----------------------------------------------------------------------------

bradford_quarterly_means <- bradford_stack %>%
  filter(quarter_year >= as.yearqtr("2019 Q1"), quarter_year <= as.yearqtr("2022 Q4")) %>%
  group_by(quarter_year, treat_group) %>%
  summarise(
    mean_injury = mean(outcome_raw),
    n_links     = n_distinct(identifier),
    total_injury = sum(outcome_raw),
    .groups = "drop"
  ) %>%
  mutate(group = if_else(treat_group == 1, "Treated", "Control"))

# Index each series to its own 2019 average, so you can see relative movement
bradford_indexed_series <- bradford_quarterly_means %>%
  group_by(group) %>%
  mutate(
    base = mean(mean_injury[quarter_year >= as.yearqtr("2019 Q1") &
                              quarter_year <= as.yearqtr("2019 Q4")]),
    index = 100 * mean_injury / base
  ) %>%
  ungroup()

ggplot(bradford_indexed_series, aes(x = as.Date(quarter_year), y = index, colour = group)) +
  geom_vline(xintercept = as.Date(as.yearqtr("2021 Q2")), linetype = "dotted") +
  geom_vline(xintercept = as.Date(as.yearqtr("2021 Q3")), linetype = "dotted") +
  geom_vline(xintercept = as.Date(as.yearqtr("2022 Q4")), linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  labs(title = "Bradford: treated vs. control, indexed to 2019 average",
       subtitle = "Dotted = event times -6/-5 (2021 Q2/Q3); dashed = actual CAZ start",
       x = NULL, y = "Index (2019 = 100)") +
  theme_minimal(base_size = 12)

# Quarter-on-quarter change for each series, to see which one moves first/more
bradford_qoq <- bradford_quarterly_means %>%
  arrange(group, quarter_year) %>%
  group_by(group) %>%
  mutate(qoq_change = mean_injury - lag(mean_injury)) %>%
  ungroup() %>%
  filter(quarter_year >= as.yearqtr("2021 Q1"), quarter_year <= as.yearqtr("2021 Q4"))

print(bradford_qoq, n = Inf)



###   treated roads recovered from COVID unusually sharply in 2021, 
# peaked around Q3, and then partially reverted. 
# That maps precisely onto the event-time coefficients: −6/−5 (2021 Q2/Q3) are the peak of the overshoot 
#and are strongly positive; by event_time 0 (2022 Q4) the series have mostly reconverged


### the recovery probleme 



plot_indexed_recovery <- function(sc, base_start = "2019 Q1", base_end = "2019 Q4") {
  d <- stacked %>%
    filter(stack_scheme == sc,
           quarter_year >= as.yearqtr("2019 Q1"), quarter_year <= as.yearqtr("2022 Q4")) %>%
    group_by(quarter_year, treat_group) %>%
    summarise(mean_injury = mean(outcome_raw), .groups = "drop") %>%
    mutate(group = if_else(treat_group == 1, "Treated", "Control")) %>%
    group_by(group) %>%
    mutate(
      base = mean(mean_injury[quarter_year >= as.yearqtr(base_start) &
                                quarter_year <= as.yearqtr(base_end)]),
      index = 100 * mean_injury / base
    ) %>%
    ungroup()
  
  ggplot(d, aes(x = as.Date(quarter_year), y = index, colour = group)) +
    geom_line(linewidth = 0.9) +
    labs(title = paste(sc, ": treated vs. control, indexed to 2019 average"),
         x = NULL, y = "Index (2019 = 100)") +
    theme_minimal(base_size = 12)
}


plot_indexed_recovery("Newcastle")
plot_indexed_recovery("Birmingham")



###   how much of the reference-period inflation is left by event_time −4 to −1

bradford_ref_vs_precovid <- bradford_stack %>%
  mutate(period = case_when(
    quarter_year >= as.yearqtr("2019 Q1") & quarter_year <= as.yearqtr("2019 Q4") ~ "pre-COVID baseline",
    event_time >= -4 & event_time <= -1 ~ "reference window",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(period, treat_group) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop") %>%
  pivot_wider(names_from = treat_group, values_from = mean_injury, names_prefix = "grp_") %>%
  mutate(treat_control_ratio = grp_1 / grp_0)

print(bradford_ref_vs_precovid)





write_csv(bradford_indexed,
          file.path(outdir, "support_bradford_indexed_pre_post.csv"))
write_csv(same_scheme_overlap,
          file.path(outdir, "support_same_scheme_oa_overlap.csv"))
write_csv(control_post_pre_index,
          file.path(outdir, "support_control_post_pre_index_by_scheme.csv"))
write_csv(reference_year_ratio,
          file.path(outdir, "support_reference_year_control_treated_ratio.csv"))



write_csv(model1_results, file.path(outdir, "model1_main_pooled_average_effect.csv"))
write_csv(model2_results, file.path(outdir, "model2_main_pooled_yearref_event_study.csv"))
write_csv(model3_results, file.path(outdir, "model3_main_scheme_average_effects.csv"))
write_csv(model4_results, file.path(outdir, "model4_main_scheme_event_studies.csv"))



cat("\n\n")
cat("###############################################################################\n")
cat("# CAZ PPML ANALYSIS CONSOLE REPORT\n")
cat("###############################################################################\n")

cat("\n\n# 0. Scheme Timing\n")
print(scheme_timing, n = Inf)

cat("\n\n# 1. Retained Analysis Sample\n")
print(sample_summary, n = Inf)

cat("\n\n# 2. Scheme Sample Sizes\n")
print(scheme_sample_summary, n = Inf)

cat("\n\n# 3. Equal-Scheme Weight Verification at Event Time -1\n")
stacked %>%
  filter(event_time == -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    n_rows = n(),
    sum_analysis_weight = sum(analysis_weight),
    sum_link_weight = sum(link_weight),
    .groups = "drop"
  ) %>%
  print(n = Inf)

cat("\n\n# 4. MAIN MODEL 1: Pooled Average Effect\n")
etable(
  m1_pooled_average,
  dict = c("post_stack" = "CAZ post-treatment")
)
cat("\nMain pooled average effect, converted to IRR and percent change:\n")
print(model1_results, n = Inf)

cat("\n\n# 5. MAIN MODEL 2: Pooled Year-Reference Event Study\n")
etable(m2_pooled_event_yearref)
cat("\nMain event-study estimates, all event times:\n")
print(model2_results, n = Inf)
cat("\nMain event-study estimates, common post-treatment window 0-5:\n")
model2_results %>%
  filter(event_time >= 0, event_time <= 5) %>%
  print(n = Inf)

cat("\n\n# 6. MAIN MODEL 3: Scheme-Specific Average Effects\n")
etable(m3_scheme_average)
cat("\nScheme-specific average effects, converted to IRR and percent change:\n")
print(model3_results, n = Inf)

cat("\n\n# 7. MAIN MODEL 4: Scheme-Specific Event Studies\n")
cat("\nScheme-specific event-study estimates, event times -8 to 8:\n")
model4_results %>%
  filter(event_time >= -8, event_time <= 8) %>%
  arrange(scheme, event_time) %>%
  print(n = Inf)

cat("\n\n# 8. SUPPORT: Model 1 Specification Sensitivity\n")
print(model1_sensitivity_results, n = Inf)

cat("\n\n# 9. SUPPORT: Quarter-Reference Event Study, Event Times -8 to 8\n")
model2_quarter_ref_results %>%
  filter(event_time >= -8, event_time <= 8) %>%
  print(n = Inf)

cat("\n\n# 10. SUPPORT: Clean Reference Periods\n")
print(clean_reference_table, n = Inf)

cat("\n\n# 11. SUPPORT: Clean-Reference Event Study\n")
cat("\nAll clean-reference event-study estimates:\n")
print(model2_cleanref_results, n = Inf)
cat("\nClean-reference estimates, common post-treatment window 0-5:\n")
model2_cleanref_results %>%
  filter(event_time >= 0, event_time <= 5) %>%
  print(n = Inf)

cat("\n\n# 12. SUPPORT: Number of Schemes Contributing by Post Event Time\n")
print(post_event_scheme_counts, n = Inf)

cat("\n\n# 13. SUPPORT: Parallel-Trends Wald Tests\n")
print(wald_summary, n = Inf)
cat("\n\nFlexible scheme-by-treatment COVID pretrend Wald test:\n")
print(pretrend_heterog_covid_wald)

cat("\n\n# 14. SUPPORT: Pretrend Scheme Composition\n")
cat("\nPre-period event-time range by scheme:\n")
print(scheme_event_range, n = Inf)
cat("\nPre-period scheme composition by event time:\n")
print(composition_wide, n = Inf)

cat("\n\n# 15. SUPPORT: Zero-Exclusion Diagnostics\n")
cat("\nObservation-level zero-exclusion summary:\n")
print(zero_obs_summary, n = Inf)
cat("\nUnit-level zero-exclusion summary:\n")
print(zero_unit_summary, n = Inf)
cat("\nZero exclusion by treatment group:\n")
print(zero_by_group, n = Inf)
cat("\nZero exclusion by scheme and treatment group:\n")
print(zero_by_scheme, n = Inf)

cat("\n\n# 16. SUPPORT: Bradford Diagnostics\n")
cat("\nBradford pre/post period summary:\n")
print(bradford_period_summary, n = Inf)
cat("\nBradford indexed pre/post comparison:\n")
print(bradford_indexed, n = Inf)
cat("\nSame-scheme OA overlap check:\n")
print(same_scheme_overlap, n = Inf)
cat("\nControl post/pre injury index by scheme:\n")
print(control_post_pre_index, n = Inf)
cat("\nControl-to-treated injury ratio in conventional reference year:\n")
print(reference_year_ratio, n = Inf)

cat("\n\n# 17. Model Objects Available in R Session\n")
cat("Main models:\n")
cat("  m1_pooled_average\n")
cat("  m2_pooled_event_yearref\n")
cat("  m3_scheme_average\n")
cat("  model4_results\n")
cat("Supporting models:\n")
cat("  m1_no_covid\n")
cat("  m1_pooled_covid\n")
cat("  m1_link_weighted\n")
cat("  m2_quarter_ref\n")
cat("  m2_clean_reference\n")
cat("  m_pretrend_heterog_covid\n")

cat("\n###############################################################################\n")
cat("# END OF CAZ PPML ANALYSIS CONSOLE REPORT\n")
cat("###############################################################################\n\n")

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
    main_models = list(
      m1_pooled_average = m1_pooled_average,
      m2_pooled_event_yearref = m2_pooled_event_yearref,
      m3_scheme_average = m3_scheme_average,
      model4_results = model4_results
    ),
    main_results = list(
      model1_results = model1_results,
      model2_results = model2_results,
      model3_results = model3_results,
      model4_results = model4_results
    ),
    support_results = list(
      model1_sensitivity_results = model1_sensitivity_results,
      model2_quarter_ref_results = model2_quarter_ref_results,
      model2_cleanref_results = model2_cleanref_results,
      clean_reference_table = clean_reference_table,
      post_event_scheme_counts = post_event_scheme_counts,
      wald_summary = wald_summary,
      pretrend_heterog_covid_wald = pretrend_heterog_covid_wald,
      bradford_indexed = bradford_indexed
    )
  ),
  here("data", "processed", "caz_ppml_reorganised_results.rds")
)

cat("\nReorganised analysis complete. Outputs saved to:", outdir, "\n")





# ==============================================================================
# MODEL 3 (REVISED): scheme-specific average effects
# Fixes vs. original: (1) comparison anchored to clean pre-COVID reference period
# only, not an unweighted pool of the entire non-stationary pre-period + controls;
# (2) common fixed post-treatment window (0-5 quarters) so all 7 schemes are
# compared on equal footing (per post_event_scheme_counts, all 7 are present
# through event_time 5).
# ==============================================================================

COMMON_POST_MAX <- 5

stacked_m3 <- stacked %>%
  filter(clean_ref_year | (event_time >= 0 & event_time <= COMMON_POST_MAX)) %>%
  mutate(
    post_common        = as.integer(treat_group == 1 &
                                      event_time >= 0 & event_time <= COMMON_POST_MAX),
    scheme_post_common  = if_else(post_common == 1, as.character(stack_scheme), "control"),
    scheme_post_common  = factor(scheme_post_common, levels = c("control", schemes_all))
  )

m3_scheme_average_clean <- feglm(
  outcome_raw ~
    i(scheme_post_common, ref = "control") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked_m3,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model3_clean_results <- extract_scheme_effects(m3_scheme_average_clean, "scheme_post_common")

cat("\nModel 3 (revised): scheme-specific effects, common 0-", COMMON_POST_MAX,
    "Q post window, clean pre-COVID reference\n", sep = "")
etable(m3_scheme_average_clean)
print(model3_clean_results, n = Inf)

# -----------------------------------------------------------------------------
# Direct before/after comparison, to document the correction in your report
# -----------------------------------------------------------------------------
model3_comparison <- bind_rows(
  model3_results %>% mutate(spec = "Original: unbalanced window, uncleaned reference"),
  model3_clean_results %>% mutate(spec = paste0("Revised: common 0-", COMMON_POST_MAX,
                                                "Q window, clean reference"))
) %>%
  select(spec, scheme, pct_change, pct_lo, pct_hi, p_value, sig) %>%
  arrange(scheme, spec)

print(model3_comparison, n = Inf)

p_model3_comparison <- ggplot(
  model3_comparison,
  aes(x = pct_change, y = fct_reorder(scheme, pct_change), colour = spec)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Scheme-specific average effects: original vs. corrected specification",
    subtitle = "Corrected: common post-treatment window + clean pre-COVID reference period",
    x = "% change in injuries", y = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "03b_model3_original_vs_clean_comparison.png"),
       p_model3_comparison, width = 10, height = 6, dpi = 300)







m1_pooled_average_clean <- feglm(
  outcome_raw ~
    post_common +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked_m3,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model1_clean_results <- tibble(
  model             = paste0("Pooled average (0-", COMMON_POST_MAX,
                             "Q, clean reference, equal-scheme weight)"),
  estimate_log_irr  = coef(m1_pooled_average_clean)["post_common"],
  se                = se(m1_pooled_average_clean)["post_common"]
) %>%
  add_irr_columns("estimate_log_irr", "se")

print(model1_clean_results)








