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

extract_scheme_effects <- function(model) {
  ct <- coeftable(model)
  
  tibble(
    term             = rownames(ct),
    estimate_log_irr = ct[, "Estimate"],
    se               = ct[, "Std. Error"],
    p_value          = ct[, "Pr(>|z|)"]
  ) %>%
    filter(str_detect(term, "^scheme_post::")) %>%
    mutate(
      scheme = str_remove(term, "scheme_post::")
    ) %>%
    add_irr_columns("estimate_log_irr", "se") %>%
    mutate(
      sig = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        p_value < 0.10  ~ ".",
        TRUE            ~ ""
      )
    ) %>%
    select(
      scheme, estimate_log_irr, se, irr, irr_lo, irr_hi,
      pct_change, pct_lo, pct_hi, p_value, sig
    ) %>%
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

stacked <- make_clean_reference(stacked)

clean_reference_table <- stacked %>%
  filter(clean_ref_year) %>%
  distinct(stack_scheme, quarter_year, event_time, covid_period) %>%
  arrange(stack_scheme, quarter_year)

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

write_csv(bradford_period_summary,
          file.path(outdir, "support_bradford_period_summary.csv"))
write_csv(bradford_indexed,
          file.path(outdir, "support_bradford_indexed_pre_post.csv"))
write_csv(same_scheme_overlap,
          file.path(outdir, "support_same_scheme_oa_overlap.csv"))
write_csv(control_post_pre_index,
          file.path(outdir, "support_control_post_pre_index_by_scheme.csv"))
write_csv(reference_year_ratio,
          file.path(outdir, "support_reference_year_control_treated_ratio.csv"))

# =============================================================================
# SAVE MAIN OUTPUTS
# =============================================================================

write_csv(model1_results, file.path(outdir, "model1_main_pooled_average_effect.csv"))
write_csv(model2_results, file.path(outdir, "model2_main_pooled_yearref_event_study.csv"))
write_csv(model3_results, file.path(outdir, "model3_main_scheme_average_effects.csv"))
write_csv(model4_results, file.path(outdir, "model4_main_scheme_event_studies.csv"))

# =============================================================================
# CONSOLE REPORT
# =============================================================================
# This section prints the main outputs and key supporting diagnostics so the
# console log can be copied into a follow-up review.

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
