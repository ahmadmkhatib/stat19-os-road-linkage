# =============================================================================
# CAZ DiD ANALYSIS - MAIN PPML ANALYSIS SCRIPT
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
# Weighting strategy:
#   The matched panel retains the original matching ratios because control
#   rows are duplicated during panel construction when matching is done with
#   replacement or at ratios greater than 1:1.
#
#   Main pooled models use one analysis_weight column. Within each scheme,
#   treated road-link units are weighted to sum to 1 and control road-link
#   units are also weighted to sum to 1:
#
#       treated row-link weight = 1 / number of treated units in scheme
#       control row-link weight = 1 / number of control units in scheme
#
#   This gives every scheme the same total treated weight and the same total
#   control weight in the pooled analysis, while retaining the original
#   matched control ratio through the duplicated matched rows.
#
#   A link-weighted sensitivity model is also estimated. It removes the
#   equal-scheme normalisation but still balances the total control weight to
#   the total treated weight within each scheme.
#
# Notes:
#   - Outcome is a sparse count (~99% zeros), so models are estimated with
#     Poisson pseudo-maximum likelihood (PPML).
#   - Standard errors are clustered at OA level to account for spatial
#     correlation and for the fact that the same control OA may appear
#     in multiple scheme stacks.
#   - Weight calculations are done after excluding road links with zero
#     injuries in every observed quarter, so the weights are outcome/sample
#     specific.
#   - Parallel-trends diagnostics should be interpreted after rerunning the
#     updated matching pipeline, because Stage 2 now balances additional
#     pre-COVID trajectory-shape variables.
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
      ci_lo      = estimate - 1.96 * se,
      ci_hi      = estimate + 1.96 * se,
      irr        = exp(estimate),
      irr_lo     = exp(ci_lo),
      irr_hi     = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo     = 100 * (irr_lo - 1),
      pct_hi     = 100 * (irr_hi - 1)
    ) %>%
    arrange(event_time)
}

extract_fixest_event_study_yearref <- function(model) {
  ct <- coeftable(model)
  tibble(
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"]
  ) %>%
    filter(str_detect(term, "event_time_ref::") &
             !str_detect(term, "ref_year")) %>%
    mutate(
      event_time = str_extract(term, "(?<=event_time_ref::)-?\\d+") %>% as.numeric(),
      ci_lo      = estimate - 1.96 * se,
      ci_hi      = estimate + 1.96 * se,
      irr        = exp(estimate),
      irr_lo     = exp(ci_lo),
      irr_hi     = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo     = 100 * (irr_lo - 1),
      pct_hi     = 100 * (irr_hi - 1)
    ) %>%
    arrange(event_time)
}

extract_fixest_event_study_year <- function(model) {
  ct <- coeftable(model)
  tibble(
    term     = rownames(ct),
    estimate = ct[, "Estimate"],
    se       = ct[, "Std. Error"]
  ) %>%
    filter(str_detect(term, "event_time_year_f::")) %>%
    mutate(
      event_time = str_extract(term, "(?<=event_time_year_f::)-?\\d+") %>% as.numeric(),
      ci_lo      = estimate - 1.96 * se,
      ci_hi      = estimate + 1.96 * se,
      irr        = exp(estimate),
      irr_lo     = exp(ci_lo),
      irr_hi     = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo     = 100 * (irr_lo - 1),
      pct_hi     = 100 * (irr_hi - 1)
    ) %>%
    arrange(event_time)
}

plot_event_study <- function(df, title, subtitle, ylab, colour = "#E74C3C") {
  ggplot(df, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = colour) +
    geom_line(linewidth = 0.8, colour = colour) +
    geom_point(size = 1.8,    colour = colour) +
    scale_x_continuous(breaks = pretty(df$event_time, n = 10)) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Quarters relative to CAZ implementation",
      y        = ylab
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# =============================================================================
# LOAD DATA
# =============================================================================
# Apply majority-quarter rule: if a CAZ starts in the second half of a
# quarter, treatment is coded as starting the following quarter.

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

# =============================================================================
# BUILD MODEL PANEL
# =============================================================================
# COVID period coding based on visual inspection of raw injury trends:
#   lockdown     : 2020 Q2 - 2021 Q1
#   recovery     : 2021 Q2 - 2021 Q4
#   non_pandemic : everything else (reference)
#
# Separating post-2021 Q4 into its own level would be collinear with the
# post-treatment event-time indicators for five of seven schemes whose
# entire post-treatment window falls after 2021 Q4.
#
# Road-link units with zero injuries across the entire study period are
# excluded before modelling - they contribute no information to PPML and
# inflate the count of singletons removed by the fixed-effect sweep.

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
    
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    ),
    
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
  # Flag and remove road-link units with zero injuries throughout study period.
  # These are dropped by PPML anyway (singleton/all-zero fixed effects) but
  # removing them explicitly here keeps the sample counts honest.
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




# =============================================================================
# DIAGNOSTIC: ALL-ZERO ROAD-LINK EXCLUSION
# =============================================================================
# IMPORTANT:
# Run this BEFORE:
#   model_panel <- model_panel %>%
#     filter(unit_total_injury_all_periods > 0) %>%
#     select(-unit_total_injury_all_periods)
#
# If you already ran that filter, rerun the script from the LOAD DATA / BUILD
# MODEL PANEL section, then run this block before the exclusion.

zero_exclusion_diag_dir <- file.path(outdir, "zero_exclusion_diagnostics")
dir.create(zero_exclusion_diag_dir, showWarnings = FALSE, recursive = TRUE)

model_panel_zero_diag <- model_panel %>%
  group_by(uid) %>%
  mutate(
    unit_total_injury_all_periods = sum(outcome_raw, na.rm = TRUE),
    zero_dropped = unit_total_injury_all_periods == 0
  ) %>%
  ungroup()

zero_exclusion_unit <- model_panel_zero_diag %>%
  mutate(
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  ) %>%
  group_by(uid, identifier, OA, scheme, treat_group, group, zero_dropped) %>%
  summarise(
    link_stack_observations = n(),
    total_inj_all_periods   = sum(outcome_raw, na.rm = TRUE),
    mean_inj_per_quarter    = mean(outcome_raw, na.rm = TRUE),
    first_quarter           = min(quarter_year, na.rm = TRUE),
    last_quarter            = max(quarter_year, na.rm = TRUE),
    n_quarters              = n_distinct(quarter_year),
    .groups = "drop"
  )

zero_obs_summary <- model_panel_zero_diag %>%
  summarise(
    link_stack_observations_total    = n(),
    link_stack_observations_dropped  = sum(zero_dropped),
    link_stack_observations_retained = sum(!zero_dropped),
    pct_observations_dropped         = 100 * mean(zero_dropped)
  )

cat("\nNumber of link-stack observations dropped:\n")
print(zero_obs_summary)

zero_unit_summary <- zero_exclusion_unit %>%
  summarise(
    link_stack_units_total       = n(),
    link_stack_units_dropped     = sum(zero_dropped),
    link_stack_units_retained    = sum(!zero_dropped),
    pct_link_stack_units_dropped = 100 * mean(zero_dropped)
  )

cat("\nNumber of unique link-stack units dropped:\n")
print(zero_unit_summary)

zero_identifier_summary <- zero_exclusion_unit %>%
  group_by(identifier) %>%
  summarise(
    appears_in_n_stacks  = n_distinct(scheme),
    dropped_in_all_stacks = all(zero_dropped),
    dropped_in_any_stack  = any(zero_dropped),
    .groups = "drop"
  ) %>%
  summarise(
    unique_road_links_total = n(),
    unique_road_links_dropped_in_all_stacks = sum(dropped_in_all_stacks),
    unique_road_links_dropped_in_any_stack  = sum(dropped_in_any_stack),
    pct_unique_links_dropped_in_all_stacks =
      100 * mean(dropped_in_all_stacks),
    pct_unique_links_dropped_in_any_stack =
      100 * mean(dropped_in_any_stack)
  )

cat("\nNumber of unique physical road links dropped:\n")
print(zero_identifier_summary)

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
  mutate(pct_within_group = 100 * link_stack_units / sum(link_stack_units)) %>%
  ungroup()

cat("\nDropped vs retained by treatment group:\n")
print(zero_by_group)

zero_by_scheme <- zero_exclusion_unit %>%
  group_by(scheme, group, zero_dropped) %>%
  summarise(
    link_stack_units = n(),
    unique_road_links = n_distinct(identifier),
    OAs = n_distinct(OA),
    .groups = "drop"
  ) %>%
  group_by(scheme, group) %>%
  mutate(pct_within_scheme_group = 100 * link_stack_units / sum(link_stack_units)) %>%
  ungroup() %>%
  arrange(scheme, group, zero_dropped)

cat("\nDropped vs retained by scheme and group:\n")
print(zero_by_scheme, n = Inf)

possible_characteristic_vars <- c(
  "prop_in_caz",
  "road_length",
  "length_m",
  "length_km",
  "mean_total_pkm",
  "trend_total_pkm",
  "recent_minus_mid_total_pkm",
  "mid_minus_early_total_pkm"
)

characteristic_vars <- intersect(possible_characteristic_vars, names(model_panel_zero_diag))

if (length(characteristic_vars) > 0) {
  road_characteristics <- model_panel_zero_diag %>%
    group_by(uid, identifier, OA, scheme, treat_group, group, zero_dropped) %>%
    summarise(
      across(
        all_of(characteristic_vars),
        ~ suppressWarnings(mean(.x, na.rm = TRUE)),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
  
  characteristic_summary <- road_characteristics %>%
    pivot_longer(
      cols = all_of(characteristic_vars),
      names_to = "characteristic",
      values_to = "value"
    ) %>%
    group_by(characteristic, zero_dropped) %>%
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      p25 = quantile(value, 0.25, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  characteristic_smd <- road_characteristics %>%
    pivot_longer(
      cols = all_of(characteristic_vars),
      names_to = "characteristic",
      values_to = "value"
    ) %>%
    group_by(characteristic) %>%
    summarise(
      mean_retained = mean(value[zero_dropped == FALSE], na.rm = TRUE),
      mean_dropped  = mean(value[zero_dropped == TRUE], na.rm = TRUE),
      sd_retained   = sd(value[zero_dropped == FALSE], na.rm = TRUE),
      sd_dropped    = sd(value[zero_dropped == TRUE], na.rm = TRUE),
      pooled_sd     = sqrt((sd_retained^2 + sd_dropped^2) / 2),
      smd_dropped_vs_retained = (mean_dropped - mean_retained) / pooled_sd,
      abs_smd = abs(smd_dropped_vs_retained),
      .groups = "drop"
    ) %>%
    arrange(desc(abs_smd))
  
  cat("\nCharacteristics of retained vs dropped roads:\n")
  print(characteristic_summary, n = Inf)
  
  cat("\nStandardised differences: dropped vs retained roads:\n")
  print(characteristic_smd, n = Inf)
  
  write_csv(characteristic_summary,
            file.path(zero_exclusion_diag_dir, "zero_exclusion_characteristic_summary.csv"))
  write_csv(characteristic_smd,
            file.path(zero_exclusion_diag_dir, "zero_exclusion_characteristic_smd.csv"))
} else {
  cat("\nNo extra road-level characteristic variables found in model_panel.\n")
  cat("This is okay if model_panel only contains identifiers, treatment, scheme, quarter, and outcome.\n")
}

write_csv(zero_obs_summary,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_observation_summary.csv"))
write_csv(zero_unit_summary,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_link_stack_unit_summary.csv"))
write_csv(zero_identifier_summary,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_unique_road_link_summary.csv"))
write_csv(zero_by_group,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_by_group.csv"))
write_csv(zero_by_scheme,
          file.path(zero_exclusion_diag_dir, "zero_exclusion_by_scheme.csv"))

cat("\nZero-exclusion diagnostics saved to:\n")
cat(zero_exclusion_diag_dir, "\n")




# remove the zeros 
model_panel <- model_panel %>%
  filter(unit_total_injury_all_periods > 0) %>%
  select(-unit_total_injury_all_periods)

cat("\nRows after excluding all-zero road-link units:", nrow(model_panel), "\n")
cat("Units after exclusion:", n_distinct(model_panel$uid), "\n")

stopifnot(all(!is.na(model_panel$covid_period)))

cat("\nCOVID period distribution:\n")
print(summary(model_panel$covid_period))

rm(road_panel)
schemes_all <- sort(unique(model_panel$scheme))

# =============================================================================
# BASIC SAMPLE CHECKS
# =============================================================================

summary_by_group <- model_panel %>%
  group_by(group) %>%
  summarise(
    units        = n_distinct(uid_int),
    observations = n(),
    mean_injury  = mean(outcome_raw),
    pct_zero     = 100 * mean(outcome_raw == 0),
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
    mean_reuse       = mean(n_rows),
    max_reuse        = max(n_rows),
    .groups = "drop"
  )

cat("\nDuplicate rows from matched-control reuse (with replacement):\n")
print(duplicate_summary)

# =============================================================================
# RAW TREND PLOT
# =============================================================================

trend_df <- model_panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop")

p_trend <- ggplot(
  trend_df,
  aes(x = as.Date(quarter_year), y = mean_injury, colour = group)
) +
  annotate("rect",
           xmin = as.Date(as.yearqtr("2020 Q2")),
           xmax = as.Date(as.yearqtr("2021 Q3") + 0.25),
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey70"
  ) +
  annotate("rect",
           xmin = as.Date(as.yearqtr("2021 Q2")),
           xmax = as.Date(as.yearqtr("2021 Q4") + 0.25),
           ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "grey70"
  ) +
  geom_line(linewidth = 0.9) +
  labs(
    title    = "Mean injuries over time",
    subtitle = "Shaded periods: COVID lockdown/disruption and recovery",
    x        = NULL,
    y        = "Mean injuries per road-link-quarter",
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_trend)
ggsave(file.path(outdir, "00_raw_trends.png"),
       p_trend, width = 10, height = 6, dpi = 300)





# =============================================================================
# BUILD STACKED DATA
# =============================================================================
# Each CAZ scheme becomes a separate matched DiD stack. The same control road
# can appear in more than one stack but receives a stack-specific unit ID
# (uid_stack) so each stack is self-contained.

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
  mutate(stack_scheme = factor(stack_scheme))

cat("\nStacked rows:", nrow(stacked), "\n")
cat("Stacked units:", n_distinct(stacked$uid_stack), "\n")

# Check which quarter is event_time = -1 for each scheme
# (important to confirm the reference quarter is not COVID-contaminated)
ref_check <- stacked %>%
  filter(event_time == -1) %>%
  distinct(stack_scheme, quarter_year, covid_period) %>%
  arrange(stack_scheme, quarter_year)

cat("\nReference quarter (event_time = -1) per scheme:\n")
print(ref_check)

# Use the year preceding implementation as reference bin rather than
# a single quarter. Averaging over four quarters is more stable and
# avoids anchoring on a potentially noisy single quarter.
stacked <- stacked %>%
  mutate(
    event_time_ref = case_when(
      event_time >= -4 & event_time <= -1 ~ "ref_year",
      TRUE ~ as.character(event_time)
    ),
    event_time_ref = relevel(factor(event_time_ref), ref = "ref_year")
  )

# =============================================================================
# CONSTRUCT FINAL WEIGHTS
# =============================================================================
#
# Analysis weights
# ----------------
# Weight calculations happen after excluding all-zero road-link units, so the
# weights are specific to the outcome/sample being modelled.
#
# Main pooled estimand: equal-scheme average effect.
# Within each scheme, treated units sum to 1 and control units sum to 1.
# This keeps the duplicated matched-control rows, so the original matching
# ratios remain represented in the regression sample.
analysis_weight_counts <- stacked %>%
  distinct(stack_scheme, treat_group, uid_stack) %>%
  count(stack_scheme, treat_group, name = "n_units")

cat("\nUnits entering analysis weights, by scheme and treatment group:\n")
print(analysis_weight_counts)

stacked <- stacked %>%
  left_join(analysis_weight_counts, by = c("stack_scheme", "treat_group")) %>%
  mutate(analysis_weight = 1 / n_units) %>%
  select(-n_units)

# Link-weighted sensitivity estimand.
# Treated rows keep weight 1. Control rows are scaled so that, within each
# scheme, total control weight equals the number of treated units. This lets
# larger schemes contribute more while still balancing treated/control totals
# inside each scheme.
link_weight_counts <- analysis_weight_counts %>%
  pivot_wider(
    names_from  = treat_group,
    values_from = n_units,
    names_prefix = "grp_"
  ) %>%
  mutate(control_link_weight = grp_1 / grp_0) %>%
  select(stack_scheme, control_link_weight)

stacked <- stacked %>%
  left_join(link_weight_counts, by = "stack_scheme") %>%
  mutate(
    link_weight = if_else(treat_group == 1, 1, control_link_weight)
  ) %>%
  select(-control_link_weight)






# --- Verification ---
# In the main analysis weights, each scheme should have total treated weight
# of 1 and total control weight of 1 in each quarter.
cat("\nAnalysis weight verification (reference quarter, event_time = -1):\n")
stacked %>%
  filter(event_time == -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    n_rows              = n(),
    sum_analysis_weight = sum(analysis_weight),
    sum_link_weight     = sum(link_weight),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = treat_group,
    values_from = c(n_rows, sum_analysis_weight, sum_link_weight),
    names_prefix = "grp_"
  ) %>%
  mutate(
    raw_control_treated_ratio      = n_rows_grp_0 / n_rows_grp_1,
    analysis_control_treated_ratio = sum_analysis_weight_grp_0 / sum_analysis_weight_grp_1,
    link_control_treated_ratio     = sum_link_weight_grp_0 / sum_link_weight_grp_1
  ) %>%
  select(
    stack_scheme,
    raw_control_treated_ratio,
    analysis_control_treated_ratio,
    link_control_treated_ratio
  ) %>%
  print()






# Test whether a quadratic pre-trend differs between treated and control
# If this is not significant, you have evidence against trend acceleration differences


# =============================================================================
# MODEL 1 - POOLED AVERAGE TREATMENT EFFECT
# =============================================================================
# Estimates one average post-treatment CAZ effect across all schemes.
# Both main and sensitivity weighted versions are fitted:
#   - analysis_weight: equal-scheme pooled estimate; treated and controls
#     each sum to 1 within every scheme
#   - link_weight: link-weighted sensitivity; larger schemes contribute more,
#     while control totals are still balanced to treated totals within scheme
# Primary headline estimate is the COVID-adjusted analysis-weighted model.

m1_no_covid <- feglm(
  outcome_raw ~ post_stack |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

m1_covid <- feglm(
  outcome_raw ~ post_stack +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

# Link-weighted version for sensitivity comparison
m1_covid_link_weighted <- feglm(
  outcome_raw ~ post_stack +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~link_weight,
  lean    = TRUE
)

cat("\nModel 1: pooled average effect\n")
etable(
  m1_no_covid,
  m1_covid,
  m1_covid_link_weighted,
  headers = c("Analysis-weighted: no COVID", "Analysis-weighted: COVID-adjusted",
              "Link-weighted: COVID-adjusted"),
  dict = c(
    "post_stack" = "CAZ post-treatment"
  )
)

model1_results <- tibble(
  model = c("Analysis-weighted: no COVID interactions",
            "Analysis-weighted: COVID-adjusted",
            "Link-weighted: COVID-adjusted"),
  estimate_log_irr = c(
    coef(m1_no_covid)["post_stack"],
    coef(m1_covid)["post_stack"],
    coef(m1_covid_link_weighted)["post_stack"]
  ),
  se = c(
    se(m1_no_covid)["post_stack"],
    se(m1_covid)["post_stack"],
    se(m1_covid_link_weighted)["post_stack"]
  )
) %>%
  mutate(
    ci_lo      = estimate_log_irr - 1.96 * se,
    ci_hi      = estimate_log_irr + 1.96 * se,
    irr        = exp(estimate_log_irr),
    irr_lo     = exp(ci_lo),
    irr_hi     = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo     = 100 * (irr_lo - 1),
    pct_hi     = 100 * (irr_hi - 1)
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
    title    = "Model 1: pooled average CAZ effect",
    subtitle = "Stacked PPML - equal-scheme main estimate and link-weighted sensitivity",
    x        = "% change in injuries",
    y        = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_model1)
ggsave(file.path(outdir, "01_model1_pooled_average_effect.png"),
       p_model1, width = 8, height = 5, dpi = 300)

# Primary headline estimate
m1_main           <- m1_covid
model1_main_result <- model1_results %>%
  filter(model == "Analysis-weighted: COVID-adjusted")

# =============================================================================
# MODEL 2 - POOLED DYNAMIC EVENT STUDY (quarter -1 reference)
# =============================================================================
# Pooled lead and lag effects relative to quarter -1.
# Weighted by analysis_weight for equal scheme contribution.

m2_event <- feglm(
  outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model2_results <- extract_fixest_event_study(m2_event)

p_model2 <- plot_event_study(
  model2_results %>% filter(event_time >= -28, event_time <= 10),
  title    = "Model 2: pooled stacked PPML event study",
  subtitle = "Analysis-weighted; reference = quarter -1",
  ylab     = "Log incidence rate ratio"
)

print(p_model2)
ggsave(file.path(outdir, "02_model2_pooled_event_study.png"),
       p_model2, width = 10, height = 7, dpi = 300)

# =============================================================================
# MODEL 2c - POOLED DYNAMIC EVENT STUDY (year -1 reference bin)
# =============================================================================
# Reference period = average of quarters -4 to -1 (year before implementation).
# Preferred presentation over single-quarter reference because:
#   (a) averaging over four quarters reduces reference period noise
#   (b) the gap before implementation is visually explicit
#   (c) avoids anchoring on a quarter that may be COVID-contaminated
#       for some schemes (particularly Bradford, where quarters -8 to
#       -5 are lockdown/recovery)

m2_event_yearref <- feglm(
  outcome_raw ~ i(event_time_ref, treat_group, ref = "ref_year") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model2c_results <- extract_fixest_event_study_yearref(m2_event_yearref)

p_model2c <- plot_event_study(
  model2c_results %>% filter(event_time >= -28, event_time <= 10),
  title    = "Model 2c: pooled stacked PPML event study",
  subtitle = "Analysis-weighted; reference = year before implementation (Q\u22124 to Q\u22121)",
  ylab     = "Log incidence rate ratio"
)

print(p_model2c)
ggsave(file.path(outdir, "02c_model2_yearref_event_study.png"),
       p_model2c, width = 10, height = 7, dpi = 300)

p_compare <- p_model2 + p_model2c +
  plot_annotation(
    title    = "Comparison of reference period choice",
    subtitle = "Left: quarter \u22121 reference  |  Right: year \u22121 reference (Q\u22124 to Q\u22121)",
    theme    = theme_minimal(base_size = 12)
  )

print(p_compare)
ggsave(file.path(outdir, "02_compare_reference_periods.png"),
       p_compare, width = 18, height = 7, dpi = 300)

write_csv(model2c_results, file.path(outdir, "model2c_yearref_event_study.csv"))

# =============================================================================
# PARALLEL-TRENDS DIAGNOSTICS
# =============================================================================
# Wald tests on pre-treatment quarters from the year-reference model.
# Quarters -4 to -1 are the reference bin and do not appear as estimated
# coefficients; tests cover earlier pre-treatment windows only.
#
# The Wald tests below provide a formal joint test of pre-treatment
# event-study coefficients. Interpret them alongside the plotted pre-period
# coefficients, because large samples can reject small residual imbalances.

wald_pretrend_yearref <- function(model, max_k) {
  pattern <- paste0(
    "event_time_ref::-(", paste(max_k:5, collapse = "|"), "):treat_group"
  )
  wald(model, keep = pattern)
}

pt_8_5  <- wald_pretrend_yearref(m2_event_yearref, 8)
pt_12_5 <- wald_pretrend_yearref(m2_event_yearref, 12)
pt_16_5 <- wald_pretrend_yearref(m2_event_yearref, 16)
pt_20_5 <- wald_pretrend_yearref(m2_event_yearref, 20)
pt_24_5 <- wald_pretrend_yearref(m2_event_yearref, 24)

cat("\nParallel-trends Wald tests (year-reference model):\n")
cat("\nWindow -8 to -5:\n");  print(pt_8_5)
cat("\nWindow -12 to -5:\n"); print(pt_12_5)
cat("\nWindow -16 to -5:\n"); print(pt_16_5)
cat("\nWindow -20 to -5:\n"); print(pt_20_5)
cat("\nWindow -24 to -5:\n"); print(pt_24_5)

wald_summary <- tibble(
  window_label = c("Quarters -8 to -5", "Quarters -12 to -5",
                   "Quarters -16 to -5", "Quarters -20 to -5",
                   "Quarters -24 to -5"),
  max_k   = c(8, 12, 16, 20, 24),
  df1     = c(pt_8_5$df1, pt_12_5$df1, pt_16_5$df1, pt_20_5$df1, pt_24_5$df1),
  stat    = c(pt_8_5$stat, pt_12_5$stat, pt_16_5$stat, pt_20_5$stat, pt_24_5$stat),
  p_value = c(pt_8_5$p, pt_12_5$p, pt_16_5$p, pt_20_5$p, pt_24_5$p)
) %>%
  mutate(conclusion = if_else(p_value < 0.05, "Reject H0", "Fail to reject H0"))

cat("\nWald test summary:\n")
print(wald_summary)

write_csv(wald_summary, file.path(outdir, "parallel_trends_wald_tests.csv"))




stacked_pre <- stacked %>%
  filter(event_time < 0) %>%
  mutate(
    event_time_c  = event_time + 1,       # reference Q-1 = 0
    event_time_sq = event_time_c^2
  )

m_pretrend_test <- feglm(
  outcome_raw ~
    treat_group:event_time_c +
    treat_group:event_time_sq |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked_pre,
  family  = "poisson",
  weights = ~analysis_weight,
  cluster = ~OA,
  lean    = TRUE
)

etable(m_pretrend_test)
wald(m_pretrend_test, keep = "treat_group:event_time")

## The flexible event-study pre-trend tests indicate some residual quarter-specific pre-treatment imbalance, but a smoother parametric test finds no evidence of systematic linear or quadratic trend differences. This suggests that the violation is not a clear monotonic divergence, but may reflect noisy quarter-to-quarter deviations












# =============================================================================
# MODEL 3 - SCHEME-SPECIFIC AVERAGE TREATMENT EFFECTS
# =============================================================================
# One average post-treatment effect per scheme estimated jointly in a
# single model. The same analysis_weight column is used, so treated and
# control units are balanced within each scheme after all-zero exclusions.

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
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

cat("\nModel 3: scheme-specific average effects\n")
etable(m3_scheme)

ct_m3 <- coeftable(m3_scheme)

model3_results <- tibble(
  term             = rownames(ct_m3),
  estimate_log_irr = ct_m3[, "Estimate"],
  se               = ct_m3[, "Std. Error"],
  p_value          = ct_m3[, "Pr(>|z|)"]
) %>%
  filter(str_detect(term, "scheme_post::")) %>%
  mutate(
    scheme     = str_remove(term, "scheme_post::"),
    ci_lo      = estimate_log_irr - 1.96 * se,
    ci_hi      = estimate_log_irr + 1.96 * se,
    irr        = exp(estimate_log_irr),
    irr_lo     = exp(ci_lo),
    irr_hi     = exp(ci_hi),
    pct_change = 100 * (irr - 1),
    pct_lo     = 100 * (irr_lo - 1),
    pct_hi     = 100 * (irr_hi - 1),
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE            ~ ""
    )
  ) %>%
  select(scheme, estimate_log_irr, se, irr, irr_lo, irr_hi,
         pct_change, pct_lo, pct_hi, p_value, sig) %>%
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
    title    = "Model 3: scheme-specific average CAZ effects",
    subtitle = "Stacked PPML; treated/control weights balanced within scheme; OA-clustered SEs",
    x        = "% change in injuries",
    y        = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_model3)
ggsave(file.path(outdir, "03_model3_scheme_average_effects.png"),
       p_model3, width = 9, height = 6, dpi = 300)

# =============================================================================
# MODEL 4 - SCHEME-SPECIFIC DYNAMIC EVENT STUDIES
# =============================================================================
# Dynamic effects estimated separately for each scheme.
# Uses analysis_weight so treated and control units are balanced within scheme.

run_scheme_event_ppml <- function(sc) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
        i(covid_period, treat_group, ref = "non_pandemic") |
        uid_stack + stack_scheme[qtr_int],
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
  
  extract_fixest_event_study(fit) %>%
    mutate(scheme = sc)
}

model4_results <- map_dfr(schemes_all, run_scheme_event_ppml)

p_model4 <- ggplot(model4_results, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~scheme, scales = "free_y", ncol = 2) +
  labs(
    title    = "Model 4: scheme-specific PPML event studies",
    subtitle = "Analysis-weighted; reference = quarter -1",
    x        = "Quarters relative to CAZ implementation",
    y        = "Log incidence rate ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

print(p_model4)
ggsave(file.path(outdir, "04_model4_scheme_event_studies.png"),
       p_model4, width = 12, height = 10, dpi = 300)

# =============================================================================
# SENSITIVITY - SPECIFICATION CHECKS
# =============================================================================

# 1. Scheme x quarter FE vs linear trend per scheme
#    Confirms linear trend is sufficient (BIC favours it; estimates stable)
m1_twoway <- feglm(
  outcome_raw ~ post_stack +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme^quarter_year,
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

cat("\nSensitivity: linear trend vs scheme x quarter FE\n")
etable(m1_covid, m1_twoway,
       headers = c("Linear trend per scheme (main)", "Scheme x Quarter FE"))

# 2. Excluding Bradford
#    Bradford has a short post-treatment window (9 quarters), COVID lockdown
#    quarters falling inside its pre-treatment event window, and a national
#    control group decline that the linear time trend partially absorbs.
#    Excluding Bradford tests whether it unduly influences the pooled estimate.
m1_excl_bradford <- feglm(
  outcome_raw ~ post_stack +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked %>% filter(stack_scheme != "Bradford"),
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

cat("\nSensitivity: pooled estimate excluding Bradford\n")
etable(m1_covid, m1_excl_bradford,
       headers = c("All schemes", "Excluding Bradford"))

# =============================================================================
# OUTPUTS
# =============================================================================

write_csv(model1_results, file.path(outdir, "model1_pooled_average_effect.csv"))
write_csv(model2_results, file.path(outdir, "model2_pooled_event_study.csv"))
write_csv(model2c_results, file.path(outdir, "model2c_yearref_event_study.csv"))
write_csv(model3_results, file.path(outdir, "model3_scheme_average_effects.csv"))
write_csv(model4_results, file.path(outdir, "model4_scheme_event_studies.csv"))
write_csv(wald_summary, file.path(outdir, "parallel_trends_wald_tests.csv"))

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
    model1_results     = model1_results,
    model1_main_result = model1_main_result,
    model2_results     = model2_results,
    model2c_results    = model2c_results,
    model3_results     = model3_results,
    model4_results     = model4_results,
    wald_summary       = wald_summary
  ),
  here("data", "processed", "caz_ppml_main_results.rds")
)

cat("\nOutputs saved to:", outdir, "\n")

# =============================================================================
# BRADFORD DIAGNOSTIC INVESTIGATION
# =============================================================================
# Bradford (CAZ opened 2022 Q4) shows an anomalous positive estimate (+43%)
# that is not consistent with the other six schemes. Diagnostics below
# establish that this reflects control group decline, not treated road
# deterioration, and that the issue is scheme-wide rather than driven by
# specific OAs or high-reuse controls.

bradford_stack <- stacked %>%
  filter(stack_scheme == "Bradford")

cat("\n--- Bradford: basic sample ---\n")
cat("Treated units:", n_distinct(bradford_stack$uid_stack[bradford_stack$treat_group == 1]), "\n")
cat("Control units:", n_distinct(bradford_stack$uid_stack[bradford_stack$treat_group == 0]), "\n")
cat("Total rows:", nrow(bradford_stack), "\n")
cat("Implementation quarter:",
    as.character(unique(bradford_stack$caz_start_q[bradford_stack$treat_group == 1])), "\n")

cat("\n--- Bradford: sparsity by group and period ---\n")
bradford_stack %>%
  mutate(period = if_else(event_time >= 0, "post", "pre")) %>%
  group_by(treat_group, period) %>%
  summarise(
    obs          = n(),
    mean_injury  = mean(outcome_raw),
    pct_zero     = 100 * mean(outcome_raw == 0),
    total_injury = sum(outcome_raw),
    .groups = "drop"
  ) %>%
  print()

cat("\n--- Bradford: raw quarterly trends ---\n")
bradford_trend <- bradford_stack %>%
  mutate(group = if_else(treat_group == 1, "Treated", "Control")) %>%
  group_by(group, quarter_year) %>%
  summarise(
    mean_injury  = mean(outcome_raw),
    total_injury = sum(outcome_raw),
    n_roads      = n_distinct(uid_stack),
    .groups = "drop"
  )

p_bradford_trend <- ggplot(
  bradford_trend,
  aes(x = as.Date(quarter_year), y = mean_injury, colour = group)
) +
  annotate("rect",
           xmin = as.Date(as.yearqtr("2020 Q2")),
           xmax = as.Date(as.yearqtr("2021 Q1") + 0.25),
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey50"
  ) +
  annotate("rect",
           xmin = as.Date(as.yearqtr("2021 Q2")),
           xmax = as.Date(as.yearqtr("2021 Q4") + 0.25),
           ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "grey70"
  ) +
  geom_vline(
    xintercept = as.Date(as.yearqtr("2022 Q4")),
    linetype = "dashed", colour = "black", linewidth = 0.8
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  labs(
    title    = "Bradford: raw injury trends",
    subtitle = "Dashed = CAZ implementation (2022 Q4); shaded = COVID periods",
    x        = NULL,
    y        = "Mean injuries per road-link-quarter",
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_bradford_trend)
ggsave(file.path(outdir, "diag_bradford_raw_trends.png"),
       p_bradford_trend, width = 10, height = 6, dpi = 300)

cat("\n--- Bradford: mean injury by event time and group ---\n")
bradford_stack %>%
  mutate(group = if_else(treat_group == 1, "Treated", "Control")) %>%
  group_by(group, event_time) %>%
  summarise(
    mean_injury  = mean(outcome_raw),
    total_injury = sum(outcome_raw),
    pct_zero     = 100 * mean(outcome_raw == 0),
    .groups = "drop"
  ) %>%
  filter(event_time >= -8, event_time <= 8) %>%
  pivot_wider(
    names_from  = group,
    values_from = c(mean_injury, total_injury, pct_zero)
  ) %>%
  arrange(event_time) %>%
  print(n = 20)

bradford_event <- model4_results %>% filter(scheme == "Bradford")

p_bradford_event <- ggplot(
  bradford_event %>% filter(event_time >= -20, event_time <= 8),
  aes(x = event_time, y = estimate)
) +
  geom_hline(yintercept = 0,    linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "#E74C3C") +
  geom_line(linewidth = 0.8, colour = "#E74C3C") +
  geom_point(size = 1.8,   colour = "#E74C3C") +
  scale_x_continuous(breaks = pretty(bradford_event$event_time, n = 10)) +
  labs(
    title    = "Bradford: scheme-specific PPML event study",
    subtitle = "Analysis-weighted; reference = quarter -1",
    x        = "Quarters relative to CAZ implementation (2022 Q4)",
    y        = "Log incidence rate ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_bradford_event)
ggsave(file.path(outdir, "diag_bradford_event_study.png"),
       p_bradford_event, width = 10, height = 6, dpi = 300)

cat("\n--- Bradford: control unit reuse ---\n")
bradford_stack %>%
  filter(treat_group == 0) %>%
  count(identifier, name = "n_appearances") %>%
  count(n_appearances, name = "n_roads") %>%
  mutate(pct = 100 * n_roads / sum(n_roads)) %>%
  print()

cat("\n--- Bradford: top 10 OAs by total post-treatment injuries (treated) ---\n")
bradford_stack %>%
  filter(treat_group == 1, event_time >= 0) %>%
  group_by(OA) %>%
  summarise(
    total_injury = sum(outcome_raw),
    n_roads      = n_distinct(identifier),
    mean_injury  = mean(outcome_raw),
    .groups = "drop"
  ) %>%
  arrange(desc(total_injury)) %>%
  slice_head(n = 10) %>%
  print()

cat("\n--- Bradford: OA overlap between treated and control ---\n")
bradford_oas_treated <- bradford_stack %>% filter(treat_group == 1) %>% distinct(OA)
bradford_oas_control <- bradford_stack %>% filter(treat_group == 0) %>% distinct(OA)
cat("Treated OAs:", nrow(bradford_oas_treated), "\n")
cat("Control OAs:", nrow(bradford_oas_control), "\n")
cat("Overlapping OAs:",
    nrow(inner_join(bradford_oas_treated, bradford_oas_control, by = "OA")), "\n")

cat("\n--- Bradford: COVID period coverage in event window ---\n")
bradford_stack %>%
  filter(treat_group == 1) %>%
  count(event_time, covid_period) %>%
  filter(event_time >= -8, event_time <= 8) %>%
  pivot_wider(names_from = covid_period, values_from = n, values_fill = 0) %>%
  arrange(event_time) %>%
  print(n = 20)

# Indexed pre/post comparison - confirms control decline drives the result,
# not treated road deterioration
bradford_indexed <- bradford_stack %>%
  mutate(
    group  = if_else(treat_group == 1, "Treated", "Control"),
    period = if_else(event_time >= 0, "post", "pre")
  ) %>%
  group_by(group, period) %>%
  summarise(mean_inj = mean(outcome_raw), .groups = "drop") %>%
  group_by(group) %>%
  mutate(index = mean_inj / mean_inj[period == "pre"] * 100) %>%
  ungroup()

cat("\n--- Bradford: pre/post indexed comparison ---\n")
print(bradford_indexed)

# Check whether same-scheme OA contamination exists
# (treated OAs appearing as controls in their own scheme)
cat("\n--- Same-scheme OA contamination check (all schemes) ---\n")
same_scheme_overlap <- stacked %>%
  group_by(stack_scheme, OA) %>%
  summarise(
    appears_treated = any(treat_group == 1),
    appears_control = any(treat_group == 0),
    .groups = "drop"
  ) %>%
  filter(appears_treated & appears_control)

print(same_scheme_overlap)
cat("Contaminated scheme-OA pairs:", nrow(same_scheme_overlap), "\n")

# Cross-scheme control contamination - Bradford control OAs that are
# treated in another scheme's stack (different issue from above)
bradford_controls_also_treated <- bradford_stack %>%
  filter(treat_group == 0) %>%
  distinct(OA) %>%
  inner_join(
    stacked %>%
      filter(treat_group == 1) %>%
      distinct(OA, stack_scheme),
    by = "OA"
  )

cat("\nBradford control OAs that are treated in another scheme:\n")
print(bradford_controls_also_treated)

# Control collapse is universal across schemes - confirms national trend,
# not Bradford-specific matching failure
cat("\n--- Control post/pre injury index by scheme ---\n")
stacked %>%
  filter(treat_group == 0) %>%
  mutate(period = if_else(event_time >= 0, "post", "pre")) %>%
  group_by(stack_scheme, period) %>%
  summarise(
    injury_per_qtr = sum(outcome_raw) / n_distinct(quarter_year),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period, values_from = injury_per_qtr) %>%
  mutate(index = 100 * post / pre) %>%
  arrange(index) %>%
  print()

# Pre-treatment control-to-treated injury ratio by scheme
# (confirms Bradford's level mismatch is not the outlier - Portsmouth is)
cat("\n--- Control-to-treated injury ratio in reference year by scheme ---\n")
stacked %>%
  filter(event_time >= -4, event_time <= -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    injury_per_qtr = sum(outcome_raw) / n_distinct(quarter_year),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = treat_group,
    values_from = injury_per_qtr,
    names_prefix = "group_"
  ) %>%
  mutate(control_to_treated_ratio = group_0 / group_1) %>%
  arrange(desc(control_to_treated_ratio)) %>%
  print()
