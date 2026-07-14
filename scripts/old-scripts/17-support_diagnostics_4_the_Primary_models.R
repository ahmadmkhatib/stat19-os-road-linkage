# =============================================================================
# CAZ INJURY DID - SUPPORTING DIAGNOSTICS AND MODEL CHOICE CHECKS
# =============================================================================
#
# Purpose
#   This script supports the four primary PPML models in
#   outputs/caz_01_primary_models.R. It keeps the accessory checks organised and
#   explains why the primary choices are preferred.
#
# How to read this script
#   The primary script answers: "What are the main estimates?"
#   This support script answers: "Why should we trust those models more than the
#   alternatives, and what caveats remain?"
#
# Sections
#   A. Run/load primary models
#   B. Basic sample, weighting, and zero-exclusion diagnostics
#   C. Raw trends
#   D. COVID/reference/weight sensitivity checks
#   E. Parallel-trends and event-time composition diagnostics
#   F. Leave-one-scheme-out influence checks
#   G. Bradford diagnostic investigation
#   H. Save compact support object and print console report
#
# =============================================================================
# =============================================================================
# A. SUPPORT HELPERS
# =============================================================================

support_heading <- function(title) {
  cat("\n\n")
  cat("-------------------------------------------------------------------------------\n")
  cat(title, "\n", sep = "")
  cat("-------------------------------------------------------------------------------\n")
}

safe_wald <- function(model, keep) {
  tryCatch(
    wald(model, keep = keep),
    error = function(e) {
      list(error = conditionMessage(e), keep = keep)
    }
  )
}

run_saturated_check <- function(sc,
                                lo_q = "2020 Q1",
                                hi_q = "2022 Q4",
                                ref_q = NULL) {
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

# =============================================================================
# B. SAMPLE, WEIGHTING, AND ZERO-EXCLUSION DIAGNOSTICS
# =============================================================================

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
  mutate(pct_within_group = 100 * link_stack_units / sum(link_stack_units)) %>%
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
  mutate(pct_within_scheme_group = 100 * link_stack_units / sum(link_stack_units)) %>%
  ungroup()

weight_check <- stacked %>%
  filter(event_time == -1) %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    n_rows = n(),
    sum_analysis_weight = sum(analysis_weight),
    sum_link_weight = sum(link_weight),
    .groups = "drop"
  )

reference_year_by_scheme <- stacked %>%
  filter(event_time == -1) %>%
  distinct(stack_scheme, quarter_year, covid_period) %>%
  arrange(stack_scheme)

# =============================================================================
# C. RAW TRENDS
# =============================================================================

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

ggsave(file.path(outdir, "support_00_raw_trends.png"),
       p_trend, width = 10, height = 6, dpi = 300)

# =============================================================================
# D. MODEL CHOICE SENSITIVITY: COVID, REFERENCE, AND WEIGHTS
# =============================================================================

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

model1_sensitivity_results <- tibble(
  model = c(
    "No COVID adjustment",
    "Pooled COVID adjustment",
    "Flexible scheme COVID adjustment",
    "Flexible scheme COVID, link-weighted"
  ),
  role = c(
    "Rejected support: ignores pandemic timing",
    "Support: partial COVID adjustment",
    "Primary: equal-scheme estimand with flexible COVID",
    "Support: changes estimand to link-weighted"
  ),
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

m2_pooled_event_yearref <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data = stacked,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
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

model2_yearref_results <- extract_event_study(m2_pooled_event_yearref, "event_time_ref")
model2_quarter_ref_results <- extract_event_study(m2_quarter_ref, "event_time_f")

reference_comparison_0_5 <- bind_rows(
  model2_yearref_results %>%
    mutate(spec = "Conventional year reference"),
  model2_results %>%
    mutate(spec = "Primary clean reference")
) %>%
  filter(event_time >= 0, event_time <= 5) %>%
  select(spec, event_time, estimate, se, pct_change, pct_lo, pct_hi)

# =============================================================================
# E. PARALLEL-TRENDS AND EVENT-TIME COMPOSITION DIAGNOSTICS
# =============================================================================

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

post_event_scheme_counts <- stacked %>%
  filter(event_time >= 0) %>%
  distinct(stack_scheme, event_time) %>%
  count(event_time, name = "n_schemes") %>%
  arrange(event_time)

p_composition <- ggplot(
  scheme_composition,
  aes(x = event_time, y = pct_of_bin, fill = stack_scheme)
) +
  geom_area(position = "stack") +
  labs(
    title = "Scheme composition of pre-treatment event-time bins",
    subtitle = "Long pre-period bins are not balanced across schemes",
    x = "Quarters relative to CAZ implementation",
    y = "% of observations in event-time bin",
    fill = "Scheme"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "support_pretrend_scheme_composition.png"),
       p_composition, width = 10, height = 6, dpi = 300)

# =============================================================================
# F. LEAVE-ONE-SCHEME-OUT INFLUENCE CHECKS
# =============================================================================

run_leave_one_out <- function(drop_scheme) {
  d <- stacked %>%
    filter(stack_scheme != drop_scheme) %>%
    droplevels()
  
  fit <- feglm(
    outcome_raw ~
      post_stack +
      i(covid_period, treat_scheme, ref = "non_pandemic") |
      uid_stack + stack_scheme[qtr_int],
    data = d,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = TRUE
  )
  
  tibble(
    dropped_scheme = drop_scheme,
    estimate_log_irr = coef(fit)["post_stack"],
    se = se(fit)["post_stack"]
  ) %>%
    add_irr_columns("estimate_log_irr", "se")
}

leave_one_out_results <- map_dfr(schemes_all, run_leave_one_out)

# =============================================================================
# G. BRADFORD DIAGNOSTIC INVESTIGATION
# =============================================================================

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

control_reuse <- stacked %>%
  filter(treat_group == 0) %>%
  distinct(identifier, stack_scheme) %>%
  count(identifier, name = "n_schemes_as_control")

bradford_shared_controls <- stacked %>%
  filter(
    stack_scheme == "Bradford",
    treat_group == 0,
    identifier %in% control_reuse$identifier[control_reuse$n_schemes_as_control > 1]
  ) %>%
  distinct(identifier) %>%
  nrow()

bradford_control_leverage <- bradford_stack %>%
  filter(treat_group == 0, event_time < 0) %>%
  group_by(identifier) %>%
  summarise(total_pre_injuries = sum(outcome_raw), .groups = "drop") %>%
  arrange(desc(total_pre_injuries)) %>%
  mutate(cum_share = cumsum(total_pre_injuries) / sum(total_pre_injuries))

top_leverage_ids <- head(bradford_control_leverage$identifier, 10)

bradford_trimmed <- bradford_stack %>%
  filter(!(treat_group == 0 & identifier %in% top_leverage_ids)) %>%
  droplevels()

m_bradford_trimmed <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + qtr_int,
  data = bradford_trimmed,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

bradford_trimmed_results <- extract_event_study(
  m_bradford_trimmed,
  "event_time_ref"
)

bradford_stack_ramp <- bradford_stack %>%
  mutate(
    qtrs_since_lockdown_start = pmax(
      0,
      qtr_int -
        (as.integer(round((as.numeric(as.yearqtr("2020 Q2")) - min_qtr) * 4)) + 1L)
    ),
    recovery_ramp = if_else(covid_period != "non_pandemic", qtrs_since_lockdown_start, 0)
  )

m_bradford_ramp <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") +
    treat_group:recovery_ramp |
    uid_stack + qtr_int,
  data = bradford_stack_ramp,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

bradford_ramp_results <- extract_event_study(m_bradford_ramp, "event_time_ref")

bradford_stack_clean <- bradford_stack %>%
  droplevels()

m_bradford_event_clean <- feglm(
  outcome_raw ~
    i(event_time_ref_clean, treat_group, ref = "ref_year") +
    i(covid_period, treat_group, ref = "non_pandemic") |
    uid_stack + qtr_int,
  data = bradford_stack_clean,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

bradford_clean_results <- extract_event_study(
  m_bradford_event_clean,
  "event_time_ref_clean"
)

bradford_clean_post_0_5_wald <- safe_wald(
  m_bradford_event_clean,
  "event_time_ref_clean::[0-5]:treat_group"
)

bradford_saturated <- bradford_stack %>%
  filter(
    quarter_year >= as.yearqtr("2020 Q1"),
    quarter_year <= as.yearqtr("2022 Q4")
  ) %>%
  mutate(qtr_f = factor(quarter_year))

m_bradford_saturated <- feglm(
  outcome_raw ~ i(qtr_f, treat_group, ref = as.character(as.yearqtr("2022 Q3"))) |
    uid_stack,
  data = bradford_saturated,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

bradford_spike_leverage <- bradford_stack %>%
  filter(event_time %in% c(-6, -5)) %>%
  group_by(identifier, treat_group) %>%
  summarise(spike_injuries = sum(outcome_raw), n_obs = n(), .groups = "drop") %>%
  arrange(desc(spike_injuries))

bradford_announce_test <- bradford_stack %>%
  mutate(
    # Replace this with the confirmed public announcement/confirmation date.
    announce_q = as.yearqtr("2021 Q3"),
    post_announce = as.integer(treat_group == 1 & quarter_year >= announce_q)
  )

m_bradford_announce <- feglm(
  outcome_raw ~ post_announce + post_stack |
    uid_stack + qtr_int,
  data = bradford_announce_test,
  family = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean = TRUE
)

saturated_check_models <- set_names(
  map(c("Bradford", "Newcastle", "Birmingham"), run_saturated_check),
  c("Bradford", "Newcastle", "Birmingham")
)

bradford_quarterly_means <- bradford_stack %>%
  filter(
    quarter_year >= as.yearqtr("2019 Q1"),
    quarter_year <= as.yearqtr("2022 Q4")
  ) %>%
  group_by(quarter_year, treat_group) %>%
  summarise(
    mean_injury = mean(outcome_raw),
    n_links = n_distinct(identifier),
    total_injury = sum(outcome_raw),
    .groups = "drop"
  ) %>%
  mutate(group = if_else(treat_group == 1, "Treated", "Control"))

bradford_indexed_series <- bradford_quarterly_means %>%
  group_by(group) %>%
  mutate(
    base = mean(mean_injury[
      quarter_year >= as.yearqtr("2019 Q1") &
        quarter_year <= as.yearqtr("2019 Q4")
    ]),
    index = 100 * mean_injury / base
  ) %>%
  ungroup()

p_bradford_indexed <- ggplot(
  bradford_indexed_series,
  aes(x = as.Date(quarter_year), y = index, colour = group)
) +
  geom_vline(xintercept = as.Date(as.yearqtr("2021 Q2")), linetype = "dotted") +
  geom_vline(xintercept = as.Date(as.yearqtr("2021 Q3")), linetype = "dotted") +
  geom_vline(xintercept = as.Date(as.yearqtr("2022 Q4")), linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Bradford: treated vs control, indexed to 2019 average",
    subtitle = "Dotted = event times -6/-5; dashed = actual CAZ start",
    x = NULL,
    y = "Index (2019 = 100)",
    colour = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(outdir, "support_bradford_indexed_series.png"),
       p_bradford_indexed, width = 10, height = 6, dpi = 300)

bradford_qoq <- bradford_quarterly_means %>%
  arrange(group, quarter_year) %>%
  group_by(group) %>%
  mutate(qoq_change = mean_injury - lag(mean_injury)) %>%
  ungroup() %>%
  filter(
    quarter_year >= as.yearqtr("2021 Q1"),
    quarter_year <= as.yearqtr("2021 Q4")
  )

bradford_ref_vs_precovid <- bradford_stack %>%
  mutate(period = case_when(
    quarter_year >= as.yearqtr("2019 Q1") &
      quarter_year <= as.yearqtr("2019 Q4") ~ "pre-COVID baseline",
    event_time >= -4 & event_time <= -1 ~ "reference window",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(period, treat_group) %>%
  summarise(mean_injury = mean(outcome_raw), .groups = "drop") %>%
  pivot_wider(names_from = treat_group, values_from = mean_injury, names_prefix = "grp_") %>%
  mutate(treat_control_ratio = grp_1 / grp_0)

# =============================================================================
# H. SAVE COMPACT SUPPORT OBJECT
# =============================================================================

support_results <- list(
  sample = list(
    zero_obs_summary = zero_obs_summary,
    zero_unit_summary = zero_unit_summary,
    zero_by_group = zero_by_group,
    zero_by_scheme = zero_by_scheme,
    weight_check = weight_check,
    reference_year_by_scheme = reference_year_by_scheme
  ),
  model_choice = list(
    model1_sensitivity_results = model1_sensitivity_results,
    model2_yearref_results = model2_yearref_results,
    model2_quarter_ref_results = model2_quarter_ref_results,
    reference_comparison_0_5 = reference_comparison_0_5,
    leave_one_out_results = leave_one_out_results
  ),
  pretrends = list(
    wald_summary = wald_summary,
    pretrend_heterog_covid_wald = pretrend_heterog_covid_wald,
    scheme_event_range = scheme_event_range,
    composition_wide = composition_wide,
    post_event_scheme_counts = post_event_scheme_counts
  ),
  bradford = list(
    bradford_period_summary = bradford_period_summary,
    bradford_indexed = bradford_indexed,
    same_scheme_overlap = same_scheme_overlap,
    control_post_pre_index = control_post_pre_index,
    reference_year_ratio = reference_year_ratio,
    bradford_control_leverage = bradford_control_leverage,
    bradford_trimmed_results = bradford_trimmed_results,
    bradford_ramp_results = bradford_ramp_results,
    bradford_clean_results = bradford_clean_results,
    bradford_clean_post_0_5_wald = bradford_clean_post_0_5_wald,
    bradford_spike_leverage = bradford_spike_leverage,
    bradford_qoq = bradford_qoq,
    bradford_ref_vs_precovid = bradford_ref_vs_precovid
  ),
  models = list(
    m1_no_covid = m1_no_covid,
    m1_pooled_covid = m1_pooled_covid,
    m1_link_weighted = m1_link_weighted,
    m2_pooled_event_yearref = m2_pooled_event_yearref,
    m2_quarter_ref = m2_quarter_ref,
    m_pretrend_heterog_covid = m_pretrend_heterog_covid,
    m_bradford_trimmed = m_bradford_trimmed,
    m_bradford_ramp = m_bradford_ramp,
    m_bradford_event_clean = m_bradford_event_clean,
    m_bradford_saturated = m_bradford_saturated,
    m_bradford_announce = m_bradford_announce,
    saturated_check_models = saturated_check_models
  )
)

saveRDS(support_results, here("data", "processed", "caz_support_diagnostics_results.rds"))

# =============================================================================
# I. CONSOLE REPORT
# =============================================================================

support_heading("SUPPORT SCRIPT: WHY THE PRIMARY MODELS WERE CHOSEN")
cat("\nPrimary model choice logic:\n")
cat("  1. No-COVID models are not credible because pandemic timing strongly affects injuries.\n")
cat("  2. Link-weighted models answer a different estimand and can be driven by larger schemes.\n")
cat("  3. Conventional year-reference event studies are useful but reference periods overlap COVID for some schemes.\n")
cat("  4. Clean-reference event studies are therefore primary for dynamics, with the distance-from-treatment caveat.\n")
cat("  5. Scheme-specific models are needed because pooled effects hide heterogeneity.\n")

support_heading("A. Sample, Weighting, And Reference Periods")
cat("\nReference quarter event_time = -1 by scheme:\n")
print(reference_year_by_scheme, n = Inf)
cat("\nClean reference periods used by the primary event-study models:\n")
print(clean_reference_table, n = Inf)
cat("\nEqual-scheme weight check at event_time = -1:\n")
print(weight_check, n = Inf)
cat("\nZero-exclusion summary:\n")
print(zero_obs_summary, n = Inf)
print(zero_unit_summary, n = Inf)
cat("\nZero exclusion by treatment group:\n")
print(zero_by_group, n = Inf)
cat("\nZero exclusion by scheme and treatment group:\n")
print(zero_by_scheme, n = Inf)

support_heading("B. Model 1 Sensitivity: COVID And Weighting")
cat("\nThis table explains why the primary pooled model uses flexible COVID terms and equal-scheme weights.\n")
print(model1_sensitivity_results, n = Inf)

support_heading("C. Event-Study Reference Sensitivity")
cat("\nConventional year-reference event-study estimates, common post window 0-5:\n")
model2_yearref_results %>%
  filter(event_time >= 0, event_time <= 5) %>%
  print(n = Inf)
cat("\nPrimary clean-reference event-study estimates, common post window 0-5:\n")
model2_results %>%
  filter(event_time >= 0, event_time <= 5) %>%
  print(n = Inf)
cat("\nSide-by-side reference comparison, event times 0-5:\n")
print(reference_comparison_0_5, n = Inf)

support_heading("D. Parallel-Trends Diagnostics")
cat("\nConventional year-reference pretrend Wald tests:\n")
print(wald_summary, n = Inf)
cat("\nFlexible scheme-by-treatment COVID pretrend Wald test:\n")
print(pretrend_heterog_covid_wald)
cat("\nInterpretation note:\n")
cat("  Conventional pretrend rejections are consistent with COVID/reference-period contamination.\n")
cat("  The flexible COVID-adjusted pretrend test is the more relevant support check.\n")

support_heading("E. Event-Time Composition")
cat("\nPre-period event-time range by scheme:\n")
print(scheme_event_range, n = Inf)
cat("\nPre-period scheme composition by event time:\n")
print(composition_wide, n = Inf)
cat("\nNumber of schemes contributing by post event time:\n")
print(post_event_scheme_counts, n = Inf)

support_heading("F. Leave-One-Scheme-Out Influence")
cat("\nPooled average model re-estimated after dropping one scheme at a time:\n")
print(leave_one_out_results, n = Inf)

support_heading("G. Bradford Diagnostics")
cat("\nBradford raw pre/post period summary:\n")
print(bradford_period_summary, n = Inf)
cat("\nBradford indexed pre/post comparison:\n")
print(bradford_indexed, n = Inf)
cat("\nSame-scheme OA overlap check:\n")
print(same_scheme_overlap, n = Inf)
cat("\nControl post/pre injury index by scheme:\n")
print(control_post_pre_index, n = Inf)
cat("\nControl-to-treated injury ratio in conventional reference year:\n")
print(reference_year_ratio, n = Inf)
cat("\nBradford control-link leverage, top 10 pre-period links:\n")
print(head(bradford_control_leverage, 10), n = Inf)
cat("\nBradford control links also used as controls elsewhere:", bradford_shared_controls, "\n")
cat("\nBradford event study excluding top-10 control leverage links, event times -8 to 8:\n")
bradford_trimmed_results %>%
  filter(event_time >= -8, event_time <= 8) %>%
  print(n = Inf)
cat("\nBradford event study with continuous COVID-recovery ramp, event times -8 to 8:\n")
bradford_ramp_results %>%
  filter(event_time >= -8, event_time <= 8) %>%
  print(n = Inf)
cat("\nBradford clean-reference event study, post event times 0-8:\n")
bradford_clean_results %>%
  filter(event_time >= 0, event_time <= 8) %>%
  print(n = Inf)
cat("\nBradford clean-reference joint post 0-5 Wald test:\n")
print(bradford_clean_post_0_5_wald)
cat("\nBradford saturated calendar-quarter x treatment model:\n")
etable(m_bradford_saturated)
cat("\nBradford event_time -6/-5 spike leverage, top 15 links:\n")
print(head(bradford_spike_leverage, 15), n = Inf)
cat("\nBradford announcement placeholder model; replace 2021 Q3 with confirmed announcement date:\n")
etable(m_bradford_announce)
cat("\nSaturated calendar-quarter x treatment checks for Bradford, Newcastle, Birmingham:\n")
for (sc in names(saturated_check_models)) {
  cat("\n--- ", sc, " ---\n", sep = "")
  print(etable(saturated_check_models[[sc]]))
}
cat("\nBradford quarter-on-quarter changes around 2021 recovery period:\n")
print(bradford_qoq, n = Inf)
cat("\nBradford reference window vs pre-COVID baseline ratio:\n")
print(bradford_ref_vs_precovid, n = Inf)

support_heading("SUPPORT SCRIPT COMPLETE")
cat("Saved support object: data/processed/caz_support_diagnostics_results.rds\n")
cat("Saved support figures in:", outdir, "\n")

