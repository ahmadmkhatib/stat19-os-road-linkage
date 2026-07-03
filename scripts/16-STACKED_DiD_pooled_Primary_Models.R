# =============================================================================
# CAZ INJURY DID - PRIMARY PPML MODELS (FINAL)
# =============================================================================
#
# Main estimand
#   Equal-weighted average effect across CAZ schemes. Within each scheme,
#   treated links sum to weight 1 and matched controls sum to weight 1. The
#   pooled estimate is an average SCHEME effect, not an average road-link
#   effect.
#
# -----------------------------------------------------------------------------
# WHY THIS IS THE FINAL PRIMARY SPECIFICATION
# -----------------------------------------------------------------------------
# The primary workflow now focuses on pooled Models 1 and 2. Scheme-specific
# Models 3 and 4 are run at the end as supplementary heterogeneity diagnostics.
#
#   FIXED EVENT-STUDY REFERENCE FOR MODEL 2
#      The primary pooled event study uses the conventional fixed year
#      reference, event times -4:-1. This keeps every dynamic coefficient
#      anchored to the same relative pre-treatment window across schemes.
#      Because some schemes' reference years overlap COVID/recovery periods,
#      the script also reports a clean-reference Model 2 sensitivity. That
#      sensitivity is important for baseline validity, especially Bradford,
#      but it trades away common event-time comparability.
#
#   BUCKETED (NOT ROW-DROPPED) IMPLEMENTATION
#      An early attempt to isolate reference-vs-post comparisons filtered the
#      data down to ONLY the reference quarters plus the post-treatment window,
#      dropping everything in between (including all of COVID). This lost
#      ~90% of rows for some schemes and forced fixest to drop >10,000 unit
#      fixed effects as singletons, inflating standard errors for no benefit.
#      The final approach instead buckets every quarter that is neither the
#      fixed reference nor the common post-treatment window into a neutral
#      "other" category, absorbed via its own scheme-by-treatment dummy
#      (`other_flag:treat_scheme`). This keeps the FULL panel (same N as the
#      original, uncorrected models) while still cleanly isolating the
#      reference-vs-post-treatment comparison.
#
# Primary models
#   Model 1: Pooled average effect
#     Headline equal-scheme average post-treatment effect, using the bucketed
#     fixed-reference average-effect construction.
#   Model 2: Pooled fixed-reference event study
#     Dynamic effect path, pooled across schemes, reference = event times -4:-1.
#
# Sensitivity analyses repeat Models 1 and 2 after excluding Bradford and after
# restricting Bradford observations to 2021 Q2 onward. Scheme-specific event
# studies are also estimated with Bradford retained.
#

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

COMMON_POST_MAX <- 5   # common post-treatment window (quarters); all 7 schemes
# have data through event_time 5 (see post_event_scheme_counts
# in the support script)

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
        term, paste0("^", var_prefix, "::(-?\\d+):treat_group$")
      )[, 2] %>% as.numeric()
    ) %>%
    filter(!is.na(event_time)) %>%
    add_irr_columns() %>%
    arrange(event_time)
}

# NOTE: var_prefix must match the actual model term exactly (e.g.
# "scheme_post_bucket", not "scheme_post") - a wrong prefix silently returns
# 0 rows rather than erroring. Always sanity-check nrow() after calling this.
extract_scheme_effects <- function(model, var_prefix) {
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
      p_value < 0.05  ~ "*",   p_value < 0.10 ~ ".", TRUE ~ ""
    )) %>%
    select(scheme, estimate_log_irr, se, irr, irr_lo, irr_hi,
           pct_change, pct_lo, pct_hi, p_value, sig) %>%
    arrange(pct_change)
}

plot_event_study <- function(df, title, subtitle, colour = "#1f77b4") {
  ggplot(df, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = colour) +
    geom_line(linewidth = 0.8, colour = colour) +
    geom_point(size = 1.8, colour = colour) +
    scale_x_continuous(breaks = pretty(df$event_time, n = 10)) +
    labs(title = title, subtitle = subtitle,
         x = "Quarters relative to CAZ implementation",
         y = "Log incidence rate ratio") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# Anchors the DiD reference period to a verified-stable pre-COVID window
# (2019 Q2 - 2020 Q1) wherever the conventional -4:-1 window overlaps COVID.
make_clean_reference <- function(stacked_data) {
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
      event_time_ref_clean = if_else(clean_ref_year, "ref_year", as.character(event_time))
    ) %>%
    ungroup() %>%
    mutate(event_time_ref_clean = relevel(factor(event_time_ref_clean), ref = "ref_year"))
}

print_heading <- function(title) {
  cat("\n\n")
  cat("===============================================================================\n")
  cat(title, "\n", sep = "")
  cat("===============================================================================\n")
}

# =============================================================================
# 2. DATA PREPARATION
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

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

model_panel_all <- road_panel %>%
  select(panel_id, identifier, OA, scheme, quarter_year,
         caz_start_q, treat_group, all_of(outcome_var)) %>%
  rename(outcome_raw = all_of(outcome_var)) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  mutate(
    uid     = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    covid_period = factor(
      case_when(
        quarter_year >= as.yearqtr("2020 Q2") & quarter_year <= as.yearqtr("2021 Q1") ~ "lockdown",
        quarter_year >= as.yearqtr("2021 Q2") & quarter_year <= as.yearqtr("2021 Q4") ~ "recovery",
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

# All-zero road-link units carry no identifying information for FE PPML and
# would be dropped as singletons regardless; excluding them explicitly keeps
# the estimation sample transparent. Full diagnostics on this exclusion live
# in the support script (model_panel_for_zero_diag is saved below for it).
model_panel_for_zero_diag <- model_panel_all

model_panel <- model_panel_all %>%
  filter(unit_total_injury_all_periods > 0) %>%
  select(-unit_total_injury_all_periods)

rm(road_panel, road_caz_props, scheme_start, model_panel_all)

schemes_all <- sort(unique(model_panel$scheme))

# Stack: one copy of each scheme's data, event_time defined relative to that
# scheme's own implementation quarter. Makes staggered adoption comparable
# across schemes and avoids "forbidden comparison" bias.
stacked <- map_dfr(schemes_all, function(sc) {
  sc_start <- scheme_timing %>% filter(scheme == sc) %>% pull(caz_start_q)
  if (length(sc_start) == 0 || is.na(sc_start)) return(NULL)
  sc_start_int <- as.integer(round((as.numeric(sc_start) - min_qtr) * 4)) + 1L
  
  model_panel %>%
    filter(scheme == sc) %>%
    mutate(
      stack_scheme = sc,
      event_time   = qtr_int - sc_start_int,
      uid_stack    = paste0(uid_int, "_", sc)
    )
}) %>%
  mutate(
    stack_scheme = factor(stack_scheme),
    treat_scheme = interaction(treat_group, stack_scheme, drop = TRUE)
  )

# Equal-scheme analysis weights: treated units sum to 1 and control units sum
# to 1 within each scheme, so no scheme dominates the pooled estimate by
# virtue of having more road links.
analysis_weight_counts <- stacked %>%
  distinct(stack_scheme, treat_group, uid_stack) %>%
  count(stack_scheme, treat_group, name = "n_units")

stacked <- stacked %>%
  left_join(analysis_weight_counts, by = c("stack_scheme", "treat_group")) %>%
  mutate(analysis_weight = 1 / n_units) %>%
  select(-n_units)

# Apply the clean-reference correction for sensitivity/diagnostic outputs.
stacked <- make_clean_reference(stacked)

stopifnot(
  "clean_ref_year" %in% names(stacked),
  "event_time_ref_clean" %in% names(stacked)
)

clean_reference_table <- stacked %>%
  filter(clean_ref_year) %>%
  distinct(stack_scheme, quarter_year, event_time, covid_period) %>%
  arrange(stack_scheme, quarter_year)

stacked <- stacked %>%
  mutate(
    fixed_ref_year = event_time >= -4 & event_time <= -1,
    event_time_ref = if_else(fixed_ref_year, "ref_year", as.character(event_time)),
    event_time_ref = relevel(factor(event_time_ref), ref = "ref_year")
  )

fixed_reference_table <- stacked %>%
  filter(fixed_ref_year) %>%
  distinct(stack_scheme, quarter_year, event_time, covid_period) %>%
  arrange(stack_scheme, quarter_year)

# Bucketed variables for Model 1: every quarter that is neither the
# fixed -4:-1 reference nor the common post-treatment window is absorbed via its
# own scheme x treatment dummy ("other_flag:treat_scheme") rather than
# dropped, preserving the full panel. See header note for why this matters.
stacked <- stacked %>%
  mutate(
    period_bucket = case_when(
      fixed_ref_year ~ "ref_year",
      treat_group == 1 & event_time >= 0 & event_time <= COMMON_POST_MAX ~ "post_common",
      TRUE ~ "other"
    ),
    post_common = as.integer(period_bucket == "post_common"),
    other_flag  = as.integer(period_bucket == "other"),
    scheme_post_bucket = if_else(post_common == 1, as.character(stack_scheme), "ref_year"),
    scheme_post_bucket = factor(scheme_post_bucket, levels = c("ref_year", schemes_all))
  )

sample_summary <- model_panel %>%
  group_by(group) %>%
  summarise(units = n_distinct(uid), observations = n(),
            total_RTI = sum(outcome_raw, na.rm = TRUE),
            mean_RTI_per_road_quarter = mean(outcome_raw, na.rm = TRUE),
            pct_zero = 100 * mean(outcome_raw == 0, na.rm = TRUE), .groups = "drop")

scheme_sample_summary <- model_panel %>%
  distinct(scheme, uid_int, treat_group) %>%
  group_by(scheme) %>%
  summarise(treated_units = sum(treat_group == 1), control_units = sum(treat_group == 0),
            .groups = "drop")

# =============================================================================
# 3. PRIMARY MODEL 1 - POOLED AVERAGE EFFECT
# =============================================================================
#
# Why this model is primary:
#   Headline equal-scheme average effect. Uses the bucketed fixed-reference
#   construction (post_common vs. event-time -4:-1 ref_year, with intervening quarters
#   absorbed via other_flag:treat_scheme) plus flexible scheme-by-treatment
#   COVID adjustment.
#
# MAIN CAVEAT: this pooled number averages scheme effects that run in
# opposite directions and are individually large and significant (see Model
# 3). A near-zero pooled effect here does NOT mean "no effect" - always
# report this alongside Model 3, never as a substitute for it.

m1_pooled_average <- feglm(
  outcome_raw ~
    post_common +
    other_flag:treat_scheme +
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
  model = "Model 1: pooled average, equal-scheme, fixed -4:-1 reference (bucketed)",
  estimate_log_irr = coef(m1_pooled_average)["post_common"],
  se = se(m1_pooled_average)["post_common"]
) %>%
  add_irr_columns("estimate_log_irr", "se")

# =============================================================================
# 4. PRIMARY MODEL 2 - POOLED FIXED-REFERENCE EVENT STUDY
# =============================================================================
#
# Why this model is primary:
#   Shows the dynamic effect path using the conventional fixed year-reference
#   window, event times -4:-1. This preserves a common event-time baseline
#   across schemes. The clean-reference version is retained immediately below
#   as a sensitivity check because some schemes' fixed reference years overlap
#   COVID/recovery periods.

m2_pooled_event_fixedref <- feglm(
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

model2_results <- extract_event_study(m2_pooled_event_fixedref, "event_time_ref")

model2_post_common_wald <- wald(
  m2_pooled_event_fixedref,
  keep = paste0("event_time_ref::[0-", COMMON_POST_MAX, "]:treat_group")
)
# NOTE: this pooled joint test dilutes any single scheme's signal with the
# other six (mostly null) schemes - it is NOT informative about any
# individual scheme. Use the within-scheme tests in Model 4 for that.

m2_pooled_event_cleanref_sensitivity <- feglm(
  outcome_raw ~
    i(event_time_ref_clean, treat_group, ref = "ref_year") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model2_cleanref_results <- extract_event_study(
  m2_pooled_event_cleanref_sensitivity,
  "event_time_ref_clean"
)

model2_cleanref_post_common_wald <- wald(
  m2_pooled_event_cleanref_sensitivity,
  keep = paste0("event_time_ref_clean::[0-", COMMON_POST_MAX, "]:treat_group")
)

# =============================================================================
# 5. SENSITIVITY ANALYSES FOR MODELS 1 AND 2
# =============================================================================

extract_model1_result <- function(model, label) {
  tibble(
    spec = label,
    estimate_log_irr = coef(model)["post_common"],
    se = se(model)["post_common"]
  ) %>%
    add_irr_columns("estimate_log_irr", "se")
}

run_models_1_2 <- function(data, label) {
  d <- data %>% droplevels()
  
  m1 <- feglm(
    outcome_raw ~
      post_common +
      other_flag:treat_scheme +
      i(covid_period, treat_scheme, ref = "non_pandemic") |
      uid_stack +
      stack_scheme[qtr_int],
    data    = d,
    family  = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean    = TRUE
  )
  
  m2 <- feglm(
    outcome_raw ~
      i(event_time_ref, treat_group, ref = "ref_year") +
      i(covid_period, treat_scheme, ref = "non_pandemic") |
      uid_stack +
      stack_scheme[qtr_int],
    data    = d,
    family  = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean    = TRUE
  )
  
  list(
    label = label,
    model1 = m1,
    model2 = m2,
    model1_results = extract_model1_result(m1, label),
    model2_results = extract_event_study(m2, "event_time_ref") %>%
      mutate(spec = label, .before = 1),
    model2_post_wald = wald(
      m2,
      keep = paste0("event_time_ref::[0-", COMMON_POST_MAX, "]:treat_group")
    )
  )
}

# Bradford's fixed -4:-1 reference begins in 2021 Q2. Change this to
# "2021 Q3" if you want the looser version of this sensitivity.
BRADFORD_RESTRICT_START_QTR <- as.yearqtr("2021 Q2")

stacked_no_bradford <- stacked %>%
  filter(stack_scheme != "Bradford") %>%
  droplevels()

stacked_bradford_post_2021q2 <- stacked %>%
  filter(stack_scheme != "Bradford" | quarter_year >= BRADFORD_RESTRICT_START_QTR) %>%
  droplevels()

sensitivity_no_bradford <- run_models_1_2(
  stacked_no_bradford,
  "Sensitivity: exclude Bradford"
)

sensitivity_bradford_post_2021q2 <- run_models_1_2(
  stacked_bradford_post_2021q2,
  "Sensitivity: Bradford pre-period restricted to >= 2021 Q2"
)

model1_comparison <- bind_rows(
  extract_model1_result(m1_pooled_average, "Primary: all schemes"),
  sensitivity_no_bradford$model1_results,
  sensitivity_bradford_post_2021q2$model1_results
) %>%
  select(spec, estimate_log_irr, se, irr, irr_lo, irr_hi,
         pct_change, pct_lo, pct_hi)

model2_post_comparison <- bind_rows(
  model2_results %>% mutate(spec = "Primary: all schemes", .before = 1),
  sensitivity_no_bradford$model2_results,
  sensitivity_bradford_post_2021q2$model2_results,
  model2_cleanref_results %>% mutate(spec = "Sensitivity: clean/flexible reference", .before = 1)
) %>%
  filter(event_time >= 0, event_time <= COMMON_POST_MAX) %>%
  select(spec, event_time, estimate, se, irr, irr_lo, irr_hi,
         pct_change, pct_lo, pct_hi) %>%
  arrange(event_time, spec)

model2_wald_comparison <- tibble(
  spec = c(
    "Primary: all schemes",
    sensitivity_no_bradford$label,
    sensitivity_bradford_post_2021q2$label,
    "Sensitivity: clean/flexible reference"
  ),
  stat = c(
    model2_post_common_wald$stat,
    sensitivity_no_bradford$model2_post_wald$stat,
    sensitivity_bradford_post_2021q2$model2_post_wald$stat,
    model2_cleanref_post_common_wald$stat
  ),
  df1 = c(
    model2_post_common_wald$df1,
    sensitivity_no_bradford$model2_post_wald$df1,
    sensitivity_bradford_post_2021q2$model2_post_wald$df1,
    model2_cleanref_post_common_wald$df1
  ),
  df2 = c(
    model2_post_common_wald$df2,
    sensitivity_no_bradford$model2_post_wald$df2,
    sensitivity_bradford_post_2021q2$model2_post_wald$df2,
    model2_cleanref_post_common_wald$df2
  ),
  p_value = c(
    model2_post_common_wald$p,
    sensitivity_no_bradford$model2_post_wald$p,
    sensitivity_bradford_post_2021q2$model2_post_wald$p,
    model2_cleanref_post_common_wald$p
  )
)

run_scheme_event_fixedref <- function(sc) {
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
      cat("Scheme fixed-reference event study failed for:", sc, "\n")
      cat("Reason:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  list(
    results = extract_event_study(fit, "event_time_ref") %>% mutate(scheme = sc),
    model = fit
  )
}

scheme_event_fixedref_fits <- map(schemes_all, run_scheme_event_fixedref)
names(scheme_event_fixedref_fits) <- schemes_all

scheme_event_fixedref_results <- map_dfr(scheme_event_fixedref_fits, "results")

scheme_event_fixedref_post_wald <- map_dfr(schemes_all, function(sc) {
  fit <- scheme_event_fixedref_fits[[sc]]$model
  if (is.null(fit)) return(NULL)
  w <- tryCatch(
    wald(fit, keep = paste0("event_time_ref::[0-", COMMON_POST_MAX, "]:treat_group")),
    error = function(e) NULL
  )
  if (is.null(w)) return(NULL)
  tibble(scheme = sc, stat = w$stat, df1 = w$df1, df2 = w$df2, p_value = w$p) %>%
    mutate(conclusion = if_else(p_value < 0.05, "Reject H0 (joint effect significant)",
                                "Fail to reject H0"))
})

# -----------------------------------------------------------------------------
# Annual (policy-facing) event study: Year 1, Year 2 post-implementation
# -----------------------------------------------------------------------------

MAX_EVENT_YEAR <- 2   # Year 1 = quarters 0-3 (fully balanced, 7 schemes);
# Year 2 = quarters 4-7 (partially balanced, see below)

stacked_annual <- stacked %>%
  mutate(
    event_year = case_when(
      treat_group == 1 & event_time >= 0 & event_time <= 3  ~ "year1",
      treat_group == 1 & event_time >= 4 & event_time <= 7  ~ "year2",
      TRUE ~ NA_character_
    ),
    period_bucket_annual = case_when(
      fixed_ref_year ~ "ref_year",
      !is.na(event_year) ~ event_year,
      TRUE ~ "other"
    ),
    scheme_year_bucket = if_else(
      period_bucket_annual %in% c("year1", "year2"),
      paste0(as.character(stack_scheme), "_", period_bucket_annual),
      "ref_year"
    ),
    scheme_year_bucket = factor(scheme_year_bucket, levels = c(
      "ref_year",
      paste0(rep(schemes_all, each = 2), "_", c("year1", "year2"))
    )),
    other_flag_annual = as.integer(period_bucket_annual == "other")
  )

# Pooled annual event study (equivalent to Model 2, but yearly)
m2_annual <- feglm(
  outcome_raw ~
    i(period_bucket_annual, treat_group, ref = "ref_year") +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked_annual %>% filter(period_bucket_annual != "other" | treat_group == 0),
  family  = "poisson", cluster = ~OA,
  weights = ~analysis_weight, lean = TRUE
)
# NOTE: control-arm rows are always kept (they define "other" for treated only
# in the pretrend/COVID sense); treated "other" quarters need the same
# other_flag absorption as before if you want full-data retention rather than
# row-dropping - swap in other_flag_annual:treat_scheme the same way as your
# COMMON_POST_MAX models if you want the fully bucketed version.

etable(m2_annual)

# =============================================================================
# 6. PRIMARY FIGURES AND SAVED MODEL OBJECTS
# =============================================================================

p_model1 <- ggplot(model1_results, aes(x = pct_change, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(title = "Model 1: pooled average CAZ effect",
       subtitle = "Equal-scheme weighted PPML; fixed -4:-1 reference (bucketed); flexible scheme-by-treatment COVID adjustment",
       x = "% change in injuries", y = NULL) +
  theme_minimal(base_size = 12)

p_model2 <- plot_event_study(
  model2_results %>% filter(event_time >= -28, event_time <= 10),
  title = "Model 2: pooled fixed-reference event study",
  subtitle = "Reference = event times -4:-1; flexible scheme-by-treatment COVID adjustment"
)

ggsave(file.path(outdir, "primary_01_pooled_average_effect.png"), p_model1, width = 9, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "primary_02_fixedref_event_study.png"), p_model2, width = 10, height = 7, dpi = 300)

primary_results <- list(
  data = list(
    scheme_timing = scheme_timing,
    sample_summary = sample_summary,
    scheme_sample_summary = scheme_sample_summary,
    fixed_reference_table = fixed_reference_table,
    clean_reference_table = clean_reference_table
  ),
  models = list(
    m1_pooled_average = m1_pooled_average,
    m2_pooled_event_fixedref = m2_pooled_event_fixedref,
    m2_pooled_event_cleanref_sensitivity = m2_pooled_event_cleanref_sensitivity,
    sensitivity_no_bradford = sensitivity_no_bradford,
    sensitivity_bradford_post_2021q2 = sensitivity_bradford_post_2021q2
  ),
  tables = list(
    model1_results = model1_results,
    model2_results = model2_results,
    model2_cleanref_results = model2_cleanref_results,
    model1_comparison = model1_comparison,
    model2_post_comparison = model2_post_comparison,
    model2_wald_comparison = model2_wald_comparison,
    scheme_event_fixedref_results = scheme_event_fixedref_results,
    scheme_event_fixedref_post_wald = scheme_event_fixedref_post_wald
  ),
  tests = list(
    model2_post_common_wald = model2_post_common_wald,
    model2_cleanref_post_common_wald = model2_cleanref_post_common_wald
  )
)

saveRDS(primary_results, here("data", "processed", "caz_primary_ppml_results.rds"))

# Objects needed by the companion support script, saved separately so it
# doesn't need to rebuild the data pipeline from scratch.
saveRDS(
  list(
    stacked = stacked,
    model_panel = model_panel,
    model_panel_for_zero_diag = model_panel_for_zero_diag,
    schemes_all = schemes_all,
    scheme_timing = scheme_timing,
    min_qtr = min_qtr,
    COMMON_POST_MAX = COMMON_POST_MAX
  ),
  here("data", "processed", "caz_support_script_inputs.rds")
)

# =============================================================================
# 8. CONSOLE REPORT
# =============================================================================

print_heading("CAZ PRIMARY PPML ANALYSIS - CONSOLE REPORT (FINAL)")

cat("\nInterpretation guide:\n")
cat("  - Log IRR estimates below are converted to IRR and percent change in tables.\n")
cat("  - The pooled estimand is equal-weighted across schemes.\n")
cat("  - Model 2 uses a fixed -4:-1 year reference as primary.\n")
cat("  - The clean-reference Model 2 is reported as a sensitivity check.\n")
cat("  - Model 1 uses the same fixed -4:-1 reference in bucketed average-effect form.\n")
cat("  - Bradford sensitivity checks are printed immediately after Models 1 and 2.\n")
cat("  - Scheme-specific Models 3 and 4 are run at the end as supplementary outputs.\n")

print_heading("0. Scheme Timing And Analysis Sample")
cat("\nScheme timing after majority-quarter adjustment:\n")
print(scheme_timing, n = Inf)
cat("\nRetained active-link analysis sample:\n")
print(sample_summary, n = Inf)
cat("\nScheme sample sizes after excluding all-zero road-link units:\n")
print(scheme_sample_summary, n = Inf)
cat("\nFixed reference periods used in primary Model 2:\n")
print(fixed_reference_table, n = Inf)
cat("\nClean reference periods used in clean-reference sensitivity models:\n")
print(clean_reference_table, n = Inf)

print_heading("1. MODEL 1 - Pooled Average Effect")
etable(m1_pooled_average, dict = c("post_common" = "CAZ post-treatment"))
cat("\nModel 1 converted to IRR and percent change:\n")
print(model1_results, n = Inf)

print_heading("2. MODEL 2 - Pooled Fixed-Reference Event Study")
etable(m2_pooled_event_fixedref)
cat("\nAll Model 2 event-study estimates:\n")
print(model2_results, n = Inf)
cat("\nCommon post-treatment window 0-", COMMON_POST_MAX, ":\n", sep = "")
model2_results %>%
  filter(event_time >= 0, event_time <= COMMON_POST_MAX) %>%
  print(n = Inf)
cat("\nPooled joint Wald test, post-treatment window (diluted across all 7 schemes\n")
cat("- NOT informative about any single scheme; see Section 4 for within-scheme tests):\n")
print(model2_post_common_wald)
cat("\nClean-reference Model 2 sensitivity, common post-treatment window 0-", COMMON_POST_MAX, ":\n", sep = "")
model2_cleanref_results %>%
  filter(event_time >= 0, event_time <= COMMON_POST_MAX) %>%
  print(n = Inf)
cat("\nClean-reference pooled joint Wald test, post-treatment window:\n")
print(model2_cleanref_post_common_wald)

print_heading("3. MODEL 1 AND 2 SENSITIVITY COMPARISON")
cat("\nModel 1 primary vs. sensitivity specifications:\n")
print(model1_comparison, n = Inf)
cat("\nModel 2 post-treatment event-time estimates, primary vs. sensitivity specifications:\n")
print(model2_post_comparison, n = Inf)
cat("\nModel 2 joint Wald tests, post-treatment window 0-", COMMON_POST_MAX, ":\n", sep = "")
print(model2_wald_comparison, n = Inf)
cat("\nScheme-specific fixed-reference event-study Wald tests, Bradford retained:\n")
print(scheme_event_fixedref_post_wald, n = Inf)

print_heading("Primary Analysis Complete")
cat("Saved primary Models 1/2 object: data/processed/caz_primary_ppml_results.rds\n")
cat("Saved support-script inputs: data/processed/caz_support_script_inputs.rds\n")
cat("Saved primary figures in:", outdir, "\n")

# =============================================================================
# 9. SUPPLEMENTARY SCHEME-SPECIFIC MODELS 3 AND 4
# =============================================================================
#
# These are intentionally kept after the primary Models 1/2 workflow. They are
# useful for heterogeneity and scheme-level diagnostics, but the primary script
# now leads with pooled Models 1 and 2 plus Bradford sensitivity checks.

# Scheme-specific annual effects (equivalent to Model 3, but yearly)
m3_annual <- feglm(
  outcome_raw ~
    i(scheme_year_bucket, ref = "ref_year") +
    other_flag_annual:treat_scheme +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack + stack_scheme[qtr_int],
  data    = stacked_annual,
  family  = "poisson", cluster = ~OA,
  weights = ~analysis_weight, lean = TRUE
)

model3_annual_results <- extract_scheme_effects(m3_annual, "scheme_year_bucket") %>%
  mutate(
    event_year = str_extract(scheme, "year\\d$"),
    scheme = str_remove(scheme, "_year\\d$")
  )

# Model 3: scheme-specific average effects under the same bucketed fixed
# reference setup as primary Model 1.
m3_scheme_average <- feglm(
  outcome_raw ~
    i(scheme_post_bucket, ref = "ref_year") +
    other_flag:treat_scheme +
    i(covid_period, treat_scheme, ref = "non_pandemic") |
    uid_stack +
    stack_scheme[qtr_int],
  data    = stacked,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

model3_results <- extract_scheme_effects(m3_scheme_average, "scheme_post_bucket")

# Model 4: scheme-specific clean-reference event studies. These are
# supplementary diagnostics because the primary pooled Model 2 uses fixed
# -4:-1 reference for event-time comparability.
run_scheme_event_cleanref <- function(sc) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~
        i(event_time_ref_clean, treat_group, ref = "ref_year") +
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
      cat("Scheme clean-reference event study failed for:", sc, "\n")
      cat("Reason:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  list(
    results = extract_event_study(fit, "event_time_ref_clean") %>% mutate(scheme = sc),
    model = fit
  )
}

scheme_event_cleanref_fits <- map(schemes_all, run_scheme_event_cleanref)
names(scheme_event_cleanref_fits) <- schemes_all

model4_results <- map_dfr(scheme_event_cleanref_fits, "results")

scheme_post_wald <- map_dfr(schemes_all, function(sc) {
  fit <- scheme_event_cleanref_fits[[sc]]$model
  if (is.null(fit)) return(NULL)
  w <- tryCatch(
    wald(fit, keep = paste0("event_time_ref_clean::[0-", COMMON_POST_MAX, "]:treat_group")),
    error = function(e) NULL
  )
  if (is.null(w)) return(NULL)
  tibble(scheme = sc, stat = w$stat, df1 = w$df1, df2 = w$df2, p_value = w$p) %>%
    mutate(conclusion = if_else(p_value < 0.05, "Reject H0 (joint effect significant)",
                                "Fail to reject H0"))
})

p_model3 <- ggplot(
  model3_results,
  aes(x = pct_change, y = fct_reorder(scheme, pct_change))
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
  geom_point(size = 3) +
  labs(title = "Supplementary Model 3: scheme-specific average CAZ effects",
       subtitle = "Fixed -4:-1 reference (bucketed); same PPML structure as Model 1",
       x = "% change in injuries", y = NULL) +
  theme_minimal(base_size = 12)

p_model4 <- ggplot(model4_results, aes(x = event_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~scheme, scales = "free_y", ncol = 2) +
  labs(title = "Supplementary Model 4: scheme-specific clean-reference event studies",
       subtitle = "See scheme_post_wald for within-scheme joint significance on the post-treatment window",
       x = "Quarters relative to CAZ implementation",
       y = "Log incidence rate ratio") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "supplementary_03_scheme_average_effects.png"), p_model3, width = 9, height = 6, dpi = 300)
ggsave(file.path(outdir, "supplementary_04_scheme_cleanref_event_studies.png"), p_model4, width = 12, height = 10, dpi = 300)

supplementary_scheme_results <- list(
  models = list(
    m3_annual = m3_annual,
    m3_scheme_average = m3_scheme_average,
    scheme_event_cleanref_fits = scheme_event_cleanref_fits
  ),
  tables = list(
    model3_annual_results = model3_annual_results,
    model3_results = model3_results,
    model4_results = model4_results,
    scheme_post_wald = scheme_post_wald
  )
)

saveRDS(
  supplementary_scheme_results,
  here("data", "processed", "caz_supplementary_scheme_models.rds")
)

print_heading("Supplementary Scheme-Specific Models Complete")
cat("\nSupplementary Model 3 converted to IRR and percent change:\n")
print(model3_results, n = Inf)
cat("\nSupplementary annual scheme-specific effects:\n")
print(model3_annual_results, n = Inf)
cat("\nSupplementary Model 4 event-study estimates, event times -8 to 8:\n")
model4_results %>%
  filter(event_time >= -8, event_time <= 8) %>%
  arrange(scheme, event_time) %>%
  print(n = Inf)
cat("\nSupplementary within-scheme joint Wald tests, post-treatment window 0-", COMMON_POST_MAX, ":\n", sep = "")
print(scheme_post_wald, n = Inf)
cat("Saved supplementary scheme model object: data/processed/caz_supplementary_scheme_models.rds\n")
