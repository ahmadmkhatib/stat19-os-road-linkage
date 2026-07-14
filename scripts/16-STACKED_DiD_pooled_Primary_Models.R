# =============================================================================
# CAZ INJURY DID -  PRIMARY PPML 
# =============================================================================
#
# 
#   1. Build a matched, stacked road-quarter panel.
#   2. Estimate the headline equal-scheme average post-CAZ effect.
#   3. Estimate the pooled event-study path.
#   4. Check whether Bradford or other schemes drive the result.
#   5. Check pre-treatment slope diagnostics after the revised matching.
#
# Main estimand
#   Equal-weighted average effect across CAZ schemes. Within each scheme,
#   treated roads sum to weight 1 and matched controls sum to weight 1, so the
#   pooled estimate is an average scheme effect, not an average road-link effect.
#
# =============================================================================

library(tidyverse)
library(arrow)
library(fixest)
library(sf)
library(here)
library(lubridate)
library(zoo)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
count  <- dplyr::count




OUTCOME_VAR <- "total_inj_adj_All"
COMMON_POST_MAX <- 5
OUTDIR <- here("output", "pooled", "All_clean")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Console Helpers
# =============================================================================

section <- function(title) {
  cat("\n\n")
  cat(strrep("=", 79), "\n", sep = "")
  cat(title, "\n", sep = "")
  cat(strrep("=", 79), "\n", sep = "")
}

subsection <- function(title) {
  cat("\n", title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
}

print_table <- function(x, n = Inf) {
  print(x, n = n, width = Inf)
}

# =============================================================================
# Model Helpers
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

add_significance <- function(df, p_col = "p_value") {
  df %>%
    mutate(sig = case_when(
      .data[[p_col]] < 0.001 ~ "***",
      .data[[p_col]] < 0.01  ~ "**",
      .data[[p_col]] < 0.05  ~ "*",
      .data[[p_col]] < 0.10  ~ ".",
      TRUE ~ ""
    ))
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
    add_significance() %>%
    select(
      scheme, estimate_log_irr, se, irr, irr_lo, irr_hi,
      pct_change, pct_lo, pct_hi, p_value, sig
    ) %>%
    arrange(pct_change)
}

average_scheme_effect <- function(model, schemes, term_prefix, label) {
  beta <- coef(model)
  V <- vcov(model)
  terms <- paste0(term_prefix, "::", schemes)
  
  missing_terms <- setdiff(terms, names(beta))
  if (length(missing_terms) > 0) {
    stop("Missing model terms: ", paste(missing_terms, collapse = ", "))
  }
  
  weights <- rep(1 / length(terms), length(terms))
  b <- beta[terms]
  V_sub <- V[terms, terms, drop = FALSE]
  
  estimate_log_irr <- as.numeric(sum(weights * b))
  se <- as.numeric(sqrt(t(weights) %*% V_sub %*% weights))
  z <- estimate_log_irr / se
  
  tibble(
    spec = label,
    n_schemes = length(schemes),
    estimate_log_irr = estimate_log_irr,
    se = se,
    z = z,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  ) %>%
    add_irr_columns("estimate_log_irr", "se")
}

run_post_wald <- function(model, var_prefix, post_max = COMMON_POST_MAX) {
  wald(
    model,
    keep = paste0(var_prefix, "::[0-", post_max, "]:treat_group")
  )
}

plot_event_study <- function(df, title, subtitle, colour = "#1f77b4") {
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
      y = "Log incidence rate ratio"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# =============================================================================
# Data Helpers
# =============================================================================

drop_geom_if_needed <- function(x) {
  if (inherits(x, "sf")) sf::st_drop_geometry(x) else x
}

adjust_scheme_start_quarter <- function(road_caz_props) {
  road_caz_props %>%
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
}

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
      event_time_ref_clean = if_else(
        clean_ref_year,
        "ref_year",
        as.character(event_time)
      )
    ) %>%
    ungroup() %>%
    mutate(
      event_time_ref_clean = relevel(
        factor(event_time_ref_clean),
        ref = "ref_year"
      )
    )
}

build_analysis_data <- function() {
  section("1. Build Matched Analysis Panel")
  
  road_panel <- arrow::read_parquet(
    here("data", "processed", "road_panel_matched_pooled.parquet")
  ) %>%
    mutate(
      quarter_year = as.yearqtr(quarter_year)
    )
  
  if (!"LAD24CD" %in% names(road_panel)) {
    oa_lad_lookup <- readRDS(here("data", "processed", "OA_matching_census.rds")) %>%
      drop_geom_if_needed() %>%
      select(OA, LAD24CD, LAD24NM) %>%
      distinct()
    
    road_panel <- road_panel %>%
      left_join(oa_lad_lookup, by = "OA")
  } else if (!"LAD24NM" %in% names(road_panel)) {
    lad_name_lookup <- readRDS(here("data", "processed", "OA_matching_census.rds")) %>%
      drop_geom_if_needed() %>%
      select(LAD24CD, LAD24NM) %>%
      distinct()
    
    road_panel <- road_panel %>%
      left_join(lad_name_lookup, by = "LAD24CD")
  }
  
  scheme_start <- readRDS(here("data", "processed", "roads_caz_props.rds")) %>%
    adjust_scheme_start_quarter() %>%
    rename(caz_start_q_from_props = caz_start_q_adj)
  
  if (!"caz_start_q" %in% names(road_panel)) {
    road_panel <- road_panel %>%
      mutate(caz_start_q = as.yearqtr(NA))
  }
  
  road_panel <- road_panel %>%
    left_join(scheme_start, by = "scheme") %>%
    mutate(
      caz_start_q = coalesce(as.yearqtr(caz_start_q), caz_start_q_from_props)
    ) %>%
    select(-caz_start_q_from_props)
  
  injuries_raw <- readRDS(here("data", "processed", "injuries_final.rds")) %>%
    drop_geom_if_needed()
  
  lad_force_lookup <- injuries_raw %>%
    filter(!is.na(LAD24CD), !is.na(police_force)) %>%
    count(LAD24CD, police_force, name = "n_records") %>%
    group_by(LAD24CD) %>%
    slice_max(n_records, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(LAD24CD, dominant_police_force = police_force)
  
  road_panel <- road_panel %>%
    left_join(lad_force_lookup, by = "LAD24CD") %>%
    mutate(police_force_fe = factor(dominant_police_force))
  
  if (anyNA(road_panel$caz_start_q)) {
    missing_start <- road_panel %>%
      filter(is.na(caz_start_q)) %>%
      distinct(scheme) %>%
      pull(scheme)
    stop("Missing CAZ start quarter for scheme(s): ", paste(missing_start, collapse = ", "))
  }
  
  if (anyNA(road_panel$police_force_fe)) {
    missing_force <- road_panel %>%
      filter(is.na(police_force_fe)) %>%
      distinct(LAD24CD, LAD24NM) %>%
      arrange(LAD24CD)
    print(missing_force, n = Inf)
    stop("Missing police-force lookup for some LADs. See printed LADs above.")
  }
  
  scheme_timing <- road_panel %>%
    filter(treat_group == 1, !is.na(caz_start_q)) %>%
    distinct(scheme, caz_start_q) %>%
    arrange(caz_start_q)
  
  min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)
  
  model_panel_all <- road_panel %>%
    select(
      panel_id, identifier, OA, LAD24CD, LAD24NM,
      dominant_police_force, police_force_fe,
      scheme, quarter_year, caz_start_q, treat_group, all_of(OUTCOME_VAR)
    ) %>%
    rename(outcome_raw = all_of(OUTCOME_VAR)) %>%
    mutate(
      uid = paste0(panel_id, "_", scheme),
      uid_int = as.integer(factor(uid)),
      qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
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
  
  model_panel_for_zero_diag <- model_panel_all
  
  model_panel <- model_panel_all %>%
    filter(unit_total_injury_all_periods > 0) %>%
    select(-unit_total_injury_all_periods)
  
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
        event_time = qtr_int - sc_start_int,
        uid_stack = paste0(uid_int, "_", sc)
      )
  }) %>%
    mutate(
      stack_scheme = factor(stack_scheme),
      treat_scheme = interaction(treat_group, stack_scheme, drop = TRUE)
    )
  
  analysis_weight_counts <- stacked %>%
    distinct(stack_scheme, treat_group, uid_stack) %>%
    count(stack_scheme, treat_group, name = "n_units")
  
  stacked <- stacked %>%
    left_join(analysis_weight_counts, by = c("stack_scheme", "treat_group")) %>%
    mutate(analysis_weight = 1 / n_units) %>%
    select(-n_units) %>%
    make_clean_reference() %>%
    mutate(
      fixed_ref_year = event_time >= -4 & event_time <= -1,
      event_time_ref = if_else(fixed_ref_year, "ref_year", as.character(event_time)),
      event_time_ref = relevel(factor(event_time_ref), ref = "ref_year"),
      period_bucket = case_when(
        fixed_ref_year ~ "ref_year",
        treat_group == 1 &
          event_time >= 0 &
          event_time <= COMMON_POST_MAX ~ "post_common",
        TRUE ~ "other"
      ),
      post_common = as.integer(period_bucket == "post_common"),
      other_flag = as.integer(period_bucket == "other"),
      scheme_post_bucket = if_else(
        post_common == 1,
        as.character(stack_scheme),
        "ref_year"
      ),
      scheme_post_bucket = factor(
        scheme_post_bucket,
        levels = c("ref_year", schemes_all)
      )
    )
  
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
  
  subsection("Scheme timing")
  print_table(scheme_timing)
  
  subsection("Analysis sample")
  print_table(sample_summary)
  
  subsection("Scheme unit counts")
  print_table(scheme_sample_summary)
  
  list(
    stacked = stacked,
    schemes_all = schemes_all,
    scheme_timing = scheme_timing,
    sample_summary = sample_summary,
    scheme_sample_summary = scheme_sample_summary,
    model_panel_for_zero_diag = model_panel_for_zero_diag
  )
}

# =============================================================================
# Primary Models
# =============================================================================

fit_model1_average <- function(data, schemes, label) {
  fit <- feglm(
    outcome_raw ~
      i(scheme_post_bucket, ref = "ref_year") +
      other_flag:treat_scheme |
      uid_stack +
      stack_scheme^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    model = fit,
    average = average_scheme_effect(
      fit,
      schemes = schemes,
      term_prefix = "scheme_post_bucket",
      label = label
    ),
    scheme_effects = extract_scheme_effects(fit, "scheme_post_bucket")
  )
}

fit_model1_bradford_anticipation <- function(data, schemes,
                                             bradford_anticipation_q = "2021 Q1",
                                             post_max = COMMON_POST_MAX) {
  if (!"Bradford" %in% schemes) {
    stop("Bradford must be present in schemes for the anticipation-adjusted model.")
  }
  
  data_ant <- data %>%
    mutate(
      bradford_anticipation_q = as.yearqtr(bradford_anticipation_q),
      bradford_anticipation_period =
        stack_scheme == "Bradford" &
        quarter_year >= bradford_anticipation_q &
        event_time < 0,
      bradford_post_period =
        stack_scheme == "Bradford" &
        event_time >= 0 &
        event_time <= post_max,
      bradford_other_period =
        stack_scheme == "Bradford" &
        !(bradford_anticipation_period | bradford_post_period),
      non_bradford_post_period =
        stack_scheme != "Bradford" &
        post_common == 1,
      non_bradford_other_period =
        stack_scheme != "Bradford" &
        other_flag == 1,
      scheme_post_ant_bucket = case_when(
        treat_group == 1 & bradford_post_period ~ "Bradford",
        treat_group == 1 & non_bradford_post_period ~ as.character(stack_scheme),
        TRUE ~ "ref_year"
      ),
      scheme_post_ant_bucket = factor(
        scheme_post_ant_bucket,
        levels = c("ref_year", schemes)
      ),
      bradford_anticipation_treat =
        as.integer(treat_group == 1 & bradford_anticipation_period),
      bradford_other_treat =
        as.integer(treat_group == 1 & bradford_other_period),
      non_bradford_other_bucket = if_else(
        treat_group == 1 & non_bradford_other_period,
        as.character(stack_scheme),
        "ref_year"
      ),
      non_bradford_other_bucket = factor(
        non_bradford_other_bucket,
        levels = c("ref_year", setdiff(schemes, "Bradford"))
      )
    )
  
  fit <- feglm(
    outcome_raw ~
      i(scheme_post_ant_bucket, ref = "ref_year") +
      bradford_anticipation_treat +
      bradford_other_treat +
      i(non_bradford_other_bucket, ref = "ref_year") |
      uid_stack +
      stack_scheme^qtr_int,
    data = data_ant,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  beta <- coef(fit)
  vc <- vcov(fit)
  terms <- paste0("scheme_post_ant_bucket::", schemes)
  missing_terms <- setdiff(terms, names(beta))
  
  if (length(missing_terms) > 0) {
    stop("Missing anticipation-adjusted post terms: ",
         paste(missing_terms, collapse = ", "))
  }
  
  weights <- rep(1 / length(terms), length(terms))
  estimate_log_irr <- as.numeric(sum(weights * beta[terms]))
  se <- as.numeric(sqrt(t(weights) %*% vc[terms, terms, drop = FALSE] %*% weights))
  z <- estimate_log_irr / se
  
  average <- tibble(
    spec = "Primary sensitivity: Bradford anticipation-adjusted",
    n_schemes = length(schemes),
    estimate_log_irr = estimate_log_irr,
    se = se,
    z = z,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  ) %>%
    add_irr_columns("estimate_log_irr", "se")
  
  bradford_anticipation <- tidy_period_terms(
    fit,
    terms = "bradford_anticipation_treat",
    labels = "Bradford anticipation: 2021 Q1 to 2022 Q3"
  )
  
  list(
    model = fit,
    average = average,
    scheme_effects = extract_scheme_effects(fit, "scheme_post_ant_bucket"),
    bradford_anticipation = bradford_anticipation,
    data = data_ant
  )
}

fit_model2_event_study <- function(data, label, ref_var = "event_time_ref") {
  formula <- as.formula(paste0(
    "outcome_raw ~ i(", ref_var, ", treat_group, ref = 'ref_year') | ",
    "uid_stack + stack_scheme^qtr_int"
  ))
  
  fit <- feglm(
    formula,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = TRUE
  )
  
  list(
    model = fit,
    results = extract_event_study(fit, ref_var) %>%
      mutate(spec = label, .before = 1),
    post_wald = run_post_wald(fit, ref_var)
  )
}

fit_model1_average_police_qtr <- function(data, schemes, label) {
  fit <- feglm(
    outcome_raw ~
      i(scheme_post_bucket, ref = "ref_year") +
      other_flag:treat_scheme |
      uid_stack +
      police_force_fe^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    model = fit,
    average = average_scheme_effect(
      fit,
      schemes = schemes,
      term_prefix = "scheme_post_bucket",
      label = label
    ),
    scheme_effects = extract_scheme_effects(fit, "scheme_post_bucket")
  )
}

fit_model2_event_study_police_qtr <- function(data, label, ref_var = "event_time_ref") {
  formula <- as.formula(paste0(
    "outcome_raw ~ i(", ref_var, ", treat_group, ref = 'ref_year') | ",
    "uid_stack + police_force_fe^qtr_int"
  ))
  
  fit <- feglm(
    formula,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = TRUE
  )
  
  list(
    model = fit,
    results = extract_event_study(fit, ref_var) %>%
      mutate(spec = label, .before = 1),
    post_wald = run_post_wald(fit, ref_var)
  )
}

run_primary_models <- function(stacked, schemes_all) {
  section("2. Primary Models")
  
  subsection("Model 1: equal-scheme average post effect")
  m1 <- fit_model1_average(
    stacked,
    schemes_all,
    "Primary: all schemes"
  )
  print_table(m1$average)
  print_table(m1$scheme_effects)
  
  subsection("Model 1b: Bradford anticipation-adjusted post effect")
  m1_bradford_anticipation <- fit_model1_bradford_anticipation(
    stacked,
    schemes_all,
    bradford_anticipation_q = "2021 Q1",
    post_max = COMMON_POST_MAX
  )
  print_table(m1_bradford_anticipation$average)
  print_table(m1_bradford_anticipation$scheme_effects)
  cat("\nBradford anticipation coefficient:\n")
  print_table(m1_bradford_anticipation$bradford_anticipation)
  
  subsection("Model 2: pooled fixed-reference event study")
  m2 <- fit_model2_event_study(
    stacked,
    "Primary: fixed reference -4:-1",
    ref_var = "event_time_ref"
  )
  print_table(m2$results %>% filter(event_time %in% c(-9, -6, -5, 0:COMMON_POST_MAX)))
  cat("\nJoint Wald test for post-treatment event times 0-", COMMON_POST_MAX, ":\n", sep = "")
  print(m2$post_wald)
  
  subsection("Model 2 sensitivity: clean pre-COVID reference where needed")
  m2_clean <- fit_model2_event_study(
    stacked,
    "Sensitivity: clean reference",
    ref_var = "event_time_ref_clean"
  )
  print_table(m2_clean$results %>% filter(event_time %in% c(-9, -6, -5, 0:COMMON_POST_MAX)))
  cat("\nClean-reference post-treatment Wald test:\n")
  print(m2_clean$post_wald)
  
  subsection("Model 1 sensitivity: police-force by quarter fixed effects")
  m1_police_qtr <- fit_model1_average_police_qtr(
    stacked,
    schemes_all,
    "Sensitivity: police force x quarter FE"
  )
  print_table(m1_police_qtr$average)
  print_table(m1_police_qtr$scheme_effects)
  
  subsection("Model 2 sensitivity: event study with police-force by quarter fixed effects")
  m2_police_qtr <- fit_model2_event_study_police_qtr(
    stacked,
    "Sensitivity: police force x quarter FE",
    ref_var = "event_time_ref"
  )
  print_table(m2_police_qtr$results %>% filter(event_time %in% c(-9, -6, -5, 0:COMMON_POST_MAX)))
  cat("\nPolice-force x quarter FE post-treatment Wald test:\n")
  print(m2_police_qtr$post_wald)
  
  list(
    m1 = m1,
    m1_bradford_anticipation = m1_bradford_anticipation,
    m2 = m2,
    m2_clean = m2_clean,
    m1_police_qtr = m1_police_qtr,
    m2_police_qtr = m2_police_qtr
  )
}

# =============================================================================
# Sensitivity Models
# =============================================================================

run_sensitivities <- function(stacked, schemes_all) {
  section("3. Sensitivity Checks")
  
  no_bradford <- stacked %>%
    filter(stack_scheme != "Bradford") %>%
    droplevels()
  
  bradford_restricted <- stacked %>%
    filter(stack_scheme != "Bradford" | quarter_year >= as.yearqtr("2021 Q2")) %>%
    droplevels()
  
  no_newcastle <- stacked %>%
    filter(stack_scheme != "Newcastle") %>%
    droplevels()
  
  m1_no_bradford <- fit_model1_average(
    no_bradford,
    setdiff(schemes_all, "Bradford"),
    "Sensitivity: exclude Bradford"
  )
  
  m1_bradford_restricted <- fit_model1_average(
    bradford_restricted,
    schemes_all,
    "Sensitivity: Bradford from 2021 Q2 onward"
  )
  
  m1_no_newcastle <- fit_model1_average(
    no_newcastle,
    setdiff(schemes_all, "Newcastle"),
    "Sensitivity: exclude Newcastle"
  )
  
  model1_comparison <- bind_rows(
    fit_model1_average(stacked, schemes_all, "Primary: all schemes")$average,
    m1_no_bradford$average,
    m1_bradford_restricted$average,
    m1_no_newcastle$average
  )
  
  subsection("Model 1 sensitivity comparison")
  print_table(model1_comparison)
  
  m2_no_bradford <- fit_model2_event_study(
    no_bradford,
    "Sensitivity: exclude Bradford",
    ref_var = "event_time_ref"
  )
  
  m2_bradford_restricted <- fit_model2_event_study(
    bradford_restricted,
    "Sensitivity: Bradford from 2021 Q2 onward",
    ref_var = "event_time_ref"
  )
  
  m2_no_newcastle <- fit_model2_event_study(
    no_newcastle,
    "Sensitivity: exclude Newcastle",
    ref_var = "event_time_ref"
  )
  
  model2_comparison <- bind_rows(
    m2_no_bradford$results,
    m2_bradford_restricted$results,
    m2_no_newcastle$results
  ) %>%
    filter(event_time %in% c(-9, -6, 0:COMMON_POST_MAX)) %>%
    select(spec, event_time, estimate, se, pct_change, pct_lo, pct_hi)
  
  subsection("Model 2 sensitivity comparison")
  print_table(model2_comparison)
  
  list(
    model1_comparison = model1_comparison,
    model2_comparison = model2_comparison,
    m1_no_bradford = m1_no_bradford,
    m1_bradford_restricted = m1_bradford_restricted,
    m1_no_newcastle = m1_no_newcastle,
    m2_no_bradford = m2_no_bradford,
    m2_bradford_restricted = m2_bradford_restricted,
    m2_no_newcastle = m2_no_newcastle
  )
}

# =============================================================================
# Scheme-Specific And Pre-Trend Diagnostics
# =============================================================================

fit_scheme_event_study <- function(stacked, sc, ref_var = "event_time_ref_clean") {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  formula <- as.formula(paste0(
    "outcome_raw ~ i(", ref_var, ", treat_group, ref = 'ref_year') | ",
    "uid_stack + qtr_int"
  ))
  
  fit <- tryCatch(
    feglm(
      formula,
      data = d,
      family = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean = TRUE
    ),
    error = function(e) {
      cat("Scheme event study failed for ", sc, ": ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  list(
    model = fit,
    results = extract_event_study(fit, ref_var) %>%
      mutate(scheme = sc, .before = 1)
  )
}

run_scheme_pretrend_test <- function(stacked, sc) {
  d <- stacked %>%
    filter(stack_scheme == sc, event_time < 0) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ event_time:treat_group | uid_stack + qtr_int,
      data = d,
      family = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean = TRUE
    ),
    error = function(e) {
      cat("Pre-trend test failed for ", sc, ": ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  ct <- coeftable(fit)
  tibble(
    scheme = sc,
    n_pre_obs = nrow(d),
    slope = ct["event_time:treat_group", "Estimate"],
    se = ct["event_time:treat_group", "Std. Error"],
    p_value = ct["event_time:treat_group", "Pr(>|z|)"]
  ) %>%
    mutate(pct_slope_per_qtr = 100 * (exp(slope) - 1)) %>%
    add_significance()
}

run_diagnostics <- function(stacked, schemes_all, primary) {
  section("4. Scheme Heterogeneity And Pre-Trend Diagnostics")
  
  scheme_event_fits <- map(schemes_all, ~fit_scheme_event_study(stacked, .x))
  names(scheme_event_fits) <- schemes_all
  
  scheme_event_results <- map_dfr(scheme_event_fits, "results")
  
  scheme_post_wald <- map_dfr(schemes_all, function(sc) {
    fit <- scheme_event_fits[[sc]]$model
    if (is.null(fit)) return(NULL)
    
    w <- tryCatch(
      run_post_wald(fit, "event_time_ref_clean"),
      error = function(e) NULL
    )
    if (is.null(w)) return(NULL)
    
    tibble(
      scheme = sc,
      stat = w$stat,
      df1 = w$df1,
      df2 = w$df2,
      p_value = w$p,
      conclusion = if_else(
        p_value < 0.05,
        "Reject H0: joint post effect",
        "Fail to reject H0"
      )
    )
  })
  
  subsection("Scheme event-study estimates around treatment")
  print_table(
    scheme_event_results %>%
      filter(event_time %in% c(-6, -5, 0:COMMON_POST_MAX)) %>%
      select(scheme, event_time, estimate, se, pct_change, pct_lo, pct_hi)
  )
  
  subsection("Within-scheme post-treatment Wald tests")
  print_table(scheme_post_wald)
  
  pretrend_slopes <- map_dfr(schemes_all, ~run_scheme_pretrend_test(stacked, .x)) %>%
    arrange(p_value)
  
  subsection("Formal pre-treatment slope tests")
  print_table(pretrend_slopes)
  
  bradford_slope <- pretrend_slopes %>%
    filter(scheme == "Bradford") %>%
    pull(slope)
  
  bradford_effect <- primary$m1$scheme_effects %>%
    filter(scheme == "Bradford") %>%
    pull(estimate_log_irr)
  
  if (length(bradford_slope) == 1 && length(bradford_effect) == 1) {
    cat(
      "\nBradford pre-trend drift by event time 2: ",
      round(bradford_slope * 2, 3),
      " log-points; Bradford Model 1 post effect: ",
      round(bradford_effect, 3),
      " log-points.\n",
      sep = ""
    )
  }
  
  suspect_calendar_map <- stacked %>%
    filter(event_time %in% c(-9, -6)) %>%
    distinct(stack_scheme, event_time, quarter_year, covid_period) %>%
    arrange(event_time, stack_scheme)
  
  subsection("Calendar mapping for suspect pre-period event times")
  print_table(suspect_calendar_map)
  
  list(
    scheme_event_fits = scheme_event_fits,
    scheme_event_results = scheme_event_results,
    scheme_post_wald = scheme_post_wald,
    pretrend_slopes = pretrend_slopes,
    suspect_calendar_map = suspect_calendar_map
  )
}

# =============================================================================
# Anticipation And Placebo Timing Checks
# =============================================================================

tidy_period_terms <- function(fit, terms, labels) {
  ct <- coeftable(fit)
  
  map_dfr(seq_along(terms), function(i) {
    term <- terms[i]
    if (!term %in% rownames(ct)) {
      return(tibble(
        period = labels[i],
        estimate_log_irr = NA_real_,
        se = NA_real_,
        p_value = NA_real_
      ))
    }
    
    tibble(
      period = labels[i],
      estimate_log_irr = ct[term, "Estimate"],
      se = ct[term, "Std. Error"],
      p_value = ct[term, "Pr(>|z|)"]
    )
  }) %>%
    add_irr_columns("estimate_log_irr", "se") %>%
    add_significance()
}

fit_anticipation_implementation_model <- function(stacked, sc, announcement_q,
                                                  post_max = COMMON_POST_MAX) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  actual_q <- d %>%
    filter(treat_group == 1, event_time == 0) %>%
    distinct(quarter_year) %>%
    pull(quarter_year) %>%
    first()
  
  q_lookup <- d %>%
    distinct(quarter_year, qtr_int)
  
  announcement_qtr_int <- q_lookup %>%
    filter(quarter_year == as.yearqtr(announcement_q)) %>%
    pull(qtr_int)
  
  implementation_qtr_int <- q_lookup %>%
    filter(quarter_year == actual_q) %>%
    pull(qtr_int)
  
  if (length(announcement_qtr_int) != 1 || length(implementation_qtr_int) != 1) {
    return(list(
      model = NULL,
      results = tibble(
        scheme = sc,
        error = paste0(
          "Announcement or implementation quarter not found: ",
          announcement_q, " / ", actual_q
        )
      )
    ))
  }
  
  d_model <- d %>%
    mutate(
      rel_to_implementation = qtr_int - implementation_qtr_int,
      pre_announcement = qtr_int < announcement_qtr_int,
      anticipation_period =
        qtr_int >= announcement_qtr_int & qtr_int < implementation_qtr_int,
      implementation_period =
        rel_to_implementation >= 0 & rel_to_implementation <= post_max,
      other_period = !(pre_announcement | anticipation_period | implementation_period),
      anticipation_treat = as.integer(anticipation_period & treat_group == 1),
      implementation_treat = as.integer(implementation_period & treat_group == 1),
      other_treat = as.integer(other_period & treat_group == 1)
    ) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ anticipation_treat + implementation_treat + other_treat |
        uid_stack + qtr_int,
      data = d_model,
      family = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      model = NULL,
      results = tibble(scheme = sc, error = conditionMessage(fit))
    ))
  }
  
  results <- tidy_period_terms(
    fit,
    terms = c("anticipation_treat", "implementation_treat"),
    labels = c(
      "Anticipation: announcement to quarter before implementation",
      paste0("Implementation: event times 0-", post_max)
    )
  ) %>%
    mutate(
      scheme = sc,
      announcement_q = as.character(as.yearqtr(announcement_q)),
      implementation_q = as.character(actual_q),
      reference_window = paste0(
        min(d_model$quarter_year[d_model$pre_announcement], na.rm = TRUE),
        " to ",
        max(d_model$quarter_year[d_model$pre_announcement], na.rm = TRUE)
      ),
      anticipation_window = paste0(
        min(d_model$quarter_year[d_model$anticipation_period], na.rm = TRUE),
        " to ",
        max(d_model$quarter_year[d_model$anticipation_period], na.rm = TRUE)
      ),
      implementation_window = paste0(
        min(d_model$quarter_year[d_model$implementation_period], na.rm = TRUE),
        " to ",
        max(d_model$quarter_year[d_model$implementation_period], na.rm = TRUE)
      ),
      error = NA_character_,
      .before = 1
    )
  
  list(model = fit, results = results)
}

average_period_effect <- function(model, schemes, term_prefix, label, period) {
  average_scheme_effect(
    model = model,
    schemes = schemes,
    term_prefix = term_prefix,
    label = label
  ) %>%
    mutate(period = period, .after = spec)
}

fit_pooled_anticipation_implementation_model <- function(stacked, policy_calendar,
                                                         post_max = COMMON_POST_MAX) {
  schemes <- policy_calendar$scheme
  
  implementation_calendar <- stacked %>%
    filter(treat_group == 1, event_time == 0) %>%
    distinct(stack_scheme, quarter_year) %>%
    transmute(
      scheme = as.character(stack_scheme),
      implementation_q = quarter_year
    )
  
  d_model <- stacked %>%
    mutate(scheme = as.character(stack_scheme)) %>%
    left_join(
      policy_calendar %>%
        transmute(
          scheme,
          announcement_q = as.yearqtr(announcement_q),
          anticipation_note = note
        ),
      by = "scheme"
    ) %>%
    left_join(implementation_calendar, by = "scheme") %>%
    filter(!is.na(announcement_q), !is.na(implementation_q)) %>%
    mutate(
      pre_announcement = quarter_year < announcement_q,
      anticipation_period =
        quarter_year >= announcement_q & quarter_year < implementation_q,
      implementation_period =
        event_time >= 0 & event_time <= post_max,
      other_period = !(pre_announcement | anticipation_period | implementation_period),
      scheme_anticipation_bucket = if_else(
        treat_group == 1 & anticipation_period,
        scheme,
        "ref_year"
      ),
      scheme_implementation_bucket = if_else(
        treat_group == 1 & implementation_period,
        scheme,
        "ref_year"
      ),
      scheme_other_bucket = if_else(
        treat_group == 1 & other_period,
        scheme,
        "ref_year"
      ),
      scheme_anticipation_bucket = factor(
        scheme_anticipation_bucket,
        levels = c("ref_year", schemes)
      ),
      scheme_implementation_bucket = factor(
        scheme_implementation_bucket,
        levels = c("ref_year", schemes)
      ),
      scheme_other_bucket = factor(
        scheme_other_bucket,
        levels = c("ref_year", schemes)
      )
    ) %>%
    droplevels()
  
  fit <- feglm(
    outcome_raw ~
      i(scheme_anticipation_bucket, ref = "ref_year") +
      i(scheme_implementation_bucket, ref = "ref_year") +
      i(scheme_other_bucket, ref = "ref_year") |
      uid_stack +
      stack_scheme^qtr_int,
    data = d_model,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  scheme_effects <- bind_rows(
    extract_scheme_effects(fit, "scheme_anticipation_bucket") %>%
      mutate(period = "Anticipation: announcement to quarter before implementation", .after = scheme),
    extract_scheme_effects(fit, "scheme_implementation_bucket") %>%
      mutate(period = paste0("Implementation: event times 0-", post_max), .after = scheme)
  )
  
  pooled_average <- bind_rows(
    average_period_effect(
      fit,
      schemes = schemes,
      term_prefix = "scheme_anticipation_bucket",
      label = "Pooled equal-scheme: anticipation",
      period = "Anticipation: announcement to quarter before implementation"
    ),
    average_period_effect(
      fit,
      schemes = schemes,
      term_prefix = "scheme_implementation_bucket",
      label = "Pooled equal-scheme: implementation",
      period = paste0("Implementation: event times 0-", post_max)
    )
  )
  
  list(
    model = fit,
    data = d_model,
    scheme_effects = scheme_effects,
    pooled_average = pooled_average
  )
}

run_one_scheme_timing <- function(stacked, sc, timing_q, timing_label,
                                  post_max = COMMON_POST_MAX) {
  d <- stacked %>%
    filter(stack_scheme == sc) %>%
    droplevels()
  
  q_lookup <- d %>%
    distinct(quarter_year, qtr_int)
  
  timing_qtr_int <- q_lookup %>%
    filter(quarter_year == as.yearqtr(timing_q)) %>%
    pull(qtr_int)
  
  if (length(timing_qtr_int) != 1) {
    return(tibble(
      scheme = sc,
      timing = timing_label,
      timing_q = timing_q,
      error = paste0("Timing quarter not found or not unique: ", timing_q)
    ))
  }
  
  d_model <- d %>%
    mutate(
      event_time_test = qtr_int - timing_qtr_int,
      ref_year = event_time_test >= -4 & event_time_test <= -1,
      post_common = event_time_test >= 0 & event_time_test <= post_max,
      other_flag = !(ref_year | post_common),
      post_treat = as.integer(post_common & treat_group == 1),
      other_treat = as.integer(other_flag & treat_group == 1)
    ) %>%
    filter(!is.na(event_time_test)) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ post_treat + other_treat | uid_stack + qtr_int,
      data = d_model,
      family = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(tibble(
      scheme = sc,
      timing = timing_label,
      timing_q = timing_q,
      error = conditionMessage(fit)
    ))
  }
  
  ct <- coeftable(fit)
  
  if (!"post_treat" %in% rownames(ct)) {
    return(tibble(
      scheme = sc,
      timing = timing_label,
      timing_q = timing_q,
      error = "post_treat coefficient not estimated"
    ))
  }
  
  est <- ct["post_treat", "Estimate"]
  se <- ct["post_treat", "Std. Error"]
  
  tibble(
    scheme = sc,
    timing = timing_label,
    timing_q = timing_q,
    reference_window = paste0(
      min(d_model$quarter_year[d_model$ref_year]), " to ",
      max(d_model$quarter_year[d_model$ref_year])
    ),
    post_window = paste0(
      min(d_model$quarter_year[d_model$post_common]), " to ",
      max(d_model$quarter_year[d_model$post_common])
    ),
    estimate_log_irr = est,
    se = se,
    p_value = ct["post_treat", "Pr(>|z|)"],
    irr = exp(est),
    pct_change = 100 * (exp(est) - 1),
    pct_lo = 100 * (exp(est - 1.96 * se) - 1),
    pct_hi = 100 * (exp(est + 1.96 * se) - 1),
    error = NA_character_
  )
}

run_all_scheme_placebos <- function(stacked, placebo_q = "2021 Q1",
                                    post_max = COMMON_POST_MAX) {
  schemes <- stacked %>%
    distinct(stack_scheme) %>%
    pull(stack_scheme) %>%
    as.character() %>%
    sort()
  
  actual_timing <- stacked %>%
    filter(treat_group == 1, event_time == 0) %>%
    distinct(stack_scheme, quarter_year) %>%
    transmute(
      scheme = as.character(stack_scheme),
      actual_q = as.character(quarter_year)
    )
  
  actual_results <- map_dfr(schemes, function(sc) {
    actual_q <- actual_timing %>%
      filter(scheme == sc) %>%
      pull(actual_q)
    
    run_one_scheme_timing(
      stacked = stacked,
      sc = sc,
      timing_q = actual_q,
      timing_label = "Actual CAZ date",
      post_max = post_max
    )
  })
  
  placebo_results <- map_dfr(schemes, function(sc) {
    run_one_scheme_timing(
      stacked = stacked,
      sc = sc,
      timing_q = placebo_q,
      timing_label = paste0("Placebo date: ", placebo_q),
      post_max = post_max
    )
  })
  
  combined <- bind_rows(actual_results, placebo_results)
  
  list(
    long = combined %>%
      select(
        scheme, timing, timing_q, reference_window, post_window,
        estimate_log_irr, se, p_value, pct_change, pct_lo, pct_hi, error
      ) %>%
      arrange(scheme, timing),
    wide = combined %>%
      filter(is.na(error)) %>%
      select(scheme, timing, estimate_log_irr, p_value, pct_change, pct_lo, pct_hi) %>%
      pivot_wider(
        names_from = timing,
        values_from = c(estimate_log_irr, p_value, pct_change, pct_lo, pct_hi)
      )
  )
}

make_raw_scheme_quarter_table <- function(stacked, scheme_timing) {
  stacked %>%
    group_by(scheme = stack_scheme, quarter_year, group) %>%
    summarise(
      n_units = n_distinct(uid_stack),
      mean_injuries = mean(outcome_raw, na.rm = TRUE),
      total_injuries = sum(outcome_raw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(scheme_timing, by = "scheme") %>%
    mutate(
      implementation_period = case_when(
        quarter_year < caz_start_q ~ "Pre implementation",
        quarter_year == caz_start_q ~ "Implementation quarter",
        quarter_year > caz_start_q ~ "Post implementation",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(scheme, quarter_year, group)
}

plot_raw_scheme_quarter <- function(raw_by_scheme_quarter) {
  ggplot(
    raw_by_scheme_quarter,
    aes(x = quarter_year, y = mean_injuries, colour = group)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    geom_vline(
      aes(xintercept = as.numeric(caz_start_q)),
      linetype = "dashed"
    ) +
    facet_wrap(~scheme, scales = "free_y", ncol = 2) +
    labs(
      title = "Raw matched mean injuries by scheme and quarter",
      subtitle = "Dashed line marks the adjusted implementation quarter",
      x = "Quarter",
      y = "Mean injuries",
      colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
}

run_anticipation_placebo_checks <- function(stacked, scheme_timing,
                                            post_max = COMMON_POST_MAX,
                                            placebo_q = "2021 Q1") {
  section("5. Anticipation And Placebo Timing Checks")
  
  cat(
    "\nInterpretation guide:\n",
    "  * Anticipation estimates use the announcement-to-implementation window.\n",
    "  * Implementation estimates are the additional effect after formal launch.\n",
    "  * Placebo estimates move the treatment date to a fake/common date.\n",
    "  * If a placebo date is after policy approval, interpret it as an early-policy timing check, not a pure placebo.\n",
    "  * Review the policy calendar dates below; they are explicit assumptions and can be edited in one place.\n",
    sep = ""
  )
  
  policy_calendar <- tibble(
    scheme = c(
      "Bath", "Birmingham", "Bradford", "Bristol",
      "Newcastle", "Portsmouth", "Sheffield"
    ),
    announcement_q = as.yearqtr(c(
      "2020 Q1", # Bath full business plan / planned 2020 launch before COVID delay
      "2019 Q1", # Birmingham charging approval / full business case period
      "2021 Q1", # Bradford council approval / planned rollout became credible
      "2022 Q3", # Bristol final start date announced
      "2022 Q3", # Newcastle/Gateshead final charging timetable becoming operationally salient
      "2021 Q2", # Portsmouth approved/confirmed before late-2021 launch
      "2022 Q3"  # Sheffield final charging timetable and support package period
    )),
    note = c(
      "Bath: full business plan and planned 2020 launch before COVID delay.",
      "Birmingham: charging approval/full business case period; long anticipation window.",
      "Bradford: council approval/planned rollout became credible; formal launch was delayed.",
      "Bristol: final CAZ start date announced before launch.",
      "Newcastle: final charging timetable became operationally salient before staged launch.",
      "Portsmouth: scheme approved/confirmed before late-2021 launch.",
      "Sheffield: final charging timetable/support package period before 2023 launch."
    )
  )
  
  policy_calendar_print <- policy_calendar %>%
    mutate(announcement_q = as.character(announcement_q))
  
  subsection("Policy calendar used for anticipation windows")
  print_table(policy_calendar_print)
  
  anticipation_models <- pmap(
    list(policy_calendar$scheme, policy_calendar$announcement_q),
    ~fit_anticipation_implementation_model(
      stacked = stacked,
      sc = ..1,
      announcement_q = ..2,
      post_max = post_max
    )
  )
  names(anticipation_models) <- policy_calendar$scheme
  
  anticipation_results <- map_dfr(anticipation_models, "results") %>%
    left_join(policy_calendar %>% select(scheme, note), by = "scheme")
  
  pooled_anticipation <- fit_pooled_anticipation_implementation_model(
    stacked = stacked,
    policy_calendar = policy_calendar,
    post_max = post_max
  )
  
  subsection("Pooled equal-scheme anticipation versus implementation effects")
  print_table(
    pooled_anticipation$pooled_average %>%
      select(
        spec, period, n_schemes, estimate_log_irr, se, pct_change,
        pct_lo, pct_hi, p_value
      )
  )
  
  subsection("Per-scheme anticipation versus implementation effects")
  print_table(
    anticipation_results %>%
      select(
        scheme, period, announcement_q, implementation_q,
        estimate_log_irr, se, pct_change, pct_lo, pct_hi, p_value, sig,
        note, error
      )
  )
  
  subsection("Per-scheme effects from pooled anticipation model")
  print_table(
    pooled_anticipation$scheme_effects %>%
      select(
        scheme, period, estimate_log_irr, se, pct_change,
        pct_lo, pct_hi, p_value, sig
      )
  )
  
  placebo_all <- run_all_scheme_placebos(
    stacked = stacked,
    placebo_q = placebo_q,
    post_max = post_max
  )
  
  subsection("Actual implementation dates versus common placebo date")
  print_table(placebo_all$long)
  
  raw_by_scheme_quarter <- make_raw_scheme_quarter_table(
    stacked = stacked,
    scheme_timing = scheme_timing
  )
  
  raw_plot <- plot_raw_scheme_quarter(raw_by_scheme_quarter)
  print(raw_plot)
  
  write_csv(
    anticipation_results,
    file.path(OUTDIR, "anticipation_implementation_by_scheme_separate_models.csv")
  )
  write_csv(
    pooled_anticipation$pooled_average,
    file.path(OUTDIR, "anticipation_implementation_pooled_equal_scheme.csv")
  )
  write_csv(
    pooled_anticipation$scheme_effects,
    file.path(OUTDIR, "anticipation_implementation_by_scheme_pooled_model.csv")
  )
  write_csv(
    policy_calendar_print,
    file.path(OUTDIR, "anticipation_policy_calendar.csv")
  )
  write_csv(
    placebo_all$long,
    file.path(OUTDIR, "placebo_actual_vs_common_date_long.csv")
  )
  write_csv(
    placebo_all$wide,
    file.path(OUTDIR, "placebo_actual_vs_common_date_wide.csv")
  )
  write_csv(
    raw_by_scheme_quarter,
    file.path(OUTDIR, "raw_mean_injuries_by_scheme_quarter.csv")
  )
  ggsave(
    file.path(OUTDIR, "raw_mean_injuries_by_scheme_quarter.png"),
    raw_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  list(
    policy_calendar = policy_calendar,
    anticipation_models = anticipation_models,
    anticipation_results = anticipation_results,
    pooled_anticipation = pooled_anticipation,
    placebo_all = placebo_all,
    raw_by_scheme_quarter = raw_by_scheme_quarter,
    raw_plot = raw_plot
  )
}

# =============================================================================
# Save Outputs And Figures
# =============================================================================

save_outputs <- function(data, primary, sensitivities, diagnostics,
                         anticipation_placebo) {
  section("6. Save Outputs")
  
  p_model1 <- ggplot(
    primary$m1$scheme_effects,
    aes(x = pct_change, y = fct_reorder(scheme, pct_change))
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = pct_lo, xmax = pct_hi), width = 0.2) +
    geom_point(size = 3) +
    labs(
      title = "Model 1: scheme-specific post-CAZ effects",
      subtitle = "Equal-scheme average reported separately in the console",
      x = "% change in injuries",
      y = NULL
    ) +
    theme_minimal(base_size = 12)
  
  p_model2 <- plot_event_study(
    primary$m2$results,
    "Model 2: pooled fixed-reference event study",
    "Reference = event times -4:-1; common post window = event times 0:5"
  )
  
  p_scheme_events <- ggplot(diagnostics$scheme_event_results, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    facet_wrap(~scheme, scales = "free_y", ncol = 2) +
    labs(
      title = "Scheme-specific clean-reference event studies",
      x = "Quarters relative to CAZ implementation",
      y = "Log incidence rate ratio"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  
  ggsave(
    file.path(OUTDIR, "primary_01_scheme_average_effects_clean.png"),
    p_model1,
    width = 9,
    height = 5,
    dpi = 300
  )
  ggsave(
    file.path(OUTDIR, "primary_02_fixedref_event_study_clean.png"),
    p_model2,
    width = 10,
    height = 7,
    dpi = 300
  )
  ggsave(
    file.path(OUTDIR, "diagnostic_scheme_cleanref_event_studies_clean.png"),
    p_scheme_events,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  results <- list(
    data = list(
      scheme_timing = data$scheme_timing,
      sample_summary = data$sample_summary,
      scheme_sample_summary = data$scheme_sample_summary
    ),
    primary = primary,
    sensitivities = sensitivities,
    diagnostics = diagnostics,
    anticipation_placebo = anticipation_placebo
  )
  
  saveRDS(results, here("data", "processed", "caz_primary_ppml_clean_results.rds"))
  
  support_objects <- list(
    stacked = data$stacked,
    schemes_all = data$schemes_all,
    model_panel_for_zero_diag = data$model_panel_for_zero_diag
  )
  saveRDS(support_objects, here("data", "processed", "caz_primary_ppml_clean_support.rds"))
  
  cat("Saved figures to: ", OUTDIR, "\n", sep = "")
  cat("Saved results: data/processed/caz_primary_ppml_clean_results.rds\n")
  cat("Saved support data: data/processed/caz_primary_ppml_clean_support.rds\n")
  
  invisible(results)
}

# =============================================================================
# Run Analysis
# =============================================================================

section("CAZ Injury DiD - Clean PPML Analysis")
cat("Outcome: ", OUTCOME_VAR, "\n", sep = "")
cat("Common post-treatment window: event times 0-", COMMON_POST_MAX, "\n", sep = "")

analysis_data <- build_analysis_data()

primary <- run_primary_models(
  stacked = analysis_data$stacked,
  schemes_all = analysis_data$schemes_all
)

sensitivities <- run_sensitivities(
  stacked = analysis_data$stacked,
  schemes_all = analysis_data$schemes_all
)

diagnostics <- run_diagnostics(
  stacked = analysis_data$stacked,
  schemes_all = analysis_data$schemes_all,
  primary = primary
)

anticipation_placebo <- run_anticipation_placebo_checks(
  stacked = analysis_data$stacked,
  scheme_timing = analysis_data$scheme_timing,
  post_max = COMMON_POST_MAX,
  placebo_q = "2021 Q1"
)

all_results <- save_outputs(
  data = analysis_data,
  primary = primary,
  sensitivities = sensitivities,
  diagnostics = diagnostics,
  anticipation_placebo = anticipation_placebo
)

section("Console Story")

cat("\nHeadline pooled average effect:\n")
print_table(primary$m1$average)

cat("\nScheme effects, ordered by percent change:\n")
print_table(primary$m1$scheme_effects)

cat("\nPooled event-study post-treatment quarters:\n")
print_table(
  primary$m2$results %>%
    filter(event_time %in% 0:COMMON_POST_MAX) %>%
    select(event_time, estimate, se, pct_change, pct_lo, pct_hi)
)

cat("\nPolice-force x quarter FE sensitivity: pooled average effect:\n")
print_table(primary$m1_police_qtr$average)

cat("\nPolice-force x quarter FE sensitivity: scheme effects:\n")
print_table(primary$m1_police_qtr$scheme_effects)

cat("\nPolice-force x quarter FE sensitivity: event-study post-treatment quarters:\n")
print_table(
  primary$m2_police_qtr$results %>%
    filter(event_time %in% 0:COMMON_POST_MAX) %>%
    select(event_time, estimate, se, pct_change, pct_lo, pct_hi)
)

cat("\nSensitivity summary:\n")
print_table(sensitivities$model1_comparison)

cat("\nPre-trend slope tests:\n")
print_table(diagnostics$pretrend_slopes)

cat("\nAnticipation versus implementation timing check:\n")
print_table(
  anticipation_placebo$pooled_anticipation$pooled_average %>%
    select(
      spec, period, n_schemes, estimate_log_irr, se,
      pct_change, pct_lo, pct_hi, p_value
    )
)

cat("\nPer-scheme anticipation versus implementation timing check:\n")
print_table(
  anticipation_placebo$anticipation_results %>%
    select(
      scheme, period, announcement_q, implementation_q,
      estimate_log_irr, se, pct_change, pct_lo, pct_hi, p_value, sig
    )
)

cat("\nActual versus placebo timing checks:\n")
print_table(anticipation_placebo$placebo_all$long)







# =============================================================================
# BRADFORD TIMING TABLE: actual date, placebo date, and actual date with
# pre-2021 data excluded from the reference (pre-intervention) window
# =============================================================================

fit_bradford_timing <- function(d, timing_label, restrict_from = NULL) {
  # d must already have event_time_test, computed relative to whichever
  # implementation date this run is testing.
  d_model <- d
  
  if (!is.null(restrict_from)) {
    d_model <- d_model %>% filter(quarter_year >= as.yearqtr(restrict_from))
  }
  
  d_model <- d_model %>%
    mutate(
      ref_year    = event_time_test >= -4 & event_time_test <= -1,
      post_common = event_time_test >= 0 & event_time_test <= 5,
      other_flag  = !(ref_year | post_common),
      post_treat  = as.integer(post_common & treat_group == 1),
      other_treat = as.integer(other_flag & treat_group == 1)
    ) %>%
    filter(!is.na(event_time_test)) %>%
    droplevels()
  
  fit <- feglm(
    outcome_raw ~ post_treat + other_treat | uid_stack + qtr_int,
    data = d_model, family = "poisson", cluster = ~OA,
    weights = ~analysis_weight, lean = FALSE
  )
  
  ct <- coeftable(fit)
  est <- ct["post_treat", "Estimate"]
  se  <- ct["post_treat", "Std. Error"]
  
  tibble(
    timing = timing_label,
    reference_window = paste0(
      min(d_model$quarter_year[d_model$ref_year]), " to ",
      max(d_model$quarter_year[d_model$ref_year])
    ),
    post_window = paste0(
      min(d_model$quarter_year[d_model$post_common]), " to ",
      max(d_model$quarter_year[d_model$post_common])
    ),
    n_obs             = nrow(d_model),
    estimate_log_irr  = est,
    se                = se,
    p_value           = ct["post_treat", "Pr(>|z|)"],
    pct_change        = 100 * (exp(est) - 1),
    pct_lo            = 100 * (exp(est - 1.96 * se) - 1),
    pct_hi            = 100 * (exp(est + 1.96 * se) - 1)
  )
}

d_bradford <- analysis_data$stacked %>%
  filter(stack_scheme == "Bradford") %>%
  droplevels()

# Look up integer quarters for the actual and placebo implementation dates
q_lookup <- d_bradford %>% distinct(quarter_year, qtr_int)

actual_qtr_int  <- q_lookup %>% filter(quarter_year == as.yearqtr("2022 Q4")) %>% pull(qtr_int)
placebo_qtr_int <- q_lookup %>% filter(quarter_year == as.yearqtr("2021 Q1")) %>% pull(qtr_int)

d_actual  <- d_bradford %>% mutate(event_time_test = qtr_int - actual_qtr_int)
d_placebo <- d_bradford %>% mutate(event_time_test = qtr_int - placebo_qtr_int)

bradford_timing_table <- bind_rows(
  fit_bradford_timing(d_actual,  "Actual CAZ date: 2022 Q4"),
  fit_bradford_timing(d_placebo, "Placebo/early date: 2021 Q1"),
  # NEW: actual implementation date, but reference window restricted to
  # exclude any quarter before 2021 Q1 -- tests whether the "early date"
  # signal depends on pre-2021 (COVID-era) data being in the comparison.
  fit_bradford_timing(d_actual,  "Actual CAZ date, pre-2021 data excluded",
                      restrict_from = "2021 Q1")
)

print(bradford_timing_table %>%
        select(timing, reference_window, post_window, n_obs,
               estimate_log_irr, pct_change, p_value),
      n = Inf)




# =============================================================================
# BRADFORD-ONLY FIXED-REFERENCE EVENT STUDY PLOT
# =============================================================================

d_bradford <- analysis_data$stacked %>%
  filter(stack_scheme == "Bradford") %>%
  droplevels()

m2_bradford <- feglm(
  outcome_raw ~
    i(event_time_ref, treat_group, ref = "ref_year") |
    uid_stack + qtr_int,
  data    = d_bradford,
  family  = "poisson",
  cluster = ~OA,
  weights = ~analysis_weight,
  lean    = TRUE
)

bradford_event_results <- extract_event_study(m2_bradford, "event_time_ref")

p_bradford_event <- plot_event_study(
  bradford_event_results,
  title    = "Bradford: fixed-reference event study",
  subtitle = "Reference = event times -4:-1; common post window = event times 0:5",
  colour   = "#c0392b"   # distinct colour so it isn't confused with the pooled plot
)

print(p_bradford_event)

ggsave(
  file.path(OUTDIR, "bradford_only_fixedref_event_study.png"),
  p_bradford_event, width = 10, height = 7, dpi = 300
)


run_scheme_pretrend_test_windowed <- function(stacked, sc, lookback = 8) {
  d <- stacked %>%
    filter(stack_scheme == sc, event_time < 0, event_time >= -lookback) %>%
    droplevels()
  
  fit <- feglm(
    outcome_raw ~ event_time:treat_group | uid_stack + qtr_int,
    data = d, family = "poisson", cluster = ~OA,
    weights = ~analysis_weight, lean = TRUE
  )
  
  ct <- coeftable(fit)
  tibble(
    scheme = sc, lookback_qtrs = lookback, n_pre_obs = nrow(d),
    slope = ct["event_time:treat_group", "Estimate"],
    se    = ct["event_time:treat_group", "Std. Error"],
    p_value = ct["event_time:treat_group", "Pr(>|z|)"]
  )
}

bradford_recent_pretrend <- run_scheme_pretrend_test_windowed(
  analysis_data$stacked, "Bradford", lookback = 8
)
print(bradford_recent_pretrend)



d_bradford_bucket <- analysis_data$stacked %>%
  filter(stack_scheme == "Bradford", event_time >= -8, event_time <= -1) %>%
  mutate(
    early_pre = as.integer(event_time >= -8 & event_time <= -5),
    ref_bucket = as.integer(event_time >= -4 & event_time <= -1)
  ) %>%
  droplevels()

m_bradford_bucket <- feglm(
  outcome_raw ~ early_pre:treat_group | uid_stack + qtr_int,
  data = d_bradford_bucket, family = "poisson", cluster = ~OA,
  weights = ~analysis_weight, lean = TRUE
)

print(coeftable(m_bradford_bucket))



# =============================================================================
# BRADFORD: RAW INJURY TRENDS, TREATED VS CONTROL, BY QUARTER
# =============================================================================

bradford_raw_trend <- analysis_data$stacked %>%
  filter(stack_scheme == "Bradford") %>%
  group_by(quarter_year, group) %>%
  summarise(
    n_units        = n_distinct(uid_stack),
    mean_injuries  = mean(outcome_raw, na.rm = TRUE),
    total_injuries = sum(outcome_raw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(quarter_year)

bradford_caz_start <- analysis_data$scheme_timing %>%
  filter(scheme == "Bradford") %>%
  pull(caz_start_q)

p_bradford_raw <- ggplot(bradford_raw_trend, aes(x = quarter_year, y = mean_injuries, colour = group)) +
  geom_vline(xintercept = as.numeric(bradford_caz_start), linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("CAZ roads" = "#e74c3c", "Matched controls" = "#2980b9")) +
  labs(
    title = "Bradford: mean injuries per road-link-quarter",
    subtitle = "CAZ roads vs. matched controls; dashed line = CAZ implementation (2022 Q4)",
    x = "Quarter",
    y = "Mean injuries per road-link-quarter",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

print(p_bradford_raw)

ggsave(
  file.path(OUTDIR, "bradford_raw_injury_trend.png"),
  p_bradford_raw, width = 10, height = 6, dpi = 300
)

p_bradford_raw_annotated <- p_bradford_raw +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey60") +
  annotate("text", x = as.numeric(as.yearqtr("2021 Q1")), y = Inf,
           label = "Announcement", angle = 90, vjust = 1.3, hjust = 1.1, size = 3, colour = "grey40") +
  annotate("text", x = as.numeric(bradford_caz_start), y = Inf,
           label = "Implementation", angle = 90, vjust = 1.3, hjust = -0.2, size = 3, colour = "grey40")

print(p_bradford_raw_annotated)
