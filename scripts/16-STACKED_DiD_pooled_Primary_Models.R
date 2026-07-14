# =============================================================================
# CAZ INJURY DIFFERENCE-IN-DIFFERENCES: PRIMARY PPML ANALYSIS
# =============================================================================
#
# Purpose
#   Estimate the average effect of seven Clean Air Zone (CAZ) schemes on road
#   injuries using the complete matched road-quarter panel.
#
# Primary estimand
#   The equal-scheme average effect during event times 0:5, relative to the
#   four quarters immediately before implementation (event times -4:-1).
#
# Primary model structure
#   - Complete panel retained.
#   - Road-by-scheme fixed effects.
#   - Scheme-by-calendar-quarter fixed effects.
#   - Treated roads and matched controls each sum to weight one per scheme.
#   - Standard errors clustered by Output Area (OA).
#
# Why police-force-by-quarter FE are not pooled across schemes
#   Bath, Bristol, Newcastle and Sheffield have no same-force control OAs;
#   Portsmouth has only three and Birmingham has only 22. A pooled police-force
#   model would therefore rely on unsupported cross-force extrapolation and can
#   become collinear. Bradford has 132 same-force control OAs, so police-force
#   adjustment is used only as a Bradford robustness analysis.
#
# Output
#   All tables and tests print to the R console. No files are saved.
#
# Run this script from the project root containing data/processed/.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(fixest)
  library(here)
  library(lubridate)
  library(sf)
  library(tidyverse)
  library(zoo)
})

# Avoid namespace conflicts in interactive sessions.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
count  <- dplyr::count

# =============================================================================
# 1. Settings
# =============================================================================

OUTCOME_VAR <- "total_inj_adj_All"

REFERENCE_MIN <- -4L
REFERENCE_MAX <- -1L
POST_MAX      <- 5L

BRADFORD_ANNOUNCEMENT_Q <- as.yearqtr("2021 Q1")

# Set to FALSE if the script is run non-interactively and plots are not wanted.
SHOW_PLOTS <- TRUE

# =============================================================================
# 2. Console and result helpers
# =============================================================================

section <- function(title) {
  cat(
    "\n\n", strrep("=", 79), "\n",
    title, "\n",
    strrep("=", 79), "\n",
    sep = ""
  )
}

subsection <- function(title) {
  cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
}

print_table <- function(x) {
  print(x, n = Inf, width = Inf)
}

drop_geometry <- function(x) {
  if (inherits(x, "sf")) sf::st_drop_geometry(x) else x
}

add_effect_scale <- function(data, estimate_col = "estimate", se_col = "se") {
  data %>%
    mutate(
      ci_lo = .data[[estimate_col]] - 1.96 * .data[[se_col]],
      ci_hi = .data[[estimate_col]] + 1.96 * .data[[se_col]],
      irr = exp(.data[[estimate_col]]),
      irr_lo = exp(ci_lo),
      irr_hi = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo = 100 * (irr_lo - 1),
      pct_hi = 100 * (irr_hi - 1)
    )
}

add_significance <- function(data, p_col = "p_value") {
  data %>%
    mutate(
      sig = case_when(
        .data[[p_col]] < 0.001 ~ "***",
        .data[[p_col]] < 0.01  ~ "**",
        .data[[p_col]] < 0.05  ~ "*",
        .data[[p_col]] < 0.10  ~ ".",
        TRUE ~ ""
      )
    )
}

tidy_coefficient <- function(model, term, label = term) {
  ct <- coeftable(model)
  available_terms <- rownames(ct)
  
  # fixest normally keeps a numeric regressor's original name, but some
  # versions append a level suffix after internal coercion. Prefer the exact
  # name; otherwise accept one unique coefficient containing the requested
  # stem. If nothing matches, report the actual model terms for diagnosis.
  resolved_term <- term
  if (!resolved_term %in% available_terms) {
    candidates <- available_terms[str_detect(available_terms, fixed(term))]
    
    if (length(candidates) == 1) {
      resolved_term <- candidates
    } else {
      stop(
        "Coefficient '", term, "' was not available. Model coefficients: ",
        paste(available_terms, collapse = ", "),
        if (length(model$collin.var) > 0) {
          paste0(
            ". Collinear variables removed: ",
            paste(model$collin.var, collapse = ", ")
          )
        } else {
          ""
        }
      )
    }
  }
  
  tibble(
    term = label,
    model_term = resolved_term,
    estimate = ct[resolved_term, "Estimate"],
    se = ct[resolved_term, "Std. Error"],
    p_value = ct[resolved_term, "Pr(>|z|)"]
  ) %>%
    add_effect_scale() %>%
    add_significance()
}

extract_scheme_effects <- function(model, prefix) {
  ct <- coeftable(model)
  
  tibble(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    p_value = ct[, "Pr(>|z|)"]
  ) %>%
    filter(str_detect(term, paste0("^", prefix, "::"))) %>%
    mutate(scheme = str_remove(term, paste0("^", prefix, "::"))) %>%
    add_effect_scale() %>%
    add_significance() %>%
    select(
      scheme, estimate_log_irr = estimate, se,
      irr, irr_lo, irr_hi,
      pct_change, pct_lo, pct_hi,
      p_value, sig
    ) %>%
    arrange(pct_change)
}

extract_event_study <- function(model, prefix = "event_time_ref") {
  ct <- coeftable(model)
  
  tibble(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    p_value = ct[, "Pr(>|z|)"]
  ) %>%
    filter(str_detect(term, paste0("^", prefix, "::"))) %>%
    mutate(
      event_time = str_match(
        term,
        paste0("^", prefix, "::(-?\\d+):treat_group$")
      )[, 2] %>% as.integer()
    ) %>%
    filter(!is.na(event_time)) %>%
    add_effect_scale() %>%
    arrange(event_time)
}

average_scheme_effect <- function(model, schemes, prefix, label) {
  beta <- coef(model)
  variance <- vcov(model)
  terms <- paste0(prefix, "::", schemes)
  
  missing_terms <- setdiff(terms, names(beta))
  if (length(missing_terms) > 0) {
    stop(
      "The primary model did not identify: ",
      paste(missing_terms, collapse = ", ")
    )
  }
  
  weights <- rep(1 / length(terms), length(terms))
  estimate <- as.numeric(sum(weights * beta[terms]))
  se <- as.numeric(
    sqrt(t(weights) %*% variance[terms, terms, drop = FALSE] %*% weights)
  )
  z <- estimate / se
  
  tibble(
    specification = label,
    n_schemes = length(schemes),
    estimate = estimate,
    se = se,
    z = z,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  ) %>%
    add_effect_scale() %>%
    add_significance()
}

joint_event_test <- function(model, period = c("recent_pre", "post"),
                             prefix = "event_time_ref") {
  period <- match.arg(period)
  
  pattern <- if (period == "recent_pre") {
    paste0(prefix, "::-[5-8]:treat_group")
  } else {
    paste0(prefix, "::[0-5]:treat_group")
  }
  
  wald(model, keep = pattern)
}

# =============================================================================
# 3. Build the complete matched road-quarter panel
# =============================================================================

adjust_scheme_start_quarter <- function(road_caz_properties) {
  road_caz_properties %>%
    distinct(scheme, startDt, caz_start_q) %>%
    filter(!is.na(startDt)) %>%
    mutate(
      start_date = dmy(startDt),
      raw_quarter = as.yearqtr(start_date),
      quarter_start = as.Date(raw_quarter),
      quarter_end = as.Date(raw_quarter + 0.25) - 1,
      quarter_midpoint = quarter_start +
        as.integer(difftime(quarter_end, quarter_start, units = "days")) / 2,
      adjusted_start_quarter = if_else(
        start_date > quarter_midpoint,
        raw_quarter + 0.25,
        raw_quarter
      )
    ) %>%
    select(scheme, adjusted_start_quarter)
}

make_force_overlap_table <- function(stacked) {
  treated_forces <- stacked %>%
    filter(treat_group == 1) %>%
    distinct(stack_scheme, police_force_fe)
  
  control_oas <- stacked %>%
    filter(treat_group == 0) %>%
    distinct(stack_scheme, OA, police_force_fe)
  
  totals <- control_oas %>%
    count(stack_scheme, name = "total_control_oas")
  
  same_force <- control_oas %>%
    inner_join(treated_forces, by = c("stack_scheme", "police_force_fe")) %>%
    count(stack_scheme, name = "same_force_control_oas")
  
  totals %>%
    left_join(same_force, by = "stack_scheme") %>%
    mutate(
      same_force_control_oas = replace_na(same_force_control_oas, 0L),
      pct_same_force = 100 * same_force_control_oas / total_control_oas,
      police_fe_assessment = case_when(
        same_force_control_oas == 0 ~ "No within-force support",
        same_force_control_oas < 10 ~ "Essentially no support",
        same_force_control_oas < 50 ~ "Thin support",
        TRUE ~ "Meaningful support"
      )
    ) %>%
    arrange(stack_scheme)
}

build_analysis_data <- function() {
  section("1. Build the matched analysis panel")
  
  road_panel <- arrow::read_parquet(
    here("data", "processed", "road_panel_matched_pooled.parquet")
  ) %>%
    mutate(quarter_year = as.yearqtr(quarter_year))
  
  oa_lookup <- readRDS(
    here("data", "processed", "OA_matching_census.rds")
  ) %>%
    drop_geometry() %>%
    select(OA, LAD24CD, LAD24NM) %>%
    distinct()
  
  if (!"LAD24CD" %in% names(road_panel)) {
    road_panel <- road_panel %>% left_join(oa_lookup, by = "OA")
  } else if (!"LAD24NM" %in% names(road_panel)) {
    road_panel <- road_panel %>%
      left_join(
        oa_lookup %>% distinct(LAD24CD, LAD24NM),
        by = "LAD24CD"
      )
  }
  
  scheme_start <- readRDS(
    here("data", "processed", "roads_caz_props.rds")
  ) %>%
    adjust_scheme_start_quarter() %>%
    rename(start_from_properties = adjusted_start_quarter)
  
  if (!"caz_start_q" %in% names(road_panel)) {
    road_panel <- road_panel %>% mutate(caz_start_q = as.yearqtr(NA))
  }
  
  road_panel <- road_panel %>%
    left_join(scheme_start, by = "scheme") %>%
    mutate(
      caz_start_q = coalesce(as.yearqtr(caz_start_q), start_from_properties)
    ) %>%
    select(-start_from_properties)
  
  injuries <- readRDS(
    here("data", "processed", "injuries_final.rds")
  ) %>%
    drop_geometry()
  
  # The dominant force is used because the injury data identify reporting force
  # at record level while the road panel is indexed by LAD and OA.
  lad_force_lookup <- injuries %>%
    filter(!is.na(LAD24CD), !is.na(police_force)) %>%
    count(LAD24CD, police_force, name = "n_records") %>%
    group_by(LAD24CD) %>%
    slice_max(n_records, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      LAD24CD,
      dominant_police_force = police_force,
      police_force_fe = factor(police_force)
    )
  
  road_panel <- road_panel %>%
    left_join(lad_force_lookup, by = "LAD24CD")
  
  if (anyNA(road_panel$caz_start_q)) {
    missing_schemes <- road_panel %>%
      filter(is.na(caz_start_q)) %>%
      distinct(scheme) %>%
      pull(scheme)
    stop("Missing implementation date for: ", paste(missing_schemes, collapse = ", "))
  }
  
  if (anyNA(road_panel$police_force_fe)) {
    missing_lads <- road_panel %>%
      filter(is.na(police_force_fe)) %>%
      distinct(LAD24CD, LAD24NM)
    print_table(missing_lads)
    stop("Police-force lookup is incomplete.")
  }
  
  scheme_timing <- road_panel %>%
    filter(treat_group == 1) %>%
    distinct(scheme, caz_start_q) %>%
    arrange(caz_start_q)
  
  if (any(count(scheme_timing, scheme)$n != 1)) {
    stop("Each scheme must have exactly one implementation quarter.")
  }
  
  first_quarter <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)
  
  panel_all <- road_panel %>%
    select(
      panel_id, identifier, OA, LAD24CD, LAD24NM,
      dominant_police_force, police_force_fe,
      scheme, quarter_year, caz_start_q, treat_group,
      all_of(OUTCOME_VAR)
    ) %>%
    rename(outcome_raw = all_of(OUTCOME_VAR)) %>%
    filter(!is.na(outcome_raw)) %>%
    mutate(
      uid_stack = paste(panel_id, scheme, sep = "__"),
      qtr_int = as.integer(
        round((as.numeric(quarter_year) - first_quarter) * 4)
      ) + 1L,
      event_time = as.integer(
        round((as.numeric(quarter_year) - as.numeric(caz_start_q)) * 4)
      ),
      group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
    ) %>%
    group_by(uid_stack) %>%
    mutate(unit_outcome_total = sum(outcome_raw, na.rm = TRUE)) %>%
    ungroup()
  
  # All-zero roads do not contribute to a conditional road-FE Poisson model.
  # Removing them before weighting matches what fixest would do automatically.
  panel <- panel_all %>%
    filter(unit_outcome_total > 0) %>%
    select(-unit_outcome_total)
  
  schemes <- sort(unique(panel$scheme))
  
  stacked <- panel %>%
    mutate(
      stack_scheme = factor(scheme, levels = schemes),
      fixed_reference = event_time >= REFERENCE_MIN &
        event_time <= REFERENCE_MAX,
      event_time_ref = if_else(
        fixed_reference,
        "ref_year",
        as.character(event_time)
      ),
      event_time_ref = relevel(factor(event_time_ref), ref = "ref_year")
    )
  
  # Equal-scheme weighting: within every scheme, all treated roads together
  # receive weight one and all matched controls together receive weight one.
  unit_counts <- stacked %>%
    distinct(stack_scheme, treat_group, uid_stack) %>%
    count(stack_scheme, treat_group, name = "n_units")
  
  stacked <- stacked %>%
    left_join(unit_counts, by = c("stack_scheme", "treat_group")) %>%
    mutate(analysis_weight = 1 / n_units) %>%
    select(-n_units)
  
  sample_summary <- stacked %>%
    distinct(stack_scheme, treat_group, uid_stack, .keep_all = TRUE) %>%
    count(stack_scheme, treat_group, name = "road_units") %>%
    arrange(stack_scheme, treat_group)
  
  # Assess police-force overlap in the original matched design, before PPML
  # removes all-zero road units. This keeps the diagnostic about the quality of
  # the matching/control pool rather than about outcome-driven model retention.
  force_overlap <- panel_all %>%
    mutate(stack_scheme = factor(scheme, levels = schemes)) %>%
    make_force_overlap_table()
  
  subsection("Scheme implementation timing")
  print_table(scheme_timing)
  
  subsection("Road units by scheme and treatment group")
  print_table(sample_summary)
  
  subsection("Police-force overlap diagnostic")
  print_table(force_overlap)
  
  list(
    stacked = stacked,
    schemes = schemes,
    scheme_timing = scheme_timing,
    sample_summary = sample_summary,
    force_overlap = force_overlap
  )
}

# =============================================================================
# 4. Full-panel static model
# =============================================================================

prepare_static_full_panel <- function(stacked, schemes) {
  stacked %>%
    filter(stack_scheme %in% schemes) %>%
    mutate(
      static_period = case_when(
        event_time < REFERENCE_MIN ~ "early_pre",
        fixed_reference ~ "reference",
        event_time >= 0 & event_time <= POST_MAX ~ "post_common",
        event_time > POST_MAX ~ "late_post",
        TRUE ~ NA_character_
      ),
      
      # The three scheme factors are defined only for treated roads. Controls
      # remain in ref_year. Omitting the reference-period treated interaction
      # makes every coefficient relative to event times -4:-1.
      scheme_early_pre = if_else(
        treat_group == 1 & static_period == "early_pre",
        as.character(stack_scheme),
        "ref_year"
      ),
      scheme_post_common = if_else(
        treat_group == 1 & static_period == "post_common",
        as.character(stack_scheme),
        "ref_year"
      ),
      scheme_late_post = if_else(
        treat_group == 1 & static_period == "late_post",
        as.character(stack_scheme),
        "ref_year"
      ),
      
      scheme_early_pre = factor(
        scheme_early_pre,
        levels = c("ref_year", schemes)
      ),
      scheme_post_common = factor(
        scheme_post_common,
        levels = c("ref_year", schemes)
      ),
      scheme_late_post = factor(
        scheme_late_post,
        levels = c("ref_year", schemes)
      )
    ) %>%
    droplevels()
}

fit_static_full_panel <- function(stacked, schemes, label) {
  data <- prepare_static_full_panel(stacked, schemes)
  
  model <- feglm(
    outcome_raw ~
      i(scheme_early_pre, ref = "ref_year") +
      i(scheme_post_common, ref = "ref_year") +
      i(scheme_late_post, ref = "ref_year") |
      uid_stack + stack_scheme^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    model = model,
    average = average_scheme_effect(
      model,
      schemes = schemes,
      prefix = "scheme_post_common",
      label = label
    ),
    scheme_effects = extract_scheme_effects(model, "scheme_post_common")
  )
}

# =============================================================================
# 5. Full-panel pooled event study
# =============================================================================

fit_event_study_full_panel <- function(stacked, schemes, label) {
  data <- stacked %>%
    filter(stack_scheme %in% schemes) %>%
    droplevels()
  
  model <- feglm(
    outcome_raw ~
      i(event_time_ref, treat_group, ref = "ref_year") |
      uid_stack + stack_scheme^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    model = model,
    results = extract_event_study(model, "event_time_ref") %>%
      mutate(specification = label, .before = 1),
    recent_pre_wald = joint_event_test(
      model,
      "recent_pre",
      "event_time_ref"
    ),
    post_wald = joint_event_test(model, "post", "event_time_ref")
  )
}

# =============================================================================
# 6. Bradford robustness analyses
# =============================================================================

prepare_bradford_periods <- function(data) {
  data %>%
    mutate(
      early_pre_treat = as.integer(
        treat_group == 1 & event_time < REFERENCE_MIN
      ),
      post_common_treat = as.integer(
        treat_group == 1 & event_time >= 0 & event_time <= POST_MAX
      ),
      late_post_treat = as.integer(
        treat_group == 1 & event_time > POST_MAX
      )
    )
}

fit_bradford_police_fe <- function(stacked) {
  # This uses Bradford's complete matched-control pool. Police-force-by-quarter
  # FE absorb reporting and calendar shocks separately for every represented
  # force. Bradford has enough controls in its own force to identify the treated
  # post coefficient. Other schemes do not, so this model is not pooled.
  data <- stacked %>%
    filter(stack_scheme == "Bradford") %>%
    prepare_bradford_periods() %>%
    droplevels()
  
  static_model <- feglm(
    outcome_raw ~
      early_pre_treat + post_common_treat + late_post_treat |
      uid_stack + police_force_fe^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  event_model <- feglm(
    outcome_raw ~
      i(event_time_ref, treat_group, ref = "ref_year") |
      uid_stack + police_force_fe^qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    static_model = static_model,
    static_result = tidy_coefficient(
      static_model,
      "post_common_treat",
      "Bradford: police-force x quarter FE"
    ),
    event_model = event_model,
    event_results = extract_event_study(event_model, "event_time_ref"),
    recent_pre_wald = joint_event_test(
      event_model,
      "recent_pre",
      "event_time_ref"
    ),
    post_wald = joint_event_test(event_model, "post", "event_time_ref")
  )
}

make_bradford_same_force_sample <- function(stacked) {
  bradford_forces <- stacked %>%
    filter(stack_scheme == "Bradford", treat_group == 1) %>%
    distinct(police_force_fe)
  
  bradford <- stacked %>%
    filter(stack_scheme == "Bradford") %>%
    semi_join(bradford_forces, by = "police_force_fe") %>%
    droplevels()
  
  # Restricting the control pool changes its size, so weights are recalculated
  # once for the complete same-force panel. They are not window-specific.
  unit_counts <- bradford %>%
    distinct(treat_group, uid_stack) %>%
    count(treat_group, name = "n_units_same_force")
  
  bradford %>%
    select(-analysis_weight) %>%
    left_join(unit_counts, by = "treat_group") %>%
    mutate(analysis_weight = 1 / n_units_same_force) %>%
    select(-n_units_same_force)
}

fit_bradford_same_force <- function(bradford) {
  data <- bradford %>%
    prepare_bradford_periods()
  
  # Because every road is now in the same police force, ordinary calendar-
  # quarter FE are equivalent to police-force-by-quarter FE.
  static_model <- feglm(
    outcome_raw ~
      early_pre_treat + post_common_treat + late_post_treat |
      uid_stack + qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  event_model <- feglm(
    outcome_raw ~
      i(event_time_ref, treat_group, ref = "ref_year") |
      uid_stack + qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  list(
    static_model = static_model,
    static_result = tidy_coefficient(
      static_model,
      "post_common_treat",
      "Bradford: same-force controls"
    ),
    event_model = event_model,
    event_results = extract_event_study(event_model, "event_time_ref"),
    recent_pre_wald = joint_event_test(
      event_model,
      "recent_pre",
      "event_time_ref"
    ),
    post_wald = joint_event_test(event_model, "post", "event_time_ref")
  )
}

fit_bradford_phase_model <- function(bradford) {
  implementation_q <- unique(bradford$caz_start_q)
  if (length(implementation_q) != 1) {
    stop("Bradford must have exactly one implementation quarter.")
  }
  
  # Use four quarters before the announcement as the phase-model reference.
  pre_announcement_start <- BRADFORD_ANNOUNCEMENT_Q - 1
  
  data <- bradford %>%
    filter(
      quarter_year >= pre_announcement_start,
      event_time <= POST_MAX
    ) %>%
    mutate(
      policy_phase = case_when(
        quarter_year < BRADFORD_ANNOUNCEMENT_Q ~ "pre_announcement",
        quarter_year < implementation_q ~ "anticipation",
        event_time >= 0 & event_time <= POST_MAX ~ "implementation",
        TRUE ~ NA_character_
      ),
      policy_phase = factor(
        policy_phase,
        levels = c("pre_announcement", "anticipation", "implementation")
      )
    ) %>%
    filter(!is.na(policy_phase)) %>%
    droplevels()
  
  model <- feglm(
    outcome_raw ~
      i(policy_phase, treat_group, ref = "pre_announcement") |
      uid_stack + qtr_int,
    data = data,
    family = "poisson",
    cluster = ~OA,
    weights = ~analysis_weight,
    lean = FALSE
  )
  
  ct <- coeftable(model)
  phase_results <- tibble(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    p_value = ct[, "Pr(>|z|)"]
  ) %>%
    filter(str_detect(term, "^policy_phase::")) %>%
    mutate(
      phase = str_match(
        term,
        "^policy_phase::([^:]+):treat_group$"
      )[, 2]
    ) %>%
    add_effect_scale() %>%
    add_significance()
  
  beta <- coef(model)
  variance <- vcov(model)
  anticipation_term <- "policy_phase::anticipation:treat_group"
  implementation_term <- "policy_phase::implementation:treat_group"
  
  if (!all(c(anticipation_term, implementation_term) %in% names(beta))) {
    stop("Could not identify both Bradford phase coefficients.")
  }
  
  contrast_estimate <- beta[implementation_term] - beta[anticipation_term]
  contrast_se <- sqrt(
    variance[implementation_term, implementation_term] +
      variance[anticipation_term, anticipation_term] -
      2 * variance[implementation_term, anticipation_term]
  )
  contrast_z <- contrast_estimate / contrast_se
  
  launch_increment <- tibble(
    term = "Implementation minus anticipation",
    estimate = as.numeric(contrast_estimate),
    se = as.numeric(contrast_se),
    p_value = 2 * pnorm(abs(contrast_z), lower.tail = FALSE)
  ) %>%
    add_effect_scale() %>%
    add_significance()
  
  list(
    model = model,
    phase_results = phase_results,
    launch_increment = launch_increment
  )
}

# =============================================================================
# 7. Optional plots (display only; never saved)
# =============================================================================

plot_scheme_effects <- function(results) {
  ggplot(
    results,
    aes(x = pct_change, y = forcats::fct_reorder(scheme, pct_change))
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.15) +
    geom_point(size = 2.2, colour = "#1f5a94") +
    labs(
      title = "Scheme-specific CAZ effects",
      subtitle = "Common post 0:5 relative to -4:-1; full-panel PPML",
      x = "Estimated percentage change in injuries",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

plot_event_study <- function(results, title, colour = "#1f5a94") {
  # The model uses the full panel. The plot focuses on the policy-relevant
  # interval so distant, noisy coefficients do not dominate the scale.
  displayed <- results %>%
    filter(event_time >= -8, event_time <= POST_MAX)
  
  ggplot(displayed, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(
      aes(ymin = ci_lo, ymax = ci_hi),
      fill = colour,
      alpha = 0.15
    ) +
    geom_line(colour = colour, linewidth = 0.8) +
    geom_point(colour = colour, size = 2) +
    scale_x_continuous(breaks = c(-8:-5, 0:POST_MAX)) +
    labs(
      title = title,
      subtitle = "Full-panel estimation; reference = event times -4:-1",
      x = "Quarters relative to implementation",
      y = "Log incidence-rate ratio"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# =============================================================================
# 8. Run analysis and print all principal output
# =============================================================================

section("CAZ injury PPML analysis")
cat("Outcome: ", OUTCOME_VAR, "\n", sep = "")
cat("Static estimand: event times 0:5 relative to -4:-1\n")
cat("Estimation sample: complete matched panel\n")
cat("Police-force x quarter FE: Bradford robustness analysis only\n")

analysis_data <- build_analysis_data()
stacked <- analysis_data$stacked
schemes_all <- analysis_data$schemes

# -----------------------------------------------------------------------------
# Primary static model
# -----------------------------------------------------------------------------

section("2. Primary full-panel static model")

primary_static <- fit_static_full_panel(
  stacked,
  schemes = schemes_all,
  label = "Primary: all seven schemes"
)

subsection("Equal-scheme average post effect")
print_table(primary_static$average)

subsection("Scheme-specific post effects")
print_table(primary_static$scheme_effects)

# -----------------------------------------------------------------------------
# Primary event study
# -----------------------------------------------------------------------------

section("3. Primary full-panel event study")

primary_event <- fit_event_study_full_panel(
  stacked,
  schemes = schemes_all,
  label = "Primary pooled event study"
)

subsection("Policy-relevant event-time coefficients")
print_table(
  primary_event$results %>%
    filter(event_time %in% c(-8:-5, 0:POST_MAX))
)

subsection("Joint recent-pre test: event times -8:-5")
print(primary_event$recent_pre_wald)

subsection("Joint post test: event times 0:5")
print(primary_event$post_wald)

# -----------------------------------------------------------------------------
# Major sensitivity: exclude Bradford
# -----------------------------------------------------------------------------

section("4. Sensitivity excluding Bradford")

schemes_no_bradford <- setdiff(schemes_all, "Bradford")

no_bradford_static <- fit_static_full_panel(
  stacked,
  schemes = schemes_no_bradford,
  label = "Sensitivity: exclude Bradford"
)

no_bradford_event <- fit_event_study_full_panel(
  stacked,
  schemes = schemes_no_bradford,
  label = "Sensitivity event study: exclude Bradford"
)

subsection("Equal-scheme average excluding Bradford")
print_table(no_bradford_static$average)

subsection("Event-study tests excluding Bradford")
print(no_bradford_event$recent_pre_wald)
print(no_bradford_event$post_wald)

# -----------------------------------------------------------------------------
# Bradford-only police-force and same-force analyses
# -----------------------------------------------------------------------------

section("5. Bradford robustness analyses")

subsection("Why police-force adjustment is restricted to Bradford")
print_table(analysis_data$force_overlap)

bradford_police <- fit_bradford_police_fe(stacked)

subsection("Bradford using full controls and police-force x quarter FE")
print_table(bradford_police$static_result)

subsection("Bradford police-FE event-study tests")
print(bradford_police$recent_pre_wald)
print(bradford_police$post_wald)

bradford_same_force_data <- make_bradford_same_force_sample(stacked)

bradford_same_force_counts <- bradford_same_force_data %>%
  distinct(treat_group, uid_stack) %>%
  count(treat_group, name = "road_units")

subsection("Bradford same-force road units")
print_table(bradford_same_force_counts)

bradford_same_force <- fit_bradford_same_force(bradford_same_force_data)

subsection("Bradford using only same-force controls")
print_table(bradford_same_force$static_result)

subsection("Bradford same-force event-study coefficients")
print_table(
  bradford_same_force$event_results %>%
    filter(event_time %in% c(-8:-5, 0:POST_MAX))
)

subsection("Bradford same-force event-study tests")
print(bradford_same_force$recent_pre_wald)
print(bradford_same_force$post_wald)

bradford_phases <- fit_bradford_phase_model(bradford_same_force_data)

subsection("Bradford announcement/anticipation and implementation phases")
print_table(bradford_phases$phase_results)

subsection("Additional implementation effect relative to anticipation")
print_table(bradford_phases$launch_increment)

# -----------------------------------------------------------------------------
# Display figures
# -----------------------------------------------------------------------------

if (SHOW_PLOTS) {
  section("6. Display figures")
  
  print(plot_scheme_effects(primary_static$scheme_effects))
  print(plot_event_study(primary_event$results, "Pooled CAZ event study"))
  print(
    plot_event_study(
      bradford_same_force$event_results,
      "Bradford event study: same-force controls",
      colour = "#b33a3a"
    )
  )
}

# Keep a compact results object in the R session for subsequent inspection.
# It is not written to disk.
results <- list(
  primary_static = primary_static,
  primary_event = primary_event,
  no_bradford_static = no_bradford_static,
  no_bradford_event = no_bradford_event,
  bradford_police = bradford_police,
  bradford_same_force = bradford_same_force,
  bradford_phases = bradford_phases,
  force_overlap = analysis_data$force_overlap
)


# =============================================================================
# FORMAL PRE-TREND TESTS FOR BATH AND BIRMINGHAM
# (self-contained -- includes both helper functions, not previously defined
# in the current script)
# =============================================================================

# -----------------------------------------------------------------------
#  Full pre-period linear slope test
#    Fits ONE linear trend in event_time, interacted with treat_group,
#    using ONLY pre-treatment quarters (event_time < 0). A significant
#    coefficient means treated and control units were on different slopes
#    before implementation.
# -----------------------------------------------------------------------

run_scheme_pretrend_test <- function(stacked, sc) {
  d <- stacked %>%
    filter(stack_scheme == sc, event_time < 0) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ event_time:treat_group | uid_stack + qtr_int,
      data    = d,
      family  = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean    = TRUE
    ),
    error = function(e) {
      cat("Pre-trend test failed for", sc, "-", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  ct <- coeftable(fit)
  tibble(
    scheme    = sc,
    n_pre_obs = nrow(d),
    slope     = ct["event_time:treat_group", "Estimate"],
    se        = ct["event_time:treat_group", "Std. Error"],
    p_value   = ct["event_time:treat_group", "Pr(>|z|)"]
  ) %>%
    mutate(
      pct_slope_per_qtr = 100 * (exp(slope) - 1),
      sig = case_when(
        p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
        p_value < 0.05  ~ "*",   p_value < 0.10 ~ ".", TRUE ~ ""
      )
    )
}

# -----------------------------------------------------------------------
# 2. Windowed (recent pre-period only) linear slope test
#    Same idea, but restricted to the most recent `lookback` pre-treatment
#    quarters -- better powered to catch a late-onset rise that a
#    full-period slope might dilute.
# -----------------------------------------------------------------------

run_scheme_pretrend_test_windowed <- function(stacked, sc, lookback = 8) {
  d <- stacked %>%
    filter(stack_scheme == sc, event_time < 0, event_time >= -lookback) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ event_time:treat_group | uid_stack + qtr_int,
      data    = d,
      family  = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean    = TRUE
    ),
    error = function(e) {
      cat("Windowed pre-trend test failed for", sc, "-", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  ct <- coeftable(fit)
  tibble(
    scheme        = sc,
    lookback_qtrs = lookback,
    n_pre_obs     = nrow(d),
    slope         = ct["event_time:treat_group", "Estimate"],
    se            = ct["event_time:treat_group", "Std. Error"],
    p_value       = ct["event_time:treat_group", "Pr(>|z|)"]
  ) %>%
    mutate(
      pct_slope_per_qtr = 100 * (exp(slope) - 1),
      sig = case_when(
        p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
        p_value < 0.05  ~ "*",   p_value < 0.10 ~ ".", TRUE ~ ""
      )
    )
}

# -----------------------------------------------------------------------
# 3. Bucketed contrast test (third, independent check used for Bradford)
#    Compares the average of event_time -8:-5 against the -4:-1 reference
#    directly, without assuming a linear relationship.
# -----------------------------------------------------------------------

run_scheme_bucket_test <- function(stacked, sc) {
  d <- stacked %>%
    filter(stack_scheme == sc, event_time >= -8, event_time <= -1) %>%
    mutate(early_pre = as.integer(event_time >= -8 & event_time <= -5)) %>%
    droplevels()
  
  fit <- tryCatch(
    feglm(
      outcome_raw ~ early_pre:treat_group | uid_stack + qtr_int,
      data    = d,
      family  = "poisson",
      cluster = ~OA,
      weights = ~analysis_weight,
      lean    = TRUE
    ),
    error = function(e) {
      cat("Bucket test failed for", sc, "-", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  ct <- coeftable(fit)
  tibble(
    scheme  = sc,
    n_obs   = nrow(d),
    estimate = ct["early_pre:treat_group", "Estimate"],
    se       = ct["early_pre:treat_group", "Std. Error"],
    p_value  = ct["early_pre:treat_group", "Pr(>|z|)"]
  ) %>%
    mutate(sig = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05  ~ "*",   p_value < 0.10 ~ ".", TRUE ~ ""
    ))
}

# -----------------------------------------------------------------------
# 4. Run all three tests for Bath and Birmingham
# -----------------------------------------------------------------------

bath_birmingham_pretrend_check <- bind_rows(
  map_dfr(c("Bath", "Birmingham"), ~run_scheme_pretrend_test(stacked, .x)) %>%
    mutate(test = "Full-period slope", .before = 1),
  map_dfr(c("Bath", "Birmingham"), ~run_scheme_pretrend_test_windowed(stacked, .x, lookback = 8)) %>%
    mutate(test = "8-qtr windowed slope", .before = 1)
)

print(bath_birmingham_pretrend_check %>%
        select(test, scheme, n_pre_obs, slope, se, p_value, sig),
      n = Inf)

bath_birmingham_bucket_check <- map_dfr(c("Bath", "Birmingham"), ~run_scheme_bucket_test(stacked, .x))
print(bath_birmingham_bucket_check, n = Inf)

# -----------------------------------------------------------------------
#  run the same three tests for ALL schemes at once, for a
#    single consolidated table (useful now that the matching logic differs
#    between schemes)
# -----------------------------------------------------------------------

all_scheme_pretrend_check <- bind_rows(
  map_dfr(schemes_all, ~run_scheme_pretrend_test(stacked, .x)) %>%
    mutate(test = "Full-period slope", .before = 1),
  map_dfr(schemes_all, ~run_scheme_pretrend_test_windowed(stacked, .x, lookback = 8)) %>%
    mutate(test = "8-qtr windowed slope", .before = 1)
)

print(all_scheme_pretrend_check %>%
        select(test, scheme, n_pre_obs, slope, p_value, sig) %>%
        arrange(scheme, test),
      n = Inf)

bath_birmingham_pretrend_check <- map_dfr(c("Bath", "Birmingham"), function(sc) {
  bind_rows(
    run_scheme_pretrend_test(stacked, sc) %>% mutate(test = "full-period slope"),
    run_scheme_pretrend_test_windowed(stacked, sc, lookback = 8) %>% mutate(test = "8-qtr windowed slope")
  )
})
print(bath_birmingham_pretrend_check)


# =============================================================================
# WITH vs. WITHOUT BRADFORD -- HEADLINE COMPARISON
# =============================================================================

headline_comparison <- bind_rows(
  primary_static$average %>% mutate(spec = "All 7 schemes (with Bradford)", .before = 1),
  no_bradford_static$average %>% mutate(spec = "6 schemes (Bradford excluded)", .before = 1)
) %>%
  select(spec, n_schemes, estimate, se, pct_change, pct_lo, pct_hi, p_value, sig)

print(headline_comparison, n = Inf)

# Event-study post-treatment window, same comparison
event_comparison <- bind_rows(
  primary_event$results %>% filter(event_time %in% 0:POST_MAX) %>%
    mutate(spec = "All 7 schemes (with Bradford)", .before = 1),
  no_bradford_event$results %>% filter(event_time %in% 0:POST_MAX) %>%
    mutate(spec = "6 schemes (Bradford excluded)", .before = 1)
) %>%
  select(spec, event_time, estimate, se, pct_change, pct_lo, pct_hi)

print(event_comparison, n = Inf)

# Joint post-treatment Wald tests, side by side
cat("\nWith Bradford:\n"); print(primary_event$post_wald)
cat("\nWithout Bradford:\n"); print(no_bradford_event$post_wald)



### ## #Bradford's individual effect collapses under police-force adjustment,
# is  explained by a West Yorkshire-wide reporting/recording shift, 
#and shows near-identical treated-vs-untreated-Bradford-OA changes
## so . . . 
 # There is no statistically  average CAZ effect on injuries across the six clean schemes.
# The headline pooled estimate, excluding the compromised scheme, is small, positive, 
#and not distinguishable from zero (+4.2%, p=0.658).
#  The one place a joint significant pattern remains is when Bradford is included 
#— and there is a strong reason to believe that's substantially an artifact of the reporting-regime 
# confound, not a genuine, broader CAZ effect being diluted by exclusion. This is the opposite conclusion from where this conversation started (where Bradford's significance looked like the strongest evidence of a real effect) — the weight of evidence has moved from "Bradford is the clearest signal" to "Bradford is the clearest artifact."
## the story






cat("\nOutcome file used by models\n")
stacked %>%
  summarise(
    rows = n(),
    road_units = n_distinct(uid_stack),
    total_injuries = sum(outcome_raw, na.rm = TRUE),
    nonzero_rows = sum(outcome_raw > 0, na.rm = TRUE)
  ) %>%
  print()

stacked %>%
  group_by(stack_scheme, treat_group) %>%
  summarise(
    road_units = n_distinct(uid_stack),
    total_injuries = sum(outcome_raw, na.rm = TRUE),
    nonzero_rows = sum(outcome_raw > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)


#no injuries anywhere in the panel

stacked %>%
  group_by(stack_scheme, treat_group, uid_stack) %>%
  summarise(
    injuries_all_quarters = sum(outcome_raw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  count(
    stack_scheme,
    treat_group,
    zero_injury_road = injuries_all_quarters == 0
  ) %>%
  print(n = Inf)



# .



