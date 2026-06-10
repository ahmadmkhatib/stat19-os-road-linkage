# =============================================================================
# POOLED MATCHING — PER-SCHEME TWO-STAGE MAHALANOBIS DISTANCE MATCHING
# =============================================================================
#
# Matches treated OAs to other-city control OAs SEPARATELY for each CAZ scheme.
#
#   Stage 1: road network, urban form, and sociodemographic covariates
#   Stage 2: POOLED total injury trend + level only (no mode-specific trends)
#
# This is the simplified pooled matching for a DiD with total injuries
# (any injury whatsoever) as the outcome. Stage 2 ratio is selected
# per scheme (1:1 to 1:10) to minimise |SMD| on the total trend variable.
#
# England only — no Scotland.
#
# OUTPUTS:
#   OA_matched_treated_pooled.rds
#   OA_matched_donors_pooled.rds
#   OA_matched_full_pooled.rds
#   OA_common_support_flags_pooled.rds
#   OA_ratio_selection_pooled.rds
#   OA_balance_tests_pooled.rds
#   OA_matching_pairs_pooled.rds
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(here)
library(MASS)
library(purrr)
library(sf)
library(ggrepel)
library(tidyverse)
library(patchwork)

set.seed(4718)

select <- dplyr::select
filter <- dplyr::filter

dir.create(here("output", "diagnostics", "pooled"),
           showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics", "pooled")

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

# Stage 1:— road, urban, business, sociodemographic
stage1_road     <- c("road_density_m_km2", "road_length_km",
                     "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban    <- c("dist_BUA_centroid", "pop_density", "area_km2")
stage1_business <- c("business_retail_per_km2")
stage1_socdem   <- c(
  "IMD",
  "cars_one_pct", "cars_twoPlus_pct",
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
  "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
  "age_45to64_pct", "age_65to84_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_business, stage1_socdem)

# Stage 2: POOLED — total injury trend + level only
stage2_trends <- c("trend_total_pkm")
stage2_levels <- c("mean_total_pkm")
stage2_vars   <- c(stage2_trends, stage2_levels)

# Log transformations
log_transform_s1        <- c("road_length_km", "pop_density", "dist_BUA_centroid",
                             "road_density_m_km2", "business_retail_per_km2")
log_nozero_s1           <- c("area_km2")
log_transform_s2_levels <- stage2_levels

log_names_s1        <- paste0("log1p_", log_transform_s1)
log_nozero_names_s1 <- paste0("log_",   log_nozero_s1)
log_names_s2        <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1, log_nozero_names_s1,
  setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1))
)
stage2_vars_log <- c(stage2_trends, log_names_s2)

# =============================================================================
# BUILD DATASET — ENGLAND ONLY
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                          ~ "Unknown"
    )
  )

treated_lads <- OA_matching_dataset %>%
  filter(treated_OA == 1, country == "England") %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

# Other-city controls only — same-city OAs excluded to prevent contamination
data_england <- OA_matching_dataset %>%
  filter(
    country           == "England",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA         == 0,
    n_roads            > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("=== ENGLAND dataset (other-city controls only) ===\n")
cat("Treated:", sum(data_england$treat_indicator == 1),
    "| Controls:", sum(data_england$treat_indicator == 0), "\n")
cat("Control-to-treated ratio:",
    round(sum(data_england$treat_indicator == 0) /
            sum(data_england$treat_indicator == 1), 1), "\n\n")

leak <- data_england %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads) %>%
  nrow()
cat("Same-city control leak (should be 0):", leak, "\n\n")

# =============================================================================
# AGE BAND AGGREGATION + WINSORISE + LOG-TRANSFORM
# =============================================================================

prep_dataset <- function(data) {
  data %>%
    mutate(
      age_under15_pct = X4under_pct  + X5to9_pct   + X10to14_pct,
      age_15to24_pct  = X15to19_pct  + X20to24_pct,
      age_25to44_pct  = X25to29_pct  + X30to34_pct + X35to39_pct + X40to44_pct,
      age_45to64_pct  = X45to49_pct  + X50to54_pct + X55to59_pct + X60to64_pct,
      age_65to84_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct
    )
}

data_england <- prep_dataset(data_england)

winsorise_and_log_s1 <- function(data, raw_vars, log_vars,
                                 log_nozero_vars = character(0)) {
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  for (v in intersect(log_nozero_vars, names(data))) {
    data[[paste0("log_", v)]] <- log(data[[v]])
  }
  data
}

winsorise_and_log_s2 <- function(data, raw_vars, log_vars) {
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  data
}

check_vars <- function(data, vars, label) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0)
    cat("WARNING —", label, "missing:", paste(missing, collapse = ", "), "\n")
  vcheck <- data %>%
    summarise(across(all_of(intersect(vars, names(data))),
                     ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")
  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)
  if (length(low) > 0)
    cat("Dropping near-zero variance (", label, "):",
        paste(low, collapse = ", "), "\n")
  setdiff(intersect(vars, names(data)), low)
}

# English schemes
english_schemes <- data_england %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("Schemes:", paste(english_schemes, collapse = ", "), "\n")
cat("N schemes:", length(english_schemes), "\n\n")

# =============================================================================
# BALANCE TEST FUNCTION
# =============================================================================

balance_test_log <- list()

run_balance_tests <- function(matchit_obj, trend_vars, label) {
  bt  <- bal.tab(matchit_obj, un = TRUE,
                 stats = c("mean.diffs", "variance.ratios"))
  bal <- bt$Balance

  mean_un  <- mean(abs(bal$Diff.Un),  na.rm = TRUE)
  mean_adj <- mean(abs(bal$Diff.Adj), na.rm = TRUE)
  test_a   <- mean_adj < mean_un
  cat(sprintf("  [TEST a] Mean |SMD|: %.3f → %.3f  %s\n",
              mean_un, mean_adj, if (test_a) "PASS" else "FAIL"))

  trend_in_bal  <- intersect(trend_vars, rownames(bal))
  max_trend_smd <- if (length(trend_in_bal) > 0)
    max(abs(bal[trend_in_bal, "Diff.Adj"]), na.rm = TRUE) else NA_real_
  test_b <- !is.na(max_trend_smd) && max_trend_smd < 0.1
  cat(sprintf("  [TEST b] Max trend |SMD|: %.4f  %s\n",
              max_trend_smd, if (test_b) "PASS" else "FAIL"))

  vr_col  <- if ("Var.Ratio.Adj" %in% names(bal)) "Var.Ratio.Adj" else NULL
  vr_fail <- character(0)
  if (!is.null(vr_col)) {
    vr      <- bal[[vr_col]]
    vr_fail <- rownames(bal)[!is.na(vr) & (vr < 0.5 | vr > 2.0)]
    test_c  <- length(vr_fail) == 0
    cat(sprintf("  [TEST c] Variance ratio [0.5, 2.0]: %d/%d pass  %s\n",
                sum(is.na(vr) | (vr >= 0.5 & vr <= 2.0), na.rm = TRUE),
                sum(!is.na(vr)),
                if (test_c) "PASS" else "FAIL"))
  } else { test_c <- NA }

  balance_test_log[[label]] <<- tibble(
    label         = label,
    mean_smd_un   = round(mean_un,       4),
    mean_smd_adj  = round(mean_adj,      4),
    max_trend_smd = round(max_trend_smd, 4),
    test_a_pass   = test_a,
    test_b_pass   = test_b,
    test_c_pass   = if (is.na(test_c)) NA else test_c,
    vr_fail_vars  = paste(vr_fail, collapse = "; ")
  )
  invisible(balance_test_log[[label]])
}

# =============================================================================
# RATIO SELECTION FUNCTIONS
# =============================================================================

prepare_s2_for_ratio <- function(data_clean, s1_vars, s2_vars, label) {
  s1v        <- check_vars(data_clean, s1_vars, paste("S1 ratio prep", label))
  formula_s1 <- reformulate(s1v, response = "treat_indicator")
  m_s1 <- matchit(formula_s1, data = data_clean, method = "nearest",
                  distance = "mahalanobis", ratio = 30, replace = TRUE)
  mm_s1       <- m_s1$match.matrix
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  s2_raw      <- bind_rows(treated_s1, controls_s1)
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  for (v in intersect(stage2_levels, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <-
      log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }
  for (v in intersect(stage2_trends, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  s2v <- check_vars(s2_raw, s2_vars, paste("S2 ratio prep", label))
  list(data = s2_raw, s2_vars = s2v)
}

select_ratio <- function(data, s2_vars, trend_vars, total_trend_var, label,
                         ratios_to_test = 1:10) {
  # trend_vars     : all trend variables (used as diagnostics)
  # total_trend_var: the single variable that drives ratio selection
  #                  (e.g. "trend_total_pkm") — parallel trends assumption
  #                  is on total injuries, not individual modes

  n_t <- sum(data$treat_indicator == 1)
  n_c <- sum(data$treat_indicator == 0)
  ratios_to_test <- ratios_to_test[ratios_to_test <= floor(n_c / n_t)]

  cat("\n--- Ratio selection:", label, "---\n")
  cat("  Treated:", n_t, "| Controls:", n_c, "\n")
  cat("  Selection criterion: total_trend_smd (|SMD| on", total_trend_var, ")\n")

  formula_s2 <- reformulate(s2_vars, response = "treat_indicator")

  ratio_results <- map_df(ratios_to_test, function(r) {
    m <- tryCatch(
      matchit(formula_s2, data = data, method = "nearest",
              distance = "mahalanobis", ratio = r, replace = TRUE),
      error = function(e) NULL
    )
    if (is.null(m)) return(tibble(ratio = r, mean_smd = NA,
                                  total_trend_smd  = NA,
                                  max_trend_smd    = NA,
                                  max_level_smd    = NA))

    bt         <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
    trend_rows <- rownames(bt)[rownames(bt) %in% trend_vars]
    level_rows <- rownames(bt)[!rownames(bt) %in% trend_vars]

    # Primary selection target: total trend only
    total_smd <- if (total_trend_var %in% rownames(bt))
      abs(bt[total_trend_var, "Diff.Adj"]) else NA_real_

    tibble(
      ratio           = r,
      mean_smd        = round(mean(abs(bt$Diff.Adj),               na.rm = TRUE), 4),
      total_trend_smd = round(total_smd,                                          4),
      max_trend_smd   = round(max(abs(bt[trend_rows, "Diff.Adj"]), na.rm = TRUE), 4),
      max_level_smd   = round(max(abs(bt[level_rows, "Diff.Adj"]), na.rm = TRUE), 4)
    )
  })

  print(ratio_results)

  # Select on total_trend_smd first; break ties with mean_smd
  best <- ratio_results %>%
    filter(!is.na(total_trend_smd)) %>%
    arrange(total_trend_smd, mean_smd) %>%
    slice(1)

  cat(sprintf("  Selected: 1:%d (total trend |SMD| = %.4f, max trend |SMD| = %.4f)\n",
              best$ratio, best$total_trend_smd, best$max_trend_smd))

  list(optimal_ratio = best$ratio, ratio_results = ratio_results, label = label)
}

# =============================================================================
# MATCHING FUNCTION
# =============================================================================

run_matching <- function(data_clean, s1_vars, s2_vars, ratio,
                         label, trend_vars) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MATCHING —", label, "| ratio 1:", ratio, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # ---- STAGE 1 ---------------------------------------------------------------
  cat("--- Stage 1:", label, "---\n")
  s1v        <- check_vars(data_clean, s1_vars, paste("S1", label))
  formula_s1 <- reformulate(s1v, response = "treat_indicator")

  m_s1 <- tryCatch(
    matchit(formula_s1, data = data_clean, method = "nearest",
            distance = "mahalanobis", ratio = 50, replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s1)) return(NULL)

  mm_s1       <- m_s1$match.matrix
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)

  cat("  Treated:", nrow(treated_s1),
      "| Unique controls in pool:", nrow(controls_s1), "\n")

  # Common support
  pool_idx <- c(as.integer(rownames(mm_s1)), ctrl_idx)
  S_s1 <- cov(data_clean[pool_idx, s1v], use = "pairwise.complete.obs")
  dist_s1 <- map_df(seq_len(nrow(mm_s1)), function(i) {
    t_idx     <- as.integer(rownames(mm_s1)[i])
    trow      <- data_clean[t_idx, , drop = FALSE]
    c_indices <- mm_s1[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    map_df(seq_along(c_indices), function(j) {
      crow <- data_clean[as.integer(c_indices[j]), , drop = FALSE]
      tibble(treated_OA = trow[["OA"]], control_OA = crow[["OA"]],
             mdist = mahalanobis(as.numeric(crow[s1v]),
                                 as.numeric(trow[s1v]), S_s1))
    })
  })

  min_dist <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")

  mad_threshold <- median(min_dist$min_dist_s1) +
    3 * mad(min_dist$min_dist_s1, constant = 1)

  isolated_OAs <- min_dist %>%
    mutate(
      structurally_isolated = min_dist_s1 > mad_threshold,
      flag_threshold        = round(mad_threshold, 4)
    )

  cat(sprintf("  Isolated OAs (median + 3*MAD threshold = %.2f): %d / %d\n",
              mad_threshold,
              sum(isolated_OAs$structurally_isolated),
              nrow(isolated_OAs)))

  cat("\n  Stage 1 balance:\n")
  run_balance_tests(m_s1, trend_vars = character(0),
                    label = paste0("S1_", label))

  # ---- STAGE 2 ---------------------------------------------------------------
  cat("\n--- Stage 2:", label, "---\n")
  s2_raw      <- bind_rows(treated_s1, controls_s1) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)

  for (v in intersect(stage2_trends, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  for (v in intersect(stage2_levels, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <-
      log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }

  s2v        <- check_vars(s2_raw, s2_vars, paste("S2", label))
  formula_s2 <- reformulate(s2v, response = "treat_indicator")

  m_s2 <- tryCatch(
    matchit(formula_s2, data = s2_raw, method = "nearest",
            distance = "mahalanobis", ratio = ratio, replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s2)) return(NULL)

  matched_data <- match.data(m_s2)

  cat("  Treated:", sum(matched_data$treat_indicator == 1),
      "| Controls:", sum(matched_data$treat_indicator == 0), "\n")

  # Distances
  mm_s2 <- m_s2$match.matrix
  S_s2  <- cov(s2_raw[as.integer(rownames(mm_s2)), s2v],
               use = "pairwise.complete.obs")
  dist_s2 <- map_df(seq_len(nrow(mm_s2)), function(i) {
    t_idx     <- as.integer(rownames(mm_s2)[i])
    trow      <- s2_raw[t_idx, , drop = FALSE]
    c_indices <- mm_s2[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- s2_raw[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(as.numeric(crow[s2v]), as.numeric(trow[s2v]), S_s2)
    })
    tibble(OA = trow[["OA"]], mdist = mean(dists))
  })
  matched_data <- matched_data %>% left_join(dist_s2, by = "OA")

  cat("\n  Stage 2 balance:\n")
  run_balance_tests(m_s2, trend_vars = trend_vars,
                    label = paste0("S2_", label))

  # Extract treated → control OA pairs
  mm_pairs    <- m_s2$match.matrix
  treated_oas <- s2_raw$OA[as.integer(rownames(mm_pairs))]
  pairs <- map_df(seq_len(nrow(mm_pairs)), function(i) {
    t_oa  <- treated_oas[i]
    c_idx <- as.integer(mm_pairs[i, ])
    c_idx <- c_idx[!is.na(c_idx)]
    if (length(c_idx) == 0) return(tibble())
    tibble(treated_OA = t_oa, control_OA = s2_raw$OA[c_idx])
  })

  list(matched_data = matched_data, isolated_OAs = isolated_OAs,
       matchit_s2 = m_s2, dist_s2 = dist_s2, pairs = pairs)
}

# =============================================================================
# PER-SCHEME MATCHING LOOP
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PER-SCHEME MATCHING — ", length(english_schemes), " schemes\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Control pool (shared across all schemes)
control_pool <- data_england %>% filter(treat_indicator == 0)
cat("\nShared control pool:", nrow(control_pool), "OAs\n")

all_results      <- list()
all_ratio_tables <- list()

for (s in english_schemes) {
  cat("\n", paste(rep("#", 60), collapse = ""), "\n")
  cat("### SCHEME:", s, "###\n")
  cat(paste(rep("#", 60), collapse = ""), "\n")

  # Build scheme-specific dataset: this scheme's treated + full control pool
  scheme_treated <- data_england %>%
    filter(treat_indicator == 1, scheme == s)

  scheme_data <- bind_rows(scheme_treated, control_pool)

  cat("Treated OAs:", nrow(scheme_treated),
      "| Control pool:", nrow(control_pool), "\n")

  # Winsorise + log-transform for this scheme's data
  scheme_clean <- winsorise_and_log_s1(scheme_data, stage1_vars,
                                       log_transform_s1, log_nozero_s1)
  scheme_clean <- winsorise_and_log_s2(scheme_clean, stage2_levels,
                                       log_transform_s2_levels)

  s1_vars_scheme <- check_vars(scheme_clean, stage1_vars_log,
                               paste("S1", s))
  s2_vars_scheme <- check_vars(scheme_clean, stage2_vars_log,
                               paste("S2", s))

  # Identify trend variables and the total trend selection target
  trend_vars_scheme <- intersect(s2_vars_scheme, stage2_trends)

  total_trend_var_scheme <- intersect(
    trend_vars_scheme,
    c("trend_total_pkm", "trend_total_ksi",
      "trend_total_slight", "trend_total")
  )
  if (length(total_trend_var_scheme) == 0) {
    stop("No total trend variable found for scheme: ", s,
         "\nAvailable trend vars: ", paste(trend_vars_scheme, collapse = ", "))
  }
  total_trend_var_scheme <- total_trend_var_scheme[1]
  cat("  Total trend variable:", total_trend_var_scheme, "\n")

  # Ratio selection for this scheme
  s2_prep <- prepare_s2_for_ratio(scheme_clean, s1_vars_scheme,
                                  s2_vars_scheme, s)

  ratio_result <- select_ratio(
    s2_prep$data, s2_prep$s2_vars,
    trend_vars      = trend_vars_scheme,
    total_trend_var = total_trend_var_scheme,
    label           = s,
    ratios_to_test  = 1:10
  )
  optimal_ratio <- ratio_result$optimal_ratio

  all_ratio_tables[[s]] <- ratio_result$ratio_results %>%
    mutate(scheme = s)

  # Run matching
  result <- run_matching(
    data_clean = scheme_clean,
    s1_vars    = s1_vars_scheme,
    s2_vars    = s2_vars_scheme,
    ratio      = optimal_ratio,
    label      = s,
    trend_vars = trend_vars_scheme
  )

  if (!is.null(result)) {
    result$matched_data <- result$matched_data %>%
      mutate(scheme = s)
    result$pairs <- result$pairs %>%
      mutate(scheme = s)
    result$isolated_OAs <- result$isolated_OAs %>%
      mutate(scheme = s)
    all_results[[s]] <- result
  } else {
    cat("WARNING: Matching FAILED for scheme:", s, "\n")
  }
}

# =============================================================================
# COMBINE ACROSS SCHEMES
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("COMBINING RESULTS ACROSS SCHEMES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

matched_full <- bind_rows(
  map(all_results, ~ .x$matched_data)
)

isolated_combined <- bind_rows(
  map(all_results, ~ .x$isolated_OAs)
)

pairs_pooled <- bind_rows(
  map(all_results, ~ .x$pairs)
)

# Per-scheme summary
scheme_summary <- matched_full %>%
  group_by(scheme) %>%
  summarise(
    n_treated  = sum(treat_indicator == 1),
    n_controls = sum(treat_indicator == 0),
    ratio      = round(n_controls / n_treated, 1),
    .groups    = "drop"
  )

cat("Per-scheme matching results:\n")
print(scheme_summary)
cat("\nTotal treated:", sum(scheme_summary$n_treated),
    "| Total control rows:", sum(scheme_summary$n_controls), "\n")
cat("Unique control OAs:", n_distinct(
  matched_full$OA[matched_full$treat_indicator == 0]), "\n\n")

# =============================================================================
# BASELINE INJURY STRATIFICATION
# =============================================================================

q_breaks <- quantile(
  matched_full %>%
    filter(treat_indicator == 1) %>%
    pull(log1p_mean_total_pkm),
  probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE
)

matched_full <- matched_full %>%
  mutate(baseline_injury_stratum = case_when(
    treat_indicator == 0                ~ NA_integer_,
    log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
    log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
    log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
    TRUE                                ~ 4L
  ))

# =============================================================================
# INTEGRITY CHECKS + SAVE
# =============================================================================

matched_treated <- matched_full %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum, scheme)

matched_controls <- matched_full %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights, scheme)

stopifnot(
  "treated weights == 1"     = all(matched_treated$weights == 1),
  "no NA control weights"    = !anyNA(matched_controls$weights),
  "no duplicate treated OAs" = !anyDuplicated(matched_treated$OA)
)
cat("All integrity checks passed.\n")

saveRDS(matched_treated,
        here("data", "processed", "OA_matched_treated_pooled.rds"))
saveRDS(matched_controls,
        here("data", "processed", "OA_matched_donors_pooled.rds"))
saveRDS(matched_full,
        here("data", "processed", "OA_matched_full_pooled.rds"))
saveRDS(isolated_combined,
        here("data", "processed", "OA_common_support_flags_pooled.rds"))
saveRDS(bind_rows(balance_test_log),
        here("data", "processed", "OA_balance_tests_pooled.rds"))
saveRDS(pairs_pooled,
        here("data", "processed", "OA_matching_pairs_pooled.rds"))

# Ratio selection table
ratio_combined <- bind_rows(all_ratio_tables)
saveRDS(ratio_combined,
        here("data", "processed", "OA_ratio_selection_pooled.rds"))

# Ratio selection plot (per scheme)
p_ratio <- ratio_combined %>%
  filter(!is.na(max_trend_smd)) %>%
  ggplot(aes(x = ratio, y = max_trend_smd,
             colour = scheme, group = scheme)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.10, linetype = "dashed", colour = "#888888") +
  geom_hline(yintercept = 0.05, linetype = "dotted", colour = "#555555") +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title    = "Ratio selection — pooled matching (per scheme)",
    subtitle = "Maximum trend |SMD| by matching ratio",
    x        = "Matching ratio (1:k)",
    y        = "Maximum trend |SMD| after matching",
    colour   = "Scheme"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "fig_ratio_selection_pooled.png"),
       p_ratio, width = 13, height = 7, dpi = 300, bg = "white")

# --- Per-scheme Stage 2 ratio summary ---
cat("\n=== STAGE 2 RATIO SELECTION PER SCHEME ===\n\n")

ratio_summary <- ratio_combined %>%
  filter(!is.na(max_trend_smd)) %>%
  group_by(scheme) %>%
  arrange(total_trend_smd, mean_smd) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(
    matched_full %>%
      group_by(scheme) %>%
      summarise(
        n_treated  = sum(treat_indicator == 1),
        n_controls = sum(treat_indicator == 0),
        .groups = "drop"
      ),
    by = "scheme"
  ) %>%
  select(scheme, n_treated, n_controls,
         selected_ratio = ratio, total_trend_smd, max_trend_smd, mean_smd) %>%
  arrange(scheme)

print(ratio_summary)

cat("\nRatio range across schemes:", min(ratio_summary$selected_ratio),
    "to", max(ratio_summary$selected_ratio), "\n")
cat("Schemes with ratio 1:1:",
    sum(ratio_summary$selected_ratio == 1), "\n\n")

# --- Overall balance summary ---
balance_summary <- bind_rows(balance_test_log) %>%
  mutate(stage = if_else(str_starts(label, "S1"), "Stage 1", "Stage 2")) %>%
  group_by(stage) %>%
  summarise(
    mean_abs_smd_unmatched = mean(mean_smd_un,  na.rm = TRUE),
    mean_abs_smd_matched   = mean(mean_smd_adj, na.rm = TRUE),
    max_trend_smd          = max(max_trend_smd, na.rm = TRUE),
    all_test_a_pass        = all(test_a_pass,   na.rm = TRUE),
    all_test_b_pass        = all(test_b_pass,   na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== OVERALL BALANCE SUMMARY (across all schemes) ===\n")
print(balance_summary)

# Per-scheme Stage 2 balance
s2_balance <- bind_rows(balance_test_log) %>%
  filter(str_starts(label, "S2_")) %>%
  mutate(scheme = str_remove(label, "^S2_"))

cat("\n=== STAGE 2 BALANCE PER SCHEME ===\n")
print(s2_balance %>% select(scheme, mean_smd_adj, max_trend_smd,
                             test_a_pass, test_b_pass))

cat("\n=== OUTPUTS SAVED ===\n")
cat("  OA_matched_full_pooled.rds\n")
cat("  OA_matched_treated_pooled.rds\n")
cat("  OA_matched_donors_pooled.rds\n")
cat("  OA_common_support_flags_pooled.rds\n")
cat("  OA_ratio_selection_pooled.rds\n")
cat("  OA_balance_tests_pooled.rds\n")
cat("  OA_matching_pairs_pooled.rds\n")
cat("  fig_ratio_selection_pooled.png\n")
