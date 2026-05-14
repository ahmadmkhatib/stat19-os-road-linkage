# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING  
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a  (DiD) 
#
#
# OUTPUTS:
#   OA_matched_treated_A.rds       — treated OA IDs + weights + stratum, A
#   OA_matched_donors_A.rds        — control OA IDs + weights, Analysis A
#   OA_matched_full_A.rds          — full matched dataset, Analysis A
#   OA_common_support_flags.rds    — structurally isolated treated OA flags
#   OA_outcome_covariates.rds      —  xformla covariate list
#   OA_balance_tests.rds           — balance improvement test results
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


select <- dplyr::select
filter <- dplyr::filter


OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))
glimpse(OA_matching_dataset)
### 

print(table(OA_matching_dataset$assignment))
### remove same city controls 

print(table(OA_matching_dataset$scheme))
cat("\n--- Zero-injury OA counts ---\n")
print(table(OA_matching_dataset$zero_injury_OA))

print(table(
  OA_matching_dataset$zero_injury_OA,
  OA_matching_dataset$treated_OA,
  dnn = c("zero_injury", "treated")
))

OA_matching_dataset %>%
  summarise(
    total_OAs  = n_distinct(OA),
    zero_roads = sum(n_roads == 0 | is.na(n_roads)),
    pct_zero   = round(100 * zero_roads / total_OAs, 2)
  ) %>% print()

weirdOAs <- OA_matching_dataset %>%
  filter((n_roads == 0 | is.na(n_roads)) & mean_total > 0)
cat("Zero-road OAs with recorded injuries:", nrow(weirdOAs),
    " (treated:", sum(weirdOAs$treated_OA == 1), ")\n")

# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

stage1_road     <- c("road_density_m_km2", "road_length_km",
                     "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban    <- c("dist_citycentre", "pop_density", "area_km2")
# Retail only: accommodation/food dropped (r = 0.80 with retail; redundant and
# worsens Mahalanobis stability without improving post-matching balance).
stage1_business <- c("business_retail_per_km2")
stage1_socdem   <- c(
  "IMD",
  # Car ownership: cars_none_pct dropped (smallest ~28%; avoids perfect collinearity)
  "cars_one_pct", "cars_twoPlus_pct",
  # Travel mode to work: Motorcycle_pct dropped (smallest share ~0.43%; reference category)
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  # Ethnicity: Other_ethnicity_pct dropped (smallest ~2%)
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
  # Age structure: 5 broad groups instead of 17 individual bands.
  # Individual bands create near-singular covariance matrices (compositional constraint);
  # broad groups reduce dimensionality while retaining age structure information.
  "age_under15_pct",   # 0–14
  "age_15to24_pct",    # 15–24
  "age_25to44_pct",    # 25–44
  "age_45to64_pct",    # 45–64
  "age_65plus_pct"     # 65–84 (85+ dropped as smallest group ~2%)
)
stage1_vars   <- c(stage1_road, stage1_urban, stage1_business, stage1_socdem)

stage2_trends <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm", "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm", "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_other_KSI_pkm", "mean_other_slight_pkm",
  "mean_total_pkm"
)

stage2_vars <- c(stage2_trends, stage2_levels)

# log-transform vars (after winsorisation)
# many vars are very right-skewed
# Percentage / bounded variables excluded — log distorts bounded scales.
# Variables to transform with log1p (may contain zeros)
log_transform_s1 <- c(
  "road_length_km",          # km: right tail 134km vs median ~0.5km; ratio ~270x
  "pop_density",             # persons/km2: right-skewed
  "dist_citycentre",         # metres: right tail 32km
  "road_density_m_km2",      # scale consistency with road length
  "business_retail_per_km2"  # density: zero-inflated right tail
)
# Variables with no zeros — use log() not log1p()
# area_km2 median is ~0.05 km2; log1p(0.05) ≈ 0.05 (no compression), log(0.05) = -3.0 (proper spread)
log_nozero_s1 <- c("area_km2")
log_transform_s2_levels <- stage2_levels
# All mean_*_pkm are injury rates — zero-inflated right-skewed counts per km.
# Trend variables: NOT transformed — already log-slopes.

cat("Stage 1 variables:", length(stage1_vars), "\n")
cat("  log1p-transformed (may have zeros):", paste(log_transform_s1, collapse = ", "), "\n")
cat("  log-transformed (no zeros):", paste(log_nozero_s1, collapse = ", "), "\n")
cat("  Untransformed (bounded/percentage):",
    paste(setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1)), collapse = ", "), "\n\n")
cat("Stage 2 variables:", length(stage2_vars), "\n")
cat("  Level vars log-transformed:", paste(log_transform_s2_levels, collapse = ", "), "\n")


# Logged variable names (used in matching formula and covariance computation)
log_names_s1       <- paste0("log1p_", log_transform_s1)
log_nozero_names_s1 <- paste0("log_", log_nozero_s1)
log_names_s2       <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1,
  log_nozero_names_s1,
  setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1))
)

# Stage 2 level vars on log scale; trend vars unchanged
stage2_vars_log <- c(stage2_trends, log_names_s2)

# =============================================================================
# — BUILD DATASET
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

# Buffer OA diagnosis
buffer_eligible <- OA_matching_dataset %>%
  filter(buffer_OA == 1, n_roads > 0) %>% nrow()
cat("--- Buffer OA pool diagnosis ---\n")
cat("  Total buffer OAs:        ", sum(OA_matching_dataset$buffer_OA == 1), "\n")
cat("  Buffer OAs with n_roads > 0:", buffer_eligible, "\n")
cat("  Excluded: contamination risk.\n\n")

# LADs that contain CAZ/LEZ treated OAs — control_group1 OAs from these cities
# are excluded to avoid contamination (same city as the intervention).
treated_lads <- OA_matching_dataset %>%
  filter(treated_OA == 1) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

cat("CAZ/LEZ cities (LADs) whose control_group1 OAs will be excluded:\n")
print(treated_lads)

# Analysis A: zero-injury treated OAs excluded upfront;
# control_group1 OAs from treated cities also excluded.
data_A <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1),
    !(control_group1_OA == 1 & LAD24CD %in% treated_lads)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("\ncontrol_group1 OAs removed (same city as CAZ/LEZ):",
    sum(OA_matching_dataset$control_group1_OA == 1 &
          OA_matching_dataset$LAD24CD %in% treated_lads), "\n")

cat("=== Analysis A (zero-injury EXCLUDED) ===\n")
cat("  Total OAs:", nrow(data_A), "| Treated:", sum(data_A$treat_indicator == 1),
    "| Controls:", sum(data_A$treat_indicator == 0), "\n\n")

table(data_A$assignment)

# Collapse 17 individual 5-year age bands into 5 broad life-stage groups before
# winsorising, to reduce near-singularity in the Mahalanobis covariance matrix.
data_A <- data_A %>%
  mutate(
    age_under15_pct = X4under_pct  + X5to9_pct   + X10to14_pct,
    age_15to24_pct  = X15to19_pct  + X20to24_pct,
    age_25to44_pct  = X25to29_pct  + X30to34_pct + X35to39_pct + X40to44_pct,
    age_45to64_pct  = X45to49_pct  + X50to54_pct + X55to59_pct + X60to64_pct,
    age_65plus_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct
  )

# =============================================================================
# — WINSORISE + LOG-TRANSFORM STAGE 1 VARIABLES
# =============================================================================

skew_fn <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x); mu <- mean(x); s <- sd(x)
  if (s == 0 || n < 3) return(NA_real_)
  sum((x - mu)^3) / (n * s^3)
}

cat("Skewness of Stage 1 variables BEFORE transformation (full dataset):\n")
map_df(c(log_transform_s1, setdiff(stage1_vars, log_transform_s1)), function(v) {
  x  <- data_A[[v]]
  q  <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
  xw <- pmin(pmax(x, q[1]), q[2])
  tibble(
    variable      = v,
    will_log      = v %in% log_transform_s1,
    skew_raw      = round(skew_fn(x), 2),
    skew_winsor   = round(skew_fn(xw), 2),
    skew_log1p    = if (v %in% log_transform_s1)
      round(skew_fn(log1p(pmax(xw, 0))), 2) else NA_real_,
    max_med_ratio = round(max(x, na.rm = TRUE) /
                            (median(x, na.rm = TRUE) + 1e-9), 1)
  )
}) %>% arrange(desc(abs(skew_winsor))) %>% print(n = Inf)

winsorise_and_log_s1 <- function(data, raw_vars, log_vars, log_nozero_vars = character(0)) {
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  # log1p for variables that contain zeros
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  # log (not log1p) for variables  have no zeros
  for (v in intersect(log_nozero_vars, names(data))) {
    data[[paste0("log_", v)]] <- log(data[[v]])
  }
  data
}

data_A_clean <- winsorise_and_log_s1(data_A, stage1_vars, log_transform_s1, log_nozero_s1)

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

data_A_clean <- winsorise_and_log_s2(data_A_clean, stage2_levels, log_transform_s2_levels)

cat("\nSkewness of log-transformed Stage 1 variables AFTER transformation:\n")
map_df(log_transform_s1, function(v) {
  raw_col <- data_A_clean[[v]]
  log_col <- data_A_clean[[paste0("log1p_", v)]]
  tibble(
    variable     = v,
    skew_winsor  = round(skew_fn(raw_col), 2),
    skew_log1p   = round(skew_fn(log_col), 2),
    improvement  = abs(skew_fn(raw_col)) > abs(skew_fn(log_col))
  )
}) %>% print()

# Near-zero variance check
check_vars <- function(data, vars, label) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0)
    cat("WARNING —", label, "missing vars:", paste(missing, collapse = ", "), "\n")
  vcheck <- data %>%
    summarise(across(all_of(intersect(vars, names(data))),
                     ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")
  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)
  if (length(low) > 0)
    cat("Dropping near-zero variance (", label, "):", paste(low, collapse = ", "), "\n")
  setdiff(intersect(vars, names(data)), low)
}

s1_vars_A <- check_vars(data_A_clean, stage1_vars_log, "Stage 1 / A")

# =============================================================================
# — STAGE 2 VARIABLE PREP: ALL TRENDS RETAINED
# =============================================================================
# No structural zero filtering — all trend variables passed to Stage 2.
# Near-zero variance check still applied as a safety net.

cat("\nStructural zero rates — Analysis A (informational only, no dropping):\n")
treated_rows_A <- data_A_clean %>% filter(treat_indicator == 1)
map_df(stage2_trends, function(v) {
  vals   <- treated_rows_A[[v]]
  n_zero <- sum(vals == 0, na.rm = TRUE)
  pct    <- round(100 * n_zero / nrow(treated_rows_A), 1)
  tibble(
    variable  = v,
    n_treated = nrow(treated_rows_A),
    n_zero    = n_zero,
    pct_zero  = pct,
    note      = case_when(
      pct > 80 ~ ">80% structural zeros",
      pct > 60 ~ ">60% structural zeros",
      pct > 40 ~ ">40% structural zeros",
      TRUE     ~ "OK"
    )
  )
}) %>% print(n = Inf)
cat("All trend variables retained.\n\n")

# Stage 2 variable sets: all trend vars + log1p level vars
s2_vars_A_raw <- c(stage2_trends, log_names_s2)
s2_vars_A     <- check_vars(data_A_clean, s2_vars_A_raw, "Stage 2 / A")

cat("Final Stage 2 trend variables — A:", paste(intersect(s2_vars_A, stage2_trends), collapse = ", "), "\n")
cat("Final Stage 2 level variables (log1p):",
    paste(intersect(s2_vars_A, log_names_s2), collapse = ", "), "\n\n")

# =============================================================================
# — BALANCE IMPROVEMENT TEST FUNCTION
# =============================================================================

balance_test_log <- list()

run_balance_tests <- function(matchit_obj, trend_vars, label) {
  
  bt  <- bal.tab(matchit_obj, un = TRUE, stats = c("mean.diffs", "variance.ratios"))
  bal <- bt$Balance
  
  smd_un  <- abs(bal$Diff.Un)
  smd_adj <- abs(bal$Diff.Adj)
  
  mean_un  <- mean(smd_un,  na.rm = TRUE)
  mean_adj <- mean(smd_adj, na.rm = TRUE)
  test_a   <- mean_adj < mean_un
  cat(sprintf("  [TEST a] SMD reduction: %.3f → %.3f  %s\n",
              mean_un, mean_adj, if (test_a) "PASS ✓" else "FAIL ✗"))
  
  trend_in_bal  <- intersect(trend_vars, rownames(bal))
  max_trend_smd <- if (length(trend_in_bal) > 0)
    max(abs(bal[trend_in_bal, "Diff.Adj"]), na.rm = TRUE) else NA_real_
  test_b <- !is.na(max_trend_smd) && max_trend_smd < 0.1
  cat(sprintf("  [TEST b] Max trend |SMD|: %.4f  %s\n",
              max_trend_smd, if (test_b) "PASS ✓" else "FAIL ✗ (>0.1 weakens parallel trends)"))
  
  vr_col <- if ("Var.Ratio.Adj" %in% names(bal)) "Var.Ratio.Adj" else NULL
  if (!is.null(vr_col)) {
    vr       <- bal[[vr_col]]
    vr_fail  <- rownames(bal)[!is.na(vr) & (vr < 0.5 | vr > 2.0)]
    test_c   <- length(vr_fail) == 0
    cat(sprintf("  [TEST c] Variance ratio [0.5, 2.0]: %d/%d vars pass  %s\n",
                sum(is.na(vr) | (vr >= 0.5 & vr <= 2.0), na.rm = TRUE),
                sum(!is.na(vr)),
                if (test_c) "PASS ✓" else "FAIL ✗"))
    if (!test_c)
      cat("    Failing vars:", paste(vr_fail, collapse = ", "), "\n")
  } else {
    test_c <- NA; vr_fail <- character(0)
    cat("  [TEST c] Variance ratio: not available in this cobalt version\n")
  }
  
  result <- tibble(
    label          = label,
    mean_smd_un    = round(mean_un,        4),
    mean_smd_adj   = round(mean_adj,       4),
    max_trend_smd  = round(max_trend_smd,  4),
    test_a_pass    = test_a,
    test_b_pass    = test_b,
    test_c_pass    = if (is.na(test_c)) NA else test_c,
    vr_fail_vars   = paste(vr_fail, collapse = "; ")
  )
  balance_test_log[[label]] <<- result
  invisible(result)
}

# =============================================================================
# — RUN MATCHING SEPARATELY BY COUNTRY THEN COMBINE
# =============================================================================

run_country_matching <- function(data, s1_vars, s2_vars, ratio, country_name, trend_vars) {
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("COUNTRY:", country_name, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  data_country <- data %>% filter(country == country_name)
  
  cat("Total OAs:", nrow(data_country), 
      "| Treated:", sum(data_country$treat_indicator == 1),
      "| Controls:", sum(data_country$treat_indicator == 0), "\n")
  
  # ---- STAGE 1 ----
  cat("\n--- Stage 1:", country_name, "---\n")
  
  s1_vars_country <- check_vars(data_country, s1_vars, paste("Stage 1", country_name))
  
  formula_s1 <- reformulate(s1_vars_country, response = "treat_indicator")
  
  m_s1 <- tryCatch(
    matchit(formula_s1, data = data_country, method = "nearest",
            distance = "mahalanobis", ratio = 10,
            replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s1)) return(NULL)
  
  mm_s1 <- m_s1$match.matrix
  
  treated_s1 <- data_country[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  control_idx_s1 <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_country[control_idx_s1, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated retained:", nrow(treated_s1), "\n")
  cat("  Unique controls in pool:", nrow(controls_s1), "\n")
  
  # Common support
  pool_idx <- c(as.integer(rownames(mm_s1)), control_idx_s1)
  S_s1 <- cov(data_country[pool_idx, s1_vars_country], use = "pairwise.complete.obs")
  
  dist_s1 <- map_df(seq_len(nrow(mm_s1)), function(i) {
    t_idx <- as.integer(rownames(mm_s1)[i])
    trow  <- data_country[t_idx, , drop = FALSE]
    c_indices <- mm_s1[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    map_df(seq_along(c_indices), function(j) {
      crow <- data_country[as.integer(c_indices[j]), , drop = FALSE]
      tibble(
        treated_OA = trow[["OA"]],
        control_OA = crow[["OA"]],
        mdist = mahalanobis(as.numeric(crow[s1_vars_country]),
                            as.numeric(trow[s1_vars_country]), S_s1)
      )
    })
  })
  
  min_dist_per_treated <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  threshold_95 <- quantile(min_dist_per_treated$min_dist_s1, 0.95)
  isolated_OAs <- min_dist_per_treated %>%
    filter(min_dist_s1 > threshold_95) %>%
    mutate(structurally_isolated = TRUE)
  
  cat("  Common support 95th pct threshold:", round(threshold_95, 3), "\n")
  cat("  Isolated OAs:", nrow(isolated_OAs), "/", nrow(treated_s1), "\n")
  
  # Stage 1 balance
  cat("\n  Stage 1 balance tests:\n")
  run_balance_tests(m_s1, trend_vars = character(0), 
                    label = paste0("S1_", country_name))
  
  s1_imbalanced <- bal.tab(m_s1, thresholds = c(m = 0.1), un = TRUE)$Balance %>%
    rownames_to_column("variable") %>%
    filter(abs(Diff.Adj) > 0.1) %>%
    pull(variable)
  cat("  Stage 1 vars |SMD| > 0.1:", paste(s1_imbalanced, collapse = ", "), "\n")
  
  # ---- PREPARE STAGE 2 DATA ----
  s2_raw <- bind_rows(treated_s1, controls_s1) %>%
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
    v_winsor <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
    s2_raw[[paste0("log1p_", v)]] <- log1p(pmax(v_winsor, 0))
  }
  
  s2_vars_country <- check_vars(s2_raw, s2_vars, paste("Stage 2", country_name))
  
  # ---- STAGE 2 ----
  cat("\n--- Stage 2:", country_name, "---\n")
  
  formula_s2 <- reformulate(s2_vars_country, response = "treat_indicator")
  
  m_s2 <- tryCatch(
    matchit(formula_s2, data = s2_raw, method = "nearest",
            distance = "mahalanobis", ratio = ratio,
            replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s2)) return(NULL)
  
  matched_data <- match.data(m_s2)
  
  cat("  Treated matched:", sum(matched_data$treat_indicator == 1), "\n")
  cat("  Controls matched:", sum(matched_data$treat_indicator == 0), "\n")
  
  # Compute distances
  mm_s2 <- m_s2$match.matrix
  S_s2  <- cov(s2_raw[as.integer(rownames(mm_s2)), s2_vars_country],
               use = "pairwise.complete.obs")
  
  dist_s2 <- map_df(seq_len(nrow(mm_s2)), function(i) {
    t_idx <- as.integer(rownames(mm_s2)[i])
    trow  <- s2_raw[t_idx, , drop = FALSE]
    c_indices <- mm_s2[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- s2_raw[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(as.numeric(crow[s2_vars_country]),
                  as.numeric(trow[s2_vars_country]), S_s2)
    })
    tibble(OA = trow[["OA"]],
           mdist_treated = mean(dists))
  })
  
  matched_data <- matched_data %>%
    left_join(dist_s2, by = "OA")
  
  # Stage 2 balance
  cat("\n  Stage 2 balance tests:\n")
  run_balance_tests(m_s2, trend_vars = trend_vars,
                    label = paste0("S2_", country_name))
  
  cat("\n  Distance summary (treated only):\n")
  print(summary(matched_data$mdist[matched_data$treat_indicator == 1]))
  
  # Weight diagnostics
  controls_out <- matched_data %>% filter(treat_indicator == 0)
  eff_n <- sum(controls_out$weights)^2 / sum(controls_out$weights^2)
  cat("\n  Effective N:", round(eff_n, 1), "| Nominal N:", nrow(controls_out),
      "| Efficiency:", round(eff_n / nrow(controls_out), 3), "\n")
  cat("  Max weight:", round(max(controls_out$weights), 3), "\n")
  
  list(
    matchit_s1    = m_s1,
    matchit_s2    = m_s2,
    matched_data  = matched_data,
    isolated_OAs  = isolated_OAs,
    s1_imbalanced = s1_imbalanced,
    s2_vars       = s2_vars_country,
    dist_s2       = dist_s2
  )
}

# =============================================================================
# — RUN FOR EACH COUNTRY
# =============================================================================

result_england  <- run_country_matching(
  data         = data_A_clean,
  s1_vars      = s1_vars_A,
  s2_vars      = s2_vars_A,
  ratio        = optimal_ratio_A,
  country_name = "England",
  trend_vars   = intersect(s2_vars_A, stage2_trends)
)

result_scotland <- run_country_matching(
  data         = data_A_clean,
  s1_vars      = s1_vars_A,
  s2_vars      = s2_vars_A,
  ratio        = optimal_ratio_A,
  country_name = "Scotland",
  trend_vars   = intersect(s2_vars_A, stage2_trends)
)

# =============================================================================
# — COMBINE RESULTS
# =============================================================================

matched_full_A <- bind_rows(
  result_england$matched_data,
  result_scotland$matched_data
)

isolated_OAs_combined <- bind_rows(
  result_england$isolated_OAs %>% mutate(country = "England"),
  result_scotland$isolated_OAs %>% mutate(country = "Scotland")
)

cat("\n=== COMBINED RESULTS ===\n")
cat("Total treated:", sum(matched_full_A$treat_indicator == 1), "\n")
cat("Total controls:", sum(matched_full_A$treat_indicator == 0), "\n")
cat("\nBy country:\n")
matched_full_A %>%
  group_by(country, treat_indicator) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

# =============================================================================
# — BASELINE INJURY STRATIFICATION ON COMBINED DATA
# =============================================================================

treated_combined <- matched_full_A %>% filter(treat_indicator == 1)
q_breaks <- quantile(treated_combined[["log1p_mean_total_pkm"]],
                     probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

matched_full_A <- matched_full_A %>%
  mutate(baseline_injury_stratum = case_when(
    treat_indicator == 0 ~ NA_integer_,
    log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
    log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
    log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
    TRUE ~ 4L
  ))

cat("\nStratum distribution (treated OAs):\n")
print(table(matched_full_A %>% 
              filter(treat_indicator == 1) %>% 
              pull(baseline_injury_stratum)))

# =============================================================================
# — EXTRACT TREATED AND CONTROLS
# =============================================================================

matched_A_treated  <- matched_full_A %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum)

matched_A_controls <- matched_full_A %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

cat("\nFinal — Treated:", nrow(matched_A_treated),
    "| Controls:", nrow(matched_A_controls), "\n")
cat("Control weight range: [",
    round(min(matched_A_controls$weights), 3), ",",
    round(max(matched_A_controls$weights), 3), "]\n")

# =============================================================================
# — INTEGRITY CHECKS AND SAVE
# =============================================================================

stopifnot(
  "A: treated weights == 1"     = all(matched_A_treated$weights == 1),
  "A: no NA control weights"    = !anyNA(matched_A_controls$weights),
  "A: no duplicate treated OAs" = !anyDuplicated(matched_A_treated$OA)
)
cat("All integrity checks passed.\n")

saveRDS(matched_A_treated,       here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A_controls,      here("data", "processed", "OA_matched_donors_A.rds"))
saveRDS(matched_full_A,          here("data", "processed", "OA_matched_full_A.rds"))
saveRDS(isolated_OAs_combined,   here("data", "processed", "OA_common_support_flags.rds"))
