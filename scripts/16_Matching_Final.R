# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING  — VERSION F
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a Difference-in-Differences (DiD)
#   analysis of a road safety intervention in Great Britain.
#
# 
#
# OUTPUTS (Section 13):
#   OA_matched_treated_A.rds       — treated OA IDs + weights + stratum, A
#   OA_matched_donors_A.rds        — control OA IDs + weights, Analysis A
#   OA_matched_treated_B.rds       — treated OA IDs + weights + stratum, B
#   OA_matched_donors_B.rds        — control OA IDs + weights, Analysis B
#   OA_matched_full_A.rds          — full matched dataset, Analysis A
#   OA_matched_full_B.rds          — full matched dataset, Analysis B
#   OA_common_support_flags.rds    — structurally isolated treated OA flags
#   OA_outcome_covariates.rds      — recommended xformla covariate list
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

# force the conflict resolution:
select <- dplyr::select
filter <- dplyr::filter


OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))
glimpse(OA_matching_dataset)

print(table(OA_matching_dataset$assignment))
cat("\n--- Zero-injury OA counts ---\n")
print(table(OA_matching_dataset$zero_injury_OA))

print(table(
  OA_matching_dataset$zero_injury_OA,
  OA_matching_dataset$treated_OA,
  dnn = c("zero_injury", "treated")
))

cat("\n--- Road network summary ---\n")
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
stage1_road   <- c("road_density_m_km2", "road_length_km",
                   "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban  <- c("dist_citycentre", "pop_density")
stage1_socdem <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                   "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct")
stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

stage2_trends <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm", "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm", "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_total_pkm"
)
stage2_vars <- c(stage2_trends, stage2_levels)

#log-transform vars (after winsorisation)
# many vars are very right-skewed  
# Percentage / bounded variables excluded — log distorts bounded scales.
log_transform_s1 <- c(
 
  "road_length_km",      # km: right tail 134km vs median ~0.5km; ratio ~270x
  "pop_density",         # persons/km2: max/median ~12x after winsorisation
  "dist_citycentre"      # metres: right tail 32km; log makes sense (distance)
)
log_transform_s2_levels <- stage2_levels
# All mean_*_pkm are injury rates — zero-inflated right-skewed counts per km.
# Trend variables: NOT transformed — already log-slopes.

cat("Stage 1 variables:", length(stage1_vars), "\n")
cat("  Log-transformed after winsorisation:", paste(log_transform_s1, collapse = ", "), "\n")
cat("  Untransformed (bounded/percentage):",
    paste(setdiff(stage1_vars, log_transform_s1), collapse = ", "), "\n\n")
cat("Stage 2 variables:", length(stage2_vars), "\n")
cat("  Level vars log-transformed:", paste(log_transform_s2_levels, collapse = ", "), "\n")


# Logged variable names (used in matching formula and covariance computation)
log_names_s1 <- paste0("log1p_", log_transform_s1)
log_names_s2 <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1,
  "road_density_m_km2",          # kept on original scale because loging made it worse skewed 
  setdiff(stage1_vars, log_transform_s1)  # untransformed percentage vars
)
# Stage 2 level vars on log scale; trend vars unchanged
stage2_vars_log <- c(stage2_trends, log_names_s2)

### remove the diplicate road_density_m 
stage1_vars_log <- c(
  log_names_s1,
  setdiff(stage1_vars, log_transform_s1)
)



# =============================================================================
# — BUILD DATASETS + BUFFER DIAGNOSIS  + ZERO-PRE FINDING
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

#  Document zero-pre-injury finding at the top of the pipeline
# - Zero-pre-injury OAs: 2.7% of OAs after excluding outcome-inactive roads
# - Contribution to total injuries: < 0.1%
#  - Scottish zero-pre-injury problematic subset: 4 roads, 4 injuries (0.00%)
#  - These OAs CANNOT materially drive ATT estimates
# - Two-stage design matches zero-pre-injury treated to zero-pre-injury controls
# - OA-level clustering absorbs any residual weight concentration
#  - Analysis B retained as robustness check; departures from A are noise

#
#Buffer OA diagnosis
buffer_eligible <- OA_matching_dataset %>%
  filter(buffer_OA == 1, n_roads > 0) %>% nrow()
cat("--- Buffer OA pool diagnosis ---\n")
cat("  Total buffer OAs:        ", sum(OA_matching_dataset$buffer_OA == 1), "\n")
cat("  Buffer OAs with n_roads > 0:", buffer_eligible, "\n")
cat("  Excluded: contamination risk.\n\n")

base_filter <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

data_A <- base_filter %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1))

 table(data_A$assignment)

data_B <- base_filter %>%
  mutate(across(
    all_of(intersect(stage2_vars, names(.))),
    ~ if_else(zero_injury_OA == 1 & treated_OA == 1, 0, .)
  ))



cat("=== Analysis A (zero-injury EXCLUDED) ===\n")
cat("  Total OAs:", nrow(data_A), "| Treated:", sum(data_A$treat_indicator == 1),
    "| Controls:", sum(data_A$treat_indicator == 0), "\n\n")

cat("=== Analysis B (zero-injury INCLUDED) ===\n")
cat("  Total OAs:", nrow(data_B), "| Treated:", sum(data_B$treat_indicator == 1),
    "  of which zero-injury:",
    sum(data_B$treat_indicator == 1 & data_B$zero_injury_OA == 1),
    "| Controls:", sum(data_B$treat_indicator == 0), "\n\n")

# =============================================================================
#  — WINSORISE + LOG-TRANSFORM STAGE 1 VARIABLES 
# =============================================================================
# (1) winsorise at p1/p99 using full-dataset quantiles,
#            (2) log1p-transform the flagged right-skewed variables,
#            (3) create new columns named log1p_<varname> so the raw values
#               are preserved 
# Percentage variables are winsorised only — no log.


# Pre-transformation skewness check
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

winsorise_and_log_s1 <- function(data, raw_vars, log_vars) {
  # Step 1: winsorise all Stage 1 vars using full-dataset quantiles
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  # Step 2: log1p-transform flagged variables; store as new columns
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  data
}

data_A_clean <- winsorise_and_log_s1(data_A, stage1_vars, log_transform_s1)
data_B_clean <- winsorise_and_log_s1(data_B, stage1_vars, log_transform_s1)

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
data_B_clean <- winsorise_and_log_s2(data_B_clean, stage2_levels, log_transform_s2_levels)

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

# Near-zero variance check — applied to the LOGGED variable names used in matching
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
s1_vars_B <- check_vars(data_B_clean, stage1_vars_log, "Stage 1 / B")

# ===========================================================
#  — STRUCTURAL ZERO DIAGNOSIS FOR STAGE 2 TRENDS
# =============================================================

diagnose_structural_zeros <- function(data, trend_vars, label) {
  treated_rows <- data %>% filter(treat_indicator == 1)
  n_treated    <- nrow(treated_rows)
  zero_diag <- map_df(trend_vars, function(v) {
    vals   <- treated_rows[[v]]
    n_zero <- sum(vals == 0, na.rm = TRUE)
    pct    <- round(100 * n_zero / n_treated, 1)
    tibble(
      variable  = v,
      n_treated = n_treated,
      n_zero    = n_zero,
      pct_zero  = pct,
      flag_drop = pct > 60,
      note      = case_when(
        pct > 80 ~ "CRITICAL: >80% structural zeros",
        pct > 60 ~ "DROP: >60% structural zeros",
        pct > 40 ~ "WARN: >40% structural zeros",
        TRUE     ~ "OK"
      )
    )
  })
  cat("\nStructural zero rates —", label, ":\n")
  print(zero_diag, n = Inf)
  dropped <- setdiff(zero_diag %>% filter(flag_drop) %>% pull(variable),
                     "trend_total_pkm")
  if (length(dropped) > 0)
    cat("\nDropping (>60% zeros):", paste(dropped, collapse = ", "), "\n")
  else
    cat("\nNo trend variables dropped.\n")
  list(diag = zero_diag, dropped = dropped)
}

zero_diag_A <- diagnose_structural_zeros(data_A_clean, stage2_trends, "Analysis A")
zero_diag_B <- diagnose_structural_zeros(data_B_clean, stage2_trends, "Analysis B")

s2_trends_A <- setdiff(stage2_trends, zero_diag_A$dropped)
s2_trends_B <- setdiff(stage2_trends, zero_diag_B$dropped)

# Stage 2 variable sets: retained trend vars (unchanged) + log1p level vars
s2_vars_A_raw <- c(s2_trends_A, log_names_s2)
s2_vars_B_raw <- c(s2_trends_B, log_names_s2)

s2_vars_A <- check_vars(data_A_clean, s2_vars_A_raw, "Stage 2 / A")
s2_vars_B <- check_vars(data_B_clean, s2_vars_B_raw, "Stage 2 / B")

cat("\nFinal Stage 2 trend variables — A:", paste(intersect(s2_vars_A, stage2_trends), collapse = ", "), "\n")
cat("Final Stage 2 trend variables — B:", paste(intersect(s2_vars_B, stage2_trends), collapse = ", "), "\n")
cat("Final Stage 2 level variables (log1p):",
    paste(intersect(s2_vars_A, log_names_s2), collapse = ", "), "\n\n")

# =============================================================================
# — BALANCE IMPROVEMENT TEST FUNCTION
# =============================================================================
# Applied after every matchit() call. Three checks:
#   (a) REDUCTION: mean |SMD| post-match < mean |SMD| pre-match
#   (b) TREND THRESHOLD: max |SMD| on trend variables < 0.1
#       (the parallel trends identifying assumption requires this)
#   (c) VARIANCE RATIO: Var(treated)/Var(control) in [0.5, 2.0] for all vars
#       (extreme variance ratios indicate matching created artificial clusters)
# All failures are warnings, not errors — pipeline continues but records fail.

# Accumulator for all balance test results across stages/analyses
balance_test_log <- list()

run_balance_tests <- function(matchit_obj, trend_vars, label) {
  
  bt  <- bal.tab(matchit_obj, un = TRUE, stats = c("mean.diffs", "variance.ratios"))
  bal <- bt$Balance
  
  smd_un  <- abs(bal$Diff.Un)
  smd_adj <- abs(bal$Diff.Adj)
  
  # (a) Reduction test
  mean_un  <- mean(smd_un,  na.rm = TRUE)
  mean_adj <- mean(smd_adj, na.rm = TRUE)
  test_a   <- mean_adj < mean_un
  cat(sprintf("  [TEST a] SMD reduction: %.3f → %.3f  %s\n",
              mean_un, mean_adj, if (test_a) "PASS ✓" else "FAIL ✗"))
  
  # (b) Trend threshold: max trend |SMD| < 0.1
  trend_in_bal  <- intersect(trend_vars, rownames(bal))
  max_trend_smd <- if (length(trend_in_bal) > 0)
    max(abs(bal[trend_in_bal, "Diff.Adj"]), na.rm = TRUE) else NA_real_
  test_b <- !is.na(max_trend_smd) && max_trend_smd < 0.1
  cat(sprintf("  [TEST b] Max trend |SMD|: %.4f  %s\n",
              max_trend_smd, if (test_b) "PASS ✓" else "FAIL ✗ (>0.1 weakens parallel trends)"))
  
  # (c) Variance ratio test: Var.Ratio.Adj in [0.5, 2.0]
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

# ===================================
# - STAGE 1 MATCHING 
# =============================

run_stage1 <- function(data, s1_vars, label, trend_vars_for_test = NULL) {
  
  cat("--- Stage 1:", label, "---\n")
  formula <- reformulate(s1_vars, response = "treat_indicator")
  
  m <- tryCatch(
    matchit(formula, data = data, method = "nearest",
            distance = "mahalanobis", ratio = 10,
            replace = TRUE, exact = ~ country),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  mm <- m$match.matrix
  cat("  Treated in match matrix:", nrow(mm), "\n")
  
  #  Pooled covariance matrix — all rows, not treated-only
  pool_idx    <- c(as.integer(rownames(mm)), unique(as.integer(mm[!is.na(mm)])))
  S_s1_pooled <- cov(data[pool_idx, s1_vars], use = "pairwise.complete.obs")
  
  dist_s1 <- map_df(seq_len(nrow(mm)), function(i) {
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    c_indices  <- mm[i, ]
    c_indices  <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    map_df(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      tibble(
        treated_OA = treated_id,
        control_OA = crow[["OA"]],
        mdist      = mahalanobis(
          x      = as.numeric(crow[s1_vars]),
          center = as.numeric(trow[s1_vars]),
          cov    = S_s1_pooled
        )
      )
    })
  })
  
  treated_matched     <- data[as.integer(rownames(mm)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  control_row_indices <- unique(as.integer(mm[!is.na(mm)]))
  controls_matched    <- data[control_row_indices, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated retained:        ", nrow(treated_matched), "\n")
  cat("  Unique controls in pool: ", nrow(controls_matched), "\n")
  
  # Common support assessment
  min_dist_per_treated <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  threshold_95 <- quantile(min_dist_per_treated$min_dist_s1, 0.95)
  isolated_OAs <- min_dist_per_treated %>%
    filter(min_dist_s1 > threshold_95) %>%
    mutate(structurally_isolated = TRUE)
  
  cat("\n  Common support:\n")
  cat("    95th pct threshold:", round(threshold_95, 3), "\n")
  cat("    Isolated OAs:", nrow(isolated_OAs), "/", nrow(treated_matched), "\n")
  cat("    Min-distance distribution:\n")
  print(summary(min_dist_per_treated$min_dist_s1))
  assign(paste0("isolated_OAs_", label), isolated_OAs, envir = .GlobalEnv)
  
  #  Balance improvement tests
  cat("\n   Balance tests:\n")
  run_balance_tests(m, trend_vars = character(0), label = paste0("S1_", label))
  
  # Balance table
  cat("\n  Stage 1 balance (SMD):\n")
  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  
  cat(names(bt$Balance), "\n")
  # Once you know the real names, e.g. "Diff.Un" and "Diff.Adj":
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    dplyr::select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  
  print(smd_df)
  
  #  Flag Stage 1 variables with |SMD| > 0.1 → outcome model covariates
  s1_imbalanced <- smd_df %>%
    filter(abs(Diff.Adj) > 0.1, variable != "country_Scotland") %>%
    pull(variable)
  if (length(s1_imbalanced) > 0) {
    cat("\n   Stage 1 vars |SMD| > 0.1 → MUST enter C&S xformla:\n")
    cat("    ", paste(s1_imbalanced, collapse = ", "), "\n")
  }
  
  cat("\n  Country distribution (treated):\n")
  print(table(treated_matched$country))
  
  # Cross-country check
  cross_check <- dist_s1 %>%
    left_join(data %>% select(OA, country) %>% rename(control_country = country),
              by = c("control_OA" = "OA")) %>%
    left_join(data %>% select(OA, country) %>% rename(treated_country = country),
              by = c("treated_OA" = "OA")) %>%
    filter(treated_country != control_country)
  cat("  Cross-country pairs:", nrow(cross_check),
      if (nrow(cross_check) == 0) "✓" else "WARNING!", "\n")
  
  lp <- love.plot(
    m, threshold = 0.1, abs = TRUE, stars = "std",
    var.order = "unadjusted",
    title     = paste("Stage 1 balance —", label),
    shapes    = c("circle filled", "triangle filled"),
    colors    = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) + theme(plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  ggsave(here("output", paste0("s1_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  list(matchit_obj = m, dist_s1 = dist_s1,
       treated = treated_matched, controls = controls_matched,
       bal = bt, s1_imbalanced = s1_imbalanced,
       min_dist_per_treated = min_dist_per_treated,
       isolated_OAs = isolated_OAs)
}

s1_A <- run_stage1(data_A_clean, s1_vars_A, "A_excl_zero")
s1_B <- run_stage1(data_B_clean, s1_vars_B, "B_incl_zero")

cat("Full balance table — Analysis A:\n")
bal.tab(s1_A$matchit_obj, un = TRUE)
cat("Full balance table — Analysis B:\n")
bal.tab(s1_B$matchit_obj, un = TRUE)

# =============================================================================
# PREPARE STAGE 2 DATA: DEDUP + WINSORISE + LOG 
# =============================================================================
# Winsorisation anchored to treated distribution (Stage 2 only).
# Log1p applied to level variables after winsorisation.
# Trend variables: winsorised only (already log-slopes; zero-spike handled
# by structural zero filter


prepare_s2_data <- function(s1_result, s2_trend_vars, s2_level_vars_raw) {
  
  s2_raw <- bind_rows(
    s1_result$treated,
    s1_result$controls
  ) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  #  Explicit loop — no cur_column() fragility
  # Winsorise trend variables (anchored to treated distribution)
  for (v in intersect(s2_trend_vars, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  
  #  Winsorise then log1p-transform level variables
  for (v in intersect(s2_level_vars_raw, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    v_winsor        <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
    s2_raw[[paste0("log1p_", v)]] <- log1p(pmax(v_winsor, 0))
  }
  
  s2_raw
}

s2_data_A <- prepare_s2_data(s1_A, s2_trends_A, stage2_levels)
s2_data_B <- prepare_s2_data(s1_B, s2_trends_B, stage2_levels)

cat("Skewness of Stage 2 level variables BEFORE and AFTER log1p:\n")
map_df(stage2_levels, function(v) {
  raw_col <- s2_data_A[[v]]
  log_col <- s2_data_A[[paste0("log1p_", v)]]
  tibble(
    variable    = v,
    skew_raw    = round(skew_fn(raw_col), 2),
    skew_log1p  = round(skew_fn(log_col), 2),
    improvement = abs(skew_fn(raw_col)) > abs(skew_fn(log_col))
  )
}) %>% print()

# =============================================================================
# — COMMON SUPPORT SENSITIVITY DATASETS 
# =============================================================================


make_restricted_data <- function(s2_data, isolated_OAs, label) {
  isolated_ids <- isolated_OAs$treated_OA
  n_removed    <- sum(s2_data$OA %in% isolated_ids & s2_data$treat_indicator == 1)
  cat(label, ": removing", n_removed, "structurally isolated treated OAs\n")
  s2_data %>% filter(!(OA %in% isolated_ids & treat_indicator == 1))
}

s2_data_A_restricted <- make_restricted_data(s2_data_A, s1_A$isolated_OAs, "Analysis A")
s2_data_B_restricted <- make_restricted_data(s2_data_B, s1_B$isolated_OAs, "Analysis B")

# =============================================================================
# — STAGE 2 MATCHING FUNCTION
# =============================================================================

run_stage2 <- function(data, s2_vars, ratio, label, trend_vars) {
  
  formula <- reformulate(s2_vars, response = "treat_indicator")
  
  m <- matchit(formula, data = data, method = "nearest",
               distance = "mahalanobis", ratio = ratio, replace = TRUE)
  
  mm <- m$match.matrix
  
  # Treated-only covariance — correct for Stage 2 ATT-focused MDM
  S_s2 <- cov(data[as.integer(rownames(mm)), s2_vars],
              use = "pairwise.complete.obs")
  
  dist_s2 <- map_df(seq_len(nrow(mm)), function(i) {
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    c_indices  <- mm[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(as.numeric(crow[s2_vars]), as.numeric(trow[s2_vars]), S_s2)
    })
    tibble(OA = treated_id,
           control_OA    = data[as.integer(c_indices), "OA", drop = FALSE] %>% pull(OA),
           mdist_control = dists,
           mdist_treated = mean(dists))
  })
  
  treated_dists <- dist_s2 %>% distinct(OA, mdist_treated) %>% rename(mdist = mdist_treated)
  control_dists <- dist_s2 %>% group_by(OA = control_OA) %>%
    summarise(mdist = mean(mdist_control), .groups = "drop")
  all_dists <- bind_rows(treated_dists, control_dists) %>%
    group_by(OA) %>% summarise(mdist = mean(mdist), .groups = "drop")
  
  matched_data <- match.data(m) %>% left_join(all_dists, by = "OA")
  cat("  Re-run with exact constraint: treated", sum(matched_data$treat_indicator==1),
      "controls", sum(matched_data$treat_indicator==0), "\n")

  cat("  Note: all_dists computed from pre-constraint run — mdist values approximate\n")
  cat("Treated matched:  ", sum(matched_data$treat_indicator == 1), "\n")
  cat("Controls matched: ", sum(matched_data$treat_indicator == 0), "\n")
  cat("NAs in mdist:     ", sum(is.na(matched_data$mdist)), "\n")
  
  # Stage 2 country integrity check
  if ("country" %in% names(matched_data)) {
    cross_s2 <- dist_s2 %>%
      left_join(data %>% select(OA, country) %>% rename(ctrl_country = country),
                by = c("control_OA" = "OA")) %>%
      left_join(data %>% select(OA, country) %>% rename(trt_country = country),
                by = c("OA")) %>%
      filter(trt_country != ctrl_country)
    cat(" Stage 2 cross-country pairs:", nrow(cross_s2),
        if (nrow(cross_s2) == 0) "✓" else "WARNING!", "\n")
    if (nrow(cross_s2) > 0) {
      cat("  Adding exact=~country to Stage 2 to resolve cross-country leakage\n")
      m <- matchit(formula, data = data, method = "nearest",
                   distance = "mahalanobis", ratio = ratio, replace = TRUE,
                   exact = ~ country)
      matched_data <- match.data(m) %>% left_join(all_dists, by = "OA")
      cat("  Re-run with exact constraint: treated", sum(matched_data$treat_indicator==1),
          "controls", sum(matched_data$treat_indicator==0), "\n")
    }
  }
  
  # Balance improvement tests
  cat("\n   Balance tests:\n")
  run_balance_tests(m, trend_vars = trend_vars, label = label)
  
  cat("Distance summary (treated only):\n")
  print(summary(matched_data$mdist[matched_data$treat_indicator == 1]))
  
  list(matchit_obj = m, primary_ratio = ratio,
       primary_data = matched_data, dist_s2 = dist_s2)
}

# =============================================================================
# — RATIO SELECTION VIA ELBOW CURVE 
# =============================================================================

compute_trend_smd_at_ratio <- function(data, s2_vars, ratio) {
  trend_vars <- intersect(s2_vars, stage2_trends)
  formula    <- reformulate(s2_vars, response = "treat_indicator")
  m <- tryCatch(
    matchit(formula, data = data, method = "nearest",
            distance = "mahalanobis", ratio = ratio, replace = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  bt  <- bal.tab(m, un = FALSE)
  smd <- bt$Balance %>% rownames_to_column("variable") %>%
    filter(variable %in% trend_vars) %>%
    summarise(max_trend_smd  = max(abs(Diff.Adj), na.rm = TRUE),
              mean_trend_smd = mean(abs(Diff.Adj), na.rm = TRUE))
  tibble(ratio = ratio,
         n_controls     = sum(match.data(m)$treat_indicator == 0),
         max_trend_smd  = round(smd$max_trend_smd, 5),
         mean_trend_smd = round(smd$mean_trend_smd, 5))
}

find_elbow <- function(curve, threshold = 0.002, fallback_ratio = 4) {
  curve <- curve %>% arrange(ratio)
  
  # Check if curve is monotonically worsening
  if (all(diff(curve$max_trend_smd) >= 0)) {
    cat("WARNING: SMD increases with ratio — control pool may be too thin.\n")
    cat("Using fallback ratio:", fallback_ratio, "\n")
    return(fallback_ratio)
  }
  
  improvements <- -diff(curve$max_trend_smd)
  elbow_pos    <- which(improvements < threshold)[1]
  
  if (is.na(elbow_pos)) {
    cat("No elbow found — using fallback ratio:", fallback_ratio, "\n")
    return(fallback_ratio)
  }
  
  elbow_ratio <- curve$ratio[elbow_pos]
  cat("Elbow at ratio", elbow_ratio, "\n")
  elbow_ratio
}


ratio_curve_A <- map_df(1:8, ~ compute_trend_smd_at_ratio(s2_data_A, s2_vars_A, .x))
print(ratio_curve_A)

ratio_curve_B <- map_df(1:8, ~ compute_trend_smd_at_ratio(s2_data_B, s2_vars_B, .x))
print(ratio_curve_B)

# Curve A: SMD stabilises around 0.053-0.055 from ratio 3 onwards
# No benefit beyond ratio 3;  ratio 5 gives a comfortable pool
optimal_ratio_A <- 4   # 1599 controls

# Curve B: SMD peaks at ratio 4 then declines — ratio 7-8 actually best
# Use ratio 7 where SMD drops back toward ratio 1 levels
optimal_ratio_B <- 7   # 2130 controls

cat("Ratio A =", optimal_ratio_A, "| trend SMD:", 
    ratio_curve_A$max_trend_smd[ratio_curve_A$ratio == optimal_ratio_A], "\n")
cat("Ratio B =", optimal_ratio_B, "| trend SMD:", 
    ratio_curve_B$max_trend_smd[ratio_curve_B$ratio == optimal_ratio_B], "\n")

p_ratio <- bind_rows(
  ratio_curve_A %>% mutate(analysis = "A: excl zero-injury"),
  ratio_curve_B %>% mutate(analysis = "B: incl zero-injury")
) %>%
  ggplot(aes(x = ratio, y = max_trend_smd, colour = analysis)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey50") +
  labs(title = "Ratio selection: max trend |SMD| vs ratio",
       x = "Ratio", y = "Max trend |SMD|", colour = "Analysis") +
  theme_minimal()
ggsave(here("output", "ratio_curve.png"), p_ratio, width = 10, height = 6, dpi = 300)

# =============================================================================
#  RUN FINAL STAGE 2 MATCHING + BALANCE REPORTING
# =============================================================================
s2_A <- run_stage2(s2_data_A, s2_vars_A, optimal_ratio_A, "A_excl_zero",
                   trend_vars = intersect(s2_vars_A, stage2_trends))
s2_B <- run_stage2(s2_data_B, s2_vars_B, optimal_ratio_B, "B_incl_zero",
                   trend_vars = intersect(s2_vars_B, stage2_trends))

s2_A_restricted <- run_stage2(s2_data_A_restricted, s2_vars_A, optimal_ratio_A,
                              "A_restricted", trend_vars = intersect(s2_vars_A, stage2_trends))
s2_B_restricted <- run_stage2(s2_data_B_restricted, s2_vars_B, optimal_ratio_B,
                              "B_restricted", trend_vars = intersect(s2_vars_B, stage2_trends))

# Balance reporting function
run_balance <- function(s2_result, label, s1_imbalanced = NULL) {
  
  cat("--- Balance (ratio:", s2_result$primary_ratio, ") —", label, "---\n")
  bt     <- bal.tab(s2_result$matchit_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    mutate(var_type = case_when(
      variable %in% stage2_trends  ~ "Trend",
      variable %in% log_names_s2   ~ "Level (log1p)",
      TRUE                         ~ "Other"
    ))
  
  cat("\n  TREND SMDs (target: all < 0.06):\n")
  trend_s <- smd_df %>% filter(var_type == "Trend") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(trend_s)
  if (any(abs(trend_s$Diff.Adj) >= 0.1, na.rm = TRUE))
    cat("  WARNING: trend SMD >= 0.1\n")
  
  cat("\n  LEVEL SMDs — log1p scale (residual imbalance → C&S outcome model):\n")
  level_s <- smd_df %>% filter(var_type == "Level (log1p)") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(level_s)
  s2_level_covs <- level_s %>% filter(abs(Diff.Adj) > 0.1) %>% pull(variable)
  
  smd_all <- abs(bt$Balance$Diff.Adj)
  cat("\n  Balanced (<0.1):", sum(smd_all < 0.1, na.rm = TRUE), "/", length(smd_all),
      "| Max |SMD|:", round(max(smd_all, na.rm = TRUE), 3),
      "| Mean |SMD|:", round(mean(smd_all, na.rm = TRUE), 3), "\n")
  
  #  Full covariate list for outcome model (Stage 1 + Stage 2 level)
  all_covs <- unique(c(s1_imbalanced, s2_level_covs))
  if (length(all_covs) > 0) {
    cat("\n   Recommended C&S xformla covariates:\n")
    cat("    ", paste(all_covs, collapse = ", "), "\n")
  }
  
  lp <- love.plot(
    s2_result$matchit_obj, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted", title = paste("Stage 2 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) + theme(axis.text.y = element_text(size = 9),
            plot.margin  = margin(10, 30, 10, 180),
            legend.position = "bottom")
  ggsave(here("output", paste0("s2_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  invisible(list(bal = bt, s2_level_covs = s2_level_covs, all_covs = all_covs))
}

bal_A <- run_balance(s2_A,            "A_excl_zero",            s1_A$s1_imbalanced)
bal_B <- run_balance(s2_B,            "B_incl_zero",            s1_B$s1_imbalanced)
run_balance(s2_A_restricted, "A_excl_zero_restricted", s1_A$s1_imbalanced)
run_balance(s2_B_restricted, "B_incl_zero_restricted", s1_B$s1_imbalanced)

# =============================================================================
#  — COMPARE A vs B + DISTANCE DIAGNOSTICS
# =============================================================================

treated_A <- s2_A$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
only_in_B <- setdiff(treated_B, treated_A)
only_in_A <- setdiff(treated_A, treated_B)
in_both   <- intersect(treated_A, treated_B)

cat("--- Sample overlap ---\n")
cat("  Treated in A:", length(treated_A), "| Treated in B:", length(treated_B), "\n")
cat("  In both:", length(in_both), "| Only in B (zero-injury):", length(only_in_B), "\n\n")

if (length(only_in_A) > 0)
  cat("WARNING: OAs in A but not B — check exclusion logic\n")

if (length(only_in_B) > 0) {
  cat("--- Structural characteristics of zero-injury OAs added in B ---\n")
  compare_grp <- s2_B$primary_data %>%
    filter(treat_indicator == 1) %>%
    mutate(grp = if_else(OA %in% only_in_B,
                         "Zero-injury (B only)", "Injury-exposed (both)"))
  s1_avail <- intersect(stage1_vars, names(compare_grp))
  cov_diff <- compare_grp %>%
    select(grp, all_of(s1_avail)) %>% pivot_longer(-grp) %>%
    group_by(name, grp) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = c(mean, sd)) %>%
    mutate(SMD = (`mean_Injury-exposed (both)` - `mean_Zero-injury (B only)`) /
             sqrt((`sd_Injury-exposed (both)`^2 + `sd_Zero-injury (B only)`^2) / 2)) %>%
    arrange(desc(abs(SMD)))
  print(cov_diff %>%
          select(name, `mean_Injury-exposed (both)`, `mean_Zero-injury (B only)`, SMD),
        n = Inf)
  
  cat("\n   Zero-pre-injury OA impact assessment:\n")
  cat("  Diagnostic confirmed these OAs contribute < 0.1% of sample injuries.\n")
  cat("  Scottish problematic subset: 4 roads, 4 injuries (0.00%).\n")
  cat("  Structural distinctiveness does not translate to outcome influence.\n\n")
  
  cat("  Country distribution of zero-injury OAs in B:\n")
  print(table(s2_B$primary_data %>% filter(OA %in% only_in_B) %>% pull(country)))
}

# Distance comparison
cat("\n--- Stage 2 distance distribution: A vs B (treated OAs only) ---\n")
bind_rows(
  s2_A$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "B: incl zero-injury")
) %>%
  group_by(analysis) %>%
  summarise(n = n(), mean = round(mean(mdist, na.rm = TRUE), 3),
            median = round(median(mdist, na.rm = TRUE), 3),
            p90 = round(quantile(mdist, 0.90, na.rm = TRUE), 3),
            p95 = round(quantile(mdist, 0.95, na.rm = TRUE), 3),
            max = round(max(mdist, na.rm = TRUE), 3), .groups = "drop") %>%
  print()

cat("\n--- Distance tail: Analysis A treated OAs ---\n")
s2_A$primary_data %>% filter(treat_indicator == 1) %>%
  summarise(n_over_20 = sum(mdist > 20), n_over_30 = sum(mdist > 30),
            n_over_50 = sum(mdist > 50), n_over_100 = sum(mdist > 100)) %>% print()

# =============================================================================
# — WEIGHT DIAGNOSTICS + CAP
# =============================================================================


weight_diagnostics <- function(s2_result, analysis_label) {
  controls     <- s2_result$primary_data %>% filter(treat_indicator == 0)
  treated_rows <- s2_result$primary_data %>% filter(treat_indicator == 1)
  trend_vars_here <- intersect(names(s2_result$primary_data), stage2_trends)
  
  cat("\n--- Weight diagnostics:", analysis_label, "---\n")
  print(quantile(controls$weights, probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)))
  
  eff_n <- sum(controls$weights)^2 / sum(controls$weights^2)
  cat("\nEffective N:", round(eff_n, 1), "vs nominal N:", nrow(controls),
      "| Efficiency:", round(eff_n / nrow(controls), 3), "\n")
  
  controls %>% group_by(country) %>%
    summarise(n = n(), mean_w = round(mean(weights), 2),
              max_w = round(max(weights), 2),
              eff_n = round(sum(weights)^2 / sum(weights^2), 1), .groups = "drop") %>%
    print()
  
  caps <- c(5, 10, 20, Inf)
  cap_results <- map_df(caps, function(cap) {
    capped    <- controls %>% mutate(w_cap = pmin(weights, cap))
    eff_n_cap <- sum(capped$w_cap)^2 / sum(capped$w_cap^2)
    trend_smds <- map_dbl(trend_vars_here, function(v) {
      t_mean    <- mean(treated_rows[[v]], na.rm = TRUE)
      c_mean    <- weighted.mean(controls[[v]], w = capped$w_cap, na.rm = TRUE)
      pooled_sd <- sqrt((var(treated_rows[[v]], na.rm = TRUE) +
                           var(controls[[v]], na.rm = TRUE)) / 2)
      if (pooled_sd == 0) return(0)
      abs(t_mean - c_mean) / pooled_sd
    })
    tibble(weight_cap    = if_else(is.infinite(cap), "No cap", paste0("Cap ", cap)),
           n_controls    = nrow(controls),
           effective_n   = round(eff_n_cap, 1),
           efficiency    = round(eff_n_cap / nrow(controls), 3),
           max_trend_smd = round(max(trend_smds), 4),
           pct_at_cap    = round(100 * mean(controls$weights >= cap), 2))
  })
  cat("\nWeight cap sensitivity:\n"); print(cap_results)
  invisible(cap_results)
}

wd_A <- weight_diagnostics(s2_A, "A_excl_zero")
wd_B <- weight_diagnostics(s2_B, "B_incl_zero")

s2_A$primary_data <- s2_A$primary_data %>% mutate(weights = pmin(weights, 5))
s2_A_restricted$primary_data <- s2_A_restricted$primary_data %>% 
  mutate(weights = pmin(weights, 5))

s2_B$primary_data <- s2_B$primary_data %>% mutate(weights = pmin(weights, 5))
s2_B_restricted$primary_data <- s2_B_restricted$primary_data %>%
  mutate(weights = pmin(weights, 5))

cat("\n--- Analysis B after cap at 5 ---\n")
s2_B$primary_data %>% filter(treat_indicator == 0) %>%
  summarise(max_weight = max(weights),
            effective_n = round(sum(weights)^2 / sum(weights^2), 1),
            efficiency = round((sum(weights)^2 / sum(weights^2)) / n(), 3)) %>% print()

# =============================================================================
# — BASELINE INJURY LEVEL STRATIFICATION 
# =============================================================================



add_baseline_stratum <- function(s2_result, label) {
  treated_rows <- s2_result$primary_data %>% filter(treat_indicator == 1)
  level_col    <- "log1p_mean_total_pkm"
  raw_col      <- "mean_total_pkm"
  use_col      <- if (level_col %in% names(treated_rows)) level_col else raw_col
  
  if (!use_col %in% names(treated_rows)) {
    cat(label, ": level column not found — skipping stratification\n")
    return(s2_result)
  }
  q_breaks <- quantile(treated_rows[[use_col]], probs = c(0, 0.25, 0.5, 0.75, 1),
                       na.rm = TRUE)
  s2_result$primary_data <- s2_result$primary_data %>%
    mutate(baseline_injury_stratum = case_when(
      treat_indicator == 0 ~ NA_integer_,
      .data[[use_col]] <= q_breaks[2] ~ 1L,
      .data[[use_col]] <= q_breaks[3] ~ 2L,
      .data[[use_col]] <= q_breaks[4] ~ 3L,
      TRUE ~ 4L
    ))
  cat(label, "— stratum distribution (treated OAs):\n")
  print(table(s2_result$primary_data %>% filter(treat_indicator == 1) %>%
                pull(baseline_injury_stratum), useNA = "ifany"))
  cat("  Quartile breaks (", use_col, "):",
      paste(round(q_breaks, 4), collapse = " | "), "\n\n")
  s2_result
}

s2_A <- add_baseline_stratum(s2_A, "Analysis A")
s2_B <- add_baseline_stratum(s2_B, "Analysis B")


#### control reuse

glimpse(
  s2_A
)
matched_data <- s2_A$primary_data


control_reuse <- matched_data %>%
  filter(treat_indicator == 0) %>%
  count(OA) %>%
  summarise(
    n_controls_used = n(),
    max_reuse = max(n),
    mean_reuse = mean(n),
    median_reuse = median(n)
  )

print(control_reuse)

matched_data %>%
  filter(treated_OA == 0) %>%
  count(weights)


matched_data %>%
  filter(treated_OA == 0) %>%
  mutate(approx_reuse = weights / min(weights)) %>%
  summarise(
    max_reuse = max(approx_reuse),
    mean_reuse = mean(approx_reuse),
    median_reuse = median(approx_reuse)
  )

### reuse within schemes 
matched_data %>%
  filter(treated_OA == 0) %>%
  group_by(scheme) %>%
  summarise(
    mean_weight = mean(weights),
    max_weight = max(weights),
    n_controls = n()
  )
### ratio 


matched_data %>%
  summarise(
    n_treated = sum(treated_OA == 1),
    n_controls = sum(treated_OA == 0),
    control_per_treated = n_controls / n_treated
  )

# =============================================================================
# — EXTRACT, CHECK, AND SAVE
# =============================================================================


matched_A_treated  <- s2_A$primary_data %>%
  filter(treat_indicator == 1) %>% select(OA, weights, baseline_injury_stratum)
matched_A_controls <- s2_A$primary_data %>%
  filter(treat_indicator == 0) %>% select(OA, weights)

matched_B_treated  <- s2_B$primary_data %>%
  filter(treat_indicator == 1) %>% select(OA, weights, baseline_injury_stratum)
matched_B_controls <- s2_B$primary_data %>%
  filter(treat_indicator == 0) %>% select(OA, weights)

cat("Analysis A — Treated:", nrow(matched_A_treated),
    "| Controls:", nrow(matched_A_controls), "\n")
cat("  Control weight range: [",
    round(min(matched_A_controls$weights), 3), ",",
    round(max(matched_A_controls$weights), 3), "]\n")
cat("Analysis B — Treated:", nrow(matched_B_treated),
    "| Controls:", nrow(matched_B_controls), "\n")
cat("  Control weight range: [",
    round(min(matched_B_controls$weights), 3), ",",
    round(max(matched_B_controls$weights), 3), "]\n\n")

# Common support flags
common_support_flags <- bind_rows(
  s1_A$isolated_OAs %>% mutate(analysis = "A"),
  s1_B$isolated_OAs %>% mutate(analysis = "B")
)

#  covariate list
outcome_covariates <- list(
  analysis_A = bal_A$all_covs,
  analysis_B = bal_B$all_covs
)
cat("Recommended C&S xformla covariates:\n")
cat("  Analysis A:", paste(outcome_covariates$analysis_A, collapse = ", "), "\n")
cat("  Analysis B:", paste(outcome_covariates$analysis_B, collapse = ", "), "\n\n")

# Consolidated balance test results 
balance_tests_summary <- bind_rows(balance_test_log)
cat("Balance test summary:\n")
print(balance_tests_summary, n = Inf)

# Integrity checks
stopifnot(
  "A: treated weights == 1"     = all(matched_A_treated$weights == 1),
  "B: treated weights == 1"     = all(matched_B_treated$weights == 1),
  "A: no NA control weights"    = !anyNA(matched_A_controls$weights),
  "B: no NA control weights"    = !anyNA(matched_B_controls$weights),
  "A: no duplicate treated OAs" = !anyDuplicated(matched_A_treated$OA),
  "B: no duplicate treated OAs" = !anyDuplicated(matched_B_treated$OA),
  "B: max control weight <= 5"  = max(matched_B_controls$weights) <= 5
)
cat("All integrity checks passed.\n\n")

# Save
saveRDS(matched_A_treated,    here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A_controls,   here("data", "processed", "OA_matched_donors_A.rds"))
saveRDS(matched_B_treated,    here("data", "processed", "OA_matched_treated_B.rds"))
saveRDS(matched_B_controls,   here("data", "processed", "OA_matched_donors_B.rds"))
saveRDS(s2_A$primary_data,    here("data", "processed", "OA_matched_full_A.rds"))
saveRDS(s2_B$primary_data,    here("data", "processed", "OA_matched_full_B.rds"))
saveRDS(common_support_flags, here("data", "processed", "OA_common_support_flags.rds"))
saveRDS(outcome_covariates,   here("data", "processed", "OA_outcome_covariates.rds"))
saveRDS(balance_tests_summary,here("data", "processed", "OA_balance_tests.rds"))

cat("Saved 9 files to data/processed/\n\n")

# =============================================================================
#### — DESIGN SUMMARY
# =============================================================================


cat("WHY MDM OVER PSM:\n")
cat("  Stage 2 has", length(s2_vars_A), "variables (low-dimensional ).\n")
cat("  PSM advantage (dimensionality relief) does not apply here.\n")
cat("  PSM collapses trajectory information into a scalar — misspecification\n")
cat("  destroys balance with no diagnostics. MDM is directly verifiable.\n\n")

cat("TRANSFORMATION PIPELINE :\n")
cat("  Stage 1: winsorise (p1/p99, pooled) → log1p for:\n")
cat("    ", paste(log_transform_s1, collapse = ", "), "\n")
cat("  Stage 2 levels: winsorise (treated-anchored) → log1p for all mean_*_pkm\n")
cat("  Stage 2 trends: winsorise only (already log-slopes; zeros handled)\n")
cat("  Percentage/bounded vars: winsorise only (no log)\n\n")

cat("STAGE 1 | MDM on", length(s1_vars_A), "vars (log-transformed where appropriate)\n")
cat("        | replace=TRUE, ratio=10, exact=~country\n")
cat("        | Pooled covariance matrix\n")
cat("        | Common support flags saved\n")
cat("        | Balance improvement tests \n\n")

cat("STAGE 2 | MDM on", length(s2_vars_A), "vars after FIX 4 + FIX 10\n")
cat("        | replace=TRUE, ratio = elbow-selected [FIX 5]\n")
cat("        | Treated-only covariance; country integrity verified \n")
cat("        | Balance improvement tests\n\n")

cat("ANALYSIS A | Primary | Zero-injury treated OAs EXCLUDED\n")
cat("           | N treated:", nrow(matched_A_treated),
    "| N controls:", nrow(matched_A_controls), "\n\n")

cat("ANALYSIS B | Robustness | Zero-injury treated OAs INCLUDED\n")
cat("           | N treated:", nrow(matched_B_treated),
    "| N controls:", nrow(matched_B_controls), "\n")
cat("           | Control weights capped at 5\n")
cat("           |  Zero-pre-injury OAs: 2.7% of OAs, <0.1% of injuries.\n")
cat("           | Scottish problematic subset: 4 roads, 4 injuries (0.00%).\n")
cat("           | Cannot drive ATT estimates.\n\n")


cat("NEXT STEPS FOR DiD:\n")
cat("  1. Join matched OA lists to panel\n")
cat("  2. att_gt(..., weightsname = 'weights'), cluster = ~OA\n")
cat("  3. xformla includes all vars in OA_outcome_covariates.rds [FIX 3]\n")
cat("     Note: covariate names are log1p_* for transformed vars\n")
cat("  4. Stratum-specific ATT by baseline_injury_stratum [FIX 8]\n")
cat("  5. Sensitivity: exclude OA_common_support_flags OAs [FIX 2]\n")
cat("  6. Check OA_balance_tests.rds for any failing tests [FIX 11]\n")
cat("  7. ATT(A) vs ATT(B): if similar → A primary, B robustness\n")
cat("     if different → heterogeneous effects by pre-period injury exposure\n")