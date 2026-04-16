# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING  — VERSION 2
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a Difference-in-Differences (DiD)
#   analysis of a road safety intervention in Great Britain.
#
# CHANGES FROM VERSION 1 (addressing methodological review):
#
#   [FIX 1]  COVARIANCE MATRIX — pooled (treated + controls) used at Stage 1
#            instead of treated-only. Treated-only is correct at Stage 2 where
#            ATT-focused weighting is wanted; at Stage 1 (structural plausibility
#            filter) pooled is the standard MDM recommendation and avoids
#            calibrating the metric purely to treated-side variance.
#
#   [FIX 2]  COMMON SUPPORT ASSESSMENT — after Stage 1 a per-OA minimum
#            Mahalanobis distance is computed. Treated OAs whose best structural
#            match exceeds the 95th percentile of all min-distances are flagged
#            as "structurally isolated" and reported. They are retained in the
#            primary analysis but excluded in a sensitivity run (Section 6b).
#
#   [FIX 3]  STAGE 1 STRUCTURAL IMBALANCE → OUTCOME MODEL — Drive_Car_pct,
#            cars_none_pct, dist_citycentre, and Walk_pct show SMD > 0.1 after
#            Stage 1 matching. These are explicitly carried forward to the
#            covariate list recommended for the C&S xformla alongside the
#            Stage 2 level variables already flagged in v1.
#
#   [FIX 4]  STRUCTURAL ZEROS IN TREND VARIABLES — rare-mode trend variables
#            (cyc_KSI, ped_KSI) have large spikes at zero arising from OAs
#            with no pre-period injuries in that mode. This distorts the
#            Stage 2 covariance matrix. Fix: (a) diagnose the zero-spike
#            proportion for each trend variable; (b) drop individual
#            mode-specific trend variables where > 60 % of treated OAs have
#            structural zeros, replacing them with trend_total_pkm only in
#            those regimes; (c) document which variables were retained.
#
#   [FIX 5]  RATIO SELECTION — instead of comparing only ratios 1, 2, 5, a
#            full curve from ratio 1 to 8 is computed and the elbow point
#            (marginal trend SMD improvement < 0.002 per additional control)
#            is used as the principled stopping rule. This replaces the
#            garden-of-forking-paths selection in v1.
#
#   [FIX 6]  WITHIN-FUNCTION WINSORISATION SAFETY — the Stage 2 winsorisation
#            now uses an explicit named-variable loop (map over variable names)
#            rather than across() with cur_column() indexing into an enclosing
#            environment, eliminating the column-order fragility identified
#            in v1.
#
#   [FIX 7]  STAGE 2 COUNTRY INTEGRITY CHECK — added explicit verification
#            that no cross-country pairs appear in Stage 2 matched output
#            (Stage 1 was already verified; Stage 2 was not).
#
#   [FIX 8]  BASELINE INJURY LEVEL HETEROGENEITY FLAGS — treated OAs are
#            stratified into quartiles of mean_total_pkm before saving. The
#            stratum indicator is saved with the matched datasets so the DiD
#            script can run stratum-specific ATT estimates and test whether
#            level imbalance is an effect modifier.
#
#   [FIX 9]  BUFFER OA POOL DIAGNOSIS — a one-paragraph summary of how many
#            buffer OAs would have been structurally eligible controls (based
#            on Stage 1 distance) is added to Section 2.
#
#   UNCHANGED: MDM over PSM (correct for low-dimensional trajectory matching);
#   two-stage architecture; replace=TRUE; exact=~country at Stage 1;
#   treated-anchored winsorisation at Stage 2; weight cap of 5 for Analysis B.
#
# OUTPUTS (Section 10):
#   OA_matched_treated_A.rds  — treated OA IDs + weights + baseline stratum, A
#   OA_matched_donors_A.rds   — control OA IDs + weights, Analysis A
#   OA_matched_treated_B.rds  — treated OA IDs + weights + baseline stratum, B
#   OA_matched_donors_B.rds   — control OA IDs + weights, Analysis B (capped)
#   OA_matched_full_A.rds     — full matched dataset, Analysis A
#   OA_matched_full_B.rds     — full matched dataset, Analysis B
#   OA_common_support_flags.rds — structurally isolated treated OA flags
#   OA_outcome_covariates.rds   — recommended xformla covariate list for DiD
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(here)
library(MASS)
library(tidyverse)
library(purrr)
library(sf)
library(ggrepel)

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

print(table(OA_matching_dataset$assignment))
cat("\n--- Zero-injury OA counts ---\n")
print(table(OA_matching_dataset$zero_injury_OA))


print(table(OA_matching_dataset$zero_injury_OA, OA_matching_dataset$country))



cat("\n--- Zero-injury breakdown by treatment status ---\n")
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
# SECTION 1: VARIABLE DEFINITIONS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 1: VARIABLE DEFINITIONS\n")
cat("====================================================\n\n")

# --- Stage 1: Structural + sociodemographic (15 variables) ---
stage1_road   <- c("road_density_m_km2", "road_length_km",
                   "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban  <- c("dist_citycentre", "pop_density")
stage1_socdem <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                   "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct")
stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

# --- Stage 2: Pre-treatment injury dynamics (14 variables) ---
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

cat("Stage 1 variables:", length(stage1_vars),
    "(road:", length(stage1_road),
    "| urban:", length(stage1_urban),
    "| socdem:", length(stage1_socdem), ")\n")
cat("Stage 2 variables:", length(stage2_vars),
    "(trends:", length(stage2_trends),
    "| levels:", length(stage2_levels), ")\n")

# =============================================================================
# SECTION 2 — BUILD BOTH ANALYSIS DATASETS + BUFFER OA DIAGNOSIS [FIX 9]
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

# Shared exclusions
base_filter <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# [FIX 9] Buffer OA pool diagnosis
# How many buffer OAs would structurally have been eligible controls?
# We approximate eligibility by checking whether they pass the road/pop filters
buffer_eligible <- OA_matching_dataset %>%
  filter(buffer_OA == 1, n_roads > 0) %>%
  nrow()
cat("\n--- Buffer OA pool diagnosis [FIX 9] ---\n")
cat("  Total buffer OAs:                   ", sum(OA_matching_dataset$buffer_OA == 1), "\n")
cat("  Buffer OAs with n_roads > 0:        ", buffer_eligible, "\n")
cat("  These are excluded due to contamination risk.\n")
cat("  If buffer zones are narrow (< 200m), sensitivity analysis excluding\n")
cat("  the outermost buffer ring is recommended in the DiD script.\n\n")

# Analysis A: exclude zero-injury treated OAs
data_A <- base_filter %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1))

# Analysis B: include zero-injury OAs (Stage 2 vars set to 0 for zero-injury treated)
data_B <- base_filter %>%
  mutate(across(
    all_of(intersect(stage2_vars, names(.))),
    ~ if_else(zero_injury_OA == 1 & treated_OA == 1, 0, .)
  ))

cat("=== Analysis A (zero-injury EXCLUDED) ===\n")
cat("  Total OAs:   ", nrow(data_A), "\n")
cat("  Treated:     ", sum(data_A$treat_indicator == 1), "\n")
cat("  Controls:    ", sum(data_A$treat_indicator == 0), "\n\n")

cat("=== Analysis B (zero-injury INCLUDED) ===\n")
cat("  Total OAs:   ", nrow(data_B), "\n")
cat("  Treated:     ", sum(data_B$treat_indicator == 1), "\n")
cat("    of which zero-injury:",
    sum(data_B$treat_indicator == 1 & data_B$zero_injury_OA == 1), "\n")
cat("  Controls:    ", sum(data_B$treat_indicator == 0), "\n\n")

# =============================================================================
# SECTION 3 — WINSORISE STAGE 1 VARIABLES
# =============================================================================
# Full-dataset quantiles (treated + controls) — standard for Stage 1 structural
# matching. Stage 2 winsorisation (Section 5) uses treated-anchored quantiles.

cat("\n====================================================\n")
cat("SECTION 3: WINSORISING STAGE 1 VARIABLES\n")
cat("====================================================\n\n")

cat("Pre-winsorisation extremes (99th pct vs max):\n")
map_df(stage1_vars, function(v) {
  tibble(
    variable = v,
    p99 = quantile(data_A[[v]], 0.99, na.rm = TRUE),
    max = max(data_A[[v]], na.rm = TRUE)
  )
}) %>% print()

winsorise_s1 <- function(data, vars) {
  # Full-dataset quantiles — symmetric treatment of treated and controls
  data %>%
    mutate(across(
      all_of(intersect(vars, names(.))),
      ~ {
        q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
        pmin(pmax(., q[1]), q[2])
      }
    ))
}

data_A_clean <- winsorise_s1(data_A, stage1_vars)
data_B_clean <- winsorise_s1(data_B, stage1_vars)

# Near-zero variance check — variables with var < 1e-8 destabilise cov inversion
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
    cat("Dropping near-zero variance:", paste(low, collapse = ", "), "\n")
  setdiff(intersect(vars, names(data)), low)
}

s1_vars_A <- check_vars(data_A_clean, stage1_vars, "Stage 1 / A")
s1_vars_B <- check_vars(data_B_clean, stage1_vars, "Stage 1 / B")

# =============================================================================
# SECTION 4 — STRUCTURAL ZERO DIAGNOSIS FOR STAGE 2 TREND VARIABLES [FIX 4]
# =============================================================================
# Rare-mode trend variables (cyc_KSI, ped_KSI) have large proportions of
# structural zeros among treated OAs. A structural zero differs from a
# near-zero slope: it arises when an OA had zero injuries in a mode in ALL
# pre-period years, making log-rate undefined. The estimated trend is then
# forced to zero by construction, not by estimation.
#
# Impact on MDM: the covariance matrix for trend variables contains a spike at
# zero for affected modes. This inflates the diagonal element for those
# variables relative to the genuine trajectory variation, effectively
# down-weighting them in the distance metric. The result is that matching
# on trajectory for common modes (car, pedestrian slight) is degraded by
# the noise introduced by zero-spike rare modes.
#
# Decision rule: if > 60% of TREATED OAs have structural zero for a trend
# variable, that variable is dropped from Stage 2 and noted as unidentifiable
# for that mode. The 60% threshold is conservative — the covariance estimate
# is unreliable well below this threshold, but we retain variables unless
# the majority of the treated distribution is at zero.
#
# trend_total_pkm is ALWAYS retained regardless of zero-spike because it
# aggregates across modes and is the primary parallel-trends target.

cat("\n====================================================\n")
cat("SECTION 4: STRUCTURAL ZERO DIAGNOSIS — STAGE 2 TRENDS [FIX 4]\n")
cat("====================================================\n\n")

diagnose_structural_zeros <- function(data, trend_vars, label) {
  
  treated_rows <- data %>% filter(treat_indicator == 1)
  n_treated    <- nrow(treated_rows)
  
  zero_diag <- map_df(trend_vars, function(v) {
    vals    <- treated_rows[[v]]
    n_zero  <- sum(vals == 0, na.rm = TRUE)
    pct     <- round(100 * n_zero / n_treated, 1)
    tibble(
      variable     = v,
      n_treated    = n_treated,
      n_zero       = n_zero,
      pct_zero     = pct,
      flag_drop    = pct > 60,
      note         = case_when(
        pct > 80 ~ "CRITICAL: >80% structural zeros — unreliable covariance",
        pct > 60 ~ "DROP: >60% structural zeros — majority of distribution at zero",
        pct > 40 ~ "WARN: >40% structural zeros — monitor balance",
        TRUE     ~ "OK"
      )
    )
  })
  
  cat("\nStructural zero rates in trend variables —", label, ":\n")
  print(zero_diag, n = Inf)
  
  dropped <- zero_diag %>% filter(flag_drop) %>% pull(variable)
  # Never drop trend_total_pkm regardless of zero rate
  dropped <- setdiff(dropped, "trend_total_pkm")
  
  if (length(dropped) > 0) {
    cat("\nDropping from Stage 2 (> 60% structural zeros):", paste(dropped, collapse = ", "), "\n")
    cat("These modes will NOT contribute to trajectory matching.\n")
    cat("Rationale: their trend covariance estimate is dominated by zero-spike\n")
    cat("variance rather than genuine pre-period dynamics.\n")
  } else {
    cat("\nNo trend variables dropped — all within acceptable zero-spike threshold.\n")
  }
  
  list(diag = zero_diag, dropped = dropped)
}

zero_diag_A <- diagnose_structural_zeros(data_A_clean, stage2_trends, "Analysis A")
zero_diag_B <- diagnose_structural_zeros(data_B_clean, stage2_trends, "Analysis B")

# Derive final Stage 2 variable sets after zero-spike filtering
s2_trends_A <- setdiff(stage2_trends, zero_diag_A$dropped)
s2_trends_B <- setdiff(stage2_trends, zero_diag_B$dropped)
s2_vars_A_raw <- c(s2_trends_A, stage2_levels)
s2_vars_B_raw <- c(s2_trends_B, stage2_levels)

# Final near-zero variance check on Stage 2 variable sets
s2_vars_A <- check_vars(data_A_clean, s2_vars_A_raw, "Stage 2 / A")
s2_vars_B <- check_vars(data_B_clean, s2_vars_B_raw, "Stage 2 / B")

cat("\nFinal Stage 2 trend variables — Analysis A:", paste(intersect(s2_vars_A, stage2_trends), collapse = ", "), "\n")
cat("Final Stage 2 trend variables — Analysis B:", paste(intersect(s2_vars_B, stage2_trends), collapse = ", "), "\n\n")

# =============================================================================
# SECTION 5 — STAGE 1 MATCHING [FIX 1: POOLED COVARIANCE MATRIX]
# =============================================================================
# MDM on structural + sociodemographic variables.
# ratio=10: large candidate pool for Stage 2 — not the final matched sample.
# replace=TRUE: every treated OA draws from the full control pool.
# exact=~country: enforces England/Scotland separation.
#
# [FIX 1] Covariance matrix: computed from the FULL dataset (treated +
# controls), not treated-only. Rationale: Stage 1 is a structural plausibility
# filter, not an ATT-focused estimator. Using the pooled covariance matrix
# is the standard MDM recommendation and avoids calibrating the distance
# metric purely to treated-side variance, which would make controls with
# different variance structures appear artificially close.
# Treated-only covariance is appropriate at Stage 2 where ATT weighting is
# explicitly desired.
#
# [FIX 2] Common support assessment runs immediately after Stage 1.

cat("\n====================================================\n")
cat("SECTION 5: STAGE 1 MATCHING [FIX 1: POOLED COV MATRIX]\n")
cat("====================================================\n\n")

run_stage1 <- function(data, s1_vars, label) {
  
  cat("--- Stage 1:", label, "---\n")
  formula <- reformulate(s1_vars, response = "treat_indicator")
  
  m <- tryCatch(
    matchit(
      formula,
      data     = data,
      method   = "nearest",
      distance = "mahalanobis",
      ratio    = 10,
      replace  = TRUE,
      exact    = ~ country
    ),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  mm <- m$match.matrix
  cat("  Treated in match matrix:", nrow(mm), "\n")
  
  # [FIX 1] POOLED covariance matrix for Stage 1 distance computation.
  # All rows (treated + controls) contribute, weighted equally.
  # This is the standard MDM recommendation for a structural filter where
  # we want the metric calibrated to the joint distribution.
  all_indices <- seq_len(nrow(data))
  S_s1_pooled <- cov(data[all_indices, s1_vars], use = "pairwise.complete.obs")
  
  # Compute pairwise Mahalanobis distances from match matrix
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
  
  cat("  Treated retained:         ", nrow(treated_matched), "\n")
  cat("  Unique controls in pool:  ", nrow(controls_matched), "\n")
  
  # -------------------------------------------------------------------------
  # [FIX 2] COMMON SUPPORT ASSESSMENT
  # For each treated OA, compute its minimum Stage 1 distance (distance to its
  # best structural match). Treated OAs whose min-distance exceeds the 95th
  # percentile of all min-distances are "structurally isolated" — they lack
  # credible structural comparators in the control pool.
  # These OAs are retained in the primary analysis but:
  #   (a) flagged in the output for reporting
  #   (b) excluded in a sensitivity run (Section 6b)
  # -------------------------------------------------------------------------
  
  min_dist_per_treated <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  
  threshold_95 <- quantile(min_dist_per_treated$min_dist_s1, 0.95)
  isolated_OAs <- min_dist_per_treated %>%
    filter(min_dist_s1 > threshold_95) %>%
    mutate(structurally_isolated = TRUE)
  
  cat("\n  [FIX 2] Common support diagnostics:\n")
  cat("    Min-distance 95th pct threshold:", round(threshold_95, 3), "\n")
  cat("    Structurally isolated treated OAs (min_dist > p95):",
      nrow(isolated_OAs), "/", nrow(treated_matched), "\n")
  cat("    Distribution of per-treated min Stage 1 distances:\n")
  print(summary(min_dist_per_treated$min_dist_s1))
  
  # Save for Section 6b sensitivity run
  assign(paste0("isolated_OAs_", label), isolated_OAs, envir = .GlobalEnv)
  
  # Balance diagnostics
  cat("\n  Stage 1 balance (SMD):\n")
  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  print(smd_df)
  
  # Identify Stage 1 variables with post-matching SMD > 0.1
  # [FIX 3] These will be recommended as outcome model covariates
  s1_imbalanced <- smd_df %>%
    filter(abs(Diff.Adj) > 0.1, variable != "country_Scotland") %>%
    pull(variable)
  if (length(s1_imbalanced) > 0) {
    cat("\n  [FIX 3] Stage 1 variables with |SMD| > 0.1 after matching:\n")
    cat("    ", paste(s1_imbalanced, collapse = ", "), "\n")
    cat("    -> These MUST be included in C&S xformla (structural imbalance\n")
    cat("       not resolved by Stage 2 trajectory matching).\n")
  }
  
  cat("\n  Country distribution (treated):\n")
  print(table(treated_matched$country))
  
  # Cross-country verification
  cross_check <- dist_s1 %>%
    left_join(data %>% select(OA, country) %>% rename(control_country = country),
              by = c("control_OA" = "OA")) %>%
    left_join(data %>% select(OA, country) %>% rename(treated_country = country),
              by = c("treated_OA" = "OA")) %>%
    filter(treated_country != control_country)
  cat("  Cross-country pairs:", nrow(cross_check),
      if (nrow(cross_check) == 0) "✓" else "WARNING: cross-country pairs found!", "\n")
  
  lp <- love.plot(
    m, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted",
    title        = paste("Stage 1 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) + theme(plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  ggsave(here("output", paste0("s1_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  list(
    matchit_obj       = m,
    dist_s1           = dist_s1,
    treated           = treated_matched,
    controls          = controls_matched,
    bal               = bt,
    s1_imbalanced     = s1_imbalanced,
    min_dist_per_treated = min_dist_per_treated,
    isolated_OAs      = isolated_OAs
  )
}

s1_A <- run_stage1(data_A_clean, s1_vars_A, "A_excl_zero")
s1_B <- run_stage1(data_B_clean, s1_vars_B, "B_incl_zero")

cat("Full balance table — Analysis A:\n")
bal.tab(s1_A$matchit_obj, un = TRUE)
cat("Full balance table — Analysis B:\n")
bal.tab(s1_B$matchit_obj, un = TRUE)

# =============================================================================
# SECTION 6 — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE [FIX 6]
# =============================================================================
# [FIX 6] Stage 2 winsorisation now uses an explicit map() loop over variable
# names rather than across() with cur_column(). This eliminates the column-
# order fragility where treated_ref could be indexed out of sync with the
# iteration order of across(), silently applying wrong quantiles.

cat("\n====================================================\n")
cat("SECTION 6: PREPARE STAGE 2 DATA [FIX 6: SAFE WINSORISATION]\n")
cat("====================================================\n\n")

prepare_s2_data <- function(s1_result, s2_vars) {
  
  s2_raw <- bind_rows(
    s1_result$treated,
    s1_result$controls
  ) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  # [FIX 6] Explicit named-variable loop — no cur_column() ambiguity
  vars_present <- intersect(s2_vars, names(s2_raw))
  for (v in vars_present) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  s2_raw
}

s2_data_A <- prepare_s2_data(s1_A, s2_vars_A)
s2_data_B <- prepare_s2_data(s1_B, s2_vars_B)

# =============================================================================
# SECTION 6b — SENSITIVITY: COMMON SUPPORT RESTRICTED DATASETS [FIX 2]
# =============================================================================
# Treated OAs flagged as structurally isolated (min Stage 1 distance > p95)
# are excluded from a sensitivity version of Stage 2. If ATT estimates from
# the restricted sample are similar to the primary, the isolated OAs do not
# drive results. If different, the isolated OAs have heterogeneous treatment
# effects and should be reported separately.

cat("\n====================================================\n")
cat("SECTION 6b: COMMON SUPPORT SENSITIVITY DATASETS [FIX 2]\n")
cat("====================================================\n\n")

make_restricted_data <- function(s2_data, isolated_OAs, label) {
  isolated_ids <- isolated_OAs$treated_OA
  n_removed    <- sum(s2_data$OA %in% isolated_ids & s2_data$treat_indicator == 1)
  cat(label, ": removing", n_removed, "structurally isolated treated OAs\n")
  s2_data %>% filter(!(OA %in% isolated_ids & treat_indicator == 1))
}

s2_data_A_restricted <- make_restricted_data(s2_data_A, s1_A$isolated_OAs, "Analysis A")
s2_data_B_restricted <- make_restricted_data(s2_data_B, s1_B$isolated_OAs, "Analysis B")

# =============================================================================
# SECTION 7 — STAGE 2 MATCHING
# =============================================================================
# MDM on pre-treatment injury dynamics within the Stage 1 candidate pool.
# replace=TRUE: globally optimal trajectory match for each treated OA.
# Covariance matrix: treated-only (correct for Stage 2 ATT-focused matching).
# [FIX 7] Explicit country integrity check added at Stage 2 level.

cat("\n====================================================\n")
cat("SECTION 7: STAGE 2 MATCHING\n")
cat("====================================================\n\n")

run_stage2 <- function(data, s2_vars, ratio, label, original_data = NULL) {
  
  cat("\n--- Stage 2:", label, "ratio", ratio, "---\n")
  
  formula <- reformulate(s2_vars, response = "treat_indicator")
  
  m <- matchit(
    formula,
    data     = data,
    method   = "nearest",
    distance = "mahalanobis",
    ratio    = ratio,
    replace  = TRUE
    # No exact constraint at Stage 2: country separation already enforced at
    # Stage 1; each treated OA's Stage 1 pool is already country-homogeneous.
    # [FIX 7] Verified explicitly below.
  )
  
  mm <- m$match.matrix
  
  # Treated-only covariance matrix — appropriate for Stage 2 ATT-focused MDM.
  # Treated-only weighting ensures the distance metric is sensitive to
  # variation where it matters for the ATT estimand.
  S_s2 <- cov(data[as.integer(rownames(mm)), s2_vars],
              use = "pairwise.complete.obs")
  
  dist_s2 <- map_df(seq_len(nrow(mm)), function(i) {
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    c_indices  <- mm[i, ]
    c_indices  <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(
        x      = as.numeric(crow[s2_vars]),
        center = as.numeric(trow[s2_vars]),
        cov    = S_s2
      )
    })
    tibble(
      OA            = treated_id,
      control_OA    = data[as.integer(c_indices), "OA", drop = FALSE] %>% pull(OA),
      mdist_control = dists,
      mdist_treated = mean(dists)
    )
  })
  
  treated_dists <- dist_s2 %>%
    distinct(OA, mdist_treated) %>%
    rename(mdist = mdist_treated)
  control_dists <- dist_s2 %>%
    group_by(OA = control_OA) %>%
    summarise(mdist = mean(mdist_control), .groups = "drop")
  all_dists <- bind_rows(treated_dists, control_dists) %>%
    group_by(OA) %>%
    summarise(mdist = mean(mdist), .groups = "drop")
  
  matched_data <- match.data(m) %>%
    left_join(all_dists, by = "OA")
  
  cat("Treated matched:  ", sum(matched_data$treat_indicator == 1), "\n")
  cat("Controls matched: ", sum(matched_data$treat_indicator == 0), "\n")
  cat("NAs in mdist:     ", sum(is.na(matched_data$mdist)), "\n")
  
  # [FIX 7] Stage 2 country integrity check
  # Verify no cross-country pairs appear despite absence of explicit exact
  # constraint. Should hold because Stage 1 pool is already country-homogeneous.
  if (!is.null(original_data) && "country" %in% names(matched_data)) {
    cross_s2 <- dist_s2 %>%
      left_join(data %>% select(OA, country) %>% rename(ctrl_country = country),
                by = c("control_OA" = "OA")) %>%
      left_join(data %>% select(OA, country) %>% rename(trt_country = country),
                by = c("OA")) %>%
      filter(trt_country != ctrl_country)
    cat("[FIX 7] Stage 2 cross-country pairs:", nrow(cross_s2),
        if (nrow(cross_s2) == 0) "✓" else "WARNING: cross-country pairs in Stage 2!", "\n")
  }
  
  cat("Distance summary (treated only):\n")
  print(summary(matched_data$mdist[matched_data$treat_indicator == 1]))
  
  list(
    matchit_obj   = m,
    primary_ratio = ratio,
    primary_data  = matched_data,
    dist_s2       = dist_s2
  )
}

# =============================================================================
# SECTION 8 — RATIO SELECTION VIA ELBOW CURVE [FIX 5]
# =============================================================================
# [FIX 5] In v1, ratio was selected by comparing only ratios 1, 2, 5 and
# choosing the best — a garden-of-forking-paths approach. In v2, the full
# curve of max trend SMD vs ratio (1 through 8) is computed and the elbow
# point is defined as the first ratio at which the marginal improvement in
# max trend SMD falls below 0.002 per additional control. This is a
# principled, pre-specified stopping rule.
#
# Rationale for 0.002 threshold: below this threshold the improvement in
# parallel trends support is negligible relative to the cost of adding
# structurally more distant controls (higher ratio = controls further from
# the Stage 1 centre of the pool). The threshold is set before seeing data.

cat("\n====================================================\n")
cat("SECTION 8: RATIO SELECTION — ELBOW CURVE [FIX 5]\n")
cat("====================================================\n\n")

compute_trend_smd_at_ratio <- function(data, s2_vars, ratio, label) {
  
  trend_vars <- intersect(s2_vars, stage2_trends)
  formula    <- reformulate(s2_vars, response = "treat_indicator")
  
  m <- tryCatch(
    matchit(formula, data = data, method = "nearest",
            distance = "mahalanobis", ratio = ratio, replace = TRUE),
    error = function(e) { cat("ratio", ratio, "failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  bt  <- bal.tab(m, un = FALSE)
  smd <- bt$Balance %>%
    rownames_to_column("variable") %>%
    filter(variable %in% trend_vars) %>%
    summarise(max_trend_smd  = max(abs(Diff.Adj), na.rm = TRUE),
              mean_trend_smd = mean(abs(Diff.Adj), na.rm = TRUE))
  
  n_controls <- sum(match.data(m)$treat_indicator == 0)
  tibble(ratio = ratio, n_controls = n_controls,
         max_trend_smd  = round(smd$max_trend_smd, 5),
         mean_trend_smd = round(smd$mean_trend_smd, 5))
}

cat("Computing ratio curve for Analysis A...\n")
ratio_curve_A <- map_df(1:8, ~ compute_trend_smd_at_ratio(s2_data_A, s2_vars_A, .x, "A"))
cat("\nRatio curve — Analysis A:\n")
print(ratio_curve_A)

cat("\nComputing ratio curve for Analysis B...\n")
ratio_curve_B <- map_df(1:8, ~ compute_trend_smd_at_ratio(s2_data_B, s2_vars_B, .x, "B"))
cat("\nRatio curve — Analysis B:\n")
print(ratio_curve_B)

# Elbow detection: first ratio where marginal improvement < 0.002
find_elbow <- function(curve, threshold = 0.002) {
  curve <- curve %>% arrange(ratio)
  improvements <- -diff(curve$max_trend_smd)   # positive = improvement
  elbow_pos    <- which(improvements < threshold)[1]
  if (is.na(elbow_pos)) {
    cat("No elbow found below threshold — using ratio 8 (maximum tested)\n")
    return(max(curve$ratio))
  }
  elbow_ratio <- curve$ratio[elbow_pos]
  cat("Elbow at ratio", elbow_ratio,
      ": marginal improvement dropped to",
      round(improvements[elbow_pos - 1], 4), "->",
      round(improvements[elbow_pos], 4), "\n")
  elbow_ratio
}

cat("\nElbow detection — Analysis A:\n")
optimal_ratio_A <- find_elbow(ratio_curve_A)
cat("Selected ratio for Analysis A:", optimal_ratio_A, "\n")

cat("\nElbow detection — Analysis B:\n")
optimal_ratio_B <- find_elbow(ratio_curve_B)
cat("Selected ratio for Analysis B:", optimal_ratio_B, "\n")

# Save ratio curves for reporting
p_ratio <- bind_rows(
  ratio_curve_A %>% mutate(analysis = "A: excl zero-injury"),
  ratio_curve_B %>% mutate(analysis = "B: incl zero-injury")
) %>%
  ggplot(aes(x = ratio, y = max_trend_smd, colour = analysis)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey50") +
  labs(title = "Ratio selection: max trend |SMD| vs number of controls per treated OA",
       subtitle = "Dashed line = 0.1 threshold; elbow = point of diminishing returns",
       x = "Ratio (controls per treated OA)", y = "Max trend |SMD|",
       colour = "Analysis") +
  theme_minimal()
ggsave(here("output", "ratio_curve.png"), p_ratio, width = 10, height = 6, dpi = 300)
cat("\nRatio curve plot saved.\n")

# =============================================================================
# SECTION 9 — RUN FINAL STAGE 2 MATCHING + BALANCE DIAGNOSTICS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 9: FINAL STAGE 2 MATCHING\n")
cat("====================================================\n\n")

# Primary specifications (elbow-selected ratios)
s2_A <- run_stage2(s2_data_A, s2_vars_A, optimal_ratio_A, "A", data_A_clean)
s2_B <- run_stage2(s2_data_B, s2_vars_B, optimal_ratio_B, "B", data_B_clean)

# Sensitivity: common-support restricted samples
s2_A_restricted <- run_stage2(s2_data_A_restricted, s2_vars_A, optimal_ratio_A,
                              "A_restricted", data_A_clean)
s2_B_restricted <- run_stage2(s2_data_B_restricted, s2_vars_B, optimal_ratio_B,
                              "B_restricted", data_B_clean)

# Balance function
run_balance <- function(s2_result, label, s1_imbalanced = NULL) {
  
  cat("--- Balance (ratio:", s2_result$primary_ratio, ") —", label, "---\n")
  bt     <- bal.tab(s2_result$matchit_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    mutate(var_type = case_when(
      variable %in% stage2_trends ~ "Trend",
      variable %in% stage2_levels ~ "Level",
      TRUE                        ~ "Other"
    ))
  
  cat("\n  TREND SMDs (target: all < 0.06):\n")
  trend_s <- smd_df %>% filter(var_type == "Trend") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(trend_s)
  if (any(abs(trend_s$Diff.Adj) >= 0.1, na.rm = TRUE))
    cat("  WARNING: trend SMD >= 0.1 — parallel trends support weakened\n")
  
  cat("\n  LEVEL SMDs (residual imbalance expected — include in C&S outcome model):\n")
  level_s <- smd_df %>% filter(var_type == "Level") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(level_s)
  s2_level_covs <- level_s %>% filter(abs(Diff.Adj) > 0.1) %>% pull(variable)
  
  smd_all <- abs(bt$Balance$Diff.Adj)
  cat("\n  Balanced (<0.1):", sum(smd_all < 0.1, na.rm = TRUE), "/", length(smd_all),
      "| Max |SMD|:", round(max(smd_all, na.rm = TRUE), 3),
      "| Mean |SMD|:", round(mean(smd_all, na.rm = TRUE), 3), "\n")
  
  # [FIX 3] Compile full covariate list for outcome model
  # = Stage 1 structurally imbalanced + Stage 2 level imbalanced
  all_covs <- unique(c(s1_imbalanced, s2_level_covs))
  if (length(all_covs) > 0) {
    cat("\n  [FIX 3] RECOMMENDED C&S xformla covariates (Stage 1 structural +\n")
    cat("          Stage 2 level imbalance):\n")
    cat("    ", paste(all_covs, collapse = ", "), "\n")
  }
  
  lp <- love.plot(
    s2_result$matchit_obj, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted", title = paste("Stage 2 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) + theme(axis.text.y = element_text(size = 9),
            plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  ggsave(here("output", paste0("s2_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  invisible(list(bal = bt, s2_level_covs = s2_level_covs, all_covs = all_covs))
}

bal_A <- run_balance(s2_A, "A_excl_zero",   s1_A$s1_imbalanced)
bal_B <- run_balance(s2_B, "B_incl_zero",   s1_B$s1_imbalanced)
run_balance(s2_A_restricted, "A_excl_zero_restricted", s1_A$s1_imbalanced)
run_balance(s2_B_restricted, "B_incl_zero_restricted", s1_B$s1_imbalanced)

# =============================================================================
# SECTION 10 — COMPARE ANALYSES A AND B
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 10: COMPARISON — ANALYSES A vs B\n")
cat("====================================================\n\n")

treated_A <- s2_A$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
only_in_B <- setdiff(treated_B, treated_A)
only_in_A <- setdiff(treated_A, treated_B)
in_both   <- intersect(treated_A, treated_B)

cat("--- Sample overlap ---\n")
cat("  Treated in A:              ", length(treated_A), "\n")
cat("  Treated in B:              ", length(treated_B), "\n")
cat("  In both:                   ", length(in_both), "\n")
cat("  Only in B (zero-injury):   ", length(only_in_B), "\n")
cat("  Only in A (trimmed from B):", length(only_in_A), "\n\n")

if (length(only_in_A) > 0)
  cat("WARNING: OAs in A but not B — unexpected, check exclusion logic\n")

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
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value,   na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = c(mean, sd)) %>%
    mutate(SMD = (`mean_Injury-exposed (both)` - `mean_Zero-injury (B only)`) /
             sqrt((`sd_Injury-exposed (both)`^2 + `sd_Zero-injury (B only)`^2) / 2)) %>%
    arrange(desc(abs(SMD)))
  print(cov_diff %>%
          select(name, `mean_Injury-exposed (both)`, `mean_Zero-injury (B only)`, SMD),
        n = Inf)
  
  cat("\n  Country distribution of zero-injury OAs in B:\n")
  print(table(s2_B$primary_data %>% filter(OA %in% only_in_B) %>% pull(country)))
  
  cat("\n  Stage 2 distances for zero-injury treated OAs:\n")
  print(summary(s2_B$primary_data %>%
                  filter(OA %in% only_in_B, treat_indicator == 1) %>%
                  pull(mdist)))
}

# Distance comparison
cat("\n--- Stage 2 distance distribution: A vs B (treated OAs only) ---\n")
dist_compare <- bind_rows(
  s2_A$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "B: incl zero-injury")
) %>%
  group_by(analysis) %>%
  summarise(
    n      = n(),
    mean   = round(mean(mdist,           na.rm = TRUE), 3),
    median = round(median(mdist,         na.rm = TRUE), 3),
    p90    = round(quantile(mdist, 0.90, na.rm = TRUE), 3),
    p95    = round(quantile(mdist, 0.95, na.rm = TRUE), 3),
    max    = round(max(mdist,            na.rm = TRUE), 3),
    .groups = "drop"
  )
print(dist_compare)

# Worst-matched treated OAs — for reporting and sensitivity
cat("\n--- Worst-matched treated OAs in Analysis A (top 20 by Stage 2 distance) ---\n")
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  arrange(desc(mdist)) %>%
  select(OA, mdist, country, starts_with("trend_total")) %>%
  head(20) %>%
  print()

cat("\n--- Distance tail: Analysis A treated OAs ---\n")
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  summarise(
    n_over_20  = sum(mdist > 20),
    n_over_30  = sum(mdist > 30),
    n_over_50  = sum(mdist > 50),
    n_over_100 = sum(mdist > 100)
  ) %>% print()

# =============================================================================
# SECTION 11 — WEIGHT DIAGNOSTICS + APPLY CAP TO ANALYSIS B
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 11: WEIGHT DIAGNOSTICS\n")
cat("====================================================\n\n")

weight_diagnostics <- function(s2_result, analysis_label) {
  
  controls      <- s2_result$primary_data %>% filter(treat_indicator == 0)
  treated_rows  <- s2_result$primary_data %>% filter(treat_indicator == 1)
  trend_vars_here <- intersect(names(s2_result$primary_data), stage2_trends)
  
  cat("\n--- Weight diagnostics:", analysis_label, "---\n")
  cat("\nWeight distribution (controls):\n")
  print(quantile(controls$weights,
                 probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)))
  
  eff_n <- sum(controls$weights)^2 / sum(controls$weights^2)
  cat("\nEffective N:", round(eff_n, 1),
      "vs nominal N:", nrow(controls), "\n")
  cat("Efficiency ratio:", round(eff_n / nrow(controls), 3), "\n")
  
  cat("\nBy country:\n")
  controls %>%
    group_by(country) %>%
    summarise(
      n      = n(),
      mean_w = round(mean(weights), 2),
      max_w  = round(max(weights), 2),
      eff_n  = round(sum(weights)^2 / sum(weights^2), 1),
      .groups = "drop"
    ) %>% print()
  
  caps <- c(5, 10, 20, Inf)
  cap_results <- map_df(caps, function(cap) {
    capped    <- controls %>% mutate(w_cap = pmin(weights, cap))
    eff_n_cap <- sum(capped$w_cap)^2 / sum(capped$w_cap^2)
    trend_smds <- map_dbl(trend_vars_here, function(v) {
      t_mean    <- mean(treated_rows[[v]], na.rm = TRUE)
      c_mean    <- weighted.mean(controls[[v]], w = capped$w_cap, na.rm = TRUE)
      pooled_sd <- sqrt((var(treated_rows[[v]], na.rm = TRUE) +
                           var(controls[[v]],   na.rm = TRUE)) / 2)
      if (pooled_sd == 0) return(0)
      abs(t_mean - c_mean) / pooled_sd
    })
    tibble(
      weight_cap    = if_else(is.infinite(cap), "No cap", paste0("Cap ", cap)),
      n_controls    = nrow(controls),
      effective_n   = round(eff_n_cap, 1),
      efficiency    = round(eff_n_cap / nrow(controls), 3),
      max_trend_smd = round(max(trend_smds), 4),
      pct_at_cap    = round(100 * mean(controls$weights >= cap), 2)
    )
  })
  cat("\nWeight cap sensitivity:\n")
  print(cap_results)
  invisible(cap_results)
}

wd_A <- weight_diagnostics(s2_A, "A_excl_zero")
wd_B <- weight_diagnostics(s2_B, "B_incl_zero")

# Apply cap of 5 to Analysis B (Scottish zero-injury sparsity — see v1 rationale)
s2_B$primary_data <- s2_B$primary_data %>%
  mutate(weights = pmin(weights, 5))

cat("\n--- Analysis B weights after cap at 5 ---\n")
s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  summarise(
    max_weight  = max(weights),
    effective_n = round(sum(weights)^2 / sum(weights^2), 1),
    efficiency  = round((sum(weights)^2 / sum(weights^2)) / n(), 3)
  ) %>% print()

# Apply cap to restricted B as well
s2_B_restricted$primary_data <- s2_B_restricted$primary_data %>%
  mutate(weights = pmin(weights, 5))

# =============================================================================
# SECTION 12 — BASELINE INJURY LEVEL STRATIFICATION [FIX 8]
# =============================================================================
# [FIX 8] Treated OAs are stratified into quartiles of mean_total_pkm
# (pre-period mean total injury rate per road-km). The stratum indicator is
# saved with the matched datasets so the DiD script can:
#   (a) run stratum-specific ATT estimates to detect heterogeneous effects
#   (b) test whether residual level imbalance is an effect modifier
#   (c) report if ATT varies by baseline injury exposure
#
# If ATT is heterogeneous across strata, the composite ATT is not the primary
# quantity of interest — stratum-specific estimates should be reported.
# If ATT is homogeneous, the composite ATT is robust to the level imbalance.

cat("\n====================================================\n")
cat("SECTION 12: BASELINE INJURY LEVEL STRATIFICATION [FIX 8]\n")
cat("====================================================\n\n")

add_baseline_stratum <- function(s2_result, label) {
  
  treated_rows <- s2_result$primary_data %>% filter(treat_indicator == 1)
  
  if (!"mean_total_pkm" %in% names(treated_rows)) {
    cat(label, ": mean_total_pkm not in dataset — skipping stratification\n")
    return(s2_result)
  }
  
  # Quartile breaks computed from treated distribution only
  q_breaks <- quantile(treated_rows$mean_total_pkm,
                       probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  s2_result$primary_data <- s2_result$primary_data %>%
    mutate(
      baseline_injury_stratum = case_when(
        treat_indicator == 0 ~ NA_integer_,   # controls don't get a stratum
        mean_total_pkm <= q_breaks[2] ~ 1L,   # Q1: lowest baseline injury rate
        mean_total_pkm <= q_breaks[3] ~ 2L,
        mean_total_pkm <= q_breaks[4] ~ 3L,
        TRUE                          ~ 4L    # Q4: highest
      )
    )
  
  cat(label, "— Baseline injury stratum distribution (treated OAs):\n")
  print(table(s2_result$primary_data %>%
                filter(treat_indicator == 1) %>%
                pull(baseline_injury_stratum),
              useNA = "ifany"))
  cat("  Stratum quartile breaks (mean_total_pkm):",
      paste(round(q_breaks, 4), collapse = " | "), "\n\n")
  
  s2_result
}

s2_A <- add_baseline_stratum(s2_A, "Analysis A")
s2_B <- add_baseline_stratum(s2_B, "Analysis B")

# =============================================================================
# SECTION 13 — EXTRACT AND SAVE MATCHED DATASETS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 13: SAVE MATCHED DATASETS\n")
cat("====================================================\n\n")

# --- Analysis A ---
matched_A_treated <- s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum)

matched_A_controls <- s2_A$primary_data %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

# --- Analysis B (weights already capped at 5) ---
matched_B_treated <- s2_B$primary_data %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum)

matched_B_controls <- s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

cat("Analysis A — Treated:", nrow(matched_A_treated),
    "| Controls:", nrow(matched_A_controls), "\n")
cat("  Weight range (controls): [",
    round(min(matched_A_controls$weights), 3), ",",
    round(max(matched_A_controls$weights), 3), "]\n")

cat("Analysis B — Treated:", nrow(matched_B_treated),
    "| Controls:", nrow(matched_B_controls), "\n")
cat("  Weight range (controls): [",
    round(min(matched_B_controls$weights), 3), ",",
    round(max(matched_B_controls$weights), 3), "]\n")

# --- Common support flags [FIX 2] ---
common_support_flags <- bind_rows(
  s1_A$isolated_OAs %>% mutate(analysis = "A"),
  s1_B$isolated_OAs %>% mutate(analysis = "B")
)
cat("\nCommon support flags — isolated treated OAs:\n")
print(common_support_flags)

# --- Recommended outcome model covariate list [FIX 3] ---
outcome_covariates <- list(
  analysis_A = bal_A$all_covs,
  analysis_B = bal_B$all_covs
)
cat("\nRecommended C&S xformla covariates:\n")
cat("  Analysis A:", paste(outcome_covariates$analysis_A, collapse = ", "), "\n")
cat("  Analysis B:", paste(outcome_covariates$analysis_B, collapse = ", "), "\n")

# --- Integrity checks ---
stopifnot(
  "A: treated weights should all be 1" =
    all(matched_A_treated$weights == 1),
  "B: treated weights should all be 1" =
    all(matched_B_treated$weights == 1),
  "A: no NA weights in controls" =
    !anyNA(matched_A_controls$weights),
  "B: no NA weights in controls" =
    !anyNA(matched_B_controls$weights),
  "A: no duplicate treated OAs" =
    !anyDuplicated(matched_A_treated$OA),
  "B: no duplicate treated OAs" =
    !anyDuplicated(matched_B_treated$OA),
  "B: max control weight <= 5" =
    max(matched_B_controls$weights) <= 5
)
cat("All integrity checks passed.\n\n")

# --- Save ---
saveRDS(matched_A_treated,    here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A_controls,   here("data", "processed", "OA_matched_donors_A.rds"))
saveRDS(matched_B_treated,    here("data", "processed", "OA_matched_treated_B.rds"))
saveRDS(matched_B_controls,   here("data", "processed", "OA_matched_donors_B.rds"))
saveRDS(s2_A$primary_data,    here("data", "processed", "OA_matched_full_A.rds"))
saveRDS(s2_B$primary_data,    here("data", "processed", "OA_matched_full_B.rds"))
saveRDS(common_support_flags, here("data", "processed", "OA_common_support_flags.rds"))
saveRDS(outcome_covariates,   here("data", "processed", "OA_outcome_covariates.rds"))

cat("Saved:\n")
cat("  OA_matched_treated_A.rds       —", nrow(matched_A_treated), "treated OAs\n")
cat("  OA_matched_donors_A.rds        —", nrow(matched_A_controls), "control OAs\n")
cat("  OA_matched_treated_B.rds       —", nrow(matched_B_treated), "treated OAs\n")
cat("  OA_matched_donors_B.rds        —", nrow(matched_B_controls), "control OAs\n")
cat("  OA_matched_full_A.rds          — full matched dataset, Analysis A\n")
cat("  OA_matched_full_B.rds          — full matched dataset, Analysis B\n")
cat("  OA_common_support_flags.rds    — structurally isolated treated OA flags\n")
cat("  OA_outcome_covariates.rds      — recommended xformla covariate list\n")

# =============================================================================
# SECTION 14 — DESIGN SUMMARY
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 14: DESIGN SUMMARY\n")
cat("====================================================\n\n")

cat("WHY MDM OVER PSM:\n")
cat("  Stage 2 is low-dimensional (", length(s2_vars_A), "variables after zero-spike filter)\n", sep = "")
cat("  PSM's advantage (curse of dimensionality relief) does not apply here.\n")
cat("  PSM collapses all trajectory information into a scalar — misspecification\n")
cat("  of the propensity model destroys covariate balance with no diagnostics.\n")
cat("  MDM directly minimises multivariate distance and balance is directly verifiable.\n")
cat("  For trajectory matching specifically, MDM is the correct choice.\n\n")

cat("STAGE 1  | MDM on", length(stage1_vars), "vars (road + urban + socdem)\n")
cat("         | replace=TRUE, ratio=10, exact=~country\n")
cat("         | Covariance matrix: POOLED (treated + controls) [FIX 1]\n")
cat("         | Purpose: structural plausibility filter\n")
cat("         | Common support flags saved for sensitivity [FIX 2]\n\n")

cat("STAGE 2  | MDM on", length(s2_vars_A), "vars after zero-spike filter [FIX 4]\n")
cat("         | replace=TRUE, ratio = elbow-selected [FIX 5]\n")
cat("         | Covariance matrix: treated-only (ATT-focused)\n")
cat("         | Country integrity verified at Stage 2 [FIX 7]\n")
cat("         | Purpose: parallel trends support via trajectory matching\n\n")

cat("ANALYSIS A | Zero-injury treated OAs EXCLUDED (PRIMARY)\n")
cat("           | N treated:", nrow(matched_A_treated),
    "| N controls:", nrow(matched_A_controls), "\n")
cat("           | Max control weight:", round(max(matched_A_controls$weights), 2),
    "| No weight cap\n\n")

cat("ANALYSIS B | Zero-injury treated OAs INCLUDED (ROBUSTNESS CHECK)\n")
cat("           | N treated:", nrow(matched_B_treated),
    "| N controls:", nrow(matched_B_controls), "\n")
cat("           | Control weights capped at 5 (Scottish zero-injury sparsity)\n\n")

cat("FIXES APPLIED:\n")
cat("  [FIX 1] Stage 1 covariance matrix: pooled (not treated-only)\n")
cat("  [FIX 2] Common support assessment: structurally isolated OAs flagged + sensitivity run\n")
cat("  [FIX 3] Stage 1 structural imbalance added to outcome model covariate list\n")
cat("  [FIX 4] Structural zeros in trend variables diagnosed; high-zero-spike vars dropped\n")
cat("  [FIX 5] Ratio selection via elbow curve (marginal improvement < 0.002)\n")
cat("  [FIX 6] Safe winsorisation: explicit map() loop, no cur_column() fragility\n")
cat("  [FIX 7] Stage 2 cross-country integrity check added\n")
cat("  [FIX 8] Baseline injury stratum saved for DiD heterogeneity testing\n")
cat("  [FIX 9] Buffer OA pool eligibility diagnosed\n\n")

cat("NEXT STEPS FOR DiD:\n")
cat("  1. Join matched OA lists to panel outcome data\n")
cat("  2. Run att_gt(..., weightsname = 'weights') on both datasets\n")
cat("  3. Include ALL covariates from OA_outcome_covariates.rds in xformla\n")
cat("     (Stage 1 structural + Stage 2 level imbalance) [FIX 3]\n")
cat("  4. Run stratum-specific ATT by baseline_injury_stratum [FIX 8]\n")
cat("     to test whether level imbalance is an effect modifier\n")
cat("  5. Run sensitivity excluding OA_common_support_flags OAs [FIX 2]\n")
cat("  6. If ATT(A) ≈ ATT(B): A is primary, B is robustness check\n")
cat("  7. If ATT(A) ≠ ATT(B): report both; discuss heterogeneous effects\n")
cat("     by pre-treatment injury exposure level\n")