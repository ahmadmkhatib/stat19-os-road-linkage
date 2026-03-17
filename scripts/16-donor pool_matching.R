# ============================================================
# OA-Level Mahalanobis Matching for Synthetic DiD
# ============================================================
#
# PURPOSE
# -------
# Match each treated OA to k = 5 control OAs using Mahalanobis
# distance on pre-treatment injury rates, injury trends, road
# characteristics, and census characteristics. Matching is
# restricted within country (England/Wales and Scotland separate)
# to avoid pairing OAs with systematically different institutional
# contexts and deprivation measures.
#
# MATCHING VARIABLES
# ------------------
# Three groups with differential weights in the distance matrix:
#   - Injury baselines & trends (per-km): weight = 3
#     These directly measure the pre-treatment outcome trajectory
#     and are the most important variables for parallel trends.
#   - Road characteristics: weight = 2
#     Strong structural confounders (road type, density).
#   - Census characteristics: weight = 1
#     Background confounders (deprivation, car ownership, travel mode).
#
# DISTANCE METRIC
# ---------------
# Mahalanobis distance with a diagonal weight matrix. The standard
# Mahalanobis distance is scaled by the covariance matrix of the
# covariates; here we additionally scale each variable group by its
# weight to give greater influence to injury-related variables.
# This follows Abadie et al. (2010) in prioritising pre-treatment
# outcome fit over covariate balance in synthetic control matching.
#
# DESIGN
# ------
# - k = 5 controls per treated OA (primary)
# - Sensitivity checks at k = 3 and k = 10
# - Matching with replacement — a control OA can be matched to
#   multiple treated OAs, which is standard in synthetic control
#   designs where the donor pool may be smaller than the treated pool
# - Zero-injury OAs retained and flagged in output
#
# OUTPUT
# ------
# data/processed/OA_matched_donors.rds
#   One row per treated-control pair. Contains:
#   - treated OA identifier and scheme
#   - matched control OA identifier
#   - Mahalanobis distance for the match
#   - zero_injury flags for both treated and control OA
#   - match rank (1 = closest match)
#   - k sensitivity indicator (k3, k5, k10)
#
# data/processed/OA_match_balance.rds
#   Pre- and post-matching covariate balance statistics (SMD)
#   for reporting in methods section.
# ============================================================

library(tidyverse)
library(here)
library(MatchIt)   # install.packages("MatchIt") if needed
library(cobalt)    # install.packages("cobalt")  — balance diagnostics

# ── Load data ─────────────────────────────────────────────────────────────────

OA_matching_census <- readRDS(
  here("data", "processed", "OA_matching_census.rds")
)

cat("Total OAs:  ", nrow(OA_matching_census), "\n")
cat("Treated:    ", sum(OA_matching_census$treated_OA == 1), "\n")
cat("Control:    ", sum(OA_matching_census$control_group2_OA == 1), "\n")

# ── Define matching variables ─────────────────────────────────────────────────
# Three groups assigned different weights in the Mahalanobis distance.
# Variables must be numeric and have no NAs in the matching sample.

# Group 1 — injury baselines per km (weight = 3)
# Pre-treatment mean quarterly casualties per road-km by mode/severity
vars_injury_baseline <- c(
  "mean_car_KSI_pkm",
  "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",
  "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",
  "mean_ped_slight_pkm",
  "mean_total_pkm"
)

# Group 2 — injury trends per km (weight = 3)
# Quasi-Poisson GLM slope (log-rate change per quarter) by mode/severity
vars_injury_trend <- c(
  "trend_car_KSI_pkm",
  "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",
  "trend_ped_slight_pkm",
  "trend_total_pkm"
)

# Group 3 — road characteristics (weight = 2)
vars_road <- c(
  "road_density_m_km2",
  "pct_A_road",
  "pct_B_road",
  "pct_minor_road",
  "road_length_km"
)

# Group 4 — census characteristics (weight = 1)
# Percentages used rather than raw counts to ensure scale comparability
vars_census <- c(
  "IMD",
  "cars_none_pct",
  "cars_one_pct",
  "Drive_Car_pct",
  "Walk_pct",
  "Bicycle_pct",
  "bus_Coach_pct",
  "Train_pct",
  "Underground_train_tram_pct",
  "White_pct",
  "Asian_pct",
  "Black_pct",
  "Mixed_pct",
  "X20to64_pct",    # working age
  "X65plus_pct",    # elderly
  "X4under_pct",    # young children
  "allInWork_pct",
  "workAthome_pct"
)

# All matching variables combined
all_match_vars <- c(
  vars_injury_baseline,
  vars_injury_trend,
  vars_road,
  vars_census
)

# ── Build weight vector ───────────────────────────────────────────────────────
# Weights are applied as multipliers on each variable before computing
# Mahalanobis distance. Higher weight = more influence on match quality.

match_weights <- c(
  rep(3, length(vars_injury_baseline)),   # injury baselines
  rep(3, length(vars_injury_trend)),      # injury trends
  rep(2, length(vars_road)),              # road characteristics
  rep(1, length(vars_census))             # census characteristics
)

names(match_weights) <- all_match_vars
cat("Matching on", length(all_match_vars), "variables\n")

# ── Matching function — runs for one country and one k ───────────────────────
# Encapsulating matching in a function allows clean iteration over
# country subsets and k sensitivity checks without code duplication.

run_matching <- function(data, k, country_label) {
  
  cat("\nMatching:", country_label, "| k =", k, "\n")
  
  # Restrict to treated + control group 2 only
  # Buffer OAs are excluded from both sides
  df <- data %>%
    filter(treated_OA == 1 | control_group2_OA == 1) %>%
    # treatment indicator for MatchIt must be 0/1
    mutate(treatment = as.integer(treated_OA == 1))
  
  cat("  Treated:", sum(df$treatment == 1),
      "| Controls:", sum(df$treatment == 0), "\n")
  
  # Check all matching variables are present and numeric
  missing_vars <- setdiff(all_match_vars, names(df))
  if (length(missing_vars) > 0) {
    cat("  WARNING — missing variables:", paste(missing_vars, collapse = ", "), "\n")
  }
  
  # Scale each variable by its weight before matching
  # Mahalanobis in MatchIt uses the scaled covariance matrix;
  # multiplying variables by their weight increases their contribution
  # to the distance proportionally
  df_scaled <- df %>%
    mutate(across(
      all_of(all_match_vars),
      ~ .x * match_weights[cur_column()]
    ))
  
  # Build formula from all matching variables
  match_formula <- reformulate(all_match_vars, response = "treatment")
  
  # Run Mahalanobis distance matching with replacement
  # ratio = k controls per treated unit
  match_out <- tryCatch(
    matchit(
      formula   = match_formula,
      data      = df_scaled,
      method    = "nearest",
      distance  = "mahalanobis",
      ratio     = k,
      replace   = TRUE     # matching with replacement — standard for synthetic control
    ),
    error = function(e) {
      cat("  ERROR in matchit:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(match_out)) return(NULL)
  
  # Extract matched pairs
  matches <- match.data(match_out) %>%
    select(OA, treatment, subclass, weights, distance) %>%
    # add back original (unscaled) data for output
    left_join(
      data %>% select(OA, scheme, treated_OA, control_group2_OA,
                      zero_injury_OA, country, all_of(all_match_vars)),
      by = "OA"
    )
  
  # Build treated → control pair table
  # MatchIt stores subclass as the treated OA's match group
  treated_oas <- matches %>%
    filter(treatment == 1) %>%
    select(treated_OA_code = OA, subclass, scheme)
  
  control_oas <- matches %>%
    filter(treatment == 0) %>%
    select(control_OA_code = OA, subclass,
           control_zero_injury = zero_injury_OA,
           mahal_distance = distance)
  
  pairs <- treated_oas %>%
    left_join(control_oas, by = "subclass") %>%
    group_by(treated_OA_code) %>%
    mutate(match_rank = row_number()) %>%   # 1 = closest match
    ungroup() %>%
    mutate(
      k             = k,
      country_group = country_label
    )
  
  cat("  Matched pairs:", nrow(pairs), "\n")
  cat("  Mean Mahalanobis distance:",
      round(mean(pairs$mahal_distance, na.rm = TRUE), 3), "\n")
  
  pairs
}

# ── Run matching separately by country ───────────────────────────────────────
# Scotland and England/Wales are matched within-country only to avoid
# pairing OAs with different deprivation indices and administrative contexts

eng_wales <- OA_matching_census %>%
  filter(country == "EnglandWales")

scotland <- OA_matching_census %>%
  filter(country == "Scotland")

cat("England/Wales — treated:", sum(eng_wales$treated_OA == 1),
    "| controls:", sum(eng_wales$control_group2_OA == 1), "\n")
cat("Scotland      — treated:", sum(scotland$treated_OA == 1),
    "| controls:", sum(scotland$control_group2_OA == 1), "\n")

# ── Primary matching: k = 5 ───────────────────────────────────────────────────

matches_eng_k5  <- run_matching(eng_wales, k = 5, "EnglandWales")
matches_scot_k5 <- run_matching(scotland,  k = 5, "Scotland")

OA_matched_k5 <- bind_rows(matches_eng_k5, matches_scot_k5)

# ── Sensitivity: k = 3 and k = 10 ────────────────────────────────────────────

matches_eng_k3  <- run_matching(eng_wales, k = 3,  "EnglandWales")
matches_scot_k3 <- run_matching(scotland,  k = 3,  "Scotland")
OA_matched_k3   <- bind_rows(matches_eng_k3, matches_scot_k3)

matches_eng_k10  <- run_matching(eng_wales, k = 10, "EnglandWales")
matches_scot_k10 <- run_matching(scotland,  k = 10, "Scotland")
OA_matched_k10   <- bind_rows(matches_eng_k10, matches_scot_k10)

# Combine all k variants into one dataset with k as identifier
OA_matched_donors <- bind_rows(
  OA_matched_k3  %>% mutate(k_label = "k3"),
  OA_matched_k5  %>% mutate(k_label = "k5"),
  OA_matched_k10 %>% mutate(k_label = "k10")
)

# ── Post-matching checks ──────────────────────────────────────────────────────

# 1. How many treated OAs were matched per scheme (k=5)
OA_matched_k5 %>%
  group_by(scheme) %>%
  summarise(
    n_treated_matched = n_distinct(treated_OA_code),
    n_pairs           = n(),
    mean_distance     = round(mean(mahal_distance, na.rm = TRUE), 3),
    max_distance      = round(max(mahal_distance,  na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# 2. How many zero-injury control OAs were selected as matches (k=5)
OA_matched_k5 %>%
  count(control_zero_injury) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# 3. Distance distribution — flag poor matches
# Large Mahalanobis distances indicate poor matches
# Threshold of > 3 SD above mean is a common rule of thumb
dist_mean <- mean(OA_matched_k5$mahal_distance, na.rm = TRUE)
dist_sd   <- sd(OA_matched_k5$mahal_distance,   na.rm = TRUE)
poor_matches <- OA_matched_k5 %>%
  filter(mahal_distance > dist_mean + 3 * dist_sd)

cat("Poor matches (distance > mean + 3SD):", nrow(poor_matches), "\n")
cat("Threshold:", round(dist_mean + 3 * dist_sd, 3), "\n")

# Distribution plot
ggplot(OA_matched_k5, aes(x = mahal_distance)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = dist_mean + 3 * dist_sd,
             colour = "red", linetype = "dashed") +
  labs(
    title = "Mahalanobis distance distribution (k = 5)",
    subtitle = "Red dashed line = mean + 3SD threshold for poor matches",
    x = "Mahalanobis distance", y = "Count"
  ) +
  theme_minimal()

# 4. Covariate balance — standardised mean differences before and after matching
# SMD < 0.1 is the standard threshold for acceptable balance (Stuart, 2010)
# Uses cobalt package which handles weighted balance automatically

bal_k5 <- bal.tab(
  reformulate(all_match_vars, response = "treatment"),
  data     = OA_matching_census %>%
    filter(treated_OA == 1 | control_group2_OA == 1) %>%
    mutate(treatment = as.integer(treated_OA == 1)),
  weights  = NULL,
  method   = "matching",
  un       = TRUE    # show pre-matching balance too
)
print(bal_k5)

# SMD plot — visual balance summary
love.plot(
  bal_k5,
  threshold = 0.1,
  abs       = TRUE,
  title     = "Covariate balance before and after matching (k = 5)"
)

# 5. How many unique control OAs used — check for over-reliance on few donors
n_unique_controls <- OA_matched_k5 %>%
  distinct(control_OA_code) %>%
  nrow()

control_usage <- OA_matched_k5 %>%
  count(control_OA_code, name = "times_used") %>%
  arrange(desc(times_used))

cat("Unique control OAs used:", n_unique_controls, "\n")
cat("Max times one control OA used:", max(control_usage$times_used), "\n")
cat("Top 10 most-used control OAs:\n")
print(head(control_usage, 10))

# ── Save ──────────────────────────────────────────────────────────────────────

saveRDS(
  OA_matched_donors,
  here("data", "processed", "OA_matched_donors.rds")
)

# Balance statistics for methods reporting
saveRDS(
  bal_k5,
  here("data", "processed", "OA_match_balance.rds")
)

cat("\nSaved: OA_matched_donors.rds\n")
cat("Saved: OA_match_balance.rds\n")

# ── Variable descriptions ─────────────────────────────────────────────────────

cat("
Output dataset: OA_matched_donors
----------------------------------
treated_OA_code   : OA code of the treated unit
control_OA_code   : OA code of the matched control unit
scheme            : CAZ scheme of the treated OA
match_rank        : 1 = closest match, up to k
mahal_distance    : Mahalanobis distance between treated and control OA
control_zero_injury: 1 = matched control OA had zero injuries in STATS19
k                 : number of controls matched per treated OA
k_label           : k3 / k5 / k10 for sensitivity analysis
country_group     : EnglandWales or Scotland
")