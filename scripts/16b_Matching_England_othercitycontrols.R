# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a DiD evaluation of CAZ/LEZ
#   interventions using two-stage Mahalanobis Distance Matching.
#
# GEOGRAPHIC SCOPE — ENGLAND ONLY (primary analysis):
#   Scotland is excluded from the matched analysis. All four major Scottish
#   cities (Glasgow, Edinburgh, Aberdeen, Dundee) implemented LEZs, which
#   means the same-city exclusion filter removes every Scottish control OA —
#   leaving 201 Scottish treated OAs with zero eligible domestic comparators.
#   Cross-country matching against English controls was considered but rejected
#   because Scotland and England differ in road safety legislation, injury
#   reporting frameworks, and baseline injury trends in ways that cannot be
#   fully addressed through covariate adjustment. Pooling cross-country pairs
#   would threaten the parallel trends assumption that underpins the DiD
#   design. Scotland is therefore acknowledged as a limitation of the
#   evaluation: near-universal LEZ adoption left no viable counterfactual
#   within Scotland. 
#
# =============================================================================

# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING — ENGLAND ONLY
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a (DiD) analysis using
#   England OAs only (Scotland excluded)
#
#
# OUTPUTS:
#   OA_matched_treated_England.rds       — treated OA IDs + weights + stratum
#   OA_matched_donors_England.rds        — control OA IDs + weights
#   OA_matched_full_England.rds          — full matched dataset
#   OA_common_support_flags_England.rds  — structurally isolated treated OA flags
#   OA_outcome_covariates_England.rds    — xformla covariate list
#   OA_balance_tests_England.rds         — balance improvement test results
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

print(table(OA_matching_dataset$assignment))
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

stage1_road <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road"
)

stage1_urban <- c(
  "dist_citycentre", "pop_density", "area_km2"
)

# Retail only: accommodation/food dropped (r = 0.80 with retail; redundant and
# worsens Mahalanobis stability without improving post-matching balance).
stage1_business <- c("business_retail_per_km2")

stage1_socdem <- c(
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
  "age_under15_pct",
  "age_15to24_pct",
  "age_25to44_pct",
  "age_45to64_pct",
  "age_65plus_pct"
)

stage1_vars <- c(
  stage1_road,
  stage1_urban,
  stage1_business,
  stage1_socdem
)

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

# =============================================================================
# LOG-TRANSFORM VARIABLES
# =============================================================================

log_transform_s1 <- c(
  "road_length_km",
  "pop_density",
  "dist_citycentre",
  "road_density_m_km2",
  "business_retail_per_km2"
)

log_nozero_s1 <- c("area_km2")

log_transform_s2_levels <- stage2_levels

cat("Stage 1 variables:", length(stage1_vars), "\n")
cat("Stage 2 variables:", length(stage2_vars), "\n")

log_names_s1        <- paste0("log1p_", log_transform_s1)
log_nozero_names_s1 <- paste0("log_", log_nozero_s1)
log_names_s2        <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1,
  log_nozero_names_s1,
  setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1))
)

stage2_vars_log <- c(stage2_trends, log_names_s2)

# =============================================================================
# — BUILD DATASET (ENGLAND same city controls ONLY)
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  ) %>%
  filter(country == "England") 

cat("England-only dataset:\n")
print(table(OA_matching_dataset$country))

# Buffer OA diagnosis
buffer_eligible <- OA_matching_dataset %>%
  filter(buffer_OA == 1, n_roads > 0) %>%
  nrow()

cat("--- Buffer OA pool diagnosis ---\n")
cat("  Total buffer OAs: ", sum(OA_matching_dataset$buffer_OA == 1), "\n")
cat("  Buffer OAs with n_roads > 0:", buffer_eligible, "\n")
cat("  Excluded: contamination risk.\n\n")

# =============================================================================
# ANALYSIS England:(othercitycontrols only (controlgroup2)) zero-injury treated OAs excluded upfront
# =============================================================================
OA_matching_dataset %>% select(control_group2_OA) %>% table()

data_England <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group2_OA == 1),
     control_group1_OA == 0,
    buffer_OA == 0,
    n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("=== Analysis England (zero-injury EXCLUDED) ===\n")
cat(
  "  Total OAs:", nrow(data_England),
  "| Treated:", sum(data_England$treat_indicator == 1),
  "| Controls:", sum(data_England$treat_indicator == 0), "\n\n"
)

table(data_England$assignment)

# Collapse age bands
data_England <- data_England %>%
  mutate(
    age_under15_pct = X4under_pct + X5to9_pct + X10to14_pct,
    age_15to24_pct  = X15to19_pct + X20to24_pct,
    age_25to44_pct  = X25to29_pct + X30to34_pct +
      X35to39_pct + X40to44_pct,
    age_45to64_pct  = X45to49_pct + X50to54_pct +
      X55to59_pct + X60to64_pct,
    age_65plus_pct  = X65to69_pct + X70to74_pct +
      X75to79_pct + X80to84_pct
  )

# =============================================================================
# WINSORISE + LOG-TRANSFORM
# =============================================================================

skew_fn <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  mu <- mean(x)
  s <- sd(x)
  
  if (s == 0 || n < 3) return(NA_real_)
  
  sum((x - mu)^3) / (n * s^3)
}

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

data_England_clean <- winsorise_and_log_s1(
  data_England,
  stage1_vars,
  log_transform_s1,
  log_nozero_s1
)

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

data_England_clean <- winsorise_and_log_s2(
  data_England_clean,
  stage2_levels,
  log_transform_s2_levels
)

# =============================================================================
# NEAR-ZERO VARIANCE CHECK
# =============================================================================

check_vars <- function(data, vars, label) {
  
  missing <- setdiff(vars, names(data))
  
  if (length(missing) > 0) {
    cat("WARNING —", label, "missing vars:",
        paste(missing, collapse = ", "), "\n")
  }
  
  vcheck <- data %>%
    summarise(across(
      all_of(intersect(vars, names(data))),
      ~ var(., na.rm = TRUE)
    )) %>%
    pivot_longer(everything(),
                 names_to = "v",
                 values_to = "var")
  
  low <- vcheck %>%
    filter(var < 1e-8) %>%
    pull(v)
  
  if (length(low) > 0) {
    cat("Dropping near-zero variance:",
        paste(low, collapse = ", "), "\n")
  }
  
  setdiff(intersect(vars, names(data)), low)
}

s1_vars_England <- check_vars(
  data_England_clean,
  stage1_vars_log,
  "Stage 1 / England"
)

s2_vars_England_raw <- c(stage2_trends, log_names_s2)

s2_vars_England <- check_vars(
  data_England_clean,
  s2_vars_England_raw,
  "Stage 2 / England"
)

# =============================================================================
# BALANCE TEST FUNCTION
# =============================================================================

balance_test_log <- list()

run_balance_tests <- function(matchit_obj, trend_vars, label) {
  
  bt  <- bal.tab(
    matchit_obj,
    un = TRUE,
    stats = c("mean.diffs", "variance.ratios")
  )
  
  bal <- bt$Balance
  
  smd_un  <- abs(bal$Diff.Un)
  smd_adj <- abs(bal$Diff.Adj)
  
  mean_un  <- mean(smd_un, na.rm = TRUE)
  mean_adj <- mean(smd_adj, na.rm = TRUE)
  
  result <- tibble(
    label         = label,
    mean_smd_un   = round(mean_un, 4),
    mean_smd_adj  = round(mean_adj, 4)
  )
  
  balance_test_log[[label]] <<- result
  
  invisible(result)
}

# =============================================================================
# STAGE 1 MATCHING
# =============================================================================

run_stage1 <- function(data, s1_vars, label) {
  
  cat("--- Stage 1:", label, "---\n")
  
  formula <- reformulate(
    s1_vars,
    response = "treat_indicator"
  )
  
  m <- matchit(
    formula,
    data = data,
    method = "nearest",
    distance = "mahalanobis",
    ratio = 10,
    replace = TRUE
  )
  
  mm <- m$match.matrix
  
  treated_matched <- data[
    as.integer(rownames(mm)),
    ,
    drop = FALSE
  ] %>%
    mutate(treat_indicator = 1L)
  
  control_row_indices <- unique(as.integer(mm[!is.na(mm)]))
  
  controls_matched <- data[
    control_row_indices,
    ,
    drop = FALSE
  ] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated retained:", nrow(treated_matched), "\n")
  cat("  Unique controls:", nrow(controls_matched), "\n")
  
  run_balance_tests(m, trend_vars = character(0),
                    label = paste0("S1_", label))
  
  bt <- bal.tab(
    m,
    thresholds = c(m = 0.1),
    un = TRUE
  )
  
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    dplyr::select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  
  print(smd_df)
  
  s1_imbalanced <- smd_df %>%
    filter(abs(Diff.Adj) > 0.1) %>%
    pull(variable)
  
  list(
    matchit_obj = m,
    treated = treated_matched,
    controls = controls_matched,
    bal = bt,
    s1_imbalanced = s1_imbalanced
  )
}

s1_England <- run_stage1(
  data_England_clean,
  s1_vars_England,
  "England_excl_zero"
)

# =============================================================================
# PREPARE STAGE 2 DATA
# =============================================================================

prepare_s2_data <- function(
    s1_result,
    s2_trend_vars,
    s2_level_vars_raw
) {
  
  s2_raw <- bind_rows(
    s1_result$treated,
    s1_result$controls
  ) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  
  treated_ref <- s2_raw %>%
    filter(treat_indicator == 1)
  
  for (v in intersect(s2_trend_vars, names(s2_raw))) {
    
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    
    s2_raw[[v]] <- pmin(
      pmax(s2_raw[[v]], q_lo),
      q_hi
    )
  }
  
  for (v in intersect(s2_level_vars_raw, names(s2_raw))) {
    
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    
    v_winsor <- pmin(
      pmax(s2_raw[[v]], q_lo),
      q_hi
    )
    
    s2_raw[[paste0("log1p_", v)]] <- log1p(pmax(v_winsor, 0))
  }
  
  s2_raw
}

s2_data_England <- prepare_s2_data(
  s1_England,
  stage2_trends,
  stage2_levels
)

# =============================================================================
# STAGE 2 MATCHING
# =============================================================================

run_stage2 <- function(
    data,
    s2_vars,
    ratio,
    label,
    trend_vars
) {
  
  formula <- reformulate(
    s2_vars,
    response = "treat_indicator"
  )
  
  m <- matchit(
    formula,
    data = data,
    method = "nearest",
    distance = "mahalanobis",
    ratio = ratio,
    replace = TRUE
  )
  
  matched_data <- match.data(m)
  
  run_balance_tests(
    m,
    trend_vars = trend_vars,
    label = label
  )
  
  list(
    matchit_obj = m,
    primary_ratio = ratio,
    primary_data = matched_data
  )
}

# =============================================================================
# RATIO SELECTION
# =============================================================================
# Test ratios 1:1 through 1:10 and select the ratio that minimises the
# maximum trend |SMD| after stage-2 matching. Trend variables drive the
# parallel trends assumption and are therefore the primary balance criterion.

trend_vars_England <- intersect(s2_vars_England, stage2_trends)

cat("\n--- Ratio selection: England ---\n")
cat("  Treated:", sum(s2_data_England$treat_indicator == 1),
    "| Controls:", sum(s2_data_England$treat_indicator == 0), "\n")

formula_s2_ratio <- reformulate(s2_vars_England, response = "treat_indicator")

ratio_results_England <- map_df(1:10, function(r) {
  m <- tryCatch(
    matchit(formula_s2_ratio, data = s2_data_England, method = "nearest",
            distance = "mahalanobis", ratio = r, replace = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(tibble(ratio = r, mean_smd = NA,
                                max_trend_smd = NA, max_level_smd = NA))
  bt         <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
  trend_rows <- rownames(bt)[rownames(bt) %in% trend_vars_England]
  level_rows <- rownames(bt)[!rownames(bt) %in% trend_vars_England]
  tibble(
    ratio         = r,
    mean_smd      = round(mean(abs(bt$Diff.Adj), na.rm = TRUE), 4),
    max_trend_smd = round(max(abs(bt[trend_rows, "Diff.Adj"]), na.rm = TRUE), 4),
    max_level_smd = round(max(abs(bt[level_rows, "Diff.Adj"]), na.rm = TRUE), 4)
  )
})

print(ratio_results_England)

best_ratio_England <- ratio_results_England %>%
  filter(!is.na(max_trend_smd)) %>%
  arrange(max_trend_smd, mean_smd) %>%
  slice(1)

cat(sprintf("  Selected: 1:%d (max trend |SMD| = %.4f)\n\n",
            best_ratio_England$ratio, best_ratio_England$max_trend_smd))

saveRDS(ratio_results_England,
        here("data", "processed", "OA_ratio_selection_England.rds"))

optimal_ratio_England <- best_ratio_England$ratio

# =============================================================================
# STAGE 2 MATCHING
# =============================================================================

s2_England <- run_stage2(
  s2_data_England,
  s2_vars_England,
  optimal_ratio_England,
  "England_excl_zero",
  trend_vars = intersect(s2_vars_England, stage2_trends)
)

# =============================================================================
# BASELINE INJURY STRATA
# =============================================================================

add_baseline_stratum <- function(s2_result, label) {
  
  treated_rows <- s2_result$primary_data %>%
    filter(treat_indicator == 1)
  
  q_breaks <- quantile(
    treated_rows$log1p_mean_total_pkm,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )
  
  s2_result$primary_data <- s2_result$primary_data %>%
    mutate(
      baseline_injury_stratum = case_when(
        treat_indicator == 0 ~ NA_integer_,
        log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
        log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
        log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
        TRUE ~ 4L
      )
    )
  
  s2_result
}

s2_England <- add_baseline_stratum(
  s2_England,
  "Analysis England"
)

# =============================================================================
# EXTRACT + SAVE
# =============================================================================

matched_England_treated <- s2_England$primary_data %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum)

matched_England_controls <- s2_England$primary_data %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

cat(
  "Analysis England — Treated:",
  nrow(matched_England_treated),
  "| Controls:",
  nrow(matched_England_controls),
  "\n"
)

# Covariates
outcome_covariates <- list(
  analysis_England = s1_England$s1_imbalanced
)

balance_tests_summary <- bind_rows(balance_test_log)

# Integrity checks
stopifnot(
  "England: treated weights == 1" =
    all(matched_England_treated$weights == 1),
  
  "England: no NA control weights" =
    !anyNA(matched_England_controls$weights),
  
  "England: no duplicate treated OAs" =
    !anyDuplicated(matched_England_treated$OA)
)

cat("All integrity checks passed.\n\n")

# =============================================================================
# SAVE
# =============================================================================

saveRDS(
  matched_England_treated,
  here("data", "processed",
       "OA_matched_treated_England.rds")
)

saveRDS(
  matched_England_controls,
  here("data", "processed",
       "OA_matched_donors_England.rds")
)

saveRDS(
  s2_England$primary_data,
  here("data", "processed",
       "OA_matched_full_England.rds")
)

saveRDS(
  outcome_covariates,
  here("data", "processed",
       "OA_outcome_covariates_England.rds")
)

saveRDS(
  balance_tests_summary,
  here("data", "processed",
       "OA_balance_tests_England.rds")
)

saveRDS(
  s2_England,
  here("data", "processed",
       "OA_s2_England.rds")
)