# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING — ENGLAND + SCOTLAND
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a DiD evaluation of CAZ/LEZ
#   interventions using two-stage Mahalanobis Distance Matching.
#
# CONTROL GROUP STRATEGY — OTHER-CITY CONTROLS ONLY:
#   Only control_group2 OAs are used as controls (other-city, i.e. OAs from
#   LADs that do not contain any treated zone). Same-city control OAs
#   (control_group1, from treated LADs) are excluded in both England and
#   Scotland to eliminate the risk of indirect treatment contamination from
#   traffic displacement and CAZ/LEZ spillover effects.
#
#   England: Other-city English controls are plentiful; the exclusion is
#     straightforward.
#
#   Scotland: Two non-LEZ Scottish LADs (S12000014, S12000050) provide
#     ~4,200 eligible other-city controls against 364 treated OAs — a
#     healthy pool. This supersedes the earlier constraint (16c) where
#     near-universal LEZ adoption left no Scottish other-city controls.
#     The updated dataset now includes those non-LEZ Scottish LADs, making
#     a symmetric same-city-exclusion strategy viable for both countries.
#
# OUTPUTS:
#   OA_matched_treated_othercity.rds       — treated OA IDs + weights + stratum
#   OA_matched_donors_othercity.rds        — control OA IDs + weights
#   OA_matched_full_othercity.rds          — combined England + Scotland
#   OA_matched_full_othercity_England.rds  — England only
#   OA_matched_full_othercity_Scotland.rds — Scotland only
#   OA_ratio_selection_othercity.rds       — ratio diagnostics by country
#   OA_balance_tests_othercity.rds         — balance improvement test results
#   OA_outcome_covariates_othercity.rds    — xformla covariate list
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

select <- dplyr::select
filter <- dplyr::filter

dir.create(here("output", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics")

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
  "trend_car_KSI_pkm",   "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",   "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",   "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm",
  "trend_total_pkm"
)

stage2_levels <- c(
  "mean_car_KSI_pkm",   "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",   "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",   "mean_ped_slight_pkm",
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
log_nozero_names_s1 <- paste0("log_",   log_nozero_s1)
log_names_s2        <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1,
  log_nozero_names_s1,
  setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1))
)

stage2_vars_log <- c(stage2_trends, log_names_s2)

# =============================================================================
# BUILD DATASETS — ENGLAND AND SCOTLAND (other-city controls only)
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

# Treated LADs — used to verify same-city exclusion is working as expected
treated_lads_england <- OA_matching_dataset %>%
  filter(treated_OA == 1, country == "England") %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

treated_lads_scotland <- OA_matching_dataset %>%
  filter(treated_OA == 1, country == "Scotland") %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

cat("English CAZ LADs (same-city controls excluded):\n")
print(treated_lads_england)
cat("\nScottish LEZ LADs (same-city controls excluded):\n")
print(treated_lads_scotland)

# --- ENGLAND: other-city controls only ---------------------------------------
# control_group2_OA are OAs from non-treated LADs. control_group1_OA == 0
# guards against any OA being double-classified; treated OAs pass through.

data_England <- OA_matching_dataset %>%
  filter(
    country          == "England",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA         == 0,
    n_roads            > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("\n=== ENGLAND dataset (other-city controls only) ===\n")
cat("  Treated:", sum(data_England$treat_indicator == 1),
    "| Controls:", sum(data_England$treat_indicator == 0), "\n\n")

# Confirm no same-city controls leaked through
eng_samecity_leak <- data_England %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_england) %>%
  nrow()
cat("Same-city control leak (should be 0):", eng_samecity_leak, "\n\n")

# --- SCOTLAND: other-city controls only --------------------------------------
# S12000014 and S12000050 are non-LEZ Scottish LADs whose OAs are classified
# as control_group2. Applying the same filter as England yields a clean
# other-city donor pool (~4,200 OAs) against 364 treated OAs.

data_Scotland <- OA_matching_dataset %>%
  filter(
    country          == "Scotland",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA         == 0,
    n_roads            > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("=== SCOTLAND dataset (other-city controls only) ===\n")
cat("  Treated:", sum(data_Scotland$treat_indicator == 1),
    "| Controls:", sum(data_Scotland$treat_indicator == 0), "\n")
cat("  Control LADs used:", paste(
    setdiff(unique(data_Scotland$LAD24CD[data_Scotland$treat_indicator == 0]),
            treated_lads_scotland),
    collapse = ", "), "\n\n")

sco_samecity_leak <- data_Scotland %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_scotland) %>%
  nrow()
cat("Same-city control leak (should be 0):", sco_samecity_leak, "\n\n")

table(data_England$assignment)
table(data_Scotland$assignment)

# =============================================================================
# AGE BAND AGGREGATION
# =============================================================================

collapse_age_bands <- function(data) {
  data %>%
    mutate(
      age_under15_pct = X4under_pct + X5to9_pct + X10to14_pct,
      age_15to24_pct  = X15to19_pct + X20to24_pct,
      age_25to44_pct  = X25to29_pct + X30to34_pct + X35to39_pct + X40to44_pct,
      age_45to64_pct  = X45to49_pct + X50to54_pct + X55to59_pct + X60to64_pct,
      age_65plus_pct  = X65to69_pct + X70to74_pct + X75to79_pct + X80to84_pct
    )
}

data_England  <- collapse_age_bands(data_England)
data_Scotland <- collapse_age_bands(data_Scotland)

# =============================================================================
# WINSORISE + LOG-TRANSFORM
# =============================================================================
# Each country is transformed independently so winsorisation thresholds
# reflect its own distribution, not the pooled one.

winsorise_and_log_s1 <- function(data, raw_vars, log_vars,
                                 log_nozero_vars = character(0)) {
  for (v in intersect(raw_vars, names(data))) {
    q         <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
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
    q         <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  data
}

data_England_clean <- winsorise_and_log_s1(data_England,  stage1_vars,
                                           log_transform_s1, log_nozero_s1)
data_England_clean <- winsorise_and_log_s2(data_England_clean, stage2_levels,
                                           log_transform_s2_levels)

data_Scotland_clean <- winsorise_and_log_s1(data_Scotland, stage1_vars,
                                            log_transform_s1, log_nozero_s1)
data_Scotland_clean <- winsorise_and_log_s2(data_Scotland_clean, stage2_levels,
                                            log_transform_s2_levels)

# =============================================================================
# NEAR-ZERO VARIANCE CHECK
# =============================================================================

check_vars <- function(data, vars, label) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0)
    cat("WARNING —", label, "missing vars:",
        paste(missing, collapse = ", "), "\n")

  vcheck <- data %>%
    summarise(across(
      all_of(intersect(vars, names(data))),
      ~ var(., na.rm = TRUE)
    )) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")

  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)

  if (length(low) > 0)
    cat("Dropping near-zero variance (", label, "):",
        paste(low, collapse = ", "), "\n")

  setdiff(intersect(vars, names(data)), low)
}

s1_vars_England  <- check_vars(data_England_clean,  stage1_vars_log, "S1 England")
s2_vars_England  <- check_vars(data_England_clean,  stage2_vars_log, "S2 England")
s1_vars_Scotland <- check_vars(data_Scotland_clean, stage1_vars_log, "S1 Scotland")
s2_vars_Scotland <- check_vars(data_Scotland_clean, stage2_vars_log, "S2 Scotland")

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
              mean_un, mean_adj, if (test_a) "PASS ✓" else "FAIL ✗"))

  trend_in_bal  <- intersect(trend_vars, rownames(bal))
  max_trend_smd <- if (length(trend_in_bal) > 0)
    max(abs(bal[trend_in_bal, "Diff.Adj"]), na.rm = TRUE) else NA_real_
  test_b <- !is.na(max_trend_smd) && max_trend_smd < 0.1
  if (!is.na(max_trend_smd))
    cat(sprintf("  [TEST b] Max trend |SMD|: %.4f  %s\n",
                max_trend_smd, if (test_b) "PASS ✓" else "FAIL ✗"))

  vr_col  <- if ("Var.Ratio.Adj" %in% names(bal)) "Var.Ratio.Adj" else NULL
  vr_fail <- character(0)
  if (!is.null(vr_col)) {
    vr      <- bal[[vr_col]]
    vr_fail <- rownames(bal)[!is.na(vr) & (vr < 0.5 | vr > 2.0)]
    test_c  <- length(vr_fail) == 0
    cat(sprintf("  [TEST c] Variance ratio [0.5, 2.0]: %d/%d pass  %s\n",
                sum(is.na(vr) | (vr >= 0.5 & vr <= 2.0), na.rm = TRUE),
                sum(!is.na(vr)),
                if (test_c) "PASS ✓" else "FAIL ✗"))
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
# STAGE 1 MATCHING
# =============================================================================

run_stage1 <- function(data, s1_vars, label) {
  cat("\n--- Stage 1:", label, "---\n")

  formula <- reformulate(s1_vars, response = "treat_indicator")

  m <- matchit(
    formula, data = data,
    method = "nearest", distance = "mahalanobis",
    ratio = 10, replace = TRUE
  )

  mm               <- m$match.matrix
  treated_matched  <- data[as.integer(rownames(mm)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx         <- unique(as.integer(mm[!is.na(mm)]))
  controls_matched <- data[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)

  cat("  Treated retained:", nrow(treated_matched), "\n")
  cat("  Unique controls:", nrow(controls_matched), "\n")

  run_balance_tests(m, trend_vars = character(0),
                    label = paste0("S1_", label))

  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    dplyr::select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  print(smd_df)

  s1_imbalanced <- smd_df %>%
    filter(abs(Diff.Adj) > 0.1) %>%
    pull(variable)

  list(
    matchit_obj   = m,
    treated       = treated_matched,
    controls      = controls_matched,
    bal           = bt,
    s1_imbalanced = s1_imbalanced
  )
}

s1_England  <- run_stage1(data_England_clean,  s1_vars_England,  "England")
s1_Scotland <- run_stage1(data_Scotland_clean, s1_vars_Scotland, "Scotland")

# =============================================================================
# PREPARE STAGE 2 DATA
# =============================================================================

prepare_s2_data <- function(s1_result, s2_trend_vars, s2_level_vars_raw) {
  s2_raw <- bind_rows(s1_result$treated, s1_result$controls) %>%
    select(-any_of(c("weights", "subclass", "distance")))

  treated_ref <- s2_raw %>% filter(treat_indicator == 1)

  for (v in intersect(s2_trend_vars, names(s2_raw))) {
    q_lo        <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi        <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }

  for (v in intersect(s2_level_vars_raw, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <-
      log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }

  s2_raw
}

s2_data_England  <- prepare_s2_data(s1_England,  stage2_trends, stage2_levels)
s2_data_Scotland <- prepare_s2_data(s1_Scotland, stage2_trends, stage2_levels)

# =============================================================================
# RATIO SELECTION
# =============================================================================
# Select the ratio that minimises maximum trend |SMD| after stage-2 matching.
# Trend variables drive the parallel trends assumption and are the primary
# balance criterion.

trend_vars_England  <- intersect(s2_vars_England,  stage2_trends)
trend_vars_Scotland <- intersect(s2_vars_Scotland, stage2_trends)

select_ratio <- function(data, s2_vars, trend_vars, label,
                         ratios_to_test = 1:10) {
  n_t <- sum(data$treat_indicator == 1)
  n_c <- sum(data$treat_indicator == 0)
  ratios_to_test <- ratios_to_test[ratios_to_test <= floor(n_c / n_t)]

  cat("\n--- Ratio selection:", label, "---\n")
  cat("  Treated:", n_t, "| Controls:", n_c, "\n")

  formula_s2 <- reformulate(s2_vars, response = "treat_indicator")

  ratio_results <- map_df(ratios_to_test, function(r) {
    m <- tryCatch(
      matchit(formula_s2, data = data, method = "nearest",
              distance = "mahalanobis", ratio = r, replace = TRUE),
      error = function(e) NULL
    )
    if (is.null(m)) return(tibble(ratio = r, mean_smd = NA,
                                  max_trend_smd = NA, max_level_smd = NA))
    bt          <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
    trend_rows  <- rownames(bt)[rownames(bt) %in% trend_vars]
    level_rows  <- rownames(bt)[!rownames(bt) %in% trend_vars]
    tibble(
      ratio         = r,
      mean_smd      = round(mean(abs(bt$Diff.Adj), na.rm = TRUE), 4),
      max_trend_smd = round(max(abs(bt[trend_rows, "Diff.Adj"]),
                                na.rm = TRUE), 4),
      max_level_smd = round(max(abs(bt[level_rows, "Diff.Adj"]),
                                na.rm = TRUE), 4)
    )
  })

  print(ratio_results)

  best <- ratio_results %>%
    filter(!is.na(max_trend_smd)) %>%
    arrange(max_trend_smd, mean_smd) %>%
    slice(1)

  cat(sprintf("  Selected: 1:%d (max trend |SMD| = %.4f)\n",
              best$ratio, best$max_trend_smd))

  list(optimal_ratio = best$ratio, ratio_results = ratio_results, label = label)
}

ratio_eng <- select_ratio(s2_data_England,  s2_vars_England,
                          trend_vars_England,  "England",  1:10)
ratio_sco <- select_ratio(s2_data_Scotland, s2_vars_Scotland,
                          trend_vars_Scotland, "Scotland", 1:10)

optimal_ratio_England  <- ratio_eng$optimal_ratio
optimal_ratio_Scotland <- ratio_sco$optimal_ratio

cat(sprintf("\nOptimal ratio — England:  1:%d\n", optimal_ratio_England))
cat(sprintf("Optimal ratio — Scotland: 1:%d\n\n", optimal_ratio_Scotland))

# Save and plot ratio selection
ratio_combined_othercity <- bind_rows(
  ratio_eng$ratio_results %>% mutate(country = "England (other-city controls)"),
  ratio_sco$ratio_results %>% mutate(country = "Scotland (other-city controls)")
)
saveRDS(ratio_combined_othercity,
        here("data", "processed", "OA_ratio_selection_othercity.rds"))

p_ratio <- ratio_combined_othercity %>%
  filter(!is.na(max_trend_smd)) %>%
  ggplot(aes(x = ratio, y = max_trend_smd,
             colour = country, group = country)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#888888", linewidth = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             colour = "#555555", linewidth = 0.5) +
  scale_colour_manual(
    values = c(
      "England (other-city controls)"  = "#2E6FAB",
      "Scotland (other-city controls)" = "#C0392B"
    )
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title    = "Ratio Selection — Maximum Trend |SMD| by Country",
    subtitle = paste0(
      "England: optimal = 1:", optimal_ratio_England,
      " | Scotland: optimal = 1:", optimal_ratio_Scotland,
      "\nBoth countries use other-city controls only (same-city exclusion applied)"
    ),
    x      = "Matching ratio (1:k)",
    y      = "Maximum trend |SMD| after matching",
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1A2E5A"),
    plot.subtitle   = element_text(colour = "#555555", size = 10),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "fig08_ratio_selection_othercity.png"),
       p_ratio, width = 13, height = 7, dpi = 300, bg = "white")
cat("Saved: fig08_ratio_selection_othercity.png\n\n")

# =============================================================================
# STAGE 2 MATCHING
# =============================================================================

run_stage2 <- function(data, s2_vars, ratio, label, trend_vars) {
  cat("\n--- Stage 2:", label, "---\n")

  formula <- reformulate(s2_vars, response = "treat_indicator")

  m <- matchit(
    formula, data = data,
    method = "nearest", distance = "mahalanobis",
    ratio = ratio, replace = TRUE
  )

  matched_data <- match.data(m)

  cat("  Treated:", sum(matched_data$treat_indicator == 1),
      "| Controls:", sum(matched_data$treat_indicator == 0), "\n")

  run_balance_tests(m, trend_vars = trend_vars, label = label)

  list(
    matchit_obj   = m,
    primary_ratio = ratio,
    primary_data  = matched_data
  )
}

s2_England  <- run_stage2(s2_data_England,  s2_vars_England,
                          optimal_ratio_England,  "England",  trend_vars_England)
s2_Scotland <- run_stage2(s2_data_Scotland, s2_vars_Scotland,
                          optimal_ratio_Scotland, "Scotland", trend_vars_Scotland)

# =============================================================================
# BASELINE INJURY STRATA
# =============================================================================

add_baseline_stratum <- function(s2_result) {
  treated_rows <- s2_result$primary_data %>% filter(treat_indicator == 1)

  q_breaks <- quantile(
    treated_rows$log1p_mean_total_pkm,
    probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE
  )

  s2_result$primary_data <- s2_result$primary_data %>%
    mutate(
      baseline_injury_stratum = case_when(
        treat_indicator == 0                ~ NA_integer_,
        log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
        log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
        log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
        TRUE                                ~ 4L
      )
    )

  s2_result
}

s2_England  <- add_baseline_stratum(s2_England)
s2_Scotland <- add_baseline_stratum(s2_Scotland)

# =============================================================================
# COMBINE + SUMMARY
# =============================================================================

matched_England  <- s2_England$primary_data  %>% mutate(country = "England")
matched_Scotland <- s2_Scotland$primary_data %>% mutate(country = "Scotland")
matched_full_othercity <- bind_rows(matched_England, matched_Scotland)

cat("\n=== FINAL COMBINED DATASET (other-city controls) ===\n")
cat("England  — Treated:", sum(matched_England$treat_indicator  == 1),
    "| Controls:", sum(matched_England$treat_indicator  == 0), "\n")
cat("Scotland — Treated:", sum(matched_Scotland$treat_indicator == 1),
    "| Controls:", sum(matched_Scotland$treat_indicator == 0), "\n")
cat("Total    — Treated:", sum(matched_full_othercity$treat_indicator == 1),
    "| Controls:", sum(matched_full_othercity$treat_indicator == 0), "\n\n")

# =============================================================================
# EXTRACT TREATED AND CONTROL TABLES
# =============================================================================

matched_othercity_treated <- matched_full_othercity %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum, country)

matched_othercity_controls <- matched_full_othercity %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights, country)

# Covariates from Stage 1 imbalanced variables
outcome_covariates_othercity <- list(
  england  = s1_England$s1_imbalanced,
  scotland = s1_Scotland$s1_imbalanced
)

balance_tests_summary <- bind_rows(balance_test_log)

# =============================================================================
# INTEGRITY CHECKS
# =============================================================================

stopifnot(
  "treated weights == 1"     = all(matched_othercity_treated$weights  == 1),
  "no NA control weights"    = !anyNA(matched_othercity_controls$weights),
  "no duplicate treated OAs" = !anyDuplicated(matched_othercity_treated$OA)
)

cat("All integrity checks passed.\n\n")

# =============================================================================
# SAVE
# =============================================================================

saveRDS(matched_othercity_treated,
        here("data", "processed", "OA_matched_treated_othercity.rds"))

saveRDS(matched_othercity_controls,
        here("data", "processed", "OA_matched_donors_othercity.rds"))

saveRDS(matched_full_othercity,
        here("data", "processed", "OA_matched_full_othercity.rds"))

saveRDS(matched_England,
        here("data", "processed", "OA_matched_full_othercity_England.rds"))

saveRDS(matched_Scotland,
        here("data", "processed", "OA_matched_full_othercity_Scotland.rds"))

saveRDS(outcome_covariates_othercity,
        here("data", "processed", "OA_outcome_covariates_othercity.rds"))

saveRDS(balance_tests_summary,
        here("data", "processed", "OA_balance_tests_othercity.rds"))

cat("=== OUTPUTS SAVED ===\n")
cat("  OA_matched_full_othercity.rds          — combined England + Scotland\n")
cat("  OA_matched_full_othercity_England.rds  — England only\n")
cat("  OA_matched_full_othercity_Scotland.rds — Scotland only\n")
cat("  OA_matched_treated_othercity.rds\n")
cat("  OA_matched_donors_othercity.rds\n")
cat("  OA_ratio_selection_othercity.rds\n")
cat("  OA_balance_tests_othercity.rds\n")
cat("  OA_outcome_covariates_othercity.rds\n")
cat("  fig08_ratio_selection_othercity.png\n")
