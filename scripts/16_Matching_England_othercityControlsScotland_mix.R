# =============================================================================
#Final matching code  OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a DiD evaluation of CAZ/LEZ
#   interventions using two-stage Mahalanobis Distance Matching.
#
# CONTROL GROUP STRATEGY (asymmetric by country):
#
#   ENGLAND — Other-city controls only (control_group2):
#     OAs from LADs that contain no treated zone. Same-city control OAs
#     (control_group1, from CAZ/LEZ LADs) are excluded to prevent indirect
#     treatment contamination from traffic displacement and CAZ spillover.
#
#   SCOTLAND — Other-city + same-city controls (control_group2 + control_group1):
#     All four major Scottish cities (Glasgow, Edinburgh, Aberdeen, Dundee)
#     implemented LEZs. Applying the same-city exclusion alone would restrict
#     Scotland to the two non-LEZ LADs (S12000014, S12000050). Retaining
#     same-city controls (outside the buffer zone) maximises pool size and
#     match quality for Scotland at the cost of potential downward bias:
#     same-city controls may have experienced indirect LEZ effects (e.g.,
#     suppressed feeder-road traffic), which would attenuate Scottish
#     estimates toward zero.
#     LIMITATION: Scottish results should be interpreted as conservative
#     lower-bound estimates and reported separately from England.
#
#   Compare 16b (other-city only for both countries) against this script to
#   assess sensitivity of Scottish estimates to same-city control inclusion.
#
# OUTPUTS:
#   OA_matched_treated_mixed.rds        — treated OA IDs + weights + stratum
#   OA_matched_donors_mixed.rds         — control OA IDs + weights
#   OA_matched_full_mixed.rds           — combined matched dataset (Eng + Scot)
#   OA_matched_full_mixed_England.rds   — England only
#   OA_matched_full_mixed_Scotland.rds  — Scotland only
#   OA_common_support_flags_mixed.rds   — structurally isolated OA flags
#   OA_ratio_selection_mixed.rds        — ratio selection diagnostics by country
#   OA_balance_tests_mixed.rds          — balance improvement test results
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


# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

stage1_road     <- c("road_density_m_km2", "road_length_km",
                     "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban    <- c("dist_citycentre", "pop_density", "area_km2")
stage1_business <- c("business_retail_per_km2")
stage1_socdem   <- c(
  "IMD",
  "cars_one_pct", "cars_twoPlus_pct",
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
  "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
  "age_45to64_pct", "age_65plus_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_business, stage1_socdem)

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

log_transform_s1        <- c("road_length_km", "pop_density", "dist_citycentre",
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
# BUILD DATASETS — ENGLAND AND SCOTLAND SEPARATELY
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                          ~ "Unknown"
    )
  )

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
cat("\nScottish LEZ LADs (same-city controls RETAINED — see header note):\n")
print(treated_lads_scotland)

# --- ENGLAND: other-city controls only (control_group2) ----------------------
# control_group1 OAs (from treated LADs) are fully excluded to prevent
# contamination from traffic displacement and CAZ spillover effects.
# control_group1_OA == 0 guards against any OA being double-classified.

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

cat("\n=== ENGLAND dataset (other-city controls only) ===\n")
cat("Treated:", sum(data_england$treat_indicator == 1),
    "| Controls:", sum(data_england$treat_indicator == 0), "\n")
cat("Control-to-treated ratio:",
    round(sum(data_england$treat_indicator == 0) /
            sum(data_england$treat_indicator == 1), 1), "\n\n")

eng_samecity_leak <- data_england %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_england) %>%
  nrow()
cat("Same-city control leak (should be 0):", eng_samecity_leak, "\n\n")

# --- SCOTLAND: other-city + same-city controls --------------------------------
# Both control_group1 (same-city, outside buffer) and control_group2
# (other-city, non-LEZ LADs) are retained to maximise pool size and match
# quality. This introduces potential downward bias if same-city controls
# were indirectly affected by the LEZ; results should be read as conservative.
# Compare with 16b (other-city only) to assess sensitivity.

data_scotland <- OA_matching_dataset %>%
  filter(
    country == "Scotland",
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1)
    # same-city exclusion deliberately omitted — see header
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("=== SCOTLAND dataset (other-city + same-city controls) ===\n")
cat("Treated:", sum(data_scotland$treat_indicator == 1),
    "| Controls:", sum(data_scotland$treat_indicator == 0), "\n")
cat("Control-to-treated ratio:",
    round(sum(data_scotland$treat_indicator == 0) /
            sum(data_scotland$treat_indicator == 1), 1), "\n")
cat("LIMITATION: same-city controls retained — potential attenuation bias.\n\n")

scotland_samecity_controls <- data_scotland %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_scotland) %>%
  nrow()
scotland_othercity_controls <- data_scotland %>%
  filter(treat_indicator == 0, !(LAD24CD %in% treated_lads_scotland)) %>%
  nrow()
cat("Scottish same-city controls:", scotland_samecity_controls, "\n")
cat("Scottish other-city controls:", scotland_othercity_controls, "\n\n")

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
      age_65plus_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct
    )
}

data_england  <- prep_dataset(data_england)
data_scotland <- prep_dataset(data_scotland)

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

# Transform each country independently so winsorisation thresholds
# reflect each country's own distribution, not the pooled distribution.
data_england_clean  <- winsorise_and_log_s1(data_england,  stage1_vars,
                                            log_transform_s1, log_nozero_s1)
data_england_clean  <- winsorise_and_log_s2(data_england_clean,  stage2_levels,
                                            log_transform_s2_levels)

data_scotland_clean <- winsorise_and_log_s1(data_scotland, stage1_vars,
                                            log_transform_s1, log_nozero_s1)
data_scotland_clean <- winsorise_and_log_s2(data_scotland_clean, stage2_levels,
                                            log_transform_s2_levels)

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

s1_vars_england  <- check_vars(data_england_clean,  stage1_vars_log, "S1 England")
s2_vars_england  <- check_vars(data_england_clean,  stage2_vars_log, "S2 England")
s1_vars_scotland <- check_vars(data_scotland_clean, stage1_vars_log, "S1 Scotland")
s2_vars_scotland <- check_vars(data_scotland_clean, stage2_vars_log, "S2 Scotland")

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
# RATIO SELECTION FUNCTION
# =============================================================================

prepare_s2_for_ratio <- function(data_clean, s1_vars, s2_vars, label) {
  s1v        <- check_vars(data_clean, s1_vars, paste("S1 ratio prep", label))
  formula_s1 <- reformulate(s1v, response = "treat_indicator")
  m_s1 <- matchit(formula_s1, data = data_clean, method = "nearest",
                  distance = "mahalanobis", ratio = 10, replace = TRUE)
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

# =============================================================================
# RATIO SELECTION — ENGLAND AND SCOTLAND SEPARATELY
# =============================================================================

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("RATIO SELECTION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

s2_prep_eng <- prepare_s2_for_ratio(data_england_clean,  s1_vars_england,
                                    s2_vars_england,  "England")
s2_prep_sco <- prepare_s2_for_ratio(data_scotland_clean, s1_vars_scotland,
                                    s2_vars_scotland, "Scotland")

trend_vars_eng <- intersect(s2_vars_england,  stage2_trends)
trend_vars_sco <- intersect(s2_vars_scotland, stage2_trends)

ratio_eng <- select_ratio(s2_prep_eng$data, s2_prep_eng$s2_vars,
                          trend_vars_eng, "England",  1:10)
# Scotland tested up to 1:10; same-city controls expand the pool so higher
# ratios are feasible, though match quality may plateau earlier.
ratio_sco <- select_ratio(s2_prep_sco$data, s2_prep_sco$s2_vars,
                          trend_vars_sco, "Scotland", 1:10)

optimal_ratio_england  <- ratio_eng$optimal_ratio
optimal_ratio_scotland <- ratio_sco$optimal_ratio

cat("\nOptimal ratio — England:  1:", optimal_ratio_england, "\n")
cat("Optimal ratio — Scotland: 1:", optimal_ratio_scotland, "\n\n")

# Save and plot ratio selection
ratio_combined <- bind_rows(
  ratio_eng$ratio_results %>% mutate(
    country = "England (other-city controls only)"),
  ratio_sco$ratio_results %>% mutate(
    country = "Scotland (other-city + same-city controls)")
)
saveRDS(ratio_combined,
        here("data", "processed", "OA_ratio_selection_mixed.rds"))

p_ratio <- ratio_combined %>%
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
      "England (other-city controls only)"            = "#2E6FAB",
      "Scotland (other-city + same-city controls)"    = "#6B3FA0"
    )
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title    = "Ratio Selection — Maximum Trend |SMD| by Country",
    subtitle = paste0(
      "England: optimal = 1:", optimal_ratio_england,
      " (other-city only)",
      " | Scotland: optimal = 1:", optimal_ratio_scotland,
      " (other-city + same-city)\n",
      "Scottish estimates may be attenuated if same-city controls were",
      " indirectly affected by the LEZ"
    ),
    x       = "Matching ratio (1:k)",
    y       = "Maximum trend |SMD| after matching",
    colour  = NULL,
    caption = paste0(
      "England: control_group2 only (same-city exclusion applied).\n",
      "Scotland: control_group1 + control_group2 (same-city exclusion not applied).\n",
      "Compare with 16b (other-city only for both) to assess sensitivity."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1A2E5A"),
    plot.subtitle   = element_text(colour = "#555555", size = 10),
    plot.caption    = element_text(colour = "#888888", hjust = 0, size = 9),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "fig08_ratio_selection_by_country.png"),
       p_ratio, width = 13, height = 7, dpi = 300, bg = "white")
cat("Saved: fig08_ratio_selection_by_country.png\n\n")

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
            distance = "mahalanobis", ratio = 10, replace = TRUE),
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

  min_dist     <- dist_s1 %>% group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  threshold_95 <- quantile(min_dist$min_dist_s1, 0.95)
  isolated_OAs <- min_dist %>% filter(min_dist_s1 > threshold_95) %>%
    mutate(structurally_isolated = TRUE)

  cat("  Isolated OAs (95th pct):", nrow(isolated_OAs), "/",
      nrow(treated_s1), "\n")

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

  # Extract treated→control OA pairs (used by scripts 17 and 18)
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
       matchit_s2 = m_s2, dist_s2 = dist_s2, country = label,
       pairs = pairs)
}

# =============================================================================
# RUN MATCHING
# =============================================================================

result_england  <- run_matching(
  data_clean  = data_england_clean,
  s1_vars     = s1_vars_england,
  s2_vars     = s2_vars_england,
  ratio       = optimal_ratio_england,
  label       = "England",
  trend_vars  = trend_vars_eng
)

result_scotland <- run_matching(
  data_clean  = data_scotland_clean,
  s1_vars     = s1_vars_scotland,
  s2_vars     = s2_vars_scotland,
  ratio       = optimal_ratio_scotland,
  label       = "Scotland",
  trend_vars  = trend_vars_sco
)

# =============================================================================
# COMBINE
# =============================================================================

matched_england  <- result_england$matched_data  %>% mutate(country = "England")
matched_scotland <- result_scotland$matched_data %>% mutate(country = "Scotland")
matched_full_mixed   <- bind_rows(matched_england, matched_scotland)

isolated_combined <- bind_rows(
  result_england$isolated_OAs  %>% mutate(country = "England"),
  result_scotland$isolated_OAs %>% mutate(
    country = "Scotland",
    note    = "Same-city controls included — potential downward bias"
  )
)

cat("\n=== FINAL COMBINED DATASET ===\n")
cat("England  — Treated:", sum(matched_england$treat_indicator  == 1),
    "| Controls:", sum(matched_england$treat_indicator  == 0),
    "(other-city only)\n")
cat("Scotland — Treated:", sum(matched_scotland$treat_indicator == 1),
    "| Controls:", sum(matched_scotland$treat_indicator == 0),
    "(other-city + same-city)\n")
cat("Total    — Treated:", sum(matched_full_mixed$treat_indicator   == 1),
    "| Controls:", sum(matched_full_mixed$treat_indicator   == 0), "\n\n")
cat("NOTE: Scottish controls include same-city OAs (outside buffer).\n")
cat("Compare Scotland results against 16b (other-city only) for sensitivity.\n")

# =============================================================================
# BASELINE INJURY STRATIFICATION
# =============================================================================

q_breaks <- quantile(
  matched_full_mixed %>%
    filter(treat_indicator == 1) %>%
    pull(log1p_mean_total_pkm),
  probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE
)

matched_full_mixed <- matched_full_mixed %>%
  mutate(baseline_injury_stratum = case_when(
    treat_indicator == 0                ~ NA_integer_,
    log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
    log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
    log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
    TRUE                                ~ 4L
  ))

# =============================================================================
# EXTRACT + INTEGRITY CHECKS + SAVE
# =============================================================================

matched_mixed_treated  <- matched_full_mixed %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum, country)

matched_mixed_controls <- matched_full_mixed %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights, country)

stopifnot(
  "treated weights == 1"     = all(matched_mixed_treated$weights == 1),
  "no NA control weights"    = !anyNA(matched_mixed_controls$weights),
  "no duplicate treated OAs" = !anyDuplicated(matched_mixed_treated$OA)
)
cat("\nAll integrity checks passed.\n")

saveRDS(matched_mixed_treated,
        here("data", "processed", "OA_matched_treated_mixed.rds"))
saveRDS(matched_mixed_controls,
        here("data", "processed", "OA_matched_donors_mixed.rds"))
saveRDS(matched_full_mixed,
        here("data", "processed", "OA_matched_full_mixed.rds"))
saveRDS(matched_england,
        here("data", "processed", "OA_matched_full_mixed_England.rds"))
saveRDS(matched_scotland,
        here("data", "processed", "OA_matched_full_mixed_Scotland.rds"))
saveRDS(isolated_combined,
        here("data", "processed", "OA_common_support_flags_mixed.rds"))
saveRDS(bind_rows(balance_test_log),
        here("data", "processed", "OA_balance_tests_mixed.rds"))

pairs_mixed <- bind_rows(
  result_england$pairs  %>% mutate(country = "England"),
  result_scotland$pairs %>% mutate(country = "Scotland")
)
saveRDS(pairs_mixed,
        here("data", "processed", "OA_matching_pairs_mixed.rds"))

cat("\n=== OUTPUTS SAVED ===\n")
cat("  OA_matched_full_mixed.rds          — combined England + Scotland\n")
cat("  OA_matched_full_mixed_England.rds  — England (other-city controls)\n")
cat("  OA_matched_full_mixed_Scotland.rds — Scotland (other-city + same-city)\n")
cat("  OA_matched_treated_mixed.rds\n")
cat("  OA_matched_donors_mixed.rds\n")
cat("  OA_common_support_flags_mixed.rds\n")
cat("  OA_ratio_selection_mixed.rds\n")
cat("  OA_balance_tests_mixed.rds\n")
cat("  OA_matching_pairs_mixed.rds\n")
cat("  fig08_ratio_selection_by_country.png\n")
