# =============================================================================
# COMPLETE SCRIPT:
# 1) Runs matching (England + Scotland)
# 2) Saves Stage1->Stage2 pools explicitly
# 3) Builds flow tables: unique and non-unique by country
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

set.seed(222)

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
  "age_45to64_pct", "age_65to84_pct"
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
# BUILD DATASETS — ENGLAND AND SCOTLAND
# =============================================================================
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  )

data_england <- OA_matching_dataset %>%
  filter(
    country == "England",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA == 0,
    n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

data_scotland <- OA_matching_dataset %>%
  filter(
    country == "Scotland",
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# =============================================================================
# PREP FUNCTIONS
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

winsorise_and_log_s1 <- function(data, raw_vars, log_vars, log_nozero_vars = character(0)) {
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
  if (length(missing) > 0) cat("WARNING —", label, "missing:", paste(missing, collapse = ", "), "\n")
  vcheck <- data %>%
    summarise(across(all_of(intersect(vars, names(data))), ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")
  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)
  setdiff(intersect(vars, names(data)), low)
}

data_england_clean  <- data_england  %>% prep_dataset() %>%
  winsorise_and_log_s1(stage1_vars, log_transform_s1, log_nozero_s1) %>%
  winsorise_and_log_s2(stage2_levels, log_transform_s2_levels)

data_scotland_clean <- data_scotland %>% prep_dataset() %>%
  winsorise_and_log_s1(stage1_vars, log_transform_s1, log_nozero_s1) %>%
  winsorise_and_log_s2(stage2_levels, log_transform_s2_levels)

s1_vars_england  <- check_vars(data_england_clean,  stage1_vars_log, "S1 England")
s2_vars_england  <- check_vars(data_england_clean,  stage2_vars_log, "S2 England")
s1_vars_scotland <- check_vars(data_scotland_clean, stage1_vars_log, "S1 Scotland")
s2_vars_scotland <- check_vars(data_scotland_clean, stage2_vars_log, "S2 Scotland")

# =============================================================================
# RATIO SELECTION
# =============================================================================
prepare_s2_for_ratio <- function(data_clean, s1_vars, s2_vars, label) {
  s1v <- check_vars(data_clean, s1_vars, paste("S1 ratio prep", label))
  m_s1 <- matchit(reformulate(s1v, response = "treat_indicator"), data = data_clean,
                  method = "nearest", distance = "mahalanobis", ratio = 10, replace = TRUE)
  mm_s1 <- m_s1$match.matrix
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>% mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>% mutate(treat_indicator = 0L)
  s2_raw <- bind_rows(treated_s1, controls_s1)
  
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  for (v in intersect(stage2_levels, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <- log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }
  for (v in intersect(stage2_trends, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  
  s2v <- check_vars(s2_raw, s2_vars, paste("S2 ratio prep", label))
  list(data = s2_raw, s2_vars = s2v)
}

select_ratio <- function(data, s2_vars, trend_vars, label, ratios_to_test = 1:10) {
  n_t <- sum(data$treat_indicator == 1)
  n_c <- sum(data$treat_indicator == 0)
  ratios_to_test <- ratios_to_test[ratios_to_test <= floor(n_c / n_t)]
  formula_s2 <- reformulate(s2_vars, response = "treat_indicator")
  
  ratio_results <- map_df(ratios_to_test, function(r) {
    m <- tryCatch(matchit(formula_s2, data = data, method = "nearest",
                          distance = "mahalanobis", ratio = r, replace = TRUE),
                  error = function(e) NULL)
    if (is.null(m)) return(tibble(ratio = r, mean_smd = NA, max_trend_smd = NA, max_level_smd = NA))
    bt <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
    trend_rows <- rownames(bt)[rownames(bt) %in% trend_vars]
    level_rows <- rownames(bt)[!rownames(bt) %in% trend_vars]
    tibble(
      ratio         = r,
      mean_smd      = round(mean(abs(bt$Diff.Adj), na.rm = TRUE), 4),
      max_trend_smd = round(max(abs(bt[trend_rows, "Diff.Adj"]), na.rm = TRUE), 4),
      max_level_smd = round(max(abs(bt[level_rows, "Diff.Adj"]), na.rm = TRUE), 4)
    )
  })
  
  best <- ratio_results %>% filter(!is.na(max_trend_smd)) %>%
    arrange(max_trend_smd, mean_smd) %>% slice(1)
  
  list(optimal_ratio = best$ratio, ratio_results = ratio_results, label = label)
}

trend_vars_eng <- intersect(s2_vars_england, stage2_trends)
trend_vars_sco <- intersect(s2_vars_scotland, stage2_trends)

s2_prep_eng <- prepare_s2_for_ratio(data_england_clean, s1_vars_england, s2_vars_england, "England")
s2_prep_sco <- prepare_s2_for_ratio(data_scotland_clean, s1_vars_scotland, s2_vars_scotland, "Scotland")

ratio_eng <- select_ratio(s2_prep_eng$data, s2_prep_eng$s2_vars, trend_vars_eng, "England", 1:10)
ratio_sco <- select_ratio(s2_prep_sco$data, s2_prep_sco$s2_vars, trend_vars_sco, "Scotland", 1:10)

optimal_ratio_england  <- ratio_eng$optimal_ratio
optimal_ratio_scotland <- ratio_sco$optimal_ratio

# =============================================================================
# MATCHING FUNCTION
# =============================================================================
run_matching <- function(data_clean, s1_vars, s2_vars, ratio, label, trend_vars) {
  
  # ---- Stage 1 ---------------------------------------------------------------
  s1v  <- check_vars(data_clean, s1_vars, paste("S1", label))
  m_s1 <- matchit(reformulate(s1v, response = "treat_indicator"), data = data_clean,
                  method = "nearest", distance = "mahalanobis", ratio = 10, replace = TRUE)
  
  mm_s1 <- m_s1$match.matrix
  
  # Unique treated and controls entering Stage 2
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  # Non-unique control slots after Stage 1:
  # each treated OA is matched to up to 10 controls (with replacement)
  # so total slots = non-NA entries in mm_s1, repeats counted
  n_treated_s1       <- nrow(mm_s1)
  n_controls_s1_slots <- sum(!is.na(mm_s1))   # with-replacement slots, not unique OAs
  
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
    s2_raw[[paste0("log1p_", v)]] <- log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }
  
  # ---- Stage 2 ---------------------------------------------------------------
  s2v  <- check_vars(s2_raw, s2_vars, paste("S2", label))
  m_s2 <- matchit(reformulate(s2v, response = "treat_indicator"), data = s2_raw,
                  method = "nearest", distance = "mahalanobis", ratio = ratio, replace = TRUE)
  
  matched_data <- match.data(m_s2)
  mm_s2        <- m_s2$match.matrix
  
  # Non-unique control slots after Stage 2
  n_treated_s2        <- nrow(mm_s2)
  n_controls_s2_slots <- sum(!is.na(mm_s2))   # with-replacement slots
  
  # Pairs
  treated_oas <- s2_raw$OA[as.integer(rownames(mm_s2))]
  pairs <- map_df(seq_len(nrow(mm_s2)), function(i) {
    c_idx <- as.integer(mm_s2[i, ])
    c_idx <- c_idx[!is.na(c_idx)]
    if (length(c_idx) == 0) return(tibble())
    tibble(treated_OA = treated_oas[i], control_OA = s2_raw$OA[c_idx])
  })
  
  list(
    matched_data         = matched_data,
    matchit_s2           = m_s2,
    pairs                = pairs,
    s2_raw               = s2_raw,
    mm_s1                = mm_s1,
    mm_s2                = mm_s2,
    # Pass slot counts through explicitly for clean table building
    n_treated_s1         = n_treated_s1,
    n_controls_s1_unique = length(ctrl_idx),
    n_controls_s1_slots  = n_controls_s1_slots,
    n_treated_s2         = n_treated_s2,
    n_controls_s2_unique = n_distinct(matched_data$OA[matched_data$treat_indicator == 0]),
    n_controls_s2_slots  = n_controls_s2_slots
  )
}

# =============================================================================
# RUN MATCHING
# =============================================================================
result_england <- run_matching(
  data_clean = data_england_clean, s1_vars = s1_vars_england,
  s2_vars = s2_vars_england, ratio = optimal_ratio_england,
  label = "England", trend_vars = trend_vars_eng
)

result_scotland <- run_matching(
  data_clean = data_scotland_clean, s1_vars = s1_vars_scotland,
  s2_vars = s2_vars_scotland, ratio = optimal_ratio_scotland,
  label = "Scotland", trend_vars = trend_vars_sco
)

matched_england  <- result_england$matched_data  %>% mutate(country = "England")
matched_scotland <- result_scotland$matched_data %>% mutate(country = "Scotland")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
saveRDS(matched_england,         here("data", "processed", "OA_matched_full_mixed_England.rds"))
saveRDS(matched_scotland,        here("data", "processed", "OA_matched_full_mixed_Scotland.rds"))
saveRDS(result_england$s2_raw,   here("data", "processed", "OA_stage1to2_pool_England.rds"))
saveRDS(result_scotland$s2_raw,  here("data", "processed", "OA_stage1to2_pool_Scotland.rds"))

pairs_mixed <- bind_rows(
  result_england$pairs  %>% mutate(country = "England"),
  result_scotland$pairs %>% mutate(country = "Scotland")
)
saveRDS(pairs_mixed, here("data", "processed", "OA_matching_pairs_mixed.rds"))

# =============================================================================
# BUILD FLOW TABLES
# =============================================================================
make_tables <- function(before_s1, result, matched, country) {
  
  # ---- UNIQUE: count distinct OA identifiers at each stage ------------------
  unique_tbl <- bind_rows(
    tibble(
      country = country,
      stage   = "Before Stage 1",
      treated = n_distinct(before_s1$OA[before_s1$treat_indicator == 1]),
      control = n_distinct(before_s1$OA[before_s1$treat_indicator == 0])
    ),
    tibble(
      country = country,
      stage   = "After Stage 1 / Before Stage 2",
      treated = result$n_treated_s1,
      control = result$n_controls_s1_unique   # unique OAs passed to Stage 2
    ),
    tibble(
      country = country,
      stage   = "After Stage 2",
      treated = result$n_treated_s2,
      control = result$n_controls_s2_unique   # unique control OAs in final match
    )
  ) %>% mutate(total = treated + control)
  
  # ---- NON-UNIQUE: count match slots (with-replacement repeats counted) -----
  # Before Stage 1: one row per OA, so same as unique
  # After Stage 1:  each treated matched to up to 10 controls — slots counted
  # After Stage 2:  each treated matched to optimal_ratio controls — slots counted
  nonunique_tbl <- bind_rows(
    tibble(
      country = country,
      stage   = "Before Stage 1",
      treated = sum(before_s1$treat_indicator == 1),
      control = sum(before_s1$treat_indicator == 0)
    ),
    tibble(
      country = country,
      stage   = "After Stage 1 / Before Stage 2",
      treated = result$n_treated_s1,
      control = result$n_controls_s1_slots    # slots in mm_s1, repeats counted
    ),
    tibble(
      country = country,
      stage   = "After Stage 2",
      treated = result$n_treated_s2,
      control = result$n_controls_s2_slots    # slots in mm_s2, repeats counted
    )
  ) %>% mutate(total = treated + control)
  
  list(unique = unique_tbl, nonunique = nonunique_tbl)
}

eng_tables <- make_tables(data_england,  result_england,  matched_england,  "England")
sco_tables <- make_tables(data_scotland, result_scotland, matched_scotland, "Scotland")

unique_england     <- eng_tables$unique
nonunique_england  <- eng_tables$nonunique
unique_scotland    <- sco_tables$unique
nonunique_scotland <- sco_tables$nonunique

cat("\n=== UNIQUE | England ===\n");     print(unique_england)
cat("\n=== NON-UNIQUE | England ===\n"); print(nonunique_england)
cat("\n=== UNIQUE | Scotland ===\n");    print(unique_scotland)
cat("\n=== NON-UNIQUE | Scotland ===\n");print(nonunique_scotland)