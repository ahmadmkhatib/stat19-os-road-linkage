# ============================================================
# Merge Pooled OA Matching Dataset with Census Characteristics
# England-only version
# ============================================================
#
# Inputs:
#   data/processed/OA_matching_data_pooled.rds
#   data/processed/outputArea_raw.csv
#   data/processed/OA_EW_businesses.csv
#   data/processed/shp_files/OA_subset.shp
#
# Outputs:
#   data/processed/OA_matching_census.rds
#   data/processed/shp_files/OA_matching_census.gpkg
#
# Notes:
#   - Uses the new pooled OA-level matching data.
#   - Uses total-injury pre-COVID variables only.
#   - No severity/mode injury variables.
#   - No per-km injury QA.
#   - England-only: keeps OAs with LAD24CD beginning with "E".
# ============================================================

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(naniar)

select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

# ============================================================
# 0. PATHS
# ============================================================

matching_path_main <- here("data", "processed", "OA_matching_data_pooled.rds")
matching_path_old  <- here("data", "processed", "OA_matching_dataset.rds")

out_rds  <- here("data", "processed", "OA_matching_census.rds")
out_gpkg <- here("data", "processed", "shp_files", "OA_matching_census.gpkg")

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD OA MATCHING DATASET
# ============================================================

if (file.exists(matching_path_main)) {
  OA_matching_dataset <- readRDS(matching_path_main)
  cat("Loaded:", matching_path_main, "\n")
} else {
  OA_matching_dataset <- readRDS(matching_path_old)
  cat("Loaded fallback:", matching_path_old, "\n")
}

# Standardise old distance variable name if needed
if ("dist_os_centre" %in% names(OA_matching_dataset) &&
    !"dist_BUA_centroid" %in% names(OA_matching_dataset)) {
  names(OA_matching_dataset)[names(OA_matching_dataset) == "dist_os_centre"] <-
    "dist_BUA_centroid"
}

# England-only filter
if (!"LAD24CD" %in% names(OA_matching_dataset)) {
  stop("LAD24CD is missing from OA_matching_dataset. Cannot apply England-only filter.")
}

OA_matching_dataset <- OA_matching_dataset %>%
  filter(substr(LAD24CD, 1, 1) == "E")

cat("\n--- OA matching dataset after England-only filter ---\n")
cat("Rows:", nrow(OA_matching_dataset), "\n")
cat("Unique OAs:", n_distinct(OA_matching_dataset$OA), "\n")

stopifnot(nrow(OA_matching_dataset) == n_distinct(OA_matching_dataset$OA))

# ============================================================
# 2. LOAD CENSUS OA DATA
# ============================================================

OA_char_raw <- read.csv(
  here("data", "processed", "outputArea_raw.csv")
)

cat("\n--- Raw census data missingness ---\n")
naniar::miss_var_summary(OA_char_raw) %>%
  print(n = Inf)

# ------------------------------------------------------------
# Correctly calculate percentage variables from raw counts
# ------------------------------------------------------------

eth_cols <- c(
  "White", "Mixed", "Asian", "Black", "Other_ethnicity"
)

age_fine <- c(
  "X4under", "X5to9", "X10to14", "X15to19", "X20to24",
  "X25to29", "X30to34", "X35to39", "X40to44", "X45to49",
  "X50to54", "X55to59", "X60to64", "X65to69", "X70to74",
  "X75to79", "X80to84", "X85plus"
)

age_agg <- c(
  "X5to19", "X20to64", "X65plus"
)

cars_cols <- c(
  "cars_none", "cars_one", "cars_two", "cars_threePlus"
)

travel_cols <- c(
  "workAthome",
  "Underground_train_tram",
  "Train",
  "bus_Coach",
  "Taxi",
  "Motorcycle",
  "Drive_Car",
  "Passenger_Car",
  "Bicycle",
  "Walk",
  "Other"
)

required_census_cols <- c(
  "OA", "Total", "totalHouseholds",
  eth_cols, age_fine, age_agg, cars_cols,
  "cars_twoPlus",
  travel_cols
)

missing_census_cols <- setdiff(required_census_cols, names(OA_char_raw))

if (length(missing_census_cols) > 0) {
  stop(
    "These required census columns are missing: ",
    paste(missing_census_cols, collapse = ", ")
  )
}

OA_char_fixed <- OA_char_raw %>%
  mutate(
    travel_base_n = rowSums(across(all_of(travel_cols)), na.rm = TRUE),
    travel_base_zero_flag = as.integer(travel_base_n == 0)
  ) %>%
  mutate(
    across(
      all_of(eth_cols),
      ~ .x / Total * 100,
      .names = "{.col}_pct"
    )
  ) %>%
  mutate(
    across(
      all_of(c(age_fine, age_agg)),
      ~ .x / Total * 100,
      .names = "{.col}_pct"
    )
  ) %>%
  mutate(
    across(
      all_of(cars_cols),
      ~ .x / totalHouseholds * 100,
      .names = "{.col}_pct"
    )
  ) %>%
  mutate(
    cars_twoPlus_pct = cars_twoPlus / totalHouseholds * 100
  ) %>%
  mutate(
    across(
      all_of(travel_cols),
      ~ if_else(travel_base_n == 0, NA_real_, .x / travel_base_n * 100),
      .names = "{.col}_pct"
    )
  )

# Sanity check: travel-mode percentages should sum to 100
travel_pct_cols <- paste0(travel_cols, "_pct")

travel_pct_check <- OA_char_fixed %>%
  filter(travel_base_zero_flag == 0) %>%
  transmute(
    travel_mode_sum = rowSums(across(all_of(travel_pct_cols)), na.rm = TRUE)
  )

stopifnot(all(abs(travel_pct_check$travel_mode_sum - 100) < 1e-6))

cat("\nTravel-mode percentage check passed.\n")

# ------------------------------------------------------------
# Rename raw count variables with _n suffix
# ------------------------------------------------------------

vars_to_rename <- setdiff(
  names(OA_char_raw),
  c("OA", "country", "Total", "IMD")
)

OA_char_raw_renamed <- OA_char_raw %>%
  rename_with(
    ~ paste0(.x, "_n"),
    all_of(vars_to_rename)
  )

OA_char_pct_fixed <- OA_char_fixed %>%
  select(
    OA,
    ends_with("_pct"),
    travel_base_n,
    travel_base_zero_flag
  )

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_fixed, by = "OA")

# Avoid IMD.x / IMD.y if matching dataset already contains IMD
if ("IMD" %in% names(OA_matching_dataset) && "IMD" %in% names(OA_census)) {
  OA_census <- OA_census %>%
    select(-IMD)
}

# ============================================================
# 3. LOAD BUSINESS DATA
# ============================================================

OA_businesses_raw <- read.csv(
  here("data", "processed", "OA_EW_businesses.csv")
)

OA_businesses <- OA_businesses_raw %>%
  rename(OA = OA21CD) %>%
  distinct() %>%
  group_by(OA) %>%
  summarise(
    business_retail_per_km2 =
      sum(business_retail_per_km2, na.rm = TRUE),
    business_accommodation_food_per_km2 =
      sum(business_accommodation_food_per_km2, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 4. LOAD OA GEOMETRY AND AREA
# ============================================================

oa_sub <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  filter(OA %in% OA_matching_dataset$OA) %>%
  mutate(
    area_m2 = as.numeric(st_area(geometry)),
    area_km2 = area_m2 / 1e6
  )

oa_area <- oa_sub %>%
  st_drop_geometry() %>%
  select(OA, area_km2)

# ============================================================
# 5. MERGE CENSUS, BUSINESS, AND AREA DATA
# ============================================================

# Remove any previous area/pop/business variables before re-joining
OA_matching_dataset <- OA_matching_dataset %>%
  select(
    -any_of(c(
      "area_km2",
      "area_m2",
      "pop_density",
      "log_pop_density",
      "log1p_pop_density",
      "business_retail_per_km2",
      "business_accommodation_food_per_km2"
    ))
  )

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by = "OA") %>%
  left_join(OA_businesses, by = "OA") %>%
  left_join(oa_area, by = "OA") %>%
  mutate(
    business_retail_per_km2 =
      replace_na(business_retail_per_km2, 0),
    business_accommodation_food_per_km2 =
      replace_na(business_accommodation_food_per_km2, 0),
    pop_density = Total / area_km2,
    log_pop_density = log1p(pop_density),
    log1p_pop_density = log1p(pop_density)
  )

# Replace missing scheme for controls
OA_matching_census <- OA_matching_census %>%
  mutate(
    scheme = replace_na(scheme, "Control")
  )

# Drop OAs with missing core census percentages
OA_matching_census <- OA_matching_census %>%
  drop_na(
    White_pct,
    Mixed_pct,
    Asian_pct,
    Black_pct,
    Other_ethnicity_pct,
    X4under_pct,
    X5to9_pct,
    X10to14_pct,
    X15to19_pct
  )

# ============================================================
# 6. BASIC POST-MERGE CHECKS
# ============================================================

cat("\n--- Post-merge checks ---\n")
cat("Rows:", nrow(OA_matching_census), "\n")
cat("Unique OAs:", n_distinct(OA_matching_census$OA), "\n")

stopifnot(nrow(OA_matching_census) == n_distinct(OA_matching_census$OA))

cat("\nPopulation density summary:\n")
summary(OA_matching_census$pop_density) %>%
  print()

cat("\nMissingness after merge:\n")
naniar::miss_var_summary(OA_matching_census) %>%
  filter(n_miss > 0) %>%
  print(n = Inf)

# ============================================================
# 7. PRE-MATCHING QA
# ============================================================

qa_pass <- TRUE

qa_fail <- function(msg) {
  cat("  FAIL ✗", msg, "\n")
  qa_pass <<- FALSE
}

qa_ok <- function(msg) {
  cat("  PASS ✓", msg, "\n")
}

cat("\n================================================\n")
cat("PRE-MATCHING QA: ENGLAND-ONLY DATASET\n")
cat("================================================\n")

# ------------------------------------------------------------
# [1] No duplicate OAs
# ------------------------------------------------------------

cat("\n[1] Duplicate OA check\n")

n_rows <- nrow(OA_matching_census)
n_oa <- n_distinct(OA_matching_census$OA)

if (n_rows == n_oa) {
  qa_ok(sprintf("No duplicate OAs: %d rows", n_rows))
} else {
  qa_fail(sprintf("%d rows but only %d distinct OAs", n_rows, n_oa))
}

# ------------------------------------------------------------
# [2] LAD prefix check
# ------------------------------------------------------------

cat("\n[2] LAD prefix check\n")

lad_check <- OA_matching_census %>%
  mutate(lad_prefix = substr(LAD24CD, 1, 1)) %>%
  count(lad_prefix, name = "n_OAs") %>%
  tibble::as_tibble()

print(lad_check, n = Inf)

n_non_english <- OA_matching_census %>%
  summarise(n = sum(substr(LAD24CD, 1, 1) != "E", na.rm = TRUE)) %>%
  pull(n)

if (n_non_english == 0) {
  qa_ok("All OAs have English LAD24CD prefix")
} else {
  qa_fail(sprintf("%d OAs do not have English LAD24CD prefix", n_non_english))
}

# ------------------------------------------------------------
# [3] Assignment counts
# ------------------------------------------------------------

cat("\n[3] Assignment counts\n")

assignment_count <- OA_matching_census %>%
  count(
    across(any_of(c(
      "treated_OA",
      "control_group1_OA",
      "control_group2_OA",
      "buffer_OA",
      "zero_injury_OA"
    ))),
    name = "n_OAs"
  ) %>%
  tibble::as_tibble()

print(assignment_count, n = Inf)

n_treated <- sum(OA_matching_census$treated_OA == 1, na.rm = TRUE)
n_ctrl1 <- sum(OA_matching_census$control_group1_OA == 1, na.rm = TRUE)
n_ctrl2 <- sum(OA_matching_census$control_group2_OA == 1, na.rm = TRUE)

if (n_treated == 0) {
  qa_fail("No treated OAs found")
} else {
  qa_ok(sprintf("Treated OAs: %d", n_treated))
}

cat(sprintf("  INFO: Control group 1 OAs: %d\n", n_ctrl1))

if (n_ctrl2 == 0) {
  qa_fail("No control group 2 OAs found")
} else {
  qa_ok(sprintf("Control group 2 OAs: %d", n_ctrl2))
}

# ------------------------------------------------------------
# [4] Injury-history profile by group
# ------------------------------------------------------------

cat("\n[4] Pre-COVID total-injury profile by group\n")

injury_required <- c(
  "mean_total_inj_pre_covid",
  "trend_total_inj_pre_covid",
  "log1p_mean_total_inj_pre_covid",
  "pct_zero_total_inj_pre_covid"
)

missing_injury_cols <- setdiff(injury_required, names(OA_matching_census))

if (length(missing_injury_cols) > 0) {
  qa_fail(sprintf(
    "Missing injury-history columns: %s",
    paste(missing_injury_cols, collapse = ", ")
  ))
} else {
  
  inj_profile <- OA_matching_census %>%
    mutate(
      group = case_when(
        treated_OA == 1 ~ "Treated",
        control_group1_OA == 1 ~ "Control_1",
        control_group2_OA == 1 ~ "Control_2",
        buffer_OA == 1 ~ "Buffer",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      pct_zero_mean_total =
        round(100 * mean(mean_total_inj_pre_covid == 0, na.rm = TRUE), 1),
      median_mean_total =
        round(median(mean_total_inj_pre_covid, na.rm = TRUE), 4),
      mean_log1p_mean_total =
        round(mean(log1p_mean_total_inj_pre_covid, na.rm = TRUE), 4),
      pct_zero_trend_total =
        round(100 * mean(trend_total_inj_pre_covid == 0, na.rm = TRUE), 1),
      median_trend_total =
        round(median(trend_total_inj_pre_covid, na.rm = TRUE), 4),
      .groups = "drop"
    ) %>%
    arrange(group) %>%
    tibble::as_tibble()
  
  print(inj_profile, n = Inf)
  qa_ok("Pre-COVID total-injury variables available")
}

# ------------------------------------------------------------
# [5] Stage 1 variable completeness: treated OAs
# ------------------------------------------------------------

cat("\n[5] Stage 1 variable completeness: treated OAs\n")

stage1_expected <- c(
  "road_density_m_km2",
  "road_length_km",
  "pct_A_road",
  "pct_B_road",
  "pct_minor_road",
  "dist_BUA_centroid",
  "pop_density",
  "log1p_pop_density",
  "IMD",
  "cars_none_pct",
  "Drive_Car_pct",
  "Walk_pct",
  "Bicycle_pct",
  "X65plus_pct",
  "X5to19_pct",
  "X20to24_pct"
)

missing_cols_s1 <- setdiff(stage1_expected, names(OA_matching_census))

if (length(missing_cols_s1) > 0) {
  qa_fail(sprintf(
    "Stage 1 columns absent: %s",
    paste(missing_cols_s1, collapse = ", ")
  ))
}

present_s1 <- intersect(stage1_expected, names(OA_matching_census))

if (length(present_s1) > 0) {
  
  na_s1 <- OA_matching_census %>%
    filter(treated_OA == 1) %>%
    summarise(
      across(
        all_of(present_s1),
        ~ sum(is.na(.))
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = "variable",
      values_to = "n_NA"
    ) %>%
    filter(n_NA > 0)
  
  if (nrow(na_s1) == 0) {
    qa_ok("No missing Stage 1 values in treated OAs")
  } else {
    cat("  Missing Stage 1 values in treated OAs:\n")
    print(na_s1, n = Inf)
    qa_fail("Missing Stage 1 values in treated OAs")
  }
}

# ------------------------------------------------------------
# [6] Stage 2 variable completeness: matching pool
# ------------------------------------------------------------

cat("\n[6] Stage 2 variable completeness: total-injury matching variables\n")

stage2_expected <- c(
  "log1p_mean_total_inj_pre_covid",
  "trend_total_inj_pre_covid",
  "mean_total_inj_pre_covid",
  "pct_zero_total_inj_pre_covid"
)

missing_cols_s2 <- setdiff(stage2_expected, names(OA_matching_census))

if (length(missing_cols_s2) > 0) {
  qa_fail(sprintf(
    "Stage 2 columns absent: %s",
    paste(missing_cols_s2, collapse = ", ")
  ))
}

present_s2 <- intersect(stage2_expected, names(OA_matching_census))

if (length(present_s2) > 0) {
  
  na_s2 <- OA_matching_census %>%
    filter(
      treated_OA == 1 | control_group2_OA == 1,
      buffer_OA == 0
    ) %>%
    summarise(
      across(
        all_of(present_s2),
        ~ sum(is.na(.))
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = "variable",
      values_to = "n_NA"
    ) %>%
    filter(n_NA > 0)
  
  if (nrow(na_s2) == 0) {
    qa_ok("No missing Stage 2 total-injury variables in matching pool")
  } else {
    cat("  Missing Stage 2 values in matching pool:\n")
    print(na_s2, n = Inf)
    qa_fail("Missing Stage 2 total-injury variables in matching pool")
  }
}

# ------------------------------------------------------------
# [7] Road variable sanity
# ------------------------------------------------------------

cat("\n[7] Road variable sanity\n")

road_sanity <- OA_matching_census %>%
  filter(treated_OA == 1 | control_group2_OA == 1) %>%
  summarise(
    n_negative_road_length =
      sum(road_length_km < 0, na.rm = TRUE),
    n_negative_density =
      sum(road_density_m_km2 < 0, na.rm = TRUE),
    n_pct_A_over_100 =
      sum(pct_A_road > 100, na.rm = TRUE),
    n_pct_B_over_100 =
      sum(pct_B_road > 100, na.rm = TRUE),
    n_pct_minor_over_100 =
      sum(pct_minor_road > 100, na.rm = TRUE),
    n_zero_roads =
      sum(n_roads == 0 | is.na(n_roads), na.rm = TRUE)
  )

print(road_sanity)

if (
  road_sanity$n_negative_road_length == 0 &&
  road_sanity$n_negative_density == 0 &&
  road_sanity$n_pct_A_over_100 == 0 &&
  road_sanity$n_pct_B_over_100 == 0 &&
  road_sanity$n_pct_minor_over_100 == 0
) {
  qa_ok("Road variables have no impossible values")
} else {
  qa_fail("Impossible values in road variables")
}

# ------------------------------------------------------------
# [8] Zero-injury flag consistency
# ------------------------------------------------------------

cat("\n[8] zero_injury_OA flag\n")

if (!"zero_injury_OA" %in% names(OA_matching_census)) {
  qa_fail("zero_injury_OA column is missing")
} else {
  
  zero_flag_check <- OA_matching_census %>%
    filter(treated_OA == 1 | control_group2_OA == 1) %>%
    count(
      treated_OA,
      control_group2_OA,
      zero_injury_OA,
      name = "n_OAs"
    ) %>%
    tibble::as_tibble()
  
  print(zero_flag_check, n = Inf)
  
  n_zero_treated <- sum(
    OA_matching_census$treated_OA == 1 &
      OA_matching_census$zero_injury_OA == 1,
    na.rm = TRUE
  )
  
  cat(sprintf(
    "  Zero-injury treated OAs: %d / %d (%.1f%%)\n",
    n_zero_treated,
    n_treated,
    100 * n_zero_treated / max(n_treated, 1)
  ))
  
  qa_ok("zero_injury_OA flag present")
}

# ------------------------------------------------------------
# [9] Census variable completeness: treated OAs
# ------------------------------------------------------------

cat("\n[9] Census variable completeness: treated OAs\n")

census_vars <- c(
  "IMD",
  "cars_none_pct",
  "Drive_Car_pct",
  "Walk_pct",
  "Bicycle_pct",
  "X65plus_pct",
  "X5to19_pct",
  "X20to24_pct",
  "pop_density",
  "log1p_pop_density"
)

missing_census_vars <- setdiff(census_vars, names(OA_matching_census))

if (length(missing_census_vars) > 0) {
  qa_fail(sprintf(
    "Census columns absent: %s",
    paste(missing_census_vars, collapse = ", ")
  ))
}

present_census_vars <- intersect(census_vars, names(OA_matching_census))

if (length(present_census_vars) > 0) {
  
  na_census <- OA_matching_census %>%
    filter(treated_OA == 1) %>%
    summarise(
      across(
        all_of(present_census_vars),
        ~ sum(is.na(.))
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = "variable",
      values_to = "n_NA"
    ) %>%
    filter(n_NA > 0)
  
  if (nrow(na_census) == 0) {
    qa_ok("No missing census variables in treated OAs")
  } else {
    cat("  Missing census values in treated OAs:\n")
    print(na_census, n = Inf)
    qa_fail("Missing census values in treated OAs")
  }
}

# ------------------------------------------------------------
# [10] Travel mode percentage consistency
# ------------------------------------------------------------

cat("\n[10] Travel mode percentage consistency\n")

travel_pct_cols <- c(
  "Drive_Car_pct",
  "Passenger_Car_pct",
  "Walk_pct",
  "Bicycle_pct",
  "bus_Coach_pct",
  "Train_pct",
  "Underground_train_tram_pct",
  "Taxi_pct",
  "workAthome_pct",
  "Other_pct",
  "Motorcycle_pct"
)

missing_travel_cols <- setdiff(travel_pct_cols, names(OA_matching_census))

if (length(missing_travel_cols) > 0) {
  qa_fail(sprintf(
    "Travel percentage columns absent: %s",
    paste(missing_travel_cols, collapse = ", ")
  ))
} else {
  
  travel_check <- OA_matching_census %>%
    filter(travel_base_zero_flag == 0) %>%
    mutate(
      travel_mode_sum =
        rowSums(across(all_of(travel_pct_cols)), na.rm = TRUE)
    )
  
  n_bad_travel <- sum(
    abs(travel_check$travel_mode_sum - 100) > 0.01,
    na.rm = TRUE
  )
  
  if (n_bad_travel == 0) {
    qa_ok(sprintf(
      "All %d OAs have travel-mode percentages summing to 100",
      nrow(travel_check)
    ))
  } else {
    qa_fail(sprintf(
      "%d OAs have travel-mode percentages not summing to 100",
      n_bad_travel
    ))
  }
}

# ------------------------------------------------------------
# [11] Final matching-variable summary
# ------------------------------------------------------------

cat("\n[11] Suggested injury-history matching variables\n")

injury_matching_vars <- c(
  "log1p_mean_total_inj_pre_covid",
  "trend_total_inj_pre_covid"
)

print(injury_matching_vars)

# ============================================================
# 8. FINAL QA RESULT
# ============================================================

cat("\n================================================\n")

if (qa_pass) {
  cat("ALL QA CHECKS PASSED — safe to save and proceed to matching\n")
} else {
  cat("ONE OR MORE QA CHECKS FAILED — fix issues above before matching\n")
}

cat("================================================\n\n")

stopifnot("Pre-matching QA failed — see output above" = qa_pass)

# ============================================================
# 9. SAVE OUTPUTS
# ============================================================

saveRDS(
  OA_matching_census,
  out_rds
)

cat("Saved RDS to:\n", out_rds, "\n")

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by = "OA") %>%
  st_as_sf()

st_write(
  OA_matching_census_sf,
  out_gpkg,
  delete_dsn = TRUE,
  quiet = TRUE
)

cat("Saved GPKG to:\n", out_gpkg, "\n")
