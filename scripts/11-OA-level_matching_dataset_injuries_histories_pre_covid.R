# ============================================================
# Create Pooled OA-Level Matching Dataset for CAZ DiD Analysis
# ============================================================
#
# Purpose:
#   Construct one pooled OA-level matching dataset for the CAZ DiD analysis.
#
# Key decisions:
#   1. Injury-history matching variables use TOTAL injuries only.
#   2. No stratification by severity or road-user mode.
#   3. Injury baseline and trend are calculated using the pre-COVID period:
#        all quarters up to and including 2019 Q4.
#   4. OA-quarter injury panel is zero-filled before calculating means/trends.
#   5. Main injury variables for Stage 2 matching are per road-km:
#        mean_total_pkm
#        log1p_mean_total_pkm
#        trend_total_pkm
#        recent_minus_mid_total_pkm
#        mid_minus_early_total_pkm
#
# Main output:
#   data/processed/OA_matching_data_pooled.rds
#
# Compatibility output:
#   data/processed/OA_matching_dataset.rds
#
# ============================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(here)
library(sf)
library(naniar)

select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

# ============================================================
# 0. SETTINGS
# ============================================================

target_crs <- 27700

# This is pre-2020 / pre-COVID
pre_covid_end_q <- as.yearqtr("2019 Q4")

# Pre-COVID trajectory windows used to add simple outcome-path shape variables.
# These complement the single whole-period slope used in Stage 2 matching.
pre_covid_recent_start_q <- as.yearqtr("2018 Q1")
pre_covid_mid_start_q    <- as.yearqtr("2016 Q1")
pre_covid_mid_end_q      <- as.yearqtr("2017 Q4")

out_rds_main   <- here("data", "processed", "OA_matching_data_pooled.rds")
out_rds_compat <- here("data", "processed", "OA_matching_dataset.rds")
out_gpkg       <- here("data", "processed", "shp_files", "OA_matching_data_pooled.gpkg")
out_plot       <- here("output", "diagnostics", "fig_pre_covid_total_pkm_trends_pooled.png")

dir.create(dirname(out_rds_main), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

# Robust conversion to year-quarter
to_yearqtr <- function(x) {
  
  if (inherits(x, "yearqtr")) return(x)
  if (inherits(x, "Date"))    return(as.yearqtr(x))
  if (inherits(x, "POSIXt"))  return(as.yearqtr(as.Date(x)))
  if (is.numeric(x))          return(as.yearqtr(x))
  
  x_chr <- as.character(x)
  
  out <- suppressWarnings(as.yearqtr(x_chr, format = "%Y Q%q"))
  
  missing <- is.na(out)
  
  if (any(missing)) {
    out[missing] <- suppressWarnings(as.yearqtr(as.Date(x_chr[missing])))
  }
  
  out
}

# ============================================================
# 1. LOAD INPUTS
# ============================================================

OA_analysis <- readRDS(
  here("data", "processed", "OA_level_from_polygons.rds")
)

OA_roads <- readRDS(
  here("data", "processed", "OA_roads.rds")
)

OA_injuries <- readRDS(
  here("data", "processed", "OA_injuries_quarterly.rds")
)

oa_sub <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(target_crs) %>%
  st_make_valid()

caz <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

# Drop geometry from OA_analysis if needed
OA_analysis_base <- OA_analysis

if (inherits(OA_analysis_base, "sf")) {
  OA_analysis_base <- st_drop_geometry(OA_analysis_base)
}

# ============================================================
# 2. BASIC INPUT CHECKS
# ============================================================

stopifnot("OA" %in% names(OA_analysis_base))
stopifnot("OA" %in% names(OA_roads))
stopifnot("OA" %in% names(OA_injuries))
stopifnot("quarter_year" %in% names(OA_injuries))
stopifnot("total_injuries" %in% names(OA_injuries))

cat("\n--- Input checks ---\n")
cat("OAs in OA_analysis:", n_distinct(OA_analysis_base$OA), "\n")
cat("OAs in OA_roads:", n_distinct(OA_roads$OA), "\n")
cat("OAs in OA_injuries:", n_distinct(OA_injuries$OA), "\n")

# ============================================================
# 3. CAZ START DATE LOOKUP
# ============================================================

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(
    caz_start_date = dmy(startDt),
    caz_start_q    = as.yearqtr(caz_start_date)
  ) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    caz_start_q    = min(caz_start_q, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- CAZ start dates ---\n")
print(caz_dates)

if (any(caz_dates$caz_start_q <= pre_covid_end_q, na.rm = TRUE)) {
  warning(
    "At least one CAZ starts before or during the pre-COVID window. ",
    "Check whether using all quarters up to 2019 Q4 is valid for every scheme."
  )
}

# ============================================================
# 4. OA SCHEME / TREATMENT LOOKUP
# ============================================================

oa_scheme_lookup <- OA_analysis_base %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE) %>%
  select(
    OA,
    any_of(c(
      "scheme",
      "treated_OA",
      "control_group1_OA",
      "control_group2_OA",
      "buffer_OA",
      "assignment"
    ))
  ) %>%
  left_join(caz_dates, by = "scheme")

cat("\n--- OA scheme lookup ---\n")
cat("Rows:", nrow(oa_scheme_lookup), "\n")
cat("Unique OAs:", n_distinct(oa_scheme_lookup$OA), "\n")

dupes_lookup <- oa_scheme_lookup %>%
  count(OA) %>%
  filter(n > 1)

cat("Duplicate OAs in lookup:", nrow(dupes_lookup), "\n")
stopifnot(nrow(dupes_lookup) == 0)

# ============================================================
# 5. CLEAN ROAD NETWORK VARIABLES
# ============================================================

OA_roads_clean <- OA_roads %>%
  group_by(OA) %>%
  summarise(
    n_roads           = sum(n_roads, na.rm = TRUE),
    total_road_length = sum(total_road_length, na.rm = TRUE),
    n_A               = sum(n_A, na.rm = TRUE),
    n_B               = sum(n_B, na.rm = TRUE),
    n_motorway        = sum(n_motorway, na.rm = TRUE),
    n_minor           = sum(n_minor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    oa_sub %>%
      st_drop_geometry() %>%
      select(OA),
    by = "OA"
  ) %>%
  mutate(
    across(
      c(n_roads, total_road_length, n_A, n_B, n_motorway, n_minor),
      ~ replace_na(.x, 0)
    ),
    road_length_km = total_road_length / 1000
  )

cat("\n--- Road data ---\n")
cat("OAs in cleaned road data:", n_distinct(OA_roads_clean$OA), "\n")
cat("OAs with no roads:", sum(OA_roads_clean$n_roads == 0), "\n")

# OA area
oa_area <- oa_sub %>%
  st_make_valid() %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(OA, area_km2)

# Road composition
# Do not include road_length_km here because it already exists in OA_roads_clean.
road_composition <- OA_roads_clean %>%
  select(OA, total_road_length, n_roads, n_A, n_B, n_motorway, n_minor) %>%
  left_join(oa_area, by = "OA") %>%
  mutate(
    pct_A_road = if_else(n_roads > 0, 100 * n_A / n_roads, 0),
    pct_B_road = if_else(n_roads > 0, 100 * n_B / n_roads, 0),
    pct_motorway_road = if_else(n_roads > 0, 100 * n_motorway / n_roads, 0),
    pct_minor_road = if_else(n_roads > 0, 100 * n_minor / n_roads, 0),
    road_density_m_km2 = if_else(
      !is.na(area_km2) & area_km2 > 0,
      total_road_length / area_km2,
      0
    )
  ) %>%
  select(
    OA,
    pct_A_road,
    pct_B_road,
    pct_motorway_road,
    pct_minor_road,
    road_density_m_km2
  )

# ============================================================
# 6. CLEAN TOTAL-INJURY DATA ONLY
# ============================================================

OA_injuries_total <- OA_injuries %>%
  mutate(
    quarter_year = to_yearqtr(quarter_year)
  ) %>%
  group_by(OA, quarter_year) %>%
  summarise(
    total_injuries = sum(total_injuries, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Injury data ---\n")
cat("Observed OA-quarter injury rows:", nrow(OA_injuries_total), "\n")
cat("Observed quarters:", n_distinct(OA_injuries_total$quarter_year), "\n")
cat("First quarter:", as.character(min(OA_injuries_total$quarter_year, na.rm = TRUE)), "\n")
cat("Last quarter:", as.character(max(OA_injuries_total$quarter_year, na.rm = TRUE)), "\n")

# ============================================================
# 7. BUILD ZERO-FILLED PRE-COVID BALANCED PANEL
# ============================================================

pre_covid_quarters <- OA_injuries_total %>%
  distinct(quarter_year) %>%
  filter(quarter_year <= pre_covid_end_q) %>%
  arrange(quarter_year) %>%
  mutate(qtr_index = row_number() - 1)

all_oas <- oa_scheme_lookup %>%
  select(OA) %>%
  distinct()

OA_injuries_pre_covid_balanced <- all_oas %>%
  cross_join(pre_covid_quarters) %>%
  left_join(
    OA_injuries_total,
    by = c("OA", "quarter_year")
  ) %>%
  mutate(
    total_injuries = replace_na(total_injuries, 0)
  ) %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  left_join(
    OA_roads_clean %>%
      select(OA, road_length_km),
    by = "OA"
  ) %>%
  mutate(
    road_length_km = replace_na(road_length_km, 0),
    log_road_km = if_else(
      road_length_km > 0,
      log(road_length_km),
      NA_real_
    )
  )

cat("\n--- Balanced pre-COVID panel ---\n")
cat("OAs:", n_distinct(OA_injuries_pre_covid_balanced$OA), "\n")
cat("Quarters:", n_distinct(OA_injuries_pre_covid_balanced$quarter_year), "\n")
cat("Rows:", nrow(OA_injuries_pre_covid_balanced), "\n")
cat(
  "Expected:",
  n_distinct(OA_injuries_pre_covid_balanced$OA) *
    n_distinct(OA_injuries_pre_covid_balanced$quarter_year),
  "\n"
)

stopifnot(
  nrow(OA_injuries_pre_covid_balanced) ==
    n_distinct(OA_injuries_pre_covid_balanced$OA) *
    n_distinct(OA_injuries_pre_covid_balanced$quarter_year)
)

# ============================================================
# 8. SLOPE FUNCTION: QUASI-POISSON TREND
# ============================================================

# Returns log-rate change per quarter.
# If offset is supplied, this estimates the trend in injury rate per road-km.
# All-zero OAs are assigned slope = 0.
# Constant-count OAs are also assigned slope = 0.

poisson_slope <- function(y, x, offset = NULL) {
  
  df <- tibble(
    y = y,
    x = x,
    off = if (is.null(offset)) rep(0, length(x)) else offset
  ) %>%
    filter(
      !is.na(y),
      !is.na(x),
      !is.na(off)
    )
  
  if (nrow(df) < 2) return(NA_real_)
  if (all(df$y == 0)) return(0)
  if (var(df$y) == 0) return(0)
  
  model <- tryCatch(
    {
      if (is.null(offset)) {
        glm(y ~ x, family = quasipoisson(link = "log"), data = df)
      } else {
        glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
      }
    },
    warning = function(w) {
      suppressWarnings(
        tryCatch(
          {
            if (is.null(offset)) {
              glm(y ~ x, family = quasipoisson(link = "log"), data = df)
            } else {
              glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
            }
          },
          error = function(e) NULL
        )
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NA_real_)
  
  slope <- coef(model)[["x"]]
  
  if (is.null(slope) || is.na(slope) || is.nan(slope) || is.infinite(slope)) {
    return(NA_real_)
  }
  
  slope
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ============================================================
# 9. TOTAL-INJURY BASELINE AND PRE-COVID TREND FEATURES
# ============================================================

injury_features_pre_covid <- OA_injuries_pre_covid_balanced %>%
  group_by(OA) %>%
  summarise(
    n_pre_covid_quarters = n(),
    
    road_length_km = first(road_length_km),
    
    total_injuries_pre_covid =
      sum(total_injuries, na.rm = TRUE),
    
    mean_total =
      mean(total_injuries, na.rm = TRUE),
    
    mean_total_inj_pre_covid =
      mean(total_injuries, na.rm = TRUE),
    
    mean_total_pkm =
      if_else(
        first(road_length_km) > 0,
        mean(total_injuries, na.rm = TRUE) / first(road_length_km),
        0
      ),
    
    pct_zero_total_inj_pre_covid =
      100 * mean(total_injuries == 0, na.rm = TRUE),
    
    trend_total =
      poisson_slope(total_injuries, qtr_index),
    
    trend_total_inj_pre_covid =
      poisson_slope(total_injuries, qtr_index),
    
    trend_total_pkm =
      poisson_slope(total_injuries, qtr_index, log_road_km),
    
    mean_total_pkm_recent_pre_covid =
      if_else(
        first(road_length_km) > 0,
        safe_mean(
          total_injuries[
            quarter_year >= pre_covid_recent_start_q &
              quarter_year <= pre_covid_end_q
          ]
        ) / first(road_length_km),
        0
      ),
    
    mean_total_pkm_mid_pre_covid =
      if_else(
        first(road_length_km) > 0,
        safe_mean(
          total_injuries[
            quarter_year >= pre_covid_mid_start_q &
              quarter_year <= pre_covid_mid_end_q
          ]
        ) / first(road_length_km),
        0
      ),
    
    mean_total_pkm_early_pre_covid =
      if_else(
        first(road_length_km) > 0,
        safe_mean(total_injuries[quarter_year < pre_covid_mid_start_q]) /
          first(road_length_km),
        0
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    trend_total =
      replace_na(trend_total, 0),
    
    trend_total_inj_pre_covid =
      replace_na(trend_total_inj_pre_covid, 0),
    
    trend_total_pkm =
      replace_na(trend_total_pkm, 0),
    
    mean_total_pkm_recent_pre_covid =
      replace_na(mean_total_pkm_recent_pre_covid, 0),
    
    mean_total_pkm_mid_pre_covid =
      replace_na(mean_total_pkm_mid_pre_covid, 0),
    
    mean_total_pkm_early_pre_covid =
      replace_na(mean_total_pkm_early_pre_covid, 0),
    
    recent_minus_mid_total_pkm =
      mean_total_pkm_recent_pre_covid - mean_total_pkm_mid_pre_covid,
    
    mid_minus_early_total_pkm =
      mean_total_pkm_mid_pre_covid - mean_total_pkm_early_pre_covid,
    
    mean_total =
      replace_na(mean_total, 0),
    
    mean_total_inj_pre_covid =
      replace_na(mean_total_inj_pre_covid, 0),
    
    mean_total_pkm =
      replace_na(mean_total_pkm, 0),
    
    log1p_total_injuries_pre_covid =
      log1p(total_injuries_pre_covid),
    
    log1p_mean_total_inj_pre_covid =
      log1p(mean_total_inj_pre_covid),
    
    log1p_mean_total_pkm =
      log1p(mean_total_pkm)
  )

cat("\n--- Pre-COVID total injury features ---\n")

summary(
  injury_features_pre_covid %>%
    select(
      mean_total,
      mean_total_pkm,
      log1p_mean_total_pkm,
      trend_total,
      trend_total_pkm,
      recent_minus_mid_total_pkm,
      mid_minus_early_total_pkm
    )
)

# ============================================================
# 10. ZERO-INJURY OA FLAG — PRE-TREATMENT PERIOD ONLY
# ============================================================

# For treated OAs: use only quarters before their scheme's CAZ start
# For control OAs: use only quarters up to the earliest CAZ start
# This prevents dropping OAs whose injuries fell to zero post-treatment

earliest_caz_q <- min(caz_dates$caz_start_q, na.rm = TRUE)

cat("\nEarliest CAZ start quarter:", as.character(earliest_caz_q), "\n")

# Join scheme start dates onto injury panel
OA_injuries_with_scheme <- OA_injuries_total %>%
  left_join(
    oa_scheme_lookup %>%
      select(OA, treated_OA, control_group1_OA,
             control_group2_OA, caz_start_q),
    by = "OA"
  )

# For each OA, determine the relevant pre-treatment cutoff
# Treated OAs: use their own scheme's CAZ start
# Control OAs: use the earliest CAZ start (most conservative)
OA_injuries_pretreat <- OA_injuries_with_scheme %>%
  mutate(
    cutoff_q = case_when(
      !is.na(caz_start_q) ~ caz_start_q,   # treated: own scheme cutoff
      TRUE                ~ earliest_caz_q  # controls: earliest CAZ
    )
  ) %>%
  filter(quarter_year < cutoff_q)           # strictly before CAZ opens

cat("\nPre-treatment injury rows:", nrow(OA_injuries_pretreat), "\n")
cat("Quarter range used:", 
    as.character(min(OA_injuries_pretreat$quarter_year)),
    "to",
    as.character(max(OA_injuries_pretreat$quarter_year)), "\n")

# Flag OAs with zero injuries in pre-treatment period only
injury_pretreat_summary <- OA_injuries_pretreat %>%
  group_by(OA) %>%
  summarise(
    total_injuries_pretreat   = sum(total_injuries, na.rm = TRUE),
    ever_injury_OA            = as.integer(total_injuries_pretreat > 0),
    .groups = "drop"
  )

# Also keep total_injuries_all_periods for reference (not used for flagging)
injury_all_periods <- OA_injuries_total %>%
  group_by(OA) %>%
  summarise(
    total_injuries_all_periods = sum(total_injuries, na.rm = TRUE),
    .groups = "drop"
  )

zero_injury_lookup <- all_oas %>%
  left_join(injury_pretreat_summary, by = "OA") %>%
  left_join(injury_all_periods,      by = "OA") %>%
  mutate(
    total_injuries_pretreat    = replace_na(total_injuries_pretreat, 0),
    total_injuries_all_periods = replace_na(total_injuries_all_periods, 0),
    ever_injury_OA             = replace_na(ever_injury_OA, 0L),
    zero_injury_OA             = if_else(ever_injury_OA == 1L, 0L, 1L)
  ) %>%
  select(
    OA,
    total_injuries_pretreat,
    total_injuries_all_periods,
    ever_injury_OA,
    zero_injury_OA
  )

cat("\n--- Zero-injury OA flag (PRE-TREATMENT ONLY) ---\n")
zero_injury_lookup %>%
  count(zero_injury_OA) %>%
  print()

# Verification: how many OAs have pre-treatment injuries = 0
# but post-treatment injuries > 0?
# These are the OAs you would have WRONGLY dropped under the old code
rescued_OAs <- zero_injury_lookup %>%
  filter(
    zero_injury_OA             == 0,   # retained (has pre-treat injuries)
    total_injuries_all_periods  > total_injuries_pretreat  # has post-treat injuries too
  ) %>%
  nrow()

wrongly_dropped <- all_oas %>%
  left_join(injury_pretreat_summary, by = "OA") %>%
  left_join(injury_all_periods,      by = "OA") %>%
  mutate(
    total_injuries_pretreat    = replace_na(total_injuries_pretreat, 0),
    total_injuries_all_periods = replace_na(total_injuries_all_periods, 0)
  ) %>%
  filter(
    total_injuries_pretreat    == 0,   # zero pre-treatment
    total_injuries_all_periods  > 0    # but has post-treatment injuries
  ) %>%
  left_join(
    oa_scheme_lookup %>% select(OA, treated_OA, scheme),
    by = "OA"
  )

cat("\nOAs with zero pre-treatment but non-zero post-treatment injuries:\n")
cat("(these would have been WRONGLY dropped by the old all-periods flag)\n")
cat("Total:", nrow(wrongly_dropped), "\n")

if (nrow(wrongly_dropped) > 0) {
  cat("\nBreakdown by treatment status:\n")
  wrongly_dropped %>%
    count(treated_OA, scheme) %>%
    print()
}
# ============================================================
# 11. ASSEMBLE POOLED OA-LEVEL MATCHING DATASET
# ============================================================

vars_to_refresh <- c(
  "n_roads",
  "total_road_length",
  "n_A",
  "n_B",
  "n_motorway",
  "n_minor",
  "road_length_km",
  "pct_A_road",
  "pct_B_road",
  "pct_motorway_road",
  "pct_minor_road",
  "road_density_m_km2",
  "log1p_road_length_km",
  "log1p_road_density_m_km2",
  
  # new total-only injury variables
  "n_pre_covid_quarters",
  "total_injuries_pre_covid",
  "mean_total",
  "mean_total_inj_pre_covid",
  "mean_total_pkm",
  "log1p_mean_total_pkm",
  "pct_zero_total_inj_pre_covid",
  "trend_total",
  "trend_total_inj_pre_covid",
  "trend_total_pkm",
  "mean_total_pkm_recent_pre_covid",
  "mean_total_pkm_mid_pre_covid",
  "mean_total_pkm_early_pre_covid",
  "recent_minus_mid_total_pkm",
  "mid_minus_early_total_pkm",
  "log1p_total_injuries_pre_covid",
  "log1p_mean_total_inj_pre_covid",
  "total_injuries_all_periods",
  "ever_injury_OA",
  "zero_injury_OA",
  
  # old names from previous versions
  "n_pre2019_quarters",
  "total_injuries_pre2019",
  "mean_total_inj_pre2019",
  "pct_zero_total_inj_pre2019",
  "trend_total_inj_pre2019",
  "log1p_total_injuries_pre2019",
  "log1p_mean_total_inj_pre2019",
  
  # old mode/severity/per-km names to remove if present
  "mean_car",
  "mean_cyc",
  "mean_ped",
  "trend_car",
  "trend_cyc",
  "trend_ped",
  "mean_car_pkm",
  "mean_cyc_pkm",
  "mean_ped_pkm",
  "trend_car_pkm",
  "trend_cyc_pkm",
  "trend_ped_pkm"
)

OA_base_unique <- OA_analysis_base %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE) %>%
  select(-any_of(vars_to_refresh))

OA_matching_data_pooled <- OA_base_unique %>%
  left_join(OA_roads_clean, by = "OA") %>%
  left_join(road_composition, by = "OA") %>%
  left_join(
    injury_features_pre_covid %>%
      select(-road_length_km),
    by = "OA"
  ) %>%
  left_join(zero_injury_lookup, by = "OA") %>%
  mutate(
    across(
      any_of(c(
        "n_roads",
        "total_road_length",
        "n_A",
        "n_B",
        "n_motorway",
        "n_minor",
        "road_length_km",
        "pct_A_road",
        "pct_B_road",
        "pct_motorway_road",
        "pct_minor_road",
        "road_density_m_km2"
      )),
      ~ replace_na(.x, 0)
    ),
    
    across(
      any_of(c(
        "n_pre_covid_quarters",
        "total_injuries_pre_covid",
        "mean_total",
        "mean_total_inj_pre_covid",
        "mean_total_pkm",
        "log1p_mean_total_pkm",
        "pct_zero_total_inj_pre_covid",
        "trend_total",
        "trend_total_inj_pre_covid",
        "trend_total_pkm",
        "mean_total_pkm_recent_pre_covid",
        "mean_total_pkm_mid_pre_covid",
        "mean_total_pkm_early_pre_covid",
        "recent_minus_mid_total_pkm",
        "mid_minus_early_total_pkm",
        "log1p_total_injuries_pre_covid",
        "log1p_mean_total_inj_pre_covid",
        "total_injuries_all_periods",
        "ever_injury_OA",
        "zero_injury_OA"
      )),
      ~ replace_na(.x, 0)
    ),
    
    log1p_road_length_km =
      log1p(road_length_km),
    
    log1p_road_density_m_km2 =
      log1p(road_density_m_km2)
  )

# Compatibility checks for matching script
stopifnot("road_length_km" %in% names(OA_matching_data_pooled))
stopifnot("mean_total_pkm" %in% names(OA_matching_data_pooled))
stopifnot("log1p_mean_total_pkm" %in% names(OA_matching_data_pooled))
stopifnot("trend_total_pkm" %in% names(OA_matching_data_pooled))
stopifnot("recent_minus_mid_total_pkm" %in% names(OA_matching_data_pooled))
stopifnot("mid_minus_early_total_pkm" %in% names(OA_matching_data_pooled))

# Matching-script reminder:
# stage2_trends <- c(
#   "trend_total_pkm",
#   "recent_minus_mid_total_pkm",
#   "mid_minus_early_total_pkm"
# )

# ============================================================
# 12. CHECKS
# ============================================================

cat("\n--- Final dataset checks ---\n")

cat("Rows:", nrow(OA_matching_data_pooled), "\n")
cat("Unique OAs:", n_distinct(OA_matching_data_pooled$OA), "\n")

dupes <- OA_matching_data_pooled %>%
  count(OA) %>%
  filter(n > 1)

cat("Duplicate OAs:", nrow(dupes), "\n")
stopifnot(nrow(dupes) == 0)

missing_from_matching <- anti_join(
  oa_sub %>%
    st_drop_geometry() %>%
    select(OA),
  OA_matching_data_pooled,
  by = "OA"
)

cat("OAs in oa_sub missing from matching data:", nrow(missing_from_matching), "\n")

if ("assignment" %in% names(OA_matching_data_pooled)) {
  
  cat("\nAssignment counts:\n")
  
  assignment_check <- OA_matching_data_pooled %>%
    count(assignment, name = "n_OAs") %>%
    tibble::as_tibble()
  
  print(assignment_check, n = Inf)
}

cat("\nTreatment/control/zero-injury counts:\n")

group_count_check <- OA_matching_data_pooled %>%
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

print(group_count_check, n = Inf)

cat("\nVariables with missing values:\n")

miss_var_summary(OA_matching_data_pooled) %>%
  filter(n_miss > 0) %>%
  print(n = Inf)

# Balance summary: treated vs other-city controls
if (all(c("treated_OA", "control_group2_OA") %in% names(OA_matching_data_pooled))) {
  
  cat("\nPre-COVID total-per-km injury features: treated vs control_group2\n")
  
  OA_matching_data_pooled %>%
    filter(treated_OA == 1 | control_group2_OA == 1) %>%
    mutate(group = if_else(treated_OA == 1, "Treated", "Control group 2")) %>%
    group_by(group) %>%
    summarise(
      n_OAs =
        n(),
      mean_total =
        mean(mean_total, na.rm = TRUE),
      mean_total_pkm =
        mean(mean_total_pkm, na.rm = TRUE),
      mean_log1p_mean_total_pkm =
        mean(log1p_mean_total_pkm, na.rm = TRUE),
      mean_trend_total_pkm =
        mean(trend_total_pkm, na.rm = TRUE),
      mean_recent_minus_mid_total_pkm =
        mean(recent_minus_mid_total_pkm, na.rm = TRUE),
      mean_mid_minus_early_total_pkm =
        mean(mid_minus_early_total_pkm, na.rm = TRUE),
      mean_pct_zero_total_inj_pre_covid =
        mean(pct_zero_total_inj_pre_covid, na.rm = TRUE),
      mean_road_length_km =
        mean(road_length_km, na.rm = TRUE),
      mean_road_density_m_km2 =
        mean(road_density_m_km2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    print()
}

# ============================================================
# 13. SIMPLE PRE-COVID PARALLEL TREND PLOT
# ============================================================

trend_plot_data <- OA_injuries_pre_covid_balanced %>%
  mutate(
    total_pkm = if_else(
      road_length_km > 0,
      total_injuries / road_length_km,
      0
    ),
    group = case_when(
      treated_OA == 1 ~ "Treated",
      control_group1_OA == 1 ~ "Control group 1",
      control_group2_OA == 1 ~ "Control group 2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  group_by(group, quarter_year) %>%
  summarise(
    mean_total_pkm = mean(total_pkm, na.rm = TRUE),
    .groups = "drop"
  )

p_pretrend <- ggplot(
  trend_plot_data,
  aes(x = as.Date(quarter_year), y = mean_total_pkm, colour = group)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  labs(
    title = "Pre-COVID injury trends per road-km",
    subtitle = "Total injuries only; zero-filled OA-quarter panel; up to 2019 Q4",
    x = NULL,
    y = "Mean total injuries per road-km per OA-quarter",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_pretrend)

ggsave(
  out_plot,
  p_pretrend,
  width = 10,
  height = 6,
  dpi = 300
)

cat("\nSaved plot to:\n", out_plot, "\n")

# ============================================================
# 14. SAVE OUTPUTS
# ============================================================

saveRDS(
  OA_matching_data_pooled,
  out_rds_main
)

# Compatibility save, in case later scripts still expect this name
saveRDS(
  OA_matching_data_pooled,
  out_rds_compat
)

cat("\nSaved main RDS to:\n", out_rds_main, "\n")
cat("Saved compatibility RDS to:\n", out_rds_compat, "\n")

OA_matching_data_pooled_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_data_pooled, by = "OA")

st_write(
  OA_matching_data_pooled_sf,
  out_gpkg,
  delete_dsn = TRUE,
  quiet = TRUE
)

cat("Saved GPKG to:\n", out_gpkg, "\n")

# ============================================================
# 15. VARIABLE DESCRIPTIONS
# ============================================================

var_description <- tibble::tribble(
  ~variable, ~description,
  "OA", "Output Area identifier. Each row is one OA.",
  "n_roads", "Total number of road segments intersecting the OA.",
  "total_road_length", "Total length of roads within the OA, in metres.",
  "road_length_km", "Total length of roads within the OA, in kilometres.",
  "n_A", "Number of A-road segments within the OA.",
  "n_B", "Number of B-road segments within the OA.",
  "n_motorway", "Number of motorway segments within the OA.",
  "n_minor", "Number of minor/local road segments within the OA.",
  "pct_A_road", "Percentage of road segments in the OA that are A-roads.",
  "pct_B_road", "Percentage of road segments in the OA that are B-roads.",
  "pct_motorway_road", "Percentage of road segments in the OA that are motorways.",
  "pct_minor_road", "Percentage of road segments in the OA that are minor roads.",
  "road_density_m_km2", "Total road length per square kilometre of OA area.",
  "log1p_road_length_km", "Log-transformed road length, log(1 + road_length_km).",
  "log1p_road_density_m_km2", "Log-transformed road density, log(1 + road_density_m_km2).",
  "n_pre_covid_quarters", "Number of zero-filled OA-quarter observations used to construct pre-COVID injury features.",
  "total_injuries_pre_covid", "Total injuries in the OA up to and including 2019 Q4.",
  "mean_total", "Mean quarterly total injuries per OA up to and including 2019 Q4. Kept for compatibility.",
  "mean_total_inj_pre_covid", "Mean quarterly total injuries per OA up to and including 2019 Q4.",
  "mean_total_pkm", "Mean quarterly total injuries per road-km up to and including 2019 Q4. Main Stage 2 level variable.",
  "log1p_mean_total_pkm", "Log-transformed mean quarterly total injuries per road-km. Main Stage 2 level variable.",
  "pct_zero_total_inj_pre_covid", "Percentage of pre-COVID quarters with zero total injuries.",
  "trend_total", "Quasi-Poisson slope for total injuries up to and including 2019 Q4. Kept for compatibility.",
  "trend_total_inj_pre_covid", "Quasi-Poisson slope for total injuries up to and including 2019 Q4.",
  "trend_total_pkm", "Quasi-Poisson slope with log road-km offset for total injury rate. Main Stage 2 trend variable.",
  "mean_total_pkm_recent_pre_covid", "Mean quarterly total injuries per road-km in the recent pre-COVID window (2018 Q1 to 2019 Q4).",
  "mean_total_pkm_mid_pre_covid", "Mean quarterly total injuries per road-km in the mid pre-COVID window (2016 Q1 to 2017 Q4).",
  "mean_total_pkm_early_pre_covid", "Mean quarterly total injuries per road-km in the early pre-COVID window (before 2016 Q1).",
  "recent_minus_mid_total_pkm", "Difference in mean quarterly total injuries per road-km between recent and mid pre-COVID windows. Stage 2 trajectory-shape variable.",
  "mid_minus_early_total_pkm", "Difference in mean quarterly total injuries per road-km between mid and early pre-COVID windows. Stage 2 trajectory-shape variable.",
  "log1p_total_injuries_pre_covid", "Log-transformed total pre-COVID injuries, log(1 + total_injuries_pre_covid).",
  "log1p_mean_total_inj_pre_covid", "Log-transformed mean quarterly pre-COVID injuries, log(1 + mean_total_inj_pre_covid).",
  "total_injuries_all_periods", "Total injuries observed in the OA across the full available injury dataset.",
  "ever_injury_OA", "Indicator equal to 1 if the OA had at least one recorded injury in the relevant pre-treatment period.",
  "zero_injury_OA", "Indicator equal to 1 if the OA had no recorded injuries in the relevant pre-treatment period."
)

cat("\nVariable descriptions:\n")
print(var_description, n = Inf)

#
