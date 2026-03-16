



# ============================================================
# Create OA-Level Matching Dataset for Synthetic Control Analysis
# ============================================================

# This script prepares a dataset of Output Areas (OAs) for use in
# constructing synthetic control units in the CAZ evaluation.
# . Restrict to pre-treatment period (up to 2020-12-31).
# . Calculate baseline injury levels per OA (mean injuries by type/severity).
# . Estimate pre-treatment trends in injuries for each OA (slope per quarter).
# . Combine OA characteristics, road characteristics, baseline injuries, and trends
#    into a single OA-level matching dataset for synthetic control construction.

## output: OA_matching_dataset.rds
#=========================================================
library(tidyverse)
library(lubridate)
library(here)
library(sf)


# OA characteristics
OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")    
)

glimpse(OA_analysis)

# road characteristics aggregated to OA
OA_roads <- readRDS(
  here("data","processed","OA_roads.rds")                 
)

glimpse(OA_roads)

# injuries aggregated by OA + quarter
OA_injuries <- readRDS(
  here("data","processed","OA_injuries_quarterly.rds")     
) 

# OA subset
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()


caz <- st_read(here("data","processed","shp_files","caz.shp"), quiet = TRUE)

# ── Build scheme → CAZ start date lookup from caz shapefile ───────────────────
# startDt can have multiple rows per scheme (multiple polygons); take the min
# date per scheme as the implementation date

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(
    caz_start_date = dmy(startDt)           # "26/09/2022" → Date
  ) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    .groups = "drop"
  )

# Quick check
print(caz_dates)

# ── Join CAZ start date onto OA_roads via scheme column ───────────────────────
# OA_roads already has a scheme column for treated OAs; controls will have NA
OA_roads <- OA_roads %>%
  left_join(caz_dates, by = "scheme")
# caz_start_date = NA for control OAs (no scheme) — this is correct

# ── Join into OA_injuries panel ───────────────────────────────────────────────
# We need per-OA: scheme, caz_start_date, treated_OA flag
oa_scheme_lookup <- OA_roads %>%
  select(OA, scheme, caz_start_date, treated_OA, control_group2_OA)

OA_injuries <- OA_injuries %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  mutate(
    quarter_date = as.Date(zoo::as.yearqtr(quarter_year))
  )

# ── Create pre/post treatment indicator ───────────────────────────────────────
# Treated OAs:  pre = before caz_start_date, post = on/after caz_start_date
# Control OAs:  pre = all quarters, post = NA
OA_injuries <- OA_injuries %>%
  mutate(
    period = case_when(
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_date <  caz_start_date ~ "pre",
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_date >= caz_start_date ~ "post",
      control_group2_OA == 1          ~ "pre",   # controls: all pre
      TRUE                            ~ NA_character_
    )
  )

# Sanity check — no treated OA should have NA period
OA_injuries %>%
  filter(treated_OA == 1) %>%
  count(period, is.na(period))

# ── [EXISTING] Baseline injury levels — computed over each OA's own pre-period ─
inj_pre <- OA_injuries %>%
  filter(period == "pre")

inj_baseline <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    mean_car_KSI    = mean(Car_Van_KSI,        na.rm = TRUE),
    mean_car_slight = mean(Car_Van_Slight,      na.rm = TRUE),
    mean_cyc_KSI    = mean(Cyclist_KSI,         na.rm = TRUE),
    mean_cyc_slight = mean(Cyclist_Slight,      na.rm = TRUE),
    mean_ped_KSI    = mean(Pedestrian_KSI,      na.rm = TRUE),
    mean_ped_slight = mean(Pedestrian_Slight,   na.rm = TRUE),
    mean_total      = mean(total_injuries,      na.rm = TRUE),
    .groups = "drop"
  )

# ── [EXISTING] Pre-treatment trends (raw counts) ──────────────────────────────
# time index = row_number() within each OA's pre-period only
inj_pre <- inj_pre %>%
  arrange(OA, quarter_date) %>%
  group_by(OA) %>%
  mutate(time = row_number()) %>%
  ungroup()

inj_trends <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI    = lm(.data[["Car_Van_KSI"]]       ~ time)$coefficients[2],
    trend_car_slight = lm(.data[["Car_Van_Slight"]]     ~ time)$coefficients[2],
    trend_cyc_KSI    = lm(.data[["Cyclist_KSI"]]        ~ time)$coefficients[2],
    trend_cyc_slight = lm(.data[["Cyclist_Slight"]]     ~ time)$coefficients[2],
    trend_ped_KSI    = lm(.data[["Pedestrian_KSI"]]     ~ time)$coefficients[2],
    trend_ped_slight = lm(.data[["Pedestrian_Slight"]]  ~ time)$coefficients[2],
    trend_total      = lm(.data[["total_injuries"]]     ~ time)$coefficients[2],
    .groups = "drop"
  )

# ── [NEW 1] Per-km baselines ───────────────────────────────────────────────────
road_lengths <- OA_roads %>%
  select(OA, total_road_length) %>%
  mutate(road_length_km = total_road_length / 1000)

inj_per_km <- inj_baseline %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    road_length_km = if_else(
      is.na(road_length_km) | road_length_km == 0, NA_real_, road_length_km
    ),
    mean_car_KSI_pkm    = mean_car_KSI    / road_length_km,
    mean_car_slight_pkm = mean_car_slight / road_length_km,
    mean_cyc_KSI_pkm    = mean_cyc_KSI    / road_length_km,
    mean_cyc_slight_pkm = mean_cyc_slight / road_length_km,
    mean_ped_KSI_pkm    = mean_ped_KSI    / road_length_km,
    mean_ped_slight_pkm = mean_ped_slight / road_length_km,
    mean_total_pkm      = mean_total      / road_length_km
  ) %>%
  select(OA, ends_with("_pkm"))

# ── [NEW 2] Per-km trend slopes ───────────────────────────────────────────────
safe_slope <- function(y, x) {
  df <- data.frame(y = y, x = x) %>% filter(!is.na(y))
  if (nrow(df) < 2) return(NA_real_)
  lm(y ~ x, data = df)$coefficients[2]
}

inj_pre_pkm <- inj_pre %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    road_length_km = if_else(
      is.na(road_length_km) | road_length_km == 0, NA_real_, road_length_km
    ),
    car_KSI_pkm    = Car_Van_KSI       / road_length_km,
    car_slight_pkm = Car_Van_Slight    / road_length_km,
    cyc_KSI_pkm    = Cyclist_KSI       / road_length_km,
    cyc_slight_pkm = Cyclist_Slight    / road_length_km,
    ped_KSI_pkm    = Pedestrian_KSI    / road_length_km,
    ped_slight_pkm = Pedestrian_Slight / road_length_km,
    total_pkm      = total_injuries    / road_length_km
  )

inj_trends_pkm <- inj_pre_pkm %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI_pkm    = safe_slope(car_KSI_pkm,    time),
    trend_car_slight_pkm = safe_slope(car_slight_pkm,  time),
    trend_cyc_KSI_pkm    = safe_slope(cyc_KSI_pkm,    time),
    trend_cyc_slight_pkm = safe_slope(cyc_slight_pkm,  time),
    trend_ped_KSI_pkm    = safe_slope(ped_KSI_pkm,    time),
    trend_ped_slight_pkm = safe_slope(ped_slight_pkm,  time),
    trend_total_pkm      = safe_slope(total_pkm,       time),
    .groups = "drop"
  )

# ── [NEW 3] Road composition & network density ────────────────────────────────
oa_area <- oa_sub %>%
  st_make_valid() %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(OA, area_km2)

road_composition <- OA_roads %>%
  select(OA, total_road_length, n_A, n_B, n_minor, n_roads) %>%
  left_join(oa_area, by = "OA") %>%
  mutate(
    total_road_length  = if_else(
      is.na(total_road_length) | total_road_length == 0, NA_real_, total_road_length
    ),
    pct_A_road         = if_else(!is.na(n_roads) & n_roads > 0,
                                 100 * n_A     / n_roads, NA_real_),
    pct_B_road         = if_else(!is.na(n_roads) & n_roads > 0,
                                 100 * n_B     / n_roads, NA_real_),
    pct_minor_road     = if_else(!is.na(n_roads) & n_roads > 0,
                                 100 * n_minor / n_roads, NA_real_),
    road_density_m_km2 = if_else(!is.na(area_km2) & area_km2 > 0,
                                 total_road_length / area_km2, NA_real_)
  ) %>%
  select(OA, pct_A_road, pct_B_road, pct_minor_road, road_density_m_km2)

# ── Assemble full matching dataset ────────────────────────────────────────────
OA_analysis <- OA_analysis %>%
  left_join(
    OA_roads %>% select(OA, scheme, caz_start_date),
    by = "OA"
  )

OA_matching_dataset <- OA_analysis %>%
  left_join(OA_roads,         by = "OA") %>%
  left_join(inj_baseline,     by = "OA") %>%
  left_join(inj_trends,       by = "OA") %>%
  left_join(inj_per_km,       by = "OA") %>%
  left_join(inj_trends_pkm,   by = "OA") %>%
  left_join(road_composition, by = "OA")

# ── Column deduplication ──────────────────────────────────────────────────────
OA_matching_dataset <- OA_matching_dataset %>%
  rename(
    LAD24CD           = LAD24CD.x,
    LAD24NM           = LAD24NM.x,
    treated_OA        = treated_OA.x,
    buffer_OA         = buffer_OA.x,
    control_group1_OA = control_group1_OA.x,
    control_group2_OA = control_group2_OA.x,
    dist_citycentre   = dist_citycentre.x,
    assignment        = assignment.x,
    scheme            = scheme.x,
    caz_start_date    = caz_start_date.x
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# ── Replace NAs with 0 for injury variables only ──────────────────────────────
# Road composition NAs are kept — OA genuinely has no roads
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(across(starts_with("mean"),  ~replace_na(.x, 0))) %>%
  mutate(across(starts_with("trend"), ~replace_na(.x, 0)))

# ── Sanity checks ─────────────────────────────────────────────────────────────
# Every treated OA should have a non-NA caz_start_date
OA_matching_dataset %>%
  filter(treated_OA == 1) %>%
  summarise(missing_caz_date = sum(is.na(caz_start_date)))

# Distribution of CAZ start dates across treated OAs
OA_matching_dataset %>%
  filter(treated_OA == 1) %>%
  count(scheme, caz_start_date) %>%
  arrange(caz_start_date)

# Check pre-period length varies by scheme (expected given staggered adoption)
inj_pre %>%
  left_join(OA_matching_dataset %>% select(OA, scheme), by = "OA") %>%
  group_by(scheme) %>%
  summarise(
    n_pre_quarters = n_distinct(quarter_date),
    pre_start      = min(quarter_date),
    pre_end        = max(quarter_date)
  )

# ── Save ───────────────────────────────────────────────────────────────────────
saveRDS(OA_matching_dataset,
        here("data","processed","OA_matching_dataset.rds"))

# ── Shapefile export ───────────────────────────────────────────────────────────
OA_matching_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_dataset, by = "OA") %>%
  st_as_sf()

st_write(
  OA_matching_sf,
  here("data","processed","shp_files","OA_matching_dataset.gpkg"),
  delete_dsn = TRUE
)





# ============================================================
# Variable descriptions for OA_matching_dataset
# ============================================================

var_description <- c(
  
  OA = "Output Area identifier (ONS 2011 OA code). Each row represents one OA.",
  
  n_roads = "Total number of road segments intersecting the OA.",
  
  total_road_length = "Total length of roads within the OA (meters).",
  
  n_A = "Number of A-road segments within the OA.",
  
  n_B = "Number of B-road segments within the OA.",
  
  n_motorway = "Number of motorway segments within the OA.",
  
  n_minor = "Number of minor/local road segments within the OA.",
  
  mean_car_KSI = "Mean quarterly number of car/van KSI (Killed or Seriously Injured) casualties in the OA during the pre-treatment period.",
  
  mean_car_slight = "Mean quarterly number of car/van slight injury casualties in the OA during the pre-treatment period.",
  
  mean_cyc_KSI = "Mean quarterly number of cyclist KSI casualties in the OA during the pre-treatment period.",
  
  mean_cyc_slight = "Mean quarterly number of cyclist slight injury casualties in the OA during the pre-treatment period.",
  
  mean_ped_KSI = "Mean quarterly number of pedestrian KSI casualties in the OA during the pre-treatment period.",
  
  mean_ped_slight = "Mean quarterly number of pedestrian slight injury casualties in the OA during the pre-treatment period.",
  
  mean_total = "Mean quarterly total number of casualties (all modes and severities) in the OA during the pre-treatment period.",
  
  trend_car_KSI = "Linear pre-treatment trend (slope per quarter) in car/van KSI casualties in the OA.",
  
  trend_car_slight = "Linear pre-treatment trend (slope per quarter) in car/van slight casualties in the OA.",
  
  trend_cyc_KSI = "Linear pre-treatment trend (slope per quarter) in cyclist KSI casualties in the OA.",
  
  trend_cyc_slight = "Linear pre-treatment trend (slope per quarter) in cyclist slight casualties in the OA.",
  
  trend_ped_KSI = "Linear pre-treatment trend (slope per quarter) in pedestrian KSI casualties in the OA.",
  
  trend_ped_slight = "Linear pre-treatment trend (slope per quarter) in pedestrian slight casualties in the OA.",
  
  trend_total = "Linear pre-treatment trend (slope per quarter) in total casualties in the OA."
)

# table 
var_description_table <- tibble(
  variable = names(var_description),
  description = unname(var_description)
)

print(var_description_table, n=30)


