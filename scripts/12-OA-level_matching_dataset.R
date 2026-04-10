# ============================================================
# Create OA-Level Matching Dataset for Controls pool for the DiD Analysis
# ============================================================
# This script prepares a dataset of Output Areas (OAs) for use in
# constructing synthetic control units in the CAZ evaluation.
#
# Pre-treatment period:
#   - Treated OAs: all quarters before the CAZ implementation date.
#   - Control OAs: all available quarters (no treatment cutoff).
#   Pre-treatment periods vary in length across schemes due to
#   staggered CAZ adoption — this is expected and handled correctly.
#
# Zero-filling:
#   OA_injuries_quarterly only contains OA-quarter rows where at least
#   one injury occurred. Zero-injury quarters are absent from the raw data.
#   A balanced panel is constructed by cross-joining all OAs × all quarters
#   and filling missing counts with zero BEFORE computing baselines or trends.
#   Without this step, trend slopes and baseline means would be estimated on
#   non-zero quarters only, biasing both upward.
#
# Baseline injury levels:
#   Mean quarterly casualties per OA over the pre-treatment period for seven
#   outcomes: car/van KSI, car/van slight, cyclist KSI, cyclist slight,
#   pedestrian KSI, pedestrian slight, and total casualties.
#
# Pre-treatment injury trends:
#   Quarterly trend slopes estimated via quasi-Poisson GLM. Raw count trends
#   use a time index as the predictor. Rate-per-km trends use the same GLM
#   with log(road_length_km) as an offset — keeping raw integer counts as the
#   outcome rather than dividing by km (which would misuse Poisson on a ratio).
#   The slope coefficient (log-rate change per quarter) is extracted.
#   Returns NA where all counts are zero or fewer than 2 observations exist.
#
# Outputs assembled into OA_matching_dataset:
#   - OA characteristics (from OA_level_from_polygons)
#   - Road network characteristics (from OA_roads)
#   - Baseline injury levels (mean per quarter, pre-period, zero-filled)
#   - Pre-treatment trend slopes (quasi-Poisson GLM, raw counts)
#   - Baseline injury rates per road-km
#   - Pre-treatment trend slopes (quasi-Poisson GLM with offset, per-km)
#   - Road composition & network density
#
# Output: data/processed/OA_matching_dataset.rds
#         data/processed/shp_files/OA_matching_dataset.gpkg
# ============================================================

library(tidyverse)
library(lubridate)
library(here)
library(sf)

# ── Load data ─────────────────────────────────────────────────────────────────

OA_analysis <- readRDS(
  here("data", "processed", "OA_level_from_polygons.rds")             # from 8 
)
names(OA_analysis)
table(OA_analysis$assignment)

# Deduplicate OA_analysis — keep first row per OA 
OA_analysis <- OA_analysis %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE)

glimpse(OA_analysis)
table(OA_analysis$scheme, useNA = "always")
table(OA_analysis$treated_OA, useNA = "always")

OA_roads <- readRDS(
  here("data", "processed", "OA_roads.rds")
)
glimpse(OA_roads)



OA_injuries <- readRDS(
  here("data", "processed", "OA_injuries_quarterly.rds")
)

oa_sub <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()

caz <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

# ── Build scheme → CAZ start date lookup ─────────────────────────────────────

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(caz_start_date = dmy(startDt)) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    .groups = "drop"
  )
print(caz_dates)

# ── Clean OA_roads — keep only road-specific columns 

OA_roads %>% count(OA) %>% filter(n > 1) %>% nrow()
# and
OA_roads %>% count(OA) %>% summarise(max(n), median(n))


# OA_roads has 69,318 — only OAs that have at least one road
# OA_roads_clean right_joins to oa_sub to recover the 1,752 road-free OAs

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
    oa_sub %>% st_drop_geometry() %>% select(OA),
    by = "OA"
  ) %>%
  mutate(across(
    c(n_roads, total_road_length, n_A, n_B, n_motorway, n_minor),
    ~ replace_na(.x, 0)
  ))


nrow(OA_roads_clean)                                          # must be 71,070
OA_roads_clean %>% count(OA) %>% filter(n > 1) %>% nrow()    # must be 0
n_distinct(OA_roads_clean %>% filter(n_roads == 0) %>% pull(OA))  # should be 1,752

anti_join(
  oa_sub %>% st_drop_geometry() %>% dplyr::select(OA),
  OA_roads_clean,
  by = "OA"
) %>% nrow()


no_road_oas <- OA_roads %>%
  filter(n_roads == 0) %>%
  pull(OA)
nrow(no_road_oas)
n_distinct(no_road_oas)  # 

has_injuries_no_roads <- OA_injuries %>%
  filter(OA %in% no_road_oas) %>%
  distinct(OA)
nrow(has_injuries_no_roads)
 
## are these genuinly no roads OAs ??
 OA_roads_clean %>%
   filter(n_roads == 0) %>%
   summarise(
     min_len = min(total_road_length),
     max_len = max(total_road_length),
     mean_len = mean(total_road_length)
   )
 
# ── Build scheme/treatment lookup ─────────────────────────────────────────────

 oa_scheme_lookup <- OA_analysis %>%
   arrange(OA, scheme) %>%
   distinct(OA, .keep_all = TRUE) %>%   # 
   dplyr::select(OA, scheme, treated_OA, control_group2_OA) %>%
   left_join(caz_dates, by = "scheme")
 
 # V
 oa_scheme_lookup %>% count(OA) %>% filter(n > 1) %>% nrow()
 
 
 

# How many injuries do these 95 OAs actually have?
OA_injuries %>%
  filter(OA %in% no_road_oas) %>%
  summarise(
    n_OAs        = n_distinct(OA),
    total_inj    = sum(total_injuries, na.rm = TRUE),
    mean_inj     = mean(total_injuries, na.rm = TRUE)
  )

# Are they concentrated in treated/buffer areas?
OA_roads_clean %>%
  filter(n_roads == 0) %>%
  left_join(OA_analysis %>% select(OA, assignment), by = "OA") %>%
  count(assignment)


nrow(oa_sub)           
nrow(OA_roads_clean)   



# ── Attach treatment info to injury panel ────────────────────────────────────
# Note: OA_injuries only has rows where injuries > 0 at this stage

OA_injuries <- OA_injuries %>%
  left_join(oa_scheme_lookup, by = "OA")


# ── Build balanced panel — fill zero-injury OA-quarters explicitly ────────────
# OA_injuries contains only quarters with at least one injury.
# Cross-joining all OAs × all quarters and filling NAs with 0 ensures:
#   (a) baseline means include zero-injury quarters in the denominator
#   (b) trend slopes are estimated over actual calendar time, not just
#       quarters where something happened

all_oas      <- oa_scheme_lookup %>%dplyr:: select(OA)   # all OAs, not just injured
all_quarters <- OA_injuries %>% distinct(quarter_year)

OA_injuries_balanced <- all_oas %>%
  cross_join(all_quarters) %>%
  left_join(OA_injuries, by = c("OA", "quarter_year")) %>%
  mutate(across(
    c(Car_Van_KSI, Car_Van_Slight, Cyclist_KSI, Cyclist_Slight,
      Pedestrian_KSI, Pedestrian_Slight, Other_KSI, Other_Slight,
      total_injuries),
    ~ replace_na(.x, 0)
  )) %>%
  # scheme/treatment info dropped by cross_join — re-attach
dplyr::   select(-any_of(c("scheme", "treated_OA", "control_group2_OA", "caz_start_date"))) %>%
  left_join(oa_scheme_lookup, by = "OA")

# Verify balance
cat("OAs:     ", n_distinct(OA_injuries_balanced$OA), "\n")
cat("Quarters:", n_distinct(OA_injuries_balanced$quarter_year), "\n")
cat("Rows:    ", nrow(OA_injuries_balanced), "\n")
cat("Expected:", n_distinct(OA_injuries_balanced$OA) *
      n_distinct(OA_injuries_balanced$quarter_year), "\n")

# ── Pre/post treatment indicator ──────────────────────────────────────────────

OA_injuries_balanced <- OA_injuries_balanced %>%
  mutate(
    period = case_when(
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_year <  caz_start_date ~ "pre",
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_year >= caz_start_date ~ "post",
      control_group2_OA == 1          ~ "pre",
      TRUE                            ~ NA_character_
    )
  )

OA_injuries_balanced %>%
  filter(treated_OA == 1) %>%
  count(period)

# ── Pre-treatment subset with calendar time index ─────────────────────────────

inj_pre <- OA_injuries_balanced %>%
  filter(period == "pre") %>%
  arrange(OA, quarter_year) %>%
  group_by(OA) %>%
  mutate(time = row_number()) %>%
  ungroup()

# ── Road lengths lookup ───────────────────────────────────────────────────────

road_lengths <- OA_roads_clean %>%
 dplyr:: select(OA, total_road_length) %>%
  mutate(
    road_length_km = if_else(
      is.na(total_road_length) | total_road_length == 0,
      NA_real_,
      total_road_length / 1000
    )
  )

# ── Quasi-Poisson GLM slope extractor ────────────────────────────────────────
# Returns log-rate change per quarter from a quasi-Poisson GLM.
# offset argument accepts log(road_length_km) for exposure-adjusted models.
# Returns NA if: fewer than 2 obs, all zeros, or constant non-zero counts.

poisson_slope <- function(y, x, offset = NULL) {
  
  df <- tibble(y = y, x = x,
               off = if (is.null(offset)) rep(0, length(x)) else offset) %>%
    filter(!is.na(y), !is.na(off))
  
  if (nrow(df) < 2)   return(NA_real_)
  if (all(df$y == 0)) return(NA_real_)
  if (var(df$y) == 0) return(NA_real_)
  
  model <- tryCatch(
    {
      if (is.null(offset)) {
        glm(y ~ x, family = quasipoisson(link = "log"), data = df)
      } else {
        glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
      }
    },
    warning = function(w) {
      if (grepl("did not converge", conditionMessage(w))) NULL else {
        withCallingHandlers(
          if (is.null(offset)) {
            glm(y ~ x, family = quasipoisson(link = "log"), data = df)
          } else {
            glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
          },
          warning = function(w2) invokeRestart("muffleWarning")
        )
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NA_real_)
  coef(model)[["x"]]
}

# ── Pre-treatment trends — raw counts ────────────────────────────────────────

inj_trends <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI    = poisson_slope(Car_Van_KSI,       time),
    trend_car_slight = poisson_slope(Car_Van_Slight,    time),
    trend_cyc_KSI    = poisson_slope(Cyclist_KSI,       time),
    trend_cyc_slight = poisson_slope(Cyclist_Slight,    time),
    trend_ped_KSI    = poisson_slope(Pedestrian_KSI,    time),
    trend_ped_slight = poisson_slope(Pedestrian_Slight, time),
    trend_total      = poisson_slope(total_injuries,    time),
    .groups = "drop"
  )

# ── Baseline injury levels (mean quarterly, pre-period, zero-filled) ──────────

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

# ── Per-km baselines ──────────────────────────────────────────────────────────

inj_per_km <- inj_baseline %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    mean_car_KSI_pkm    = mean_car_KSI    / road_length_km,
    mean_car_slight_pkm = mean_car_slight / road_length_km,
    mean_cyc_KSI_pkm    = mean_cyc_KSI    / road_length_km,
    mean_cyc_slight_pkm = mean_cyc_slight / road_length_km,
    mean_ped_KSI_pkm    = mean_ped_KSI    / road_length_km,
    mean_ped_slight_pkm = mean_ped_slight / road_length_km,
    mean_total_pkm      = mean_total      / road_length_km
  ) %>%
dplyr::  select(OA, ends_with("_pkm"))

# ── Per-km trend slopes — quasi-Poisson GLM with offset ──────────────────────

inj_pre_offset <- inj_pre %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    log_road_km = if_else(
      !is.na(road_length_km) & road_length_km > 0,
      log(road_length_km),
      NA_real_
    )
  )

inj_trends_pkm <- inj_pre_offset %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI_pkm    = poisson_slope(Car_Van_KSI,       time, log_road_km),
    trend_car_slight_pkm = poisson_slope(Car_Van_Slight,    time, log_road_km),
    trend_cyc_KSI_pkm    = poisson_slope(Cyclist_KSI,       time, log_road_km),
    trend_cyc_slight_pkm = poisson_slope(Cyclist_Slight,    time, log_road_km),
    trend_ped_KSI_pkm    = poisson_slope(Pedestrian_KSI,    time, log_road_km),
    trend_ped_slight_pkm = poisson_slope(Pedestrian_Slight, time, log_road_km),
    trend_total_pkm      = poisson_slope(total_injuries,    time, log_road_km),
    .groups = "drop"
  )

# ── Road composition & network density ───────────────────────────────────────

oa_area <- oa_sub %>%
  st_make_valid() %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(OA, area_km2)

road_composition <- OA_roads_clean %>%
  select(OA, total_road_length, n_A, n_B, n_minor, n_roads) %>%
  left_join(oa_area, by = "OA") %>%
  mutate(
    total_road_length  = if_else(
      is.na(total_road_length) | total_road_length == 0,
      NA_real_, total_road_length),
    pct_A_road         = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_A     / n_roads, NA_real_),
    pct_B_road         = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_B     / n_roads, NA_real_),
    pct_minor_road     = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_minor / n_roads, NA_real_),
    road_density_m_km2 = if_else(
      !is.na(area_km2) & area_km2 > 0,
      total_road_length / area_km2, NA_real_)
  ) %>%
dplyr::  select(OA, pct_A_road, pct_B_road, pct_minor_road, road_density_m_km2)

# ── Assemble OA matching dataset ─────────────────────────────────────────────
# OA_analysis is the base — all other tables are one-row-per-OA so no
# row inflation can occur. OA_roads_clean has been stripped of any columns
# already in OA_analysis to prevent .x/.y suffix collisions.

OA_matching_dataset <- OA_analysis %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE) %>%   # ← add this
  left_join(OA_roads_clean,   by = "OA") %>%
  left_join(inj_baseline,     by = "OA") %>%
  left_join(inj_trends,       by = "OA") %>%
  left_join(inj_per_km,       by = "OA") %>%
  left_join(inj_trends_pkm,   by = "OA") %>%
  left_join(road_composition, by = "OA") %>%
  mutate(road_length_km = total_road_length / 1000)

# Replace NAs with 0 for injury count variables only
# _pkm NAs kept — undefined where OA has no roads
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(across(starts_with("mean_"),  ~ replace_na(.x, 0))) %>%
  mutate(across(starts_with("trend_"), ~ replace_na(.x, 0)))

# ── Checks ────────────────────────────────────────────────────────────────────

# 1. Unit counts by assignment group
OA_matching_dataset %>%
  count(treated_OA, control_group2_OA, buffer_OA) %>%
  print()

# 2. No OAs missing from oa_sub
anti_join(
  oa_sub %>% st_drop_geometry() %>% select(OA),
  OA_matching_dataset, by = "OA"
) %>% nrow()

# 3. Pre-period length by scheme — now using balanced panel so all OAs
#    should have the same number of quarters within each scheme
inj_pre %>%
  group_by(scheme, treated_OA) %>%
  summarise(
    n_OAs           = n_distinct(OA),
    min_quarters    = min(table(OA)),
    max_quarters    = max(table(OA)),
    median_quarters = median(as.numeric(table(OA))),
    .groups = "drop"
  ) %>%
  print()

# 4. Short pre-period flag (should be 0 or very few after zero-filling)
short_pre <- inj_pre %>%
  filter(treated_OA == 1) %>%
  count(OA, name = "n_pre_quarters") %>%
  filter(n_pre_quarters < 4)
cat("OAs with fewer than 4 pre-treatment quarters:", nrow(short_pre), "\n")

# 5. Zero inflation in baseline injury variables
inj_baseline %>%
  left_join(OA_analysis %>% select(OA, treated_OA, control_group2_OA), by = "OA") %>%
  group_by(treated_OA, control_group2_OA) %>%
  summarise(
    pct_zero_car_KSI = mean(mean_car_KSI  == 0) * 100,
    pct_zero_cyc_KSI = mean(mean_cyc_KSI  == 0) * 100,
    pct_zero_ped_KSI = mean(mean_ped_KSI  == 0) * 100,
    pct_zero_total   = mean(mean_total    == 0) * 100,
    .groups = "drop"
  ) %>%
  print()

# 6. NAs remaining in matching variables for treated OAs
OA_matching_dataset %>%
  filter(treated_OA == 1) %>%
  summarise(across(
    c(starts_with("trend_"), starts_with("mean_"), ends_with("_pkm")),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  filter(n_NA > 0) %>%
  arrange(desc(n_NA)) %>%
  print()

# 7. Parallel trends plot — uses inj_pre directly (already has treated_OA etc.)
trend_plot <- inj_pre %>%
  mutate(group = case_when(
    treated_OA        == 1 ~ "Treated",
    control_group2_OA == 1 ~ "Control",
    TRUE                   ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(group, quarter_year) %>%
  summarise(mean_total = mean(total_injuries, na.rm = TRUE), .groups = "drop")

ggplot(trend_plot, aes(x = quarter_year, y = mean_total, colour = group)) +
  geom_line() +
  geom_point(size = 1) +
  labs(
    title = "Pre-treatment parallel trends check",
    x = "Quarter", y = "Mean total casualties per OA", colour = NULL
  ) +
  theme_minimal()

# 8. Baseline balance — treated vs control
OA_matching_dataset %>%
  filter(treated_OA == 1 | control_group2_OA == 1) %>%
  mutate(group = if_else(treated_OA == 1, "Treated", "Control")) %>%
  group_by(group) %>%
  summarise(
    n_OAs             = n(),
    mean_total_inj    = mean(mean_total,         na.rm = TRUE),
    mean_road_km      = mean(road_length_km,     na.rm = TRUE),
    mean_road_density = mean(road_density_m_km2, na.rm = TRUE),
    mean_pct_A        = mean(pct_A_road,         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# 9. Duplicate OA check — must be zero
dupes <- OA_matching_dataset %>% count(OA) %>% filter(n > 1)
cat("Duplicate OAs:", nrow(dupes), "\n")
stopifnot(nrow(dupes) == 0)

# 10. Balanced panel verification
cat("Balanced panel rows:    ", nrow(inj_pre), "\n")
cat("Expected (n_OA × n_qtr):",
    n_distinct(inj_pre$OA) * n_distinct(inj_pre$quarter_year), "\n")



# Which OAs are in oa_scheme_lookup but never appeared in OA_injuries?
missing_oas <- anti_join(
  oa_scheme_lookup %>% select(OA),
  OA_injuries %>% select(OA),
  by = "OA"
)
nrow(missing_oas)
## flag no injuries OAs
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(zero_injury_OA = if_else(
    OA %in% (missing_oas %>% pull(OA)), 1L, 0L
  ))

cat("OAs never in OA_injuries:", nrow(missing_oas), "\n")
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  count(treated_OA, control_group2_OA)

# Are any of these treated?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  count(treated_OA, control_group2_OA)

# Where are these zero-injury treated OAs geographically?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  count(scheme) %>%
  arrange(desc(n))

# Do the OA codes look the same format?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  pull(OA) %>%
  head(10)

OA_injuries %>% pull(OA) %>% head(10)


missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  left_join(OA_roads_clean, by = "OA") %>%
  summarise(
    n_no_roads       = sum(is.na(n_roads) | n_roads == 0),
    n_has_roads      = sum(!is.na(n_roads) & n_roads > 0),
    mean_road_length = mean(total_road_length, na.rm = TRUE)
  )

## — these are CAZ OAs that happen to have no recorded STATS19 casualties across the whole study period.

#flag to excluded them from the matching pool — they cannot serve as useful donors or treated units 

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(zero_injury_OA = if_else(
    OA %in% (missing_oas %>% pull(OA)), 1L, 0L
  ))

# Summary of what gets flagged
OA_matching_dataset %>%
  count(treated_OA, control_group2_OA, zero_injury_OA) %>%
  arrange(treated_OA, control_group2_OA, zero_injury_OA)

# Effective analysis sample after flagging
OA_matching_dataset %>%
  filter(zero_injury_OA == 0) %>%
  count(treated_OA, control_group2_OA, buffer_OA)



miss_var_summary(OA_matching_dataset)


# replace the no roads OA with zeros in the 1752 NAs

OA_matching_dataset <- OA_matching_dataset    %>%   
  mutate(across(
    c(pct_A_road, pct_B_road, pct_minor_road, road_density_m_km2),
    ~ replace_na(.x, 0)
  ))



# ── Save ──────────────────────────────────────────────────────────────────────

saveRDS(
  OA_matching_dataset,
  here("data", "processed", "OA_matching_dataset.rds")
)
OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_dataset.rds"))



######################################################################################
#######################################################################################
####################################################################################
##################################################################################
##################################################################################
###########################################################################
# ── Variable descriptions ─────────────────────────────────────────────────────

var_description <- c(
  OA                   = "Output Area identifier (ONS 2011 OA code). Each row is one OA.",
  n_roads              = "Total number of road segments intersecting the OA.",
  total_road_length    = "Total length of roads within the OA (metres).",
  road_length_km       = "Total length of roads within the OA (kilometres).",
  n_A                  = "Number of A-road segments within the OA.",
  n_B                  = "Number of B-road segments within the OA.",
  n_motorway           = "Number of motorway segments within the OA.",
  n_minor              = "Number of minor/local road segments within the OA.",
  mean_car_KSI         = "Mean quarterly car/van KSI casualties, pre-treatment period (zero-filled).",
  mean_car_slight      = "Mean quarterly car/van slight casualties, pre-treatment period (zero-filled).",
  mean_cyc_KSI         = "Mean quarterly cyclist KSI casualties, pre-treatment period (zero-filled).",
  mean_cyc_slight      = "Mean quarterly cyclist slight casualties, pre-treatment period (zero-filled).",
  mean_ped_KSI         = "Mean quarterly pedestrian KSI casualties, pre-treatment period (zero-filled).",
  mean_ped_slight      = "Mean quarterly pedestrian slight casualties, pre-treatment period (zero-filled).",
  mean_total           = "Mean quarterly total casualties (all modes/severities), pre-treatment (zero-filled).",
  trend_car_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for car/van KSI, pre-treatment.",
  trend_car_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for car/van slight, pre-treatment.",
  trend_cyc_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for cyclist KSI, pre-treatment.",
  trend_cyc_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for cyclist slight, pre-treatment.",
  trend_ped_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for pedestrian KSI, pre-treatment.",
  trend_ped_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for pedestrian slight, pre-treatment.",
  trend_total          = "Quasi-Poisson GLM slope (log-rate/quarter) for total casualties, pre-treatment.",
  mean_car_KSI_pkm     = "Mean quarterly car/van KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_car_slight_pkm  = "Mean quarterly car/van slight casualties per road-km, pre-treatment (zero-filled).",
  mean_cyc_KSI_pkm     = "Mean quarterly cyclist KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_cyc_slight_pkm  = "Mean quarterly cyclist slight casualties per road-km, pre-treatment (zero-filled).",
  mean_ped_KSI_pkm     = "Mean quarterly pedestrian KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_ped_slight_pkm  = "Mean quarterly pedestrian slight casualties per road-km, pre-treatment (zero-filled).",
  mean_total_pkm       = "Mean quarterly total casualties per road-km, pre-treatment (zero-filled).",
  trend_car_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for car/van KSI rate.",
  trend_car_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for car/van slight rate.",
  trend_cyc_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for cyclist KSI rate.",
  trend_cyc_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for cyclist slight rate.",
  trend_ped_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for pedestrian KSI rate.",
  trend_ped_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for pedestrian slight rate.",
  trend_total_pkm      = "Quasi-Poisson GLM slope with log(road_km) offset for total casualty rate.",
  pct_A_road           = "Percentage of road segments in the OA that are A-roads.",
  pct_B_road           = "Percentage of road segments in the OA that are B-roads.",
  pct_minor_road       = "Percentage of road segments in the OA that are minor roads.",
  road_density_m_km2   = "Total road length per km² of OA area (metres per km²).",
  zero_injury_OA       = "OA with no injury"
  
)

tibble(
  variable    = names(var_description),
  description = unname(var_description)
) %>%
  print(n = Inf)