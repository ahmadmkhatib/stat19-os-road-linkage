# ============================================================
# OA-Level Analysis from Injuries
# ============================================================

library(tidyverse)
library(sf)
library(here)


# Injuries matched to roads (points)
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))

# Output Areas
oa_2011 <- st_read("../stat19-os-road-linkage-data/infuse_oa_lyr_2011_clipped.shp") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  rename(OA_CODE = geo_code)

# LADs subset
lads_sub <- readRDS(here("data", "processed", "LADs_sub.rds"))

# CAZ polygons
caz_boundaries <- st_read(here("data", "processed", "shp_files", "CAZ_areas.shp")) %>%
  st_transform(st_crs(oa_2011)) %>%
  st_make_valid()

# -----------------------------
# Attach LADs to OAs
# -----------------------------
oa_2011 <- st_join(
  oa_2011,
  lads_sub %>% select(LAD24CD),
  join = st_intersects,
  largest = TRUE
) %>%
  filter(!is.na(LAD24CD))

# -----------------------------
#Map injuries → OAs
# -----------------------------
injuries_sf <- injuries %>%
  st_as_sf(coords = c("x_coord", "y_coord"), crs = 27700)

injuries_oa <- st_join(
  injuries_sf,
  oa_2011 %>% select(OA_CODE, LAD24CD),
  join = st_intersects,
  left = FALSE
)

# -----------------------------
# Aggregate injuries to OA level
# -----------------------------
oa_injuries <- injuries_oa %>%
  st_drop_geometry() %>%
  group_by(OA_CODE.x, LAD24CD.x) %>%
  summarise(
    KSI_adj = sum(KSI_adj, na.rm = TRUE),
    Slight_adj = sum(Slight_adj, na.rm = TRUE),
    KSI_unadj = sum(KSI_unadj, na.rm = TRUE),
    Slight_unadj = sum(Slight_unadj, na.rm = TRUE),
    total_inj_adj = sum(KSI_adj + Slight_adj, na.rm = TRUE),
    total_inj_unadj = sum(KSI_unadj + Slight_unadj, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Assign OA Groups (same logic as original OA_analysis)
# -----------------------------

# Treated OAs
treated_OAs <- st_join(
  oa_2011,
  caz_boundaries,
  join = st_intersects,
  left = FALSE
) %>%
  pull(OA_CODE) %>% unique()

# Spillover buffer OAs (1 km outside CAZ)
caz_buffer <- st_buffer(caz_boundaries, 1000) %>%
  st_difference(caz_boundaries)

buffer_OAs <- oa_2011 %>%
  filter(!OA_CODE %in% treated_OAs) %>%
  st_join(caz_buffer, join = st_intersects, left = FALSE) %>%
  pull(OA_CODE) %>% unique()

# Treated LADs (for same-city controls)
treated_LADs <- oa_2011 %>%
  filter(OA_CODE %in% treated_OAs) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

# -----------------------------
# Initialize OA classification
# -----------------------------
OA_classification <- oa_injuries %>%
  mutate(
    treated_OA = if_else(OA_CODE.x %in% treated_OAs, 1, 0),
    buffer_OA = if_else(OA_CODE.x %in% buffer_OAs, 1, 0),
    control_group1_OA = if_else(
      LAD24CD.x %in% treated_LADs & treated_OA == 0 & buffer_OA == 0, 1, 0
    )
  )

# -----------------------------
# Control Group 2 (non-CAZ city centres)
# -----------------------------
# OAs not treated or buffer
OA_sf_for_centres <- oa_2011 %>%
  filter(!OA_CODE %in% treated_OAs & !OA_CODE %in% buffer_OAs) %>%
  select(OA_CODE, LAD24CD, geometry)

# City centroids & 1 km buffer
city_centroids <- OA_sf_for_centres %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid()

city_buffers <- st_buffer(city_centroids, 1000)

# OAs intersecting city-centres
city_centre_OAs <- st_join(OA_sf_for_centres, city_buffers, join = st_intersects, left = FALSE) %>%
  pull(OA_CODE)

# Assign control_group2
OA_classification <- OA_classification %>%
  mutate(
    control_group2_OA = if_else(
      OA_CODE.x %in% city_centre_OAs &
        treated_OA == 0 &
        buffer_OA == 0 &
        control_group1_OA == 0,
      1, 0
    )
  )

# -----------------------------
# Restrict to relevant OAs
# -----------------------------
OA_analysis_injuries <- OA_classification %>%
  filter(treated_OA == 1 | buffer_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1)

# Check mutual exclusivity
OA_analysis_injuries %>%
  mutate(group_count = treated_OA + buffer_OA + control_group1_OA + control_group2_OA) %>%
  count(group_count)

##############
# -----------------------------
saveRDS(OA_analysis_injuries, here("data", "processed", "OA_level_from_injuries.rds"))

#  overlap with original OA_analysis
common_OAs <- intersect(OA_analysis$OA_CODE, OA_analysis_injuries$OA_CODE.x)
length(common_OAs)
