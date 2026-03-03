# ============================================================
# OA-Level Analysis from Injuries
# ============================================================

library(tidyverse)
library(sf)
library(here)

# -----------------------------
# 1️⃣ Load Data
# -----------------------------

# Injuries matched to roads
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))

# Road attributes (with geometry)
road_attributes <- readRDS(here("data", "processed", "road_attributes.rds"))

# Road → CAZ mapping
road_caz_prop <- readRDS(here("data", "processed", "roads_caz_props.rds")) %>%
  mutate(ever_treated = 1)

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
# 2️⃣ Attach LADs to OAs
# -----------------------------
oa_2011 <- st_join(
  oa_2011,
  lads_sub %>% select(LAD24CD),
  join = st_intersects,
  largest = TRUE
) %>%
  filter(!is.na(LAD24CD))

# -----------------------------
# 3️⃣ Map injuries → roads → OAs
# -----------------------------
# Make sure injuries has geometry (point or road line)
injuries_sf <- injuries %>%
  st_as_sf(coords = c("x_coord", "y_coord"), crs = 27700)  # replace with your point columns

# Spatial join injuries → OAs
injuries_oa <- st_join(
  injuries_sf,
  oa_2011 %>% select(OA_CODE, LAD24CD),
  join = st_intersects,
  left = FALSE
)

# Map injuries to OAs via spatial join
injuries_oa <- st_join(
  injuries,
  oa_2011 %>% select(OA_CODE, LAD24CD),
  join = st_intersects,
  left = FALSE
)

injuries_oa$OA
# -----------------------------
# 4️⃣ Aggregate injuries to OA level
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
# Assign OA Groups (as in OA_analysis)
# -----------------------------
# Treated OAs (intersect CAZ)
treated_OAs <- st_join(
  oa_2011,
  caz_boundaries,
  join = st_intersects,
  left = FALSE
) %>%
  pull(OA_CODE) %>% unique()

# Buffer OAs (1 km outside CAZ)
caz_buffer <- st_buffer(caz_boundaries, 1000) %>%
  st_difference(caz_boundaries)

buffer_OAs <- oa_2011 %>%
  filter(!OA_CODE %in% treated_OAs) %>%
  st_join(caz_buffer, join = st_intersects, left = FALSE) %>%
  pull(OA_CODE) %>% unique()

# Treated LADs for control_group1
treated_LADs <- oa_2011 %>%
  filter(OA_CODE %in% treated_OAs) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

# Initialize OA classification
OA_classification1 <- oa_injuries %>%
  mutate(
    treated_OA = if_else(OA_CODE.x %in% treated_OAs, 1, 0),
    buffer_OA = if_else(OA_CODE.x %in% buffer_OAs, 1, 0),
    control_group1_OA = if_else(
      LAD24CD.x %in% treated_LADs & treated_OA == 0 & buffer_OA == 0, 1, 0
    )
  )

# -----------------------------
# 6️⃣ Control group 2 (non-CAZ city centres)
# -----------------------------
OA_classification_sf <- oa_2011 %>%
  filter(!OA_CODE %in% treated_OAs & !OA_CODE %in% buffer_OAs) %>%
  select(OA_CODE, LAD24CD, geometry)

# City centroids
city_centroids <- OA_classification_sf %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_centroid()

# 1 km buffer
city_buffers <- st_buffer(city_centroids, 1000)

# OAs intersecting city-centres
city_centre_OAs <- st_join(OA_classification_sf, city_buffers, join = st_intersects, left = FALSE) %>%
  pull(OA_CODE)

# Assign control_group2
OA_classification1 <- OA_classification1 %>%
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
# 7️⃣ Restrict to relevant OAs
# -----------------------------
OA_analysis2 <- OA_classification1 %>%
  filter(treated_OA == 1 | buffer_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1)

# Check mutual exclusivity
OA_analysis2 %>%
  mutate(group_count = treated_OA + buffer_OA + control_group1_OA + control_group2_OA) %>%
  count(group_count)

# -----------------------------
# 8️⃣ Save OA-level dataset
# -----------------------------
saveRDS(OA_analysis2, here("data", "processed", "OA_level_from_injuries.rds"))





common_OAs <- intersect(OA_analysis$OA_CODE, OA_analysis2$OA_CODE.x)
length(common_OAs)




