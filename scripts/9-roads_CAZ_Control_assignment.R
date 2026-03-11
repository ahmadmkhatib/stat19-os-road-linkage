# ============================================================
# Assign Roads to Output Areas (OA)
# ============================================================
# This script prepares the road-level spatial data and assigns
# each road link to a single Output Area (OA).
#
# Assignment rule:
# A road is assigned to the OA containing the largest share
# of its length.
#
# also calculates the proportion of each road
# inside CAZ polygons,  stored for reference.
# + OA-level treatment classification 
# ============================================================

library(tidyverse)
library(sf)
library(here)
library(zoo)

# ============================================================
# Load Data
# ============================================================

# CAZ polygons
caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_make_valid()

# LAD subset
lads_sub <- readRDS(
  here("data","processed","LADs_sub.rds")
)

# Output Areas
oa <- st_read(
  here("data","processed","shp_files","OAs_comb.shp")
) %>%
  st_transform(27700) %>%
  st_make_valid()

# Roads
road_attributes <- readRDS(
  here("data","processed","roads_filtered.rds")
) %>%
  st_make_valid()

# ============================================================
# CRS consistency
# ============================================================

if (st_crs(road_attributes) != st_crs(oa)) {
  road_attributes <- st_transform(road_attributes, st_crs(oa))
}

if (st_crs(caz) != st_crs(oa)) {
  caz <- st_transform(caz, st_crs(oa))
}

# ============================================================
# Prepare CAZ polygons
# ============================================================

caz <- caz %>%
  mutate(
    start_date = as.Date(startDt, "%d/%m/%Y")
  )

# dissolve polygons by CAZ scheme
caz <- caz %>%
  group_by(scheme, start_date, LocAuth1) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ============================================================
# Calculate road length inside CAZ
# 
road_caz <- st_intersection(
  road_attributes,
  caz %>% select(scheme, start_date)
) %>%
  mutate(int_length = st_length(st_geometry(.)))%>%
  st_drop_geometry()

# total road length
road_length <- road_attributes %>%
  mutate(total_length = st_length(.)) %>%
  st_drop_geometry() %>%
  select(identifier, total_length)

road_caz_prop <- road_caz %>%
  left_join(road_length, by = "identifier") %>%
  mutate(
    prop_inside = as.numeric(int_length / total_length),
    caz_start_q = zoo::as.yearqtr(start_date)
  )

# Save CAZ exposure proportions
saveRDS(
  road_caz_prop,
  here("data","processed","roads_caz_props.rds")
)

# ============================================================
# Restrict roads to selected LADs
# ============================================================

roads_lad <- st_join(
  road_attributes,
  lads_sub %>% select(LAD24CD, LAD24NM),
  join = st_intersects,
  left = FALSE
)

# ============================================================
# Assign each road to ONE OA
# ============================================================

roads_oa <- st_intersection(
  roads_lad,
  oa_sub %>% select(OA)
) %>%
  mutate(
    int_length = st_length(.)
  )

# Select OA with largest overlap
road_attributes_OA <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup()

# ============================================================
# Save final road dataset
# ============================================================

st_write(
  road_attributes_OA,
  here("data","processed","road_attributes_OA.gpkg"),
  delete_dsn = TRUE
)
