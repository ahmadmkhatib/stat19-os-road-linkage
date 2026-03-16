# ============================================================
# Assign Roads to Output Areas (OA)
# ============================================================
# Purpose:
# 1. Assign each road link to ONE Output Area
#    based on the largest share of road length.
# 2. Attach OA treatment assignment
# 3. Create OA-level road characteristics dataset
# # # output : OA_roads_dataset.rds
# ============================================================

library(sf)
library(tidyverse)
library(here)


# Roads
roads <- readRDS(
  here("data","processed","roads_filtered.rds")
) %>%
  st_make_valid()

glimpse(roads)

# OA shpfile 
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()

# OA classification dataset
OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")         # from 8
)


if (st_crs(roads) != st_crs(oa_sub)) {
  roads <- st_transform(roads, st_crs(oa_sub))
}

# ============================================================
# 3 Intersect Roads with OAs
# ============================================================

roads_oa <- st_intersection(
  roads,
  oa_sub %>% select(OA)
)

# length of road segment inside OA
roads_oa <- roads_oa %>%
  mutate(
    int_length = st_length(.)
  )

# ============================================================
# 4 Assign each road to ONE OA (largest overlap)
# ============================================================

roads_oa <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup()

# ============================================================
# 5 Attach OA treatment classification
# ============================================================

road_attributes_OA <- roads_oa %>%
  left_join(
    OA_analysis %>%
      select(OA, assignment),
    by = "OA"
  )

glimpse(OA_analysis)
# ============================================================
# 6 Save road-level dataset
# ============================================================
road_attributes_OA<-st_read(here("data","processed","road_attributes_OA.gpkg"))

st_write(
  road_attributes_OA,
  here("data","processed","road_attributes_OA.gpkg"),
  delete_dsn = TRUE
)

# ============================================================
# 7 Aggregate road characteristics to OA
# ============================================================

OA_roads <- road_attributes_OA %>%
  st_drop_geometry() %>%
  group_by(OA) %>%
  summarise(
    scheme = first(na.omit(scheme)),   # scheme affecting that OA
    n_roads = n(),
    total_road_length = sum(as.numeric(int_length)),
    n_A = sum(road_class == "A"),
    n_B = sum(road_class == "B"),
    n_motorway = sum(road_class == "Motorway"),
    n_minor = sum(road_class == "minor"),
    .groups = "drop"
  )

# attach OA classification
OA_roads <- OA_roads %>%
  left_join(OA_analysis, by = "OA")




glimpse(OA_roads)

OA_missing_roads <- OA_analysis %>%
  filter(!OA %in% OA_roads$OA)

nrow(OA_missing_roads)

OA_missing_roads %>%
  count(assignment)

# ============================================================
# Save OA road dataset
# ============================================================

saveRDS(
  OA_roads,
  here("data","processed","OA_roads.rds")
)






