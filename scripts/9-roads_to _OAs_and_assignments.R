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

table(OA_analysis$assignment)

if (st_crs(roads) != st_crs(oa_sub)) {
  roads <- st_transform(roads, st_crs(oa_sub))
}

# ============================================================
#  Intersect Roads with OAs
# ============================================================

roads_oa <- st_intersection(
  roads,
  oa_sub %>%dplyr::  select(OA)
)

# length of road segment inside OA
roads_oa <- roads_oa %>%
  mutate(
    int_length = st_length(.)
  )

# ============================================================
#  Assign each road to ONE OA (largest overlap)
# ============================================================

roads_oa <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup()

# ============================================================
#Attach OA treatment classification
# ============================================================

road_attributes_OA <- roads_oa %>%
  left_join(
    OA_analysis %>%
    dplyr::  select(OA, assignment,scheme),
    by = "OA"
  )


n_distinct(roads$identifier)
n_distinct(road_attributes_OA$identifier)

missing_roads <- roads %>%
  filter(!identifier %in% road_attributes_OA$identifier)

nrow(missing_roads)

# ============================================================
# road-level dataset
# ============================================================

st_write(
  road_attributes_OA,
  here("data","processed","road_attributes_OA.gpkg"),
  delete_dsn = TRUE
)

# ============================================================
# Aggregate road characteristics to OA
# ============================================================
# ── Keep full intersection (one row per road-OA pair) ─────────────────────
## Each OA gets credit for every road that touches it

OA_roads <- roads_oa %>%          # 
  st_drop_geometry() %>%
  group_by(OA) %>%
  summarise(
    n_roads           = n_distinct(identifier),
    total_road_length = sum(as.numeric(int_length), na.rm = TRUE),
    n_A               = sum(road_class == "A",        na.rm = TRUE),
    n_B               = sum(road_class == "B",        na.rm = TRUE),
    n_motorway        = sum(road_class == "Motorway", na.rm = TRUE),
    n_minor           = sum(road_class == "minor",    na.rm = TRUE),
    .groups = "drop"
  )

# ── One-OA-per-road assignment (largest overlap) ──────────────────────────
# 
road_attributes_OA <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(
    OA_analysis %>%
      arrange(OA) %>%
      distinct(OA, .keep_all = TRUE) %>%    
      select(OA, assignment),
    by = "OA"
  )

# attach OA classification
OA_roads <- OA_roads %>%
  left_join(
    OA_analysis,  by = "OA"  )

glimpse(OA_roads)
OA_roads %>% count(OA) %>% filter(n > 1) %>% nrow()  # no road belong to more than 1 OA
nrow(OA_roads)  

OA_missing_roads <- OA_analysis %>%
  filter(!OA %in% OA_roads$OA)

nrow(OA_missing_roads)

OA_missing_roads %>%
  count(assignment)


table(OA_roads$assignment)


OA_analysis %>%
  count(OA) %>%
  filter(n > 1)
n_distinct(OA_analysis$OA)
nrow(OA_analysis)


# ============================================================
# Save OA road dataset
# ============================================================

saveRDS(
  OA_roads,
  here("data","processed","OA_roads.rds")
)




