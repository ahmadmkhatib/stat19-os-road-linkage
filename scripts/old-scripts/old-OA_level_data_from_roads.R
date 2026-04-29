# ============================================================
# OA-Level Analysis from Roads
# ============================================================
# Purpose:
# This script assigns treatment status to roads based on their
# overlap with Clean Air Zone (CAZ) boundaries.
#
# Steps:
# 1. Calculate proportion of each road inside CAZ
# 3. Create treatment indicators
# 4. Identify roads within a 1km CAZ spillover buffer
# 5. Compare road-based vs OA-based classification
# 6. Map disagreements between methods
# ============================================================

library(sf)
library(tidyverse)
library(here)
library(tmap)

# ============================================================
# 1. Load spatial data
# ============================================================

# Road network with attributes
road_attributes <- st_read(
  here("data","processed","road_attributes.gpkg")
)

# CAZ boundaries
caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp")
) %>%
  st_transform(st_crs(road_attributes)) %>%
  st_make_valid()

# Output Areas
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp")
) %>%
  st_transform(st_crs(road_attributes)) %>%
  st_make_valid()


# LAD subset
lads_sub <- readRDS(
  here("data","processed","LADs_sub.rds")
)

oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid() 


oa_sub_from_roads <-oa %>% filter(OA %in%  unique(road_attributes$OA))

# Attach LAD to each OA
oa_sub_roads <-st_join(
  oa_sub_from_roads,
  lads_sub %>% select(LAD24CD),
  join = st_within
) %>%  as.data.frame()

oa_sub_roads %>% filter(LAD24CD %in% lads_sub$LAD24CD) %>% nrow()
road_attributes %>% filter(LAD24CD %in% lads_sub$LAD24CD) %>% nrow()
# ============================================================
# 2. Calculate road length and CAZ overlap
# ============================================================

# Total length of each road
road_total_len <- road_attributes %>%
  mutate(total_length = st_length(st_geometry(.))) %>%
  st_drop_geometry() %>%
  select(identifier, total_length)

# Intersections between roads and CAZ
road_caz_intersection <- road_attributes %>%
  st_intersection(
    caz %>% select(scheme)
  ) %>%
  mutate(int_length = st_length(st_geometry(.))) %>%
  st_drop_geometry() %>%
  group_by(identifier) %>%
  summarise(
    int_length = sum(int_length),
    .groups = "drop"
  )

# ============================================================
# 3. Calculate proportion of each road inside CAZ
# ============================================================

road_caz_prop <- road_total_len %>%
  left_join(road_caz_intersection, by = "identifier") %>%
  mutate(
    int_length = replace_na(int_length, units::set_units(0,"m")),
    prop_inside = as.numeric(int_length / total_length)
  )

# ============================================================
# 4. Create road treatment indicators
# ============================================================

road_treatment <- road_caz_prop %>%
  mutate(
    treated_any = if_else(prop_inside > 0,1,0),
    treated_50pct = if_else(prop_inside >= 0.5,1,0)
  )

# ============================================================
# 5. Create 1km CAZ spillover buffer
# ============================================================

caz_buffer <- st_buffer(caz, 1000) %>%
  st_difference(caz)

# Identify roads intersecting the buffer
roads_buffer <- road_attributes %>%
  st_join(
    caz_buffer,
    join = st_intersects,
    left = FALSE
  ) %>%
  pull(identifier) %>%
  unique()

# Add buffer indicator
road_treatment <- road_treatment %>%
  mutate(
    buffer_road = if_else(identifier %in% roads_buffer,1,0)
  )

# ============================================================
# 6. Load OA-level classification
# ============================================================

OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")
)

# ============================================================
# 7. Compare road vs OA treatment assignment
# ============================================================

road_compare <- road_attributes %>%
  st_drop_geometry() %>%
  left_join(road_treatment, by="identifier") %>%
  left_join(OA_analysis, by="OA")

# Any overlap comparison
table(
  road_method = road_compare$treated_any,
  oa_method   = road_compare$treated_OA_any
)

# 50% threshold comparison
table(
  road50 = road_compare$treated_50pct,
  oa50   = road_compare$treated_OA_50pct
)

# ============================================================
# 8. Calculate disagreement rates
# ============================================================

road_compare %>%
  mutate(
    disagree_any = treated_any != treated_OA_any,
    disagree_50  = treated_50pct != treated_OA_50pct
  ) %>%
  summarise(
    pct_disagree_any = mean(disagree_any, na.rm = TRUE) * 100,
    pct_disagree_50  = mean(disagree_50,  na.rm = TRUE) * 100
  )

# ============================================================
# 9. Identify misclassified roads
# ============================================================

misclassified_roads <- road_compare %>%
  filter(treated_50pct != treated_OA_50pct)

nrow(misclassified_roads)

# ============================================================
# 10. Map treatment classifications
# ============================================================

roads_map <- road_attributes %>%
  left_join(road_compare, by="identifier")

# Road-based treatment
tm_shape(roads_map) +
  tm_lines("treated_50pct") +
  tm_layout(title="Road-based treatment (≥50% inside CAZ)")

# OA-based treatment
tm_shape(roads_map) +
  tm_lines("treated_OA_50pct") +
  tm_layout(title="OA-based treatment (≥50% OA inside CAZ)")

# ============================================================
# 11. Map disagreement between methods
# ============================================================

roads_map$misclassified <-
  roads_map$treated_50pct != roads_map$treated_OA_50pct

tm_shape(roads_map) +
  tm_lines(
    col = "misclassified",
    palette = c("grey","red")
  ) +
  tm_layout(title="Disagreement between road and OA methods")

# ============================================================
# 12. Final analysis dataset
# ============================================================

road_caz_prop %>% filter(prop_inside > 0 &  prop_inside < 0.5) %>% nrow()
#partially overlapping roads should be removed from control.

names(analysis_roads)


analysis_roads <- analysis_roads %>%
  mutate(
    road_group = case_when(
      prop_inside >= 0.5 ~ "treated",
      prop_inside == 0 & buffer_road == 1 ~ "buffer",
      prop_inside == 0 & buffer_road == 0 ~ "control",
      prop_inside > 0 & prop_inside < 0.5 ~ "excluded_partial"
    )
  )
analysis_roads %>%
  count(road_group)
analysis_roads %>%
  summarise(
    na_prop_inside = sum(is.na(prop_inside)),
    na_buffer = sum(is.na(buffer_road))
  )
analysis_roads %>%
  filter(is.na(prop_inside)) %>%
  select(identifier, total_length, int_length) %>%
  head()
