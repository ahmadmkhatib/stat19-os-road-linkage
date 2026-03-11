#============================================================
# OA-Level Data Construction
#============================================================
# Purpose:
# Construct OA dataset defining treatment, buffer, and control
# areas.
#
# Classification uses majority-area (>=50%) spatial rules.
#============================================================

library(tidyverse)
library(sf)
library(here)

#============================================================
# 1. Load Data
#============================================================

lads_sub <- readRDS(here("data","processed","LADs_sub.rds"))

oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid()

caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_transform(st_crs(oa)) %>%
  st_make_valid()

#============================================================
# 2. Keep OAs within study LADs
#============================================================

oa_sub <- st_join(
  oa,
  lads_sub %>% select(LAD24CD, LAD24NM),
  join = st_within
) %>%
  filter(!is.na(LAD24CD))

#============================================================
# 3. Calculate OA area
#============================================================

oa_area <- oa_sub %>%
  mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  select(OA, total_area)

#============================================================
# 4. OA proportion inside CAZ
#============================================================

oa_caz_intersection <- st_intersection(
  oa_sub %>% select(OA),
  caz %>% select(scheme)
) %>%
  mutate(int_area = st_area(geometry)) %>%
  st_drop_geometry()

oa_caz_prop <- oa_caz_intersection %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_caz = as.numeric(int_area / total_area))

treated_OAs <- oa_caz_prop %>%
  filter(prop_caz >= 0.5) %>%
  pull(OA)

#============================================================
# 5. Create 1km CAZ buffer
#============================================================

caz_buffer <- st_buffer(caz,1000) %>%
  st_difference(caz)

oa_buffer_intersection <- st_intersection(
  oa_sub %>% select(OA),
  caz_buffer
) %>%
  mutate(int_area = st_area(geometry)) %>%
  st_drop_geometry()

oa_buffer_prop <- oa_buffer_intersection %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_buffer = as.numeric(int_area / total_area))

buffer_OAs <- oa_buffer_prop %>%
  filter(prop_buffer >= 0.5) %>%
  pull(OA)

#============================================================
# 6. Identify CAZ LADs
#============================================================

treated_LADs <- oa_sub %>%
  filter(OA %in% treated_OAs) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

#============================================================
# 7. Distance to city centre
#============================================================

city_centroids <- oa_sub %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry), .groups="drop") %>%
  st_centroid()

oa_centroids <- st_centroid(oa_sub)

dist_citycentre <- st_distance(
  st_geometry(oa_centroids),
  city_centroids[match(oa_centroids$LAD24CD,
                       city_centroids$LAD24CD),]$geometry,
  by_element=TRUE
) %>%
  as.numeric()

#============================================================
# 8. OA classification
#============================================================

OA_analysis <- oa_sub %>%
  st_drop_geometry() %>%
  select(OA, LAD24CD, LAD24NM) %>%
  mutate(
    
    treated_OA = if_else(OA %in% treated_OAs,1,0),
    
    buffer_OA = if_else(
      OA %in% buffer_OAs &
        treated_OA == 0,
      1,0
    ),
    
    control_group1_OA = if_else(
      LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    ),
    
    control_group2_OA = if_else(
      !LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    )
  )

OA_analysis$dist_citycentre <- dist_citycentre


OA_analysis %>%  select(dist_citycentre) %>%   summary()
OA_analysis  %>% filter(treated_OA==1)%>% select(dist_citycentre) %>%   summary()
## # is this resonable 10K 



OA_analysis <- OA_analysis %>%
  mutate(
    assignment = case_when(
      treated_OA == 1 ~ "treated",
      buffer_OA == 1 ~ "buffer",
      control_group1_OA == 1 ~ "control_same_city",
      control_group2_OA == 1 ~ "control_other_city",
      TRUE ~ NA_character_
    )
  )


OA_analysis %>%
  group_by(assignment) %>%
  summarise(
    max_dist = max(dist_citycentre),
    median_dist = median(dist_citycentre)
  )

table(OA_analysis$assignment)

#distributions 
OA_analysis %>%
  count(assignment) %>%
  mutate(pct = n / sum(n))


#============================================================
# 10. Save outputs
#============================================================

saveRDS(
  OA_analysis,
  here("data","processed","OA_level_from_polygons.rds")
)

st_write(
  oa_sub,
  here("data","processed","shp_files","OA_subset.shp"),
  delete_dsn=TRUE
)
