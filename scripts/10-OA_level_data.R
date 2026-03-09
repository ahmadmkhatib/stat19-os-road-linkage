#============================================================
  # OA-Level Data
  # ============================================================

# This script constructs Output Area (OA)
# dataset used to define treatment and control groups for
# the CAZ policy evaluation.
#
# Treatment and control status are defined  using
# spatial relationships between OAs and CAZ boundaries.
#
# No road-level data are used in the classification.
#
# ------------------------------------------------------------
# GROUP DEFINITIONS
# ------------------------------------------------------------
#
# Treatment OAs
# ------------------------------------------------------------
# Output Areas that intersect CAZ boundaries.
#
#
# Buffer / Spillover OAs
# ------------------------------------------------------------
# Output Areas within 1 km outside CAZ boundaries.
#
# These areas may experience behavioural spillovers
# due to traffic rerouting or displacement.
#
#
# Control Group 1: Same-City Controls
# ------------------------------------------------------------
# OAs located in the same cities as CAZ zones but
# outside both the CAZ boundary and the 1 km buffer.
#
# These areas share the same macro environment but
# are unlikely to experience direct policy effects.
#
#
# Control Group 2: # OAs  of cities without CAZ.
#
# with distance to city centre 
#
#
# ------------------------------------------------------------
# OUTPUT DATASET
# ------------------------------------------------------------
#
# OA_CODE
# LAD24CD
# LAD24NM
#
# treated_OA
# buffer_OA
# control_group1_OA
# control_group2_OA
# distance 
#
# ============================================================

library(tidyverse)
library(sf)
library(here)

# ============================================================
# Load Data
# ============================================================
### lads subset
lads_sub<-readRDS (here("data", "processed", "LADs_sub.rds"))


# Output Area geometries
oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid()

# CAZ boundaries
caz_boundaries <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

st_crs(caz_boundaries)

caz_boundaries <-  st_transform(caz_boundaries, st_crs(oa)) %>% st_make_valid()
oa<-st_make_valid(oa)


# Attach LAD to each OA
oa_sub <- st_join(
  oa,
  lads_sub %>% select(LAD24CD),
  join = st_within
)

sum(duplicated(oa_sub$OA_CODE))
oa_sub <-oa_sub %>% filter(!is.na(LAD24CD))    # crop OA that are not in the samepl 

# Treatment OAs (inside CAZ)

treated_OAs_any <- oa_sub %>%
  st_join(caz_boundaries, join = st_intersects, left = FALSE) %>%
  pull(OA) %>% unique()

# ------------------------------
# Majority (>=50% area inside CAZ)
# ------------------------------

oa_intersection <- st_intersection(
  oa_sub %>% select(OA),
  caz_boundaries %>% select(scheme)
)

oa_intersection <- oa_intersection %>%
  mutate(intersection_area = st_area(geometry))

oa_area <- oa_sub %>%
  mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  select(OA, total_area)

oa_prop <- oa_intersection %>%
  st_drop_geometry() %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_inside =as.numeric (intersection_area / total_area))


summary(oa_prop$prop_inside)

treated_OAs_50pct <- oa_prop %>%
  filter(prop_inside >= 0.5) %>%
  pull(OA)


table(
  any = oa_prop$prop_inside > 0,
  pct50 = oa_prop$prop_inside >= 0.5
)


# Spillover Buffer OAs (1 km outside CAZ)
# ============================================================
caz_buffer <- st_buffer(caz_boundaries, 1000) %>%
  st_difference(caz_boundaries)


buffer_OAs <- oa_sub %>%
  filter(!OA %in% treated_OAs) %>%  # exclude treated
  st_join(caz_buffer, join = st_intersects, left = FALSE) %>%
 pull(OA) %>% unique()

# Treated LADs 
#Treated LADs 
treated_LADs <- oa_sub %>%
  filter(OA %in% treated_OAs) %>%
  distinct(LAD24CD) %>% pull(LAD24CD)

n_distinct(treated_LADs)
n_distinct(caz_boundaries$scheme)

unique(treated_LADs)
unique(caz_boundaries$scheme)

lads_sub %>% filter(LAD24CD %in% treated_LADs) %>% View() 


# city centroids for non-CAZ cities

city_centroids <- oa_sub %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid()


# ============================================================
# Distance from OA to City Centre
# ============================================================

oa_centroids <- st_centroid(oa_sub)

oa_with_centroid <- oa_centroids %>%
  left_join(
    city_centroids %>%
      st_drop_geometry() %>%
      mutate(city_geometry = city_centroids$geometry),
    by = "LAD24CD"
  )

dist_citycentre <- st_distance(
  st_geometry(oa_centroids),
  city_centroids[match(oa_centroids$LAD24CD, city_centroids$LAD24CD), ]$geometry,
  by_element = TRUE
) %>%
  as.numeric()


OA_classification <- oa_sub %>%
  st_drop_geometry() %>%
  select(OA, LAD24CD) %>%
  mutate(
    treated_OA_any = if_else(OA %in% treated_OAs_any, 1, 0),
    
    treated_OA_50pct = if_else(OA %in% treated_OAs_50pct, 1, 0),
    
    buffer_OA = if_else(OA %in% buffer_OAs, 1, 0),
    
    control_group1_OA = if_else(
      LAD24CD %in% treated_LADs &
        treated_OA_any == 0 &
        buffer_OA == 0,
      1,0
    ),
    
    control_group2_OA = if_else(
      !LAD24CD %in% treated_LADs &
        treated_OA_any == 0 &
        buffer_OA == 0,
      1,0
    )
  )

OA_classification$dist_citycentre <- dist_citycentre

# ============================================================
# Restrict to relevant OAs
# ============================================================

OA_analysis <- OA_classification %>%
  filter(
    treated_OA_any == 1 |
      buffer_OA == 1 |
      control_group1_OA == 1 |
      control_group2_OA == 1
  )


###checks
OA_analysis %>%
  mutate(group_count =
           treated_OA_any +
           buffer_OA +
           control_group1_OA +
           control_group2_OA) %>%
  count(group_count)

summary(OA_analysis$dist_citycentre)
##  ?!?!?! 

OA_analysis %>%
  filter(dist_citycentre > 20000) %>%
  count(LAD24CD, sort = TRUE)

saveRDS(OA_analysis, here("data", "processed", "OA_level_from_polygons.rds"))


                      