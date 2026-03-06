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


oa<- oa %>%
  rename(OA_CODE = OA) 


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

treated_OAs <- oa_sub %>%
  st_join(caz_boundaries, join = st_intersects, left = FALSE) %>%
  pull(OA_CODE) %>% unique()

# Spillover Buffer (1 km outside CAZ)
# ============================================================
caz_buffer <- st_buffer(caz_boundaries, 1000) %>%
  st_difference(caz_boundaries)


buffer_OAs <- oa_sub %>%
  filter(!OA_CODE %in% treated_OAs) %>%  # exclude treated
  st_join(caz_buffer, join = st_intersects, left = FALSE) %>%
 pull(OA_CODE) %>% unique()

# Treated LADs 
#Treated LADs 
treated_LADs <- oa_sub %>%
  filter(OA_CODE %in% treated_OAs) %>%
  distinct(LAD24CD) %>% pull(LAD24CD)

n_distinct(treated_LADs)
n_distinct(caz_boundaries$scheme)

unique(treated_LADs)
unique(caz_boundaries$scheme)

lads_sub %>% filter(LAD24CD %in% treated_LADs) %>% View() 
### there are 1 scheam for 2 LADs 



# city centroids for non-CAZ cities
city_centroids <- oa_s %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid()
# ============================================================
# Distance from OA to City Centre
# ============================================================

oa_centroids <- st_centroid(oa_2011)

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


OA_classification <- oa_2011 %>%
  st_drop_geometry() %>%
  select(OA_CODE, LAD24CD) %>%
  mutate(
    
    treated_OA = if_else(OA_CODE %in% treated_OAs, 1, 0),
    
    buffer_OA = if_else(OA_CODE %in% buffer_OAs, 1, 0),
    
    control_group1_OA = if_else(
      LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1, 0
    ),
    
    control_group2_OA = if_else(
      !LAD24CD %in% treated_LADs  &
        treated_OA == 0 &
        buffer_OA == 0,
      1, 0
    )
    
  )

# add dist
OA_classification$dist_citycentre <- dist_citycentre


# relevant OAs
OA_analysis <- OA_classification %>%
  filter(treated_OA == 1 | buffer_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1)

OA_analysis %>%
  mutate(group_count =
           treated_OA +
           buffer_OA +
           control_group1_OA +
           control_group2_OA) %>%
  count(group_count)



summary(OA_analysis$dist_citycentre)

saveRDS(OA_analysis, here("data", "processed", "OA_level_from_polygons.rds"))
OA_analysi_old<-read_rds(here("data", "processed", "OA_level_from_polygons.rds"))

                      