# ============================================================
# OA-Level Analysis from roads
# ============================================================


# assigns each road to the OA it mostly overlaps
# computes road–CAZ overlap proportions
# creates treatment indicators
# compares classification methods
# maps disagreements





library(sf)
library(tidyverse)
library(here)

# roads
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
  here("data","processed","shp_files","OAs_subset.shp")
) 
oa_sub <-oa_sub%>%
  st_transform(st_crs(road_attributes)) %>%
  st_make_valid()

# LAD subset
lads_sub <- readRDS(
  here("data","processed","LADs_sub.rds")
)

#Roads intersecting CAZ with proportions 
road_total_len <- road_attributes %>%
  mutate(total_length = st_length(geom)) %>%
  st_drop_geometry() %>%
  select(identifier, total_length)

road_caz_intersection <- road_attributes %>%
  st_intersection(
    caz %>% select(scheme)
  ) %>%
  mutate(int_length = st_length(geometry)) %>%
  st_drop_geometry() %>%
  group_by(identifier) %>%
  summarise(int_length = sum(int_length), .groups="drop")

road_caz_prop <- road_total_len %>%
  left_join(road_caz_intersection, by="identifier") %>%
  mutate(
    int_length = replace_na(int_length, units::set_units(0,"m")),
    prop_inside = as.numeric(int_length / total_length)
  )


road_treatment <- road_caz_prop %>%
  mutate(
    treated_any = if_else(prop_inside > 0,1,0),
    treated_50pct = if_else(prop_inside >= 0.5,1,0)
  )

##1KM buffer
caz_buffer <- st_buffer(caz,1000) %>%
  st_difference(caz)


### roads in the buffer 
roads_buffer <- road_attributes %>%
  filter(!identifier %in% roads_treated) %>%
  st_join(
    caz_buffer,
    join = st_intersects,
    left = FALSE
  ) %>%
  pull(identifier) %>%
  unique()


road_treatment <- road_treatment %>%
  mutate(
    buffer_road = if_else(identifier %in% roads_buffer,1,0)
  )





####   compare with the OA method 
OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")
)


road_OA_method <- roads_with_OA %>%
  st_drop_geometry() %>%
  left_join(OA_classification, by="OA")


road_compare <- road_attributes %>%
  st_drop_geometry() %>%
  left_join(road_treatment, by="identifier") %>%
  left_join(OA_analysis, by="OA")

#Road any vs OA any
table(
  road_method = road_compare$treated_any,
  oa_method   = road_compare$treated_OA_any
)

#Road 50% vs OA 50%
table(
  road50 = road_compare$treated_50pct,
  oa50   = road_compare$treated_OA_50pct
)

## % disgreement 

road_compare %>%
  mutate(
    disagree_any = treated_any != treated_OA_any,
    disagree_50  = treated_50pct != treated_OA_50pct
  ) %>%
  summarise(
    pct_disagree_any = mean(disagree_any)*100,
    pct_disagree_50  = mean(disagree_50)*100
  )

table(road_compare$treated_road, road_compare$treated_OA_any)

road_compare %>%
  mutate(disagree = treated_road != treated_OA_any) %>%
  summarise(pct_disagree = mean(disagree)*100)


## missclassified ones 
misclassified_roads <- road_compare %>%
  filter(treated_50pct != treated_OA_50pct)

nrow(misclassified_roads)



roads_map <- road_attributes %>%
  left_join(road_compare, by="identifier")

#road based treatment 
tm_shape(roads_map) +
  tm_lines("treated_50pct") +
  tm_layout(title="Road-based treatment (≥50% inside CAZ)")

## OA_based treatment 
tm_shape(roads_map) +
  tm_lines("treated_OA_50pct") +
  tm_layout(title="OA-based treatment (≥50% OA inside CAZ)")

##

### misclassification map 
roads_map$misclassified <-
  roads_map$treated_50pct != roads_map$treated_OA_50pct

tm_shape(roads_map) +
  tm_lines(
    col="misclassified",
    palette=c("grey","red")
  ) +
  tm_layout(title="Disagreement between road and OA methods")




analysis_roads <- road_compare %>%
  mutate(
    treated_group_road_any = treated_any,
    treated_group_road_50pct = treated_50pct
  )