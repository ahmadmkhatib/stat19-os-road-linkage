# ============================================================
#  Assign CAZ Treatment/control  to Road Links
#Road links considered exposed if it intersects CAZ polygon at any length.
# ============================================================

library(tidyverse)
library(lubridate)
library(sf)
library(here)
library(zoo)

# -----------
#  data 
# ------

caz_polygons <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)


lads_sub<-readRDS(here("data", "processed", "LADs_sub.rds"))


oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid() 


###  raods subset and thier attributes 
road_attributes<-readRDS(here("data", "processed", "roads_filtered.rds"))
### save it as gpkg for spatial handeling 
st_write(
  road_attributes,
  here("data","processed","road_attributes.gpkg"),
  delete_dsn = TRUE
)                         

#read road shp file                         
road_attributes <- st_read(
  here("data","processed","road_attributes.gpkg"),
  quiet = TRUE
) %>%  select(identifier, road_class, length)

st_crs(road_attributes)
st_crs(lads_sub)
st_crs(oa)

summary(road_attributes$length)
# ------------------------------------------------------------
# prepare CAZ data
# ------------------------------------------------------------

glimpse(caz_polygons)
table(caz_polygons$scheme)
table(caz_polygons$LocAuth1)
class(caz_polygons$startDt)

caz <- caz_polygons %>%
  mutate(
    start_date = as.Date(startDt, "%d/%m/%Y"),
    caz_id = row_number()
  )   %>%    st_make_valid()


# recode caz class
table(caz$Class)
caz <- caz %>%
  mutate(
    Class = case_when(
      LocAuth1 == "Bradford" ~ "C",
      LocAuth1 %in% c("Aberdeen City", "Dundee City", "City of Edinburgh", "Glasgow City") ~ "LEZ",
      LocAuth1 == "Oxford" ~ "ZEZ",
      TRUE ~ Class
    )
  )


table(caz$start_date)

# dissolve polygons by scheme
caz <- caz %>%
   group_by(scheme, start_date, LocAuth1, Class) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")


# same crs

if (st_crs(caz) != st_crs(road_attributes)) {
  caz <- st_transform(caz, st_crs(road_attributes))
}

# ------------------------------------------------------------
#   intersection (road link inside CAZ (any dist))
# ------------------------------------------------------------
road_attributes <- st_make_valid(road_attributes)
caz <- st_make_valid(caz)


road_caz <- st_intersection(
  road_attributes,
  caz %>% select(start_date)
) %>%
  st_drop_geometry() %>%
  group_by(identifier) 


## Instead of "any overlap", compute proportion of road length inside CAZ:
 road_caz_prop <- st_intersection(
    road_attributes,
    caz %>% select(scheme, start_date)
  ) %>%
   mutate(int_length = st_length(st_geometry(.))) %>%
  st_drop_geometry()
 
sum(duplicated(road_caz_prop$identifier))
 
 summary(road_caz_prop$length)
 summary(road_caz_prop$int_length)
 
 road_lengths <- road_attributes %>%
   mutate(total_length = st_length(st_geometry(.))) %>%
   st_drop_geometry() %>%
   select(identifier, total_length)
 
 road_caz_prop  <- road_caz_prop %>%
   left_join(road_lengths, by = "identifier") 
 
 ###Some roads intersect minimally with CAZ polygons and also with OA polygons- #
 # so I created two scenarios, any or 50% intersections and they differ by ~300 road links- 

 road_caz_prop <- road_caz_prop%>%
   mutate(
     prop_inside = as.numeric(int_length / total_length),  # proportion inside
     ever_treated_any = 1,                                  # any overlap = treated
     ever_treated_50pct = if_else(prop_inside >= 0.5, 1, 0) # 50%+ treated
   )
 
 # prop_inside = proportion of link inside CAZ

summary(road_caz_prop$prop_inside)

table(road_caz_prop$ever_treated_any)
table(road_caz_prop$ever_treated_50pct)
## 

table(road_caz_prop$scheme)

# -----------------#-----------------#
# Convert CAZ start date → quarter
# ----------------###------------------#
glimpse(road_caz_prop)
road_caz_prop <- road_caz_prop%>%
  mutate(
    caz_start_q = as.yearqtr(start_date)
  )


#all raod links that intersect with CAZ at any length (we may filter more in difineing treatment areas) 

saveRDS(road_caz_prop,  here("data", "processed", "roads_caz_props.rds"))

names(road_attributes)

## create a subset of OAs
# Pre-crop OAs to roads extent
oa_sub <- st_filter(
  oa %>% select(OA),
  road_attributes,
  .predicate = st_intersects
)

roads_lad <- st_join(
  road_attributes,
  lads_sub %>% select(LAD24CD, LAD24NM),
  join = st_intersects,
  left = FALSE
)

# Intersect with OAs
roads_oa <- st_intersection(
  roads_lad,
  oa_sub %>% select(OA)
  ) %>% 
  mutate(int_length= st_length(geom))

# Keep OA with largest overlap per road
road_attributes_OA <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(identifier, road_class, LAD24CD, LAD24NM, OA, geom)


### save

st_write(
  road_attributes_OA,
  dsn = here("data", "processed", "road_attributes.gpkg"),
  layer = "road_attributes",
  delete_layer = TRUE
)
