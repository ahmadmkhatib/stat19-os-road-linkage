# ============================================================
#  Assign CAZ Treatment to Road Link × Quarter Panel
#Road link is treated if it intersects CAZ polygon at any length.
# ============================================================

library(tidyverse)
library(lubridate)
library(sf)
library(here)
library(zoo)

# ------------------------------------------------------------
#  Load data 
# ------------------------------------------------------------

caz_polygons <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)


### all raods and thier attributes 
road_attributes<- read_rds(here("data", "processed", "roads_filtered.rds")) %>%
  select(identifier, road_class, geom) %>%
  distinct(identifier, .keep_all = TRUE) %>%
  st_as_sf()


## geo data 
lads_sub<-readRDS (here("data", "processed", "LADs_sub.rds"))
oa_2011 <- st_read("../stat19-os-road-linkage-data/infuse_oa_lyr_2011_clipped.shp") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  rename(OA_CODE = geo_code) 

st_crs(road_attributes)
st_crs(lads_sub)
st_crs(oa_2011)





# ------------------------------------------------------------
# prepare CAZ data
# ------------------------------------------------------------

glimpse(caz_polygons)
class(caz_polygons$startDt)

caz <- caz_polygons %>%
  mutate(
    start_date = as.Date(startDt, "%d/%m/%Y"),
    caz_id = row_number()
  )   %>%    st_make_valid()



View(caz)

# recode caz  

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
road_caz <- st_intersection(
  road_attributes,
  caz %>% select(start_date)
) %>%
  st_drop_geometry() %>%
  group_by(identifier) 


road_caz %>%
  count(identifier) %>%
  filter(n > 1)




## Instead of "any overlap", compute proportion of road length inside CAZ:
 road_caz_prop <- st_intersection(
    road_attributes,
    caz %>% select(scheme, start_date)
  ) %>%
   mutate(int_length = st_length(st_geometry(.))) %>%
  st_drop_geometry()

 road_lengths <- road_attributes %>%
   mutate(total_length = st_length(st_geometry(.))) %>%
   st_drop_geometry() %>%
   select(identifier, total_length)
 
 road_caz_prop <- road_caz_prop %>%
   left_join(road_lengths, by = "identifier") %>%
   mutate(prop_inside = int_length / total_length)
 # prop_inside = proportion of link inside CAZ

summary(road_caz_prop$prop_inside)
table(road_caz_prop$scheme)
# ------------------------------------------------------------
# Convert CAZ start date → quarter
# ------------------------------------------------------------

road_caz_prop <- road_caz_prop%>%
  mutate(
    caz_start_q = as.yearqtr(start_date)
  )



saveRDS(road_caz_prop,  here("data", "processed", "roads_caz_props.rds"))

names(road_attributes)
attr(road_attributes, "sf_column")

### create roads attributes dataset 
road_attributes <- road_attributes %>%
  
  # Attach LAD
  st_join(
    lads_sub %>% select(LAD24CD, LAD24NM),
    join = st_intersects,
    left = TRUE
  ) %>%
  
  #  Compute total road length
  mutate(total_length = st_length(.)) %>%
  
  # Intersect with OA polygons
  st_intersection(
    oa_2011 %>% select(OA_CODE)
  ) %>%
  
  # Compute share inside each OA
  mutate(
    segment_length = st_length(.),
    share = as.numeric(segment_length / total_length)
  ) %>%
  # Keep OA with largest share
  group_by(identifier) %>%
  slice_max(order_by = share, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    identifier,
    road_class,
    LAD24CD,
    LAD24NM,
    OA_CODE
  )
  )

saveRDS(road_attributes, here("data", "processed", "road_attributes.rds"))
