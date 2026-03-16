
#=========================================================
library(tidyverse)
library(lubridate)
library(here)
library(sf)




OA_matching_dataset<- readRDS(here("data","processed","OA_matching_dataset.rds"))
glimpse(OA_matching_dataset)

####### as a shapefile   # # # # 



# OA shpfile 
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()



# Join geometry to matching dataset
OA_matching_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_dataset, by = "OA") %>%
  st_as_sf()

# Quick check
print(OA_matching_sf)


table(st_is_valid(OA_matching_sf))
st_geometry_type(OA_matching_sf) %>% table()


# Save as GeoPackage
st_write(
  OA_matching_sf,
  here("data","processed","shp_files","OA_matching_dataset.gpkg"),
  delete_dsn = TRUE
)


OA_matching_sf<- st_read(here("data","processed","shp_files","OA_matching_dataset.gpkg"))

OA_matching_sf %>%
  ggplot() +
  geom_sf(aes(fill = treated_OA), colour = NA) +
  theme_minimal()
## check the CAZ shapfile 


