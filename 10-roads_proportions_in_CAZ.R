
library(sf)
library(tidyverse)
library(here)
library(units)



# # # - - - - load the roads cleaned data from data 9
roads<-st_read(here("data","processed","road_attributes_OA.gpkg"))
glimpse(roads)

# ── Load CAZ boundaries ────────────────────────────────────
caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp")
) %>%
  st_make_valid()

# ── Ensure same CRS ────────────────────────────────────────
roads <- st_transform(roads, st_crs(caz))

# ── Total road length ──────────────────────────────────────
roads <- roads %>%
  mutate(
    road_length = st_length(.)
  )

# ── Intersect roads with CAZ polygons ──────────────────────
roads_caz <- st_intersection(roads, caz)

# ── Length of road inside CAZ ──────────────────────────────
roads_caz <- roads_caz %>%
  mutate(
    length_in_caz = st_length(.)
  )

# ── Calculate proportion of road inside CAZ ────────────────

roads_caz_props <- roads_caz %>%
  st_drop_geometry() %>%
  select(identifier, scheme, length_in_caz) %>%
  left_join(
    roads %>%
      st_drop_geometry() %>%
      select(identifier, road_length),
    by = "identifier"
  ) %>%
  mutate(
    prop_in_caz = as.numeric(length_in_caz / road_length)
  )

# ── Save dataset ───────────────────────────────────────────
saveRDS(
  roads_caz_props,
  here("data","processed","roads_caz_props.rds")
)

