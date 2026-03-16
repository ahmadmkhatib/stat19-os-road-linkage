
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

glimpse(caz)
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
  group_by(identifier, scheme) %>%
  summarise(length_in_caz = sum(length_in_caz), .groups = "drop") %>%
  left_join(
    roads %>% st_drop_geometry() %>% select(identifier, road_length),
    by = "identifier"
  ) %>%
  mutate(prop_in_caz = as.numeric(length_in_caz / road_length))


# ── Save dataset ───────────────────────────────────────────

# ── Join CAZ proportion back to roads ───────────────────────
roads <- roads %>%
  left_join(
    roads_caz_props %>% select(identifier, scheme, prop_in_caz),
    by = "identifier"
  )

# ── Save final roads dataset ─────────────────────────────────
st_write(
  roads,
  here("data","processed","roads_final.gpkg"),
  delete_dsn = TRUE
)


####
#
#


saveRDS(
  roads_caz_props,
  here("data","processed","roads_caz_props.rds")
)

