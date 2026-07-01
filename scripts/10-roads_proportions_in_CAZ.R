
library(sf)
library(tidyverse)
library(lubridate)
library(zoo)
library(here)
library(units)
library(lubridate)



# # # - - - - load the roads cleaned data from data 9
roads<-st_read(here("data","processed","road_attributes_OA.gpkg"))
glimpse(roads)

# ── Load CAZ boundaries ────────────────────────────────────
caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp")
) %>%
  st_make_valid()

table(caz$scheme, caz$startDt)

# Replace incorrect start date for Portsmouth
caz <- caz %>%
  mutate(
    startDt = as.character(startDt),
    startDt = if_else(scheme == "Portsmouth", "29/11/2021", startDt)
  )


caz %>%
  st_drop_geometry() %>%
  filter(scheme == "Portsmouth") %>%
  distinct(scheme, startDt) %>%
  print()


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
    roads %>% st_drop_geometry() %>% dplyr::select(identifier, road_length),
    by = "identifier"
  ) %>%
  mutate(prop_in_caz = as.numeric(length_in_caz / road_length)) %>%
  left_join(
    caz %>%
      st_drop_geometry() %>%
      dplyr:: select(scheme, startDt) %>%
      distinct(scheme, .keep_all = TRUE) %>%
      mutate(
        start_date  = dmy(startDt),
        raw_q       = as.yearqtr(start_date),
        q_start     = as.Date(raw_q),
        q_end       = as.Date(raw_q + 0.25) - 1,
        q_mid       = q_start + as.integer(difftime(q_end, q_start, units = "days")) / 2,
        # Majority rule: if treatment started past the quarter midpoint,
        # most of that quarter was pre-treatment → assign to next quarter
        caz_start_q = if_else(start_date > q_mid, raw_q + 0.25, raw_q)
      ) %>%
      select(scheme, startDt, caz_start_q),
    by = "scheme"
  )


# ── Join CAZ proportion back to roads ───────────────────────
roads <- roads %>%
  left_join(
    roads_caz_props %>% dplyr::select(identifier, scheme, prop_in_caz),
    by = "identifier"
  )


st_write(
  roads,
  here("data","processed","roads_final.gpkg"),
  delete_dsn = TRUE
)



saveRDS(
  roads_caz_props,
  here("data","processed","roads_caz_props.rds")
)

# because portsmoth was corrected
st_write(
  caz,
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  delete_dsn = TRUE
)

