# ============================================================
# Merge OA Matching Dataset with Census Characteristics
# ============================================================

library(tidyverse)
library(lubridate)
library(here)
library(sf)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

OA_char_raw <- read.csv(here("data","processed","outputArea_raw.csv"))
OA_char_percent <- read.csv(here("data","processed","outputArea_percent.csv"))

OA_matching_dataset <- readRDS(
  here("data","processed","OA_matching_dataset.rds")
)

# ------------------------------------------------------------
# Load OA shapefile and compute area
# ------------------------------------------------------------

oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  mutate(
    area_m2 = as.numeric(st_area(.)),
    area_km2 = area_m2 / 1e6
  )

# ------------------------------------------------------------
# Pre-merge checks
# ------------------------------------------------------------

cat("OA_matching_dataset rows:", nrow(OA_matching_dataset), "\n")
cat("OA_char_raw rows:", nrow(OA_char_raw), "\n")
cat("OA_char_percent rows:", nrow(OA_char_percent), "\n")

cat("OA_char_raw dups:",
    OA_char_raw %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

cat("OA_char_percent dups:",
    OA_char_percent %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

# ------------------------------------------------------------
# Rename census variables
# ------------------------------------------------------------

vars_to_rename <- setdiff(
  names(OA_char_raw),
  c("OA","country","Total","IMD")
)

OA_char_raw_renamed <- OA_char_raw %>%
  rename_with(~ paste0(.x,"_n"), all_of(vars_to_rename))

names(OA_char_percent)
OA_char_pct_renamed <- OA_char_percent %>%
  rename_with(~ paste0(.x,"_pct"), all_of(vars_to_rename)) %>%
  select(-country,-Total,-IMD)

# ------------------------------------------------------------
# Merge census tables
# ------------------------------------------------------------

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_renamed, by="OA")

# ------------------------------------------------------------
# Merge census onto OA matching dataset
# ------------------------------------------------------------

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by="OA")

# ------------------------------------------------------------
# Add OA area + population density
# ------------------------------------------------------------

OA_matching_census <- OA_matching_census %>%
  left_join(
    oa_sub %>%
      st_drop_geometry() %>%
      select(OA, area_km2),
    by = "OA"
  ) %>%
  mutate(
    pop_density = Total / area_km2,
    log_pop_density = log1p(pop_density)
  )

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

stopifnot(nrow(OA_matching_census) == nrow(OA_matching_dataset))

stopifnot(
  OA_matching_census %>%
    count(OA) %>%
    filter(n > 1) %>%
    nrow() == 0
)

# ------------------------------------------------------------
# Save non-spatial dataset
# ------------------------------------------------------------

saveRDS(
  OA_matching_census,
  here("data","processed","OA_matching_census.rds")
)

# ------------------------------------------------------------
# Create spatial version
# ------------------------------------------------------------

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by="OA") %>%
  st_as_sf()

# ------------------------------------------------------------
# Save spatial dataset
# ------------------------------------------------------------

st_write(
  OA_matching_census_sf,
  here("data","processed","shp_files","OA_matching_census.gpkg"),
  delete_dsn = TRUE
)

cat("Saved: OA_matching_census.rds and OA_matching_census.gpkg\n")

# ------------------------------------------------------------
# Quick visual check
# ------------------------------------------------------------

OA_matching_census_sf %>%
  ggplot() +
  geom_sf(aes(fill = treated_OA), colour = NA) +
  theme_minimal()