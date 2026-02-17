# ==========================================================
# Script: 00_create_large_city_LAD_subset.R
# Purpose: Create a subset of large UK cities (>=100k population)
#          and assign corresponding Local Authority District (LAD)
#          codes using spatial joins.
# Output:  data/processed/big_cities_with_LADs.rds
# ==========================================================

library(sf)
library(tidyverse)
library(readxl)
library(stringr)
library(here)

# ----------------------------------------------------------
# Paths
# ----------------------------------------------------------

lads_path          <- here("data", "raw", "LAD_DEC_24_UK_BGC.shp")
buas_path          <- here("data", "raw", "BUAs.xlsx")
scotland_pop_path  <- here("data", "raw", "scot_pop.xlsx")
cities_path        <- here("data", "raw", "ukcities.csv")



dir.create(here("data", "processed"), showWarnings = FALSE)
output_path <- here("data", "processed", "big_cities_with_LADs.rds")

# ----------------------------------------------------------
# Load spatial boundaries
# ----------------------------------------------------------

LADs <- st_read(lads_path, quiet = TRUE)

# ----------------------------------------------------------
# Load Built-Up Areas (England + Wales)
# ----------------------------------------------------------

england <- read_excel(buas_path, sheet = 1) %>%
  mutate(country = "England")

wales <- read_excel(buas_path, sheet = 2) %>%
  mutate(country = "Wales")


buas <- bind_rows(england, wales) %>%
  filter(Counts >= 100000)

# Remove large conurbations handled separately if required
buas <- buas %>%
  filter(!`BUA name` %in% c(
    "Bath", "Birmingham", "Bradford", "Bristol",
    "Newcastle upon Tyne", "Portsmouth", "Sheffield"
  ))

# ----------------------------------------------------------
# Clean city names for matching
# ----------------------------------------------------------

buas <- buas %>%
  mutate(
    `BUA name` = str_trim(`BUA name`),
    BUA_clean = str_replace(`BUA name`, " \\(.*\\)", ""),
    BUA_clean = case_when(
      `BUA name` == "Southend-on-Sea" ~ "Southend",
      `BUA name` == "St Helens (St. Helens)" ~ "Saint Helens",
      `BUA name` == "Swansea" ~ "Abertawe",
      `BUA name` == "Kingswood and Fishponds" ~ "Kingswood",
      TRUE ~ BUA_clean
    )
  )

# ----------------------------------------------------------
# Load city coordinates
# ----------------------------------------------------------

cities <- read_csv(cities_path, show_col_types = FALSE) %>%
  mutate(
    city = str_trim(city),
    city = case_when(
      city == "Caerdydd" ~ "Cardiff",
      city == "Brighton" ~ "Brighton and Hove",
      TRUE ~ city
    )
  )

# ----------------------------------------------------------
# Join coordinates to BUAs
# ----------------------------------------------------------

cities_joined <- buas %>%
  left_join(cities, by = c("BUA_clean" = "city")) %>%
  select(
    `BUA name`,
    population = Counts,
    lat,
    lng,
    country 
  )

# ----------------------------------------------------------
# Convert to sf and match to LAD boundaries
# ----------------------------------------------------------

cities_sf <- st_as_sf(
  cities_joined,
  coords = c("lng", "lat"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(st_crs(LADs))

cities_with_LAD <- st_join(
  cities_sf,
  LADs,
  join = st_within,
  left = TRUE
) %>%
  st_drop_geometry() %>%
  select(
    `BUA name`,
    population,
    lat,
    lng,
    country,
    LAD24CD,
    LAD24NM
  )

# ----------------------------------------------------------
# Add Scotland settlements (>=100k)   
# form https://www.nrscotland.gov.uk/publications/population-estimates-for-settlements-and-localities-in-scotland-mid-2020/
# ----------------------------------------------------------

scotland <- read_excel(scotland_pop_path) %>%
  filter(Sex == "All") %>%
    select(
    `Settlement name`,
    `Settlement code`,
    population = `All ages`
      ) %>%
  rename(
    `BUA name` = `Settlement name`,
    LAD24NM    = `Settlement name`,
    LAD24CD    = `Settlement code`
  ) %>%
  filter(population >= 100000) %>% 
  mutate(country = "Scotland")
# ----------------------------------------------------------
# Combine all countries
# ----------------------------------------------------------

final_dataset <- bind_rows(
  cities_with_LAD,
  scotland
)

# ----------------------------------------------------------
# Save output
# ----------------------------------------------------------

saveRDS(final_dataset, output_path)


