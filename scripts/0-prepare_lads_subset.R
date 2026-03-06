# ==========================================================
# Script: 00_create_large_city_LAD_subset.R
# Purpose: Create a subset of large UK cities (>=100k population) without the CAZ
#          and assign the corresponding Local Authority District (LAD24)
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

lads_path          <- here("../stat19-os-road-linkage-data", "LAD_DEC_24_UK_BGC.shp")
buas_path          <- here("../stat19-os-road-linkage-data", "BUAs.xlsx")
scotland_pop_path  <- here("../stat19-os-road-linkage-data", "scot_pop.xlsx")
cities_path        <- here("../stat19-os-road-linkage-data", "ukcities.csv")

dir.create(here("data", "processed"), showWarnings = FALSE)

# ----------------------------------------------------------
# Load spatial boundaries
# ----------------------------------------------------------

LADs <- st_read(lads_path, quiet = TRUE)

# ----------------------------------------------------------
# Load Built-Up Areas (England + Wales)
# ----------------------------------------------------------

england <- read_excel(buas_path, sheet = 1) %>%
  mutate(country = "England")

### wales <- read_excel(buas_path, sheet = 2) %>%   mutate(country = "Wales")    ### decided to remove 


buas <- bind_rows(england) %>%
  filter(Counts >= 100000)
buas$`BUA name`

# Remove cities with CAZ -- handled separately in the next script
buas <- buas %>%
  filter(!`BUA name` %in% c(
    "Bath", "Birmingham", "Bradford", "Bristol",
    "Newcastle upon Tyne", "Portsmouth", "Sheffield", "Oxford"
  ))

# -------------------------------------------------------------
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
    BUA_name=`BUA name`,
    population = Counts,
    lat,
    lng,
    country, id 
  )

# --------------------
# clean cities joined 
#--------------------
cities_joined %>%  View()
cities_joined_clean <- cities_joined %>%
  
  # Remove Kingswood 
  filter(BUA_name != "Kingswood and Fishponds") %>%
  
  #  Keep only the decided correct city IDs where duplicates exist
  filter(
    !(BUA_name == "Blackburn (Blackburn with Darwen)" & id != 1826802533),
    !(BUA_name == "Gillingham (Medway)" & id != 1826642243),
    !(BUA_name == "St Helens (St. Helens)" & id != 1826775771),
    !(BUA_name == "Swindon (Swindon)" & id != 1826498106),
    !(BUA_name == "Southend-on-Sea" & id != 1826524208)
  )


cities_joined_clean %>% select(BUA_name) %>%  n_distinct()

# ----------------------------------------------------------
# Convert to sf and match to LAD boundaries
# ----------------------------------------------------------

cities_sf <- st_as_sf(
  cities_joined_clean,
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
    `BUA_name`,
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
    `BUA_name` = `Settlement name`,
    LAD24NM    = `Settlement name`,
    LAD24CD    = `Settlement code`
  ) %>%
  filter(population >= 100000) %>% 
  mutate(country = "Scotland")

# Remove cities with LEZ -- handled separately in the next script
scotland$LAD24NM

scotland<- scotland %>%
  filter(!LAD24NM %in% c(
    "Aberdeen, Milltimber, and Peterculter", "Edinburgh", "Greater Glasgow"
  ))
# ---------------
# Combined
# ----------------

big_cities_with_LADs <- bind_rows(
  cities_with_LAD,
  scotland
)

# ----------------------------------------------------------
# Save output
# ----------------------------------------------------------

saveRDS(big_cities_with_LADs, here("data", "processed", "big_cities_with_LADs.rds"))



