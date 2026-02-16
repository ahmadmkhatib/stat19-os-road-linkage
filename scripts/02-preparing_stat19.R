# ======================================================================
# SCRIPT: PREPARE STATS19 DATA FOR ROAD LINK MATCHING
# ======================================================================
# PURPOSE
# -------
# - Import STATS19 collision, vehicle, and casualty datasets
# - Clean and harmonise variables
# - Filter to 2015+
# - Attach LAD geography
# - Remove London
# - Classify injuries by road hierarchy
# - Output RTIs_final.rds (EPSG:27700)
# ======================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)

# ----------------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------------

collisions_path  <- here("data", "raw",
                         "dft-road-casualty-statistics-collision-1979-latest-published-year.csv")

vehicles_path    <- here("data", "raw",
                         "dft-road-casualty-statistics-vehicle-1979-latest-published-year.csv")

casualties_path  <- here("data", "raw",
                         "dft-road-casualty-statistics-casualty-1979-latest-published-year.csv")

lad_path <- here("data", "raw", "LAD_DEC_24_UK_BGC.shp")

# ----------------------------------------------------------------------
# 2. Load data
# ----------------------------------------------------------------------

collisions <- read_csv(collisions_path, show_col_types = FALSE)
vehicles   <- read_csv(vehicles_path, show_col_types = FALSE)
casualties <- read_csv(casualties_path, show_col_types = FALSE)

names(collisions) <- tolower(names(collisions))
names(vehicles)   <- tolower(names(vehicles))
names(casualties) <- tolower(names(casualties))

# ----------------------------------------------------------------------
# Parse date and filter 2015+
# ----------------------------------------------------------------------

if (!"date" %in% names(collisions)) {
  stop("Date column not found in collisions dataset.")
}

collisions <- collisions %>%
  mutate(date = dmy(date)) %>%
  filter(year(date) >= 2015)

# ----------------------------------------------------------------------
# Join datasets (collision -> vehicle -> casualty)
# ----------------------------------------------------------------------

vehicles_trim <- vehicles %>%
  select(-any_of(names(collisions)[names(collisions) != "collision_index"]))

joined_cv <- collisions %>%
  left_join(vehicles_trim, by = "collision_index")

casualties_trim <- casualties %>%
  select(-any_of(names(joined_cv)[names(joined_cv) != "collision_index"]))

stats19 <- joined_cv %>%
  left_join(casualties_trim, by = "collision_index")

# ----------------------------------------------------------------------
# Road class recoding
# ----------------------------------------------------------------------

safe_int <- function(x) suppressWarnings(as.integer(x))

stats19 <- stats19 %>%
  mutate(
    first_road_class_num = safe_int(first_road_class),
    
    first_road_class_label = case_when(
      first_road_class_num == 1 ~ "Motorway",
      first_road_class_num %in% c(2,3) ~ "A",
      first_road_class_num == 4 ~ "B",
      first_road_class_num == 5 ~ "C",
      first_road_class_num == 6 ~ "Unclassified",
      TRUE ~ NA_character_
    ),
    
    injury_class = case_when(
      first_road_class_label %in% c("Motorway","A","B") ~ first_road_class_label,
      first_road_class_label %in% c("C","Unclassified") ~ "minor",
      TRUE ~ NA_character_
    )
  )

# ----------------------------------------------------------------------
# Convert to sf (British National Grid)
# ----------------------------------------------------------------------

if (all(c("location_easting_osgr","location_northing_osgr") %in% names(stats19))) {
  
  stats19_sf <- stats19 %>%
    filter(!is.na(location_easting_osgr),
           !is.na(location_northing_osgr)) %>%
    st_as_sf(
      coords = c("location_easting_osgr","location_northing_osgr"),
      crs = 27700
    )
  
} else {
  stop("OSGR coordinates not found.")
}

# ----------------------------------------------------------------------
#
# ----------------------------------------------------------------------
# LADs_sub from the previous script
stats19_sf <- st_join(stats19_sf, LADs_sub, join = st_within)

# filter on the LADs of intrest
stats19_sf <- stats19_sf %>%
  filter(LAD24CD %in% LADs_sub$LAD24CD)

# ----------------------------------------------------------------------
# Create injury ID
# ----------------------------------------------------------------------

stats19_sf <- stats19_sf %>%
  mutate(injury_id = row_number())

summary(stats19$date)

# ----------------------------------------------------------------------
# Save 
# ----------------------------------------------------------------------
saveRDS(
  stats19_sf,
  here("data","processed","RTIs_final.rds")
)
