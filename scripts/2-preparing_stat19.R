# ======================================================================
# PREPARE STATS19 INJURY-LEVEL DATASET (NO DUPLICATION)
# ======================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)

# ----------------------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------------------

collisions  <- read_csv(here("data","raw",
                             "dft-road-casualty-statistics-collision-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

vehicles    <- read_csv(here("data","raw",
                             "dft-road-casualty-statistics-vehicle-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

casualties  <- read_csv(here("data","raw",
                             "dft-road-casualty-statistics-casualty-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

names(collisions)  <- tolower(names(collisions))
names(vehicles)    <- tolower(names(vehicles))
names(casualties)  <- tolower(names(casualties))

# ----------------------------------------------------------------------
# 2. Filter collisions to 2015+
# ----------------------------------------------------------------------

collisions <- collisions %>%
  mutate(date = dmy(date)) %>%
  filter(year(date) >= 2015)

# ----------------------------------------------------------------------
# 3. Select required columns (memory efficient)
# ----------------------------------------------------------------------

collisions_sel <- collisions %>%
  select(
    collision_index,
    date,
    location_easting_osgr,
    location_northing_osgr,
    first_road_class,
    second_road_class,
    speed_limit,
    urban_or_rural_area
  )

vehicles_sel <- vehicles %>%
  select(
    collision_index,
    vehicle_reference,
    vehicle_type,
    propulsion_code,
    car_passenger
  ) %>%
  distinct(collision_index, vehicle_reference, .keep_all = TRUE)

casualties_sel <- casualties %>%
  select(
    collision_index,
    vehicle_reference,
    casualty_reference,
    casualty_severity,
    casualty_type
  ) %>%
  distinct(collision_index, casualty_reference, .keep_all = TRUE)

# ----------------------------------------------------------------------
# 4. Build TRUE injury-level dataset
# ----------------------------------------------------------------------

injuries <- casualties_sel %>%
  left_join(collisions_sel, by = "collision_index") %>%
  left_join(vehicles_sel,
            by = c("collision_index","vehicle_reference")) %>%
  mutate(injury_id = paste0(collision_index,"_",casualty_reference))

stopifnot(sum(duplicated(injuries$injury_id)) == 0)

# ----------------------------------------------------------------------
# 5. Recode variables
# ----------------------------------------------------------------------

safe_int <- function(x) suppressWarnings(as.integer(x))

injuries <- injuries %>%
  mutate(
    # Road class
    first_road_class_num = safe_int(first_road_class),
    first_road_class_label = case_when(
      first_road_class_num == 1 ~ "Motorway",
      first_road_class_num %in% c(2,3) ~ "A",
      first_road_class_num == 4 ~ "B",
      first_road_class_num == 5 ~ "C",
      first_road_class_num == 6 ~ "Unclassified",
      TRUE ~ NA_character_
    ),
    first_road_class_major_minor = case_when(
      first_road_class_label %in% c("Motorway","A","B") ~ "Major",
      first_road_class_label %in% c("C","Unclassified") ~ "Minor",
      TRUE ~ NA_character_
    ),
    
    # Severity
    casualty_severity_num = safe_int(casualty_severity),
    casualty_severity_label = case_when(
      casualty_severity_num %in% c(1,2) ~ "KSI",
      casualty_severity_num == 3 ~ "Slight",
      TRUE ~ NA_character_
    ),
    
    # Propulsion
    propulsion_code_num = safe_int(propulsion_code),
    propulsion_label = case_when(
      propulsion_code_num == 1 ~ "Petrol",
      propulsion_code_num == 2 ~ "Diesel",
      propulsion_code_num %in% c(3,11) ~ "Electric",
      propulsion_code_num %in% c(8,12) ~ "Hybrid",
      propulsion_code_num == -1 ~ "Unknown",
      TRUE ~ "Other"
    ),
    
    # Vehicle type simplified
    vehicle_type_num = safe_int(vehicle_type),
    vehicle_type_label = case_when(
      vehicle_type_num == 1 ~ "Cycle",
      vehicle_type_num %in% c(2,3,4,5) ~ "Motorcycle",
      vehicle_type_num %in% c(8,9,108,109) ~ "Car/Taxi",
      vehicle_type_num %in% c(20,21) ~ "HGV",
      TRUE ~ "Other"
    )
  ) %>%
  select(-ends_with("_num"))

# ----------------------------------------------------------------------
# 6. Convert to sf (BNG)
# ----------------------------------------------------------------------

injuries_sf <- injuries %>%
  filter(!is.na(location_easting_osgr),
         !is.na(location_northing_osgr)) %>%
  st_as_sf(coords = c("location_easting_osgr",
                      "location_northing_osgr"),
           crs = 27700)

# ----------------------------------------------------------------------
# 7. Attach LAD safely (NO duplication)
# ----------------------------------------------------------------------

LADs <- st_read(here("data","raw","LAD_DEC_24_UK_BGC.shp"), quiet = TRUE) %>%
  st_transform(27700)

LADs_sub <- readRDS(here("data","processed","LADs_sub.rds"))

LADs_filtered <- LADs %>%
  filter(LAD24CD %in% LADs_sub$LAD24CD)

injuries_sf <- st_join(injuries_sf,
                       LADs_filtered,
                       join = st_within,
                       left = TRUE)

# ----------------------------------------------------------------------
# 8. Final integrity check
# ----------------------------------------------------------------------

cat("Rows:", nrow(injuries_sf), "\n")
cat("Unique injury_id:", length(unique(injuries_sf$injury_id)), "\n")
cat("Duplicates:", sum(duplicated(injuries_sf$injury_id)), "\n")

stopifnot(sum(duplicated(injuries_sf$injury_id)) == 0)

# ----------------------------------------------------------------------
# 9. Save
# ----------------------------------------------------------------------

saveRDS(injuries_sf,
        here("data","processed","injuries_final.rds"))
