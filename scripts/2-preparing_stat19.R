# ======================================================================
# PREPARE STATS19 INJURY-LEVEL DATASET
# ======================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)

# ----------------------------------------------------------------------
# data
# ----------------------------------------------------------------------

collisions  <- read_csv(here("../stat19-os-road-linkage-data",
                             "dft-road-casualty-statistics-collision-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

vehicles    <- read_csv(here("../stat19-os-road-linkage-data",
                             "dft-road-casualty-statistics-vehicle-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

casualties  <- read_csv(here("../stat19-os-road-linkage-data",
                             "dft-road-casualty-statistics-casualty-1979-latest-published-year.csv"),
                        show_col_types = FALSE)

names(collisions)  <- tolower(names(collisions))
names(vehicles)    <- tolower(names(vehicles))
names(casualties)  <- tolower(names(casualties))

# -------------------------------------------------------
# Filter to 2015+
# --------------------------------------------------------

collisions <- collisions %>%
  mutate(date = dmy(date)) %>%
  filter(year(date) >= 2015)

### keep injuries for thos collisions 
casualties <- casualties %>%
  filter(collision_index %in% collisions$collision_index)

### and vehiculs 

vehicles <- vehicles %>%
  filter(collision_index %in% collisions$collision_index)


### Ensure keys unique BEFORE join
# ----------------------------------------------------------------------

nrow(casualties) - (casualties %>%
  distinct(collision_index, casualty_reference, .keep_all = TRUE) %>%  nrow()
)
casualties <- casualties %>%
  distinct(collision_index, casualty_reference, .keep_all = TRUE)

# ----------------------------------------------------------------------
# Build injury-level dataset
# ----------------------------------------------------------------------

stats19 <- casualties %>%
  left_join(collisions, by = "collision_index") %>%
  left_join(vehicles,
            by = c("collision_index", "vehicle_reference")) %>%
  mutate(injury_id = paste0(collision_index, "_", casualty_reference))


### ----------------------------------------------------------------------
# Recode variables
# ----------------------------------------------------------------------
# safe integer conversion
safe_int <- function(x) {
  if (is.null(x)) return(NULL)
  suppressWarnings(as.integer(x))
}
stats19 <- stats19 %>%
  mutate(
    # =========================
    # FIRST ROAD CLASS
    # =========================
    first_road_class_num = safe_int(first_road_class),
    
    first_road_class_label1 = case_when(
      first_road_class_num %in% c(2, 3) ~ "A",
      first_road_class_num == 1 ~ "Motorway",
      first_road_class_num == 4 ~ "B",
      first_road_class_num == 5 ~ "C",
      first_road_class_num == 6 ~ "Unclassified",
      first_road_class_num == -1 ~ "Data missing or out of range",
      TRUE ~ NA_character_
    ),
    
    first_road_class1 = case_when(
      first_road_class_label1 %in% c("Motorway", "A", "B") ~ "Major",
      first_road_class_label1 %in% c("C", "Unclassified") ~ "Minor",
      first_road_class_label1 == "Data missing or out of range" ~ "Data missing",
      TRUE ~ NA_character_
    ),
    
    # =========================
    # SECOND ROAD CLASS
    # =========================
    second_road_class_num = safe_int(second_road_class),
    
    second_road_class_label1 = case_when(
      second_road_class_num == 0 ~ "Not at junction or within 20 metres",
      second_road_class_num == 1 ~ "Motorway",
      second_road_class_num %in% c(2, 3) ~ "A",
      second_road_class_num == 4 ~ "B",
      second_road_class_num == 5 ~ "C",
      second_road_class_num == 6 ~ "Unclassified",
      second_road_class_num == 9 ~ "Unknown (self rep only)",
      second_road_class_num == -1 ~ "Data missing or out of range",
      TRUE ~ NA_character_
    ),
    
    second_road_class1 = case_when(
      second_road_class_label1 %in% c("Motorway", "A", "B") ~ "Major",
      second_road_class_label1 %in% c("C", "Unclassified") ~ "Minor",
      second_road_class_label1 %in% c("Unknown (self rep only)", 
                                      "Data missing or out of range") ~ "Data missing",
      TRUE ~ NA_character_
    ),
    
    # Junction indicator
    atJunction = case_when(
      is.na(second_road_class_num) ~ NA_integer_,
      second_road_class_num == 0 ~ 0L,
      TRUE ~ 1L
    ),
    
    # =========================
    # CASUALTY SEVERITY
    # =========================
    casualty_severity_num = safe_int(casualty_severity),
    
    casualty_severity1 = case_when(
      casualty_severity_num %in% c(1, 2) ~ "KSI",
      casualty_severity_num == 3 ~ "Slight",
      TRUE ~ NA_character_
    ),
    
    # =========================
    # CASUALTY TYPE
    # =========================
    casualty_type_num = safe_int(casualty_type),
    
    casualty_type1 = case_when(
      casualty_type_num == 1 ~ "Cyclist",
      casualty_type_num == 0 ~ "Pedestrian",
      casualty_type_num %in% c(9, 8, 19, 10) ~ "Car or van driver or occupant",
      car_passenger %in% c(1, 2) ~ "Car or van driver or occupant",
      TRUE ~ "Other"
    ),
    
    # =========================
    # PROPULSION CODE
    # =========================
    propulsion_code_num = safe_int(propulsion_code),
    
    propulsion_code1 = case_when(
      propulsion_code_num == 1 ~ "Petrol",
      propulsion_code_num == 2 ~ "Diesel",
      propulsion_code_num %in% c(3, 11) ~ "Electric",
      propulsion_code_num %in% c(8, 12) ~ "Hybrid",
      propulsion_code_num == -1 ~ "Unknown",
      TRUE ~ "Other"
    ),
    
    # =========================
    # VEHICLE TYPE
    # =========================
    vehicle_type_num = safe_int(vehicle_type),
    
    vehicle_type1 = case_when(
      vehicle_type_num == 1 ~ "Pedal cycle",
      vehicle_type_num %in% c(2,3,4,5,23,97,103,104,105,106) ~ "Motorcycle",
      vehicle_type_num %in% c(8,9,108,109) ~ "Car/Taxi",
      vehicle_type_num %in% c(10,110) ~ "Minibus",
      vehicle_type_num == 11 ~ "Bus/Coach",
      vehicle_type_num == 19 ~ "Van",
      vehicle_type_num %in% c(20,21,98,113) ~ "HGV",
      vehicle_type_num == 16 ~ "Ridden horse",
      vehicle_type_num == 17 ~ "Agricultural vehicle",
      vehicle_type_num == 18 ~ "Tram",
      vehicle_type_num == 22 ~ "Mobility scooter",
      vehicle_type_num %in% c(90,99) ~ "Other vehicle",
      vehicle_type_num == -1 ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    vehicle_type1_simplified = case_when(
      vehicle_type1 == "Pedal cycle" ~ "Cycle",
      vehicle_type1 == "Motorcycle" ~ "Motorcycle",
      vehicle_type1 == "Car/Taxi" ~ "Car/Taxi",
      vehicle_type1 %in% c("Minibus", "Bus/Coach") ~ "Minibus/Bus/Coach",
      vehicle_type1 == "Van" ~ "Van",
      vehicle_type1 == "HGV" ~ "HGV",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Remove only helper numeric columns
  select(
    -first_road_class_num,
    -second_road_class_num,
    -casualty_severity_num,
    -casualty_type_num,
    -propulsion_code_num,
    -vehicle_type_num
  )
# ----------------------------------------------------------------------
# Convert to sf (BNG)
# ----------------------------------------------------------------------

injuries_sf <- stats19 %>%
  filter(!is.na(location_easting_osgr),
         !is.na(location_northing_osgr)) %>%
  st_as_sf(coords = c("location_easting_osgr",
                      "location_northing_osgr"),
           crs = 27700)
nrow(injuries_sf)
nrow(stats19)


### miswsing coords
collisions_missing <- collisions %>%
  filter(is.na(location_easting_osgr) | is.na(location_northing_osgr))

nrow(collisions_missing)

# ----------------------------------------------------------------------
# Attach LAD 
# ----------------------------------------------------------------------

LADs <- st_read(here("data","raw","LAD_DEC_24_UK_BGC.shp"), quiet = TRUE) %>%
  st_transform(27700)

LADs_sub <- readRDS(here("data","processed","LADs_sub.rds"))

LADs_filtered <- LADs %>%
  filter(LAD24CD %in% LADs_sub$LAD24CD)


injuries_sf <- st_join(injuries_sf,
                       LADs_filtered,
                       join = st_intersects,
                       left = F)


st_write(LADs,
         here("data","processed","shp_files","LADs.shp"))

st_write(LADs_filtered,
         here("data","processed","shp_files","LADs_filtered.shp"))



#filter injuries on the LADs subset 
injuries_sf <- injuries_sf  %>%
  filter(LAD24CD %in% LADs_sub$LAD24CD)


# ----------------------------------------------------------------------
#  checks
# ----------------------------------------------------------------------

cat("Rows in the LADs subset:", nrow(injuries_sf), "\n")
cat("Rows in all LADs:", nrow(stats19), "\n")
cat("Unique injury_id:", length(unique(injuries_sf$injury_id)), "\n")
cat("Duplicates:", sum(duplicated(injuries_sf$injury_id)), "\n")

stopifnot(sum(duplicated(injuries_sf$injury_id)) == 0)


# Save


saveRDS(injuries_sf,
        here("data","processed","injuries_final.rds"))
