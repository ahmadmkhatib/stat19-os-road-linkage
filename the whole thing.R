
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

### CAN YOU SEE THIS CHANGE? ####

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
# Load Built-Up Areas (England + Scotland)
# ----------------------------------------------------------

england <- read_excel(buas_path, sheet = 1) %>%
  mutate(country = "England")


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
  
  # first create new columns
  mutate(
    BUA_name = `Settlement name`,
    LAD24NM  = `Settlement name`,
    LAD24CD  = `Settlement code`
  ) %>%
  
  filter(population >= 100000) %>%
  
  mutate(
    country = "Scotland",
    
    # harmonise names for coordinate matching
    BUA_name = case_when(
      BUA_name == "Greater Glasgow" ~ "Glasgow",
      BUA_name == "Motherwell and Wishaw" ~ "Motherwell",
      TRUE ~ BUA_name
    )
  ) %>%
  
  left_join(
    cities %>%
      select(city, lat, lng),
    by = c("BUA_name" = "city")
  )

# Remove cities with LEZ -- handled separately in the next script
scotland$LAD24NM

scotland <- scotland %>%
  filter(!LAD24NM %in% c(
    "Aberdeen, Milltimber, and Peterculter", "Edinburgh", "Greater Glasgow",
    "Dundee"   # Dundee implemented an LEZ — exclude from control pool
  )) %>%
  # Replace NRS settlement codes (S200xxxxx) with proper LAD24CD codes (S12xxxxx)
  # so that downstream joins against LAD boundary files work correctly
  mutate(LAD24CD = case_when(
    LAD24NM == "Motherwell and Wishaw" ~ "S12000050",  # North Lanarkshire
    LAD24NM == "Falkirk"              ~ "S12000014",   # Falkirk
    TRUE                              ~ LAD24CD
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






# ==========================================================
## Purpose: Prepare OS Open Roads data within selected LADs
#          defined in Script 0 (large city subset)
# Output:  data/processed/roads_filtered.rds
# ==========================================================

library(sf)
library(tidyverse)
library(stringr)
library(here)

# ----------------------------------------------------------
#Paths
# ----------------------------------------------------------

#  must download:
# - OS Open Roads
# - LAD boundaries (ONS)

roads_path <- here("../stat19-os-road-linkage-data", "OS highways all.shp")
lads_path  <- here("../stat19-os-road-linkage-data", "LAD_DEC_24_UK_BGC.shp")
cities_path  <- here("data", "processed", "big_cities_with_LADs.rds")

# ----------------------------------------------------------
# Load LAD Boundaries
# ----------------------------------------------------------
## LADs with CAZ 

CAZ_LADs <- c(
  # CAZ
  "E06000022","E08000025","E08000032","E06000023",
  "E06000044","E08000019","E08000021","E08000018", 
  "S12000049","S12000036","S12000033", "S12000042","E07000178")

LADs <- st_read(lads_path, quiet = TRUE)


cities_with_LADs <- readRDS(cities_path)

# ----------------------------------------------------------
# Define LAD Selection 
# ----------------------------------------------------------

control_LADs <- cities_with_LADs %>%
  pull(LAD24CD)  %>%
  setdiff(CAZ_LADs) 


selected_lads <- c(
  #  Intervention
  CAZ_LADs,
  #  Controls
  control_LADs )


LADs_sub <- LADs %>%
  filter(LAD24CD %in% selected_lads)

LADs_sub <- LADs_sub %>% filter(!grepl("^W", LAD24CD))


### save the LADs sub 

saveRDS(LADs_sub, here("data", "processed", "LADs_sub.rds"))

lads_union <- st_union(LADs_sub)

# ----------------------------------------------------------
# Load OS Open Roads (only for the lads subset)
# ----------------------------------------------------------
#roads <- st_read(roads_path, quiet = TRUE)
# for speed ## dealing with .gpkg is faster than .shp 

#roads <- roads %>% st_zm( drop = TRUE, what = "ZM") %>% select(-fid)
#st_write(roads, "data/processed/roads.gpkg", delete_dsn = TRUE)
#rm(roads)

roads <- st_read(
  "data/processed/roads.gpkg",
  wkt_filter = st_as_text(lads_union),
  quiet = TRUE
)    # # #  this filters roads to only those in the LADs subset at any length 

summary(roads$length)
roads %>% filter(length>1000) %>%select(length) %>%  nrow()
roads %>% filter(length>3000) %>%select(length) %>%  nrow()
###!!!
# ----------------------------------------------------------
# Recode Road Classification
# ----------------------------------------------------------
table(roads$roadClassi)


roads <- roads %>%
  mutate(
    road_class = case_when(
      roadClassi == "A Road" ~ "A",
      roadClassi == "B Road" ~ "B",
      roadClassi == "Motorway" ~ "Motorway",
      roadClassi == "Classified Unnumbered" ~ "minor",
      roadClassi %in% c("Not Classified", "Unclassified", "Unknown") ~ "minor",
      TRUE ~ NA_character_
    )
  )




saveRDS(roads, here("data", "processed", "roads_filtered.rds"))



# ======================================================================
# PREPARE STATS19 INJURY-LEVEL DATASET
# ======================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)



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




collisions <- collisions %>%
  mutate(date = dmy(date)) %>%
  filter(year(date) >= 2015)

### keep injuries for thos collisions 
casualties <- casualties %>%
  filter(collision_index %in% collisions$collision_index)

### and vehiculs 

vehicles <- vehicles %>%
  filter(collision_index %in% collisions$collision_index)


### Ensure keys unique
# ----------------------------------------

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
      casualty_type_num %in% c(9, 8, 10, 16, 19) ~ "Car or van driver or occupant",
      car_passenger %in% c(1, 2) ~ "Car or van driver or occupant",
      TRUE ~ "Other"
    ),
    
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


LADs <- st_read( here("../stat19-os-road-linkage-data", "LAD_DEC_24_UK_BGC.shp"), quiet = TRUE) %>%
  st_transform(27700)

LADs_sub <- readRDS(here("data","processed","LADs_sub.rds"))

LADs_filtered <- LADs %>%
  filter(LAD24CD %in% LADs_sub$LAD24CD)


injuries_sf <- st_join(injuries_sf,
                       LADs_filtered,
                       join = st_intersects,
                       left = F)

st_write(LADs_filtered,
         here("data","processed","shp_files","LADs_filtered.shp"), append=FALSE)

sum(duplicated(injuries_sf$injury_id))

saveRDS(injuries_sf,
        here("data","processed","injuries_final.rds"))

# ==========================================================
# STAT19 – OS Open Roads Linkage Framework
# Script: 4_match_injuries_to_roads_by_type.R
# Purpose: Match STATS19 injuries subset to nearest plausible road link
# CRS: British National Grid (EPSG:27700)
# ==========================================================

library(sf)
library(tidyverse)
library(here)

# ----------------------------------------------------------
# Load Processed Data
# ----------------------------------------------------------
injuries <- readRDS(here("data", "processed", "injuries_final.rds"))
roads    <- readRDS( here("data", "processed", "roads_filtered.rds"))
names(roads)

table(injuries$LAD24NM)
injuries %>%
  filter(str_starts(LAD24CD, "S")) %>%
  distinct(LAD24CD, LAD24NM) %>%
  arrange(LAD24NM)
# ----------------------------------------------------------
# Split Roads by Class
# --------------------------------------------------------
table(roads$road_class)
table(injuries$first_road_class_label1)


#### recode class
injuries <- injuries %>% 
  mutate (road_class = 
            case_when( first_road_class_label1 %in% c("Motorway","A","B") ~ first_road_class_label1, 
                       first_road_class_label1 %in% c("C","Unclassified") ~ "minor", TRUE ~ NA_character_))

table(injuries$road_class)

roads_A        <- roads %>% filter(road_class == "A")
roads_B        <- roads %>% filter(road_class == "B")
roads_motorway <- roads %>% filter(road_class == "Motorway")
roads_minor    <- roads %>% filter(road_class == "minor")

# -----------------
# Matching Function
# ----

match_one_class <- function(injuries, roads_same, roads_any,
                            max_same = 50, max_any = 100) {
  
  # nearest ANY class
  idx_any <- st_nearest_feature(injuries, roads_any)
  dist_any <- as.numeric(
    st_distance(injuries, roads_any[idx_any, ], by_element = TRUE)
  )
  
  # nearest SAME class
  idx_same <- st_nearest_feature(injuries, roads_same)
  dist_same <- as.numeric(
    st_distance(injuries, roads_same[idx_same, ], by_element = TRUE)
  )
  
  injuries %>%
    mutate(
      matched_roadID = case_when(
        dist_same <= max_same ~ roads_same$identifier[idx_same],
        dist_any  <= max_any  ~ roads_any$identifier[idx_any],
        TRUE ~ NA_character_
      ),
      match_type = case_when(
        dist_same <= max_same ~ "same_class",
        dist_any  <= max_any  ~ "fallback_any",
        TRUE ~ "dropped"
      ),
      dist_same = dist_same,
      dist_any  = dist_any
    ) %>%
    filter(match_type != "dropped")
}
# --------------------------------
# Perform Matching by Injury Class
# --------------------------------

matched_A <- match_one_class(
  injuries %>% filter(road_class == "A"),
  roads_A,
  roads
)

matched_B <- match_one_class(
  injuries %>% filter(road_class == "B"),
  roads_B,
  roads
)

matched_motorway <- match_one_class(
  injuries %>% filter(road_class == "Motorway"),
  roads_motorway,
  roads
)

matched_minor <- match_one_class(
  injuries %>% filter(road_class == "minor"),
  roads_minor,
  roads
)

matched <- bind_rows(
  matched_A,
  matched_B,
  matched_motorway,
  matched_minor
)

# -----------------
#   Diagnostics 
# -------------

# total unique injuries entering matching
n_total <- injuries %>%
  distinct(injury_id) %>%
  nrow()

# unique injuries  matched
n_matched <- matched %>%
  distinct(injury_id) %>%
  nrow()

n_unmatched <- n_total - n_matched

summary_table <- tibble(
  injuries_total    = n_total,
  injuries_matched  = n_matched,
  injuries_unmatched = n_unmatched,
  pct_unmatched     = 100 * n_unmatched / n_total
)

print(summary_table)

table(matched$match_type)
matched %>%
  summarise(
    mean_dist = mean(dist_any),
    median_dist = median(dist_any),
    p90 = quantile(dist_any, 0.9),
    max_dist = max(dist_any)
  )

table(matched$LAD24NM)

matched %>%
  filter(str_starts(LAD24CD, "S")) %>%
  distinct(LAD24CD, LAD24NM) %>%
  arrange(LAD24NM)

saveRDS(matched, here("data", "processed", "injuries_matched.rds"))



####     KSI and Sevirity adjustment from stats method

library(tidyverse)
library(lubridate)
library(here)


# -------------------------------
injuries <- read_rds(here("data", "processed", "injuries_with_oa.rds"))
# -------------------------------
#  severity adjustment 
# --------------------------------

# STATS19 coding:
# 1 = Fatal
# 2 = Serious
# 3 = Slight
# KSI = Fatal + Serious

injuries <- injuries %>%
  mutate(
    # ----------------------------
    # Unadjusted indicators
    # ----------------------------
    KSI_unadj = as.numeric(casualty_severity %in% c(1, 2)),
    Slight_unadj = as.numeric(casualty_severity == 3),
    
    # ----------------------------
    # Adjusted (proper KSI)
    # ----------------------------
    KSI_adj = case_when(
      casualty_severity == 1 ~ 1,  # Fatal
      casualty_severity == 2 ~ casualty_adjusted_severity_serious,
      casualty_severity == 3 ~ 0
    ),
    
    Slight_adj = case_when(
      casualty_severity == 3 ~ casualty_adjusted_severity_slight,
      TRUE ~ 0
    )
  )
# -----------
#checks

stopifnot(
  all(injuries$KSI_adj >= 0 & injuries$KSI_adj <= 1),
  all(injuries$Slight_adj >= 0 & injuries$Slight_adj <= 1)
)

# :  comparison summary
national_totals <- injuries %>%
  summarise(
    KSI_unadj_total = sum(KSI_unadj, na.rm = TRUE),
    KSI_adj_total   = sum(KSI_adj, na.rm = TRUE),
    Slight_unadj_total = sum(Slight_unadj, na.rm = TRUE),
    Slight_adj_total   = sum(Slight_adj, na.rm = TRUE)
  )

print(national_totals)


injuries %>%
  filter(casualty_severity == 1) %>%
  summarise(mean_KSI_adj = mean(KSI_adj))

# # # #Adjusted KSI is very close to unadjusted KSI (difference = 20). 
##Slight decreases a lot after adjustment - - 


write_rds(
  injuries,
  here("data", "processed", "injuries_matched_final.rds")
)



#============================================================
# OA-Level Data 
#============================================================
# Purpose:
# Construct OA dataset defining treatment, buffer, and control
# areas.
## Classification uses majority-area (>=50%) spatial rules.
### output :  OA_level_from_polygons.rds
###           OA_subset.shp
#============================================================

library(tidyverse)
library(sf)
library(here)

lads_sub <- readRDS(here("data","processed","LADs_sub.rds"))

oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid()

caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_transform(st_crs(27700)) %>%
  st_make_valid()

st_crs(caz)
#============================================================
# Keep OAs within study LADs
#============================================================

##  keeping the LAD with the largest overlap
oa_sub <- st_join(
  oa,
  lads_sub %>% dplyr::select(LAD24CD, LAD24NM),
  join = st_intersects,
  largest = TRUE          # keeps only the LAD with greatest overlap area
) %>%
  filter(!is.na(LAD24CD))



st_crs(oa_sub)$epsg
st_crs(caz)$epsg
#============================================================
# Calculate OA area
#============================================================

oa_area <- oa_sub %>%
  mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  dplyr:: select(OA, total_area)

#============================================================
#  OA proportion inside CAZ
#============================================================
stopifnot(st_crs(oa_sub) == st_crs(caz))
oa_caz_intersection <- st_intersection(
  oa_sub %>% dplyr:: select(OA),
  caz %>% dplyr::select(scheme)
) %>%
  mutate(int_area = st_area(geometry)) %>%
  st_drop_geometry()

oa_caz_prop <- oa_caz_intersection %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_caz = as.numeric(int_area / total_area))

treated_OAs <- oa_caz_prop %>%
  filter(prop_caz >= 0.5) %>%
  arrange(OA, desc(prop_caz)) %>%   #
  distinct(OA, .keep_all = TRUE) %>%
  select(OA, scheme)

#============================================================
# Create 1km CAZ buffer
#============================================================
caz_buffer <- st_buffer(caz,1000) %>%
  st_difference(caz)

oa_buffer_intersection <- st_intersection(
  oa_sub %>%dplyr::  select(OA),
  caz_buffer
) %>%
  mutate(int_area = st_area(geometry)) %>%
  st_drop_geometry()

oa_buffer_prop <- oa_buffer_intersection %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_buffer = as.numeric(int_area / total_area))

buffer_OAs <- oa_buffer_prop %>%
  filter(prop_buffer >= 0.5) %>%
  pull(OA)

#============================================================
# Identify CAZ LADs
#============================================================

treated_LADs <- oa_sub %>%
  filter(OA %in% treated_OAs$OA) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

#============================================================
# Distance to city centre
#============================================================

city_centroids <- oa_sub %>%
  group_by(LAD24CD) %>%
  summarise(geometry = st_union(geometry), .groups="drop") %>%
  st_centroid()

oa_centroids <- st_centroid(oa_sub)

dist_citycentre <- st_distance(
  st_geometry(oa_centroids),
  city_centroids[match(oa_centroids$LAD24CD,
                       city_centroids$LAD24CD),]$geometry,
  by_element=TRUE
) %>%
  as.numeric()

#============================================================
# OA classification
####
#dups
oa_sub %>% st_drop_geometry() %>% count(OA) %>% filter(n > 1) %>% nrow()





OA_analysis <- oa_sub %>%
  st_drop_geometry() %>%
  dplyr::select(OA, LAD24CD, LAD24NM) %>%
  mutate(
    treated_OA = if_else(OA %in% treated_OAs$OA,1,0),
    
    buffer_OA = if_else(
      OA %in% buffer_OAs &
        treated_OA == 0,
      1,0
    ),
    
    control_group1_OA = if_else(
      LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    ),
    
    control_group2_OA = if_else(
      !LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    )
  )

OA_analysis$dist_citycentre <- dist_citycentre

nrow(OA_analysis)
OA_analysis %>% distinct(OA) %>% nrow()

# Are all oa_sub OAs actually in OA_analysis?
anti_join(
  oa_sub %>% st_drop_geometry() %>% dplyr::select(OA),
  OA_analysis, by = "OA"
) %>% nrow()


OA_analysis %>%  dplyr:: select(dist_citycentre) %>%   summary()
table(OA_analysis$treated_OA)

OA_analysis  %>% filter(treated_OA==1)%>% dplyr::select(dist_citycentre) %>%   summary()

## # is this resonable 10K 

OA_analysis <- OA_analysis %>%
  mutate(
    assignment = case_when(
      treated_OA == 1 ~ "treated",
      buffer_OA == 1 ~ "buffer",
      control_group1_OA == 1 ~ "control_same_city",
      control_group2_OA == 1 ~ "control_other_city",
      TRUE ~ NA_character_
    )
  )
nrow(treated_OAs)

OA_analysis <- OA_analysis %>%
  left_join(
    treated_OAs,
    by = "OA"
  )


OA_analysis %>%
  group_by(assignment) %>%
  summarise(
    max_dist = max(dist_citycentre),
    median_dist = median(dist_citycentre)
  )

table(OA_analysis$assignment)

#distri
OA_analysis %>%
  count(assignment) %>%
  mutate(pct = n / sum(n))

glimpse(OA_analysis)


saveRDS(
  OA_analysis,
  here("data","processed","OA_level_from_polygons.rds")
)

st_write(
  oa_sub,
  here("data","processed","shp_files","OA_subset.shp"),
  delete_dsn=TRUE
)


# ============================================================
# Assign Roads to Output Areas (OA)
# ============================================================
# Purpose:
# 1. Assign each road link to ONE Output Area
#    based on the largest share of road length.
# 2. Attach OA treatment assignment
# create road level data 
####  road_attributes_OA 
# 3. Create OA-level road characteristics dataset
# # # output : OA_roads_dataset.rds
# ============================================================

library(sf)
library(tidyverse)
library(here)


# Roads
roads <- readRDS(
  here("data","processed","roads_filtered.rds")
) %>%
  st_make_valid()

glimpse(roads)
mean(roads$length)

# OA shpfile 
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()

# OA classification dataset
OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")         # from 8
)

table(OA_analysis$assignment)

if (st_crs(roads) != st_crs(oa_sub)) {
  roads <- st_transform(roads, st_crs(oa_sub))
}

# ==================================############
#  Intersect Roads with OAs
# =================================

roads_oa <- st_intersection(
  roads,
  oa_sub %>%dplyr::  select(OA)
)

# length of road segment inside OA
roads_oa <- roads_oa %>%
  mutate(
    int_length = st_length(.)
  )

# ============================================================
#  Assign each road to ONE OA (largest overlap)
# ============================================================

roads_oa <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup()

# ===================================================
#Attach OA treatment classification
# =======================================

road_attributes_OA <- roads_oa %>%
  left_join(
    OA_analysis %>%
      dplyr::  select(OA, assignment,scheme),
    by = "OA"
  )


n_distinct(roads$identifier)
n_distinct(road_attributes_OA$identifier)

missing_roads <- roads %>%
  filter(!identifier %in% road_attributes_OA$identifier)

nrow(missing_roads)

# ==================
# road-level dataset
#========================

st_write(
  road_attributes_OA,
  here("data","processed","road_attributes_OA.gpkg"),
  delete_dsn = TRUE
)

# =====================================
# Aggregate road characteristics to OA
# =====================================
# ──# Each road contributes only to its dominant OA (largest overlap)

OA_roads <- roads_oa %>%          # 
  st_drop_geometry() %>%
  group_by(OA) %>%
  summarise(
    n_roads           = n_distinct(identifier),
    total_road_length = sum(as.numeric(int_length), na.rm = TRUE),
    n_A               = sum(road_class == "A",        na.rm = TRUE),
    n_B               = sum(road_class == "B",        na.rm = TRUE),
    n_motorway        = sum(road_class == "Motorway", na.rm = TRUE),
    n_minor           = sum(road_class == "minor",    na.rm = TRUE),
    .groups = "drop"
  )

# ── One-OA-per-road assignment (largest overlap) ──────────────────────────
# 
road_attributes_OA <- roads_oa %>%
  group_by(identifier) %>%
  slice_max(int_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(
    OA_analysis %>%
      arrange(OA) %>%
      distinct(OA, .keep_all = TRUE) %>%    
      select(OA, assignment),
    by = "OA"
  )

# attach OA classification
OA_roads <- OA_roads %>%
  left_join(
    OA_analysis,  by = "OA"  )

glimpse(OA_roads)
OA_roads %>% count(OA) %>% filter(n > 1) %>% nrow()  # no road belong to more than 1 OA
nrow(OA_roads)  

OA_missing_roads <- OA_analysis %>%
  filter(!OA %in% OA_roads$OA)

nrow(OA_missing_roads)

OA_missing_roads %>%
  count(assignment)


table(OA_roads$assignment)


OA_analysis %>%
  count(OA) %>%
  filter(n > 1)
n_distinct(OA_analysis$OA)
nrow(OA_analysis)



saveRDS(
  OA_roads,
  here("data","processed","OA_roads.rds")
)




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
  mutate(prop_in_caz = as.numeric(length_in_caz / road_length)) %>%
  left_join(
    caz %>%
      st_drop_geometry() %>%
      select(scheme, startDt) %>%
      distinct(scheme, .keep_all = TRUE) %>%
      mutate(caz_start_q =as.yearqtr(dmy(startDt))),
    by = "scheme"
  )


# ── Join CAZ proportion back to roads ───────────────────────
roads <- roads %>%
  left_join(
    roads_caz_props %>% select(identifier, scheme, prop_in_caz),
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

# ============================================================
# OA-Level Injury Dataset Construction
# ============================================================
#
# Purpose
# ------------------------------------------------------------
# aggregates road traffic injuries from the
# matched STATS19–road dataset to the Output Area (OA) level.
#
# Injuries are aggregated by:
#   • OA
#   • Quarter
#   • Casualty type
#   • Injury severity (adjusted)
#
# The script produces an OA × quarter panel dataset of injury
# counts that is later used to construct:
#
#   • Pre-treatment injury baselines
#   • Pre-treatment injury trends
#   • Matching variables for synthetic control construction
#
# Steps
# ------------------------------------------------------------
#
# Aggregate injuries by OA × quarter × casualty type.
# Use adjusted counts (KSI_adj and Slight_adj).
# Pivot the dataset to wide format.
# Get total injuries per OA × quarter across all
# casualty types and severities.
#
# Output
# ------------------------------------------------------------
# OA_injuries_quarterly.rds
#
# Dataset structure:
#   OA                : Output Area identifier
#   quarter_year      : Quarter of observation
#   *_KSI             : Killed or Seriously Injured counts
#   *_Slight          : Slight injury counts
#   total_injuries    : Total injuries across all types
#
# This dataset feeds into the construction of the
# OA_matching_dataset used for synthetic control matching.
#
# ============================================================

library(tidyverse)
library(lubridate)
library(here)


injuries <- readRDS(
  here("data", "processed","injuries_matched_final.rds")
)

# ------------------------------------------------------------
# Inspect casualty types
# ------------------------------------------------------------
table(injuries$casualty_type1)

# remove "/" from Car/Van for cleaner column names
injuries <- injuries %>%
  mutate(
    casualty_type1 = str_replace_all(casualty_type1, "/", "_")
  )

summary(injuries$KSI_adj)
summary(injuries$Slight_adj)

# ------------------------------------------------------------
# Aggregate injuries by OA × quarter × casualty type
# using adjusted injury counts
# ------------------------------------------------------------
OA_injuries_q <- injuries %>%
  group_by(OA, quarter_year, casualty_type1) %>%
  summarise(
    KSI = sum(KSI_adj, na.rm = TRUE),
    Slight = sum(Slight_adj, na.rm = TRUE),
    .groups = "drop"
  )

glimpse(OA_injuries_q)

# ------------------------------------------------------------
# Pivot to wide format: type + severity columns
# ------------------------------------------------------------
OA_injuries_wide <- OA_injuries_q %>%
  pivot_longer(
    cols = c(KSI, Slight),
    names_to = "severity",
    values_to = "n_injuries"
  ) %>%
  unite(type_severity, casualty_type1, severity) %>%
  pivot_wider(
    names_from = type_severity,
    values_from = n_injuries,
    values_fill = 0
  )

# ------------------------------------------------------------
# Compute total injuries per OA × quarter
# ------------------------------------------------------------
inj_cols <- c(
  "Car_Van_KSI", "Car_Van_Slight",
  "Cyclist_KSI", "Cyclist_Slight",
  "Pedestrian_KSI", "Pedestrian_Slight",
  "Other_KSI", "Other_Slight"
)

OA_injuries_wide <- OA_injuries_wide %>%
  mutate(
    total_injuries = rowSums(across(any_of(inj_cols)), na.rm = TRUE)
  )

glimpse(OA_injuries_wide)


saveRDS(
  OA_injuries_wide,
  here("data","processed","OA_injuries_quarterly.rds")
)



# ============================================================
# Create OA-Level Matching Dataset for Controls pool for the DiD Analysis
# ============================================================
# This script prepares a dataset of Output Areas (OAs) for use in
# constructing synthetic control units in the CAZ evaluation.
#
# Pre-treatment period:
#   - Treated OAs: all quarters before the CAZ implementation date.
#   - Control OAs: all available quarters (no treatment cutoff).
#   Pre-treatment periods vary in length across schemes due to
#   staggered CAZ adoption — 
#
# Zero-filling:
#   OA_injuries_quarterly only contains OA-quarter rows where at least
#   one injury occurred. Zero-injury quarters are absent from the raw data.
#   A balanced panel is constructed by cross-joining all OAs × all quarters
#   and filling missing counts with zero BEFORE computing baselines or trends.
#   Without this step, trend slopes and baseline means would be estimated on
#   non-zero quarters only, biasing both upward.
#
# Baseline injury levels:
#   Mean quarterly casualties per OA over the pre-treatment period for seven
#   outcomes: car/van KSI, car/van slight, cyclist KSI, cyclist slight,
#   pedestrian KSI, pedestrian slight, and total casualties.
#
# Pre-treatment injury trends:
#   Quarterly trend slopes estimated via quasi-Poisson GLM. Raw count trends
#   use a time index as the predictor. Rate-per-km trends use the same GLM
#   with log(road_length_km) as an offset — keeping raw integer counts as the
#   outcome rather than dividing by km (which would misuse Poisson on a ratio).
#   The slope coefficient (log-rate change per quarter) is extracted.
#   Returns NA where all counts are zero or fewer than 2 observations exist.
#
# Outputs assembled into OA_matching_dataset:
#   - OA characteristics (from OA_level_from_polygons)
#   - Road network characteristics (from OA_roads)
#   - Baseline injury levels (mean per quarter, pre-period, zero-filled)
#   - Pre-treatment trend slopes (quasi-Poisson GLM, raw counts)
#   - Baseline injury rates per road-km
#   - Pre-treatment trend slopes (quasi-Poisson GLM with offset, per-km)
#   - Road composition & network density
#
# Output: data/processed/OA_matching_dataset.rds
#         data/processed/shp_files/OA_matching_dataset.gpkg
# ============================================================

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(naniar)



OA_analysis <- readRDS(
  here("data", "processed", "OA_level_from_polygons.rds")             # from 8 
)
names(OA_analysis)
table(OA_analysis$assignment)



glimpse(OA_analysis)
table(OA_analysis$scheme, useNA = "always")
table(OA_analysis$treated_OA, useNA = "always")

OA_roads <- readRDS(
  here("data", "processed", "OA_roads.rds")  
)
glimpse(OA_roads)



OA_injuries <- readRDS(
  here("data", "processed", "OA_injuries_quarterly.rds")
)
glimpse(OA_injuries)


OA_injuries %>%
  mutate(
    computed_total = Car_Van_KSI + Car_Van_Slight + 
      Cyclist_KSI + Cyclist_Slight +
      Pedestrian_KSI + Pedestrian_Slight +
      Other_KSI + Other_Slight
  ) %>%
  summarise(
    total_matches_computed = mean(total_injuries == computed_total, na.rm = TRUE),
    n_discrepant           = sum(total_injuries != computed_total, na.rm = TRUE),
    max_diff               = max(abs(total_injuries - computed_total), na.rm = TRUE)
  )


OA_injuries %>%
  summarise(
    total_casualties    = sum(total_injuries, na.rm = TRUE),
    other_casualties    = sum(Other_KSI + Other_Slight, na.rm = TRUE),
    other_pct_of_total  = 100 * other_casualties / total_casualties
  )



oa_sub <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()

caz <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

# ── Build scheme → CAZ start date lookup ─────────────────────────────────────

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(caz_start_date = dmy(startDt)) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    .groups = "drop"
  )
print(caz_dates)



## OA_roads is OA with roads 
# OA_roads_clean right_joins to oa_sub to recover the 1,752 road-free OAs

OA_roads_clean <- OA_roads %>%
  group_by(OA) %>%
  summarise(
    n_roads           = sum(n_roads, na.rm = TRUE),
    total_road_length = sum(total_road_length, na.rm = TRUE),
    n_A               = sum(n_A, na.rm = TRUE),
    n_B               = sum(n_B, na.rm = TRUE),
    n_motorway        = sum(n_motorway, na.rm = TRUE),
    n_minor           = sum(n_minor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    oa_sub %>% st_drop_geometry() %>% select(OA),
    by = "OA"
  ) %>%
  mutate(across(
    c(n_roads, total_road_length, n_A, n_B, n_motorway, n_minor),
    ~ replace_na(.x, 0)
  ))


nrow(OA_roads_clean)                                          
OA_roads_clean %>% count(OA) %>% filter(n > 1) %>% nrow()    
n_distinct(OA_roads_clean %>% filter(n_roads == 0) %>% pull(OA))  # Oa without roads 

anti_join(
  oa_sub %>% st_drop_geometry() %>% dplyr::select(OA),
  OA_roads_clean,
  by = "OA"
) %>% nrow()


no_road_oas <- OA_roads %>%
  filter(n_roads == 0) %>%
  pull(OA)
nrow(no_road_oas)
n_distinct(no_road_oas)  # 

has_injuries_no_roads <- OA_injuries %>%
  filter(OA %in% no_road_oas) %>%
  distinct(OA)
nrow(has_injuries_no_roads)


# ── Build scheme/treatment lookup ─────────────────────────────────────────────

oa_scheme_lookup <- OA_analysis %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE) %>%
  dplyr::select(OA, scheme, treated_OA, control_group1_OA, control_group2_OA) %>%
  left_join(caz_dates, by = "scheme")



# where Are they concentrated 
OA_roads_clean %>%
  filter(n_roads == 0) %>%
  left_join(OA_analysis %>% select(OA, assignment), by = "OA") %>%
  count(assignment)


# ── Attach treatment info to injury panel ────────────────────────────────────
# Note: OA_injuries only has rows where injuries > 0  

OA_injuries <- OA_injuries %>%
  left_join(oa_scheme_lookup, by = "OA")


# ── Build balanced panel — fill zero-injury OA-quarters explicitly ────────────
# OA_injuries contains only quarters with at least one injury.
# Cross-joining all OAs × all quarters and filling NAs with 0 ensures:
#   (a) baseline means include zero-injury quarters in the denominator
#   (b) trend slopes are estimated over actual calendar time, not just
#       quarters where something happened

all_oas      <- oa_scheme_lookup %>%dplyr:: select(OA)   # all OAs, not just injured
all_quarters <- OA_injuries %>% distinct(quarter_year)

OA_injuries_balanced <- all_oas %>%
  cross_join(all_quarters) %>%
  left_join(OA_injuries, by = c("OA", "quarter_year")) %>%
  mutate(across(
    c(Car_Van_KSI, Car_Van_Slight, Cyclist_KSI, Cyclist_Slight,
      Pedestrian_KSI, Pedestrian_Slight, Other_KSI, Other_Slight,
      total_injuries),
    ~ replace_na(.x, 0)
  )) %>%
  # scheme/treatment info dropped by cross_join — re-attach
  dplyr::   select(-any_of(c("scheme", "treated_OA", "control_group1_OA", "control_group2_OA", "caz_start_date"))) %>%
  left_join(oa_scheme_lookup, by = "OA")

# Verify balance
cat("OAs:     ", n_distinct(OA_injuries_balanced$OA), "\n")
cat("Quarters:", n_distinct(OA_injuries_balanced$quarter_year), "\n")
cat("Rows:    ", nrow(OA_injuries_balanced), "\n")
cat("Expected:", n_distinct(OA_injuries_balanced$OA) *
      n_distinct(OA_injuries_balanced$quarter_year), "\n")

# ── Pre/post treatment indicator ──────────────────────────────────────────────

OA_injuries_balanced <- OA_injuries_balanced %>%
  mutate(
    period = case_when(
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_year <  caz_start_date ~ "pre",
      treated_OA == 1 & !is.na(caz_start_date) &
        quarter_year >= caz_start_date ~ "post",
      control_group1_OA == 1          ~ "pre",
      control_group2_OA == 1          ~ "pre",
      TRUE                            ~ NA_character_
    )
  )

OA_injuries_balanced %>%
  filter(treated_OA == 1) %>%
  count(period)

# ── Pre-treatment subset with calendar time index ─────────────────────────────

inj_pre <- OA_injuries_balanced %>%
  filter(period == "pre") %>%
  arrange(OA, quarter_year) %>%
  group_by(OA) %>%
  mutate(time = row_number()) %>%
  ungroup()

# ── Road lengths lookup ───────────────────────────────────────────────────────

road_lengths <- OA_roads_clean %>%
  dplyr:: select(OA, total_road_length) %>%
  mutate(
    road_length_km = if_else(
      is.na(total_road_length) | total_road_length == 0,
      NA_real_,
      total_road_length / 1000
    )
  )

# ── Quasi-Poisson GLM slope extractor ────────────────────────────────────────
# Returns log-rate change per quarter from a quasi-Poisson GLM.
# offset argument accepts log(road_length_km) for exposure-adjusted models.
# Returns NA if: fewer than 2 obs, all zeros, or constant non-zero counts.

poisson_slope <- function(y, x, offset = NULL) {
  
  df <- tibble(y = y, x = x,
               off = if (is.null(offset)) rep(0, length(x)) else offset) %>%
    filter(!is.na(y), !is.na(off))
  
  if (nrow(df) < 2)   return(NA_real_)
  if (all(df$y == 0)) return(NA_real_)
  if (var(df$y) == 0) return(NA_real_)
  
  model <- tryCatch(
    {
      if (is.null(offset)) {
        glm(y ~ x, family = quasipoisson(link = "log"), data = df)
      } else {
        glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
      }
    },
    warning = function(w) {
      if (grepl("did not converge", conditionMessage(w))) NULL else {
        withCallingHandlers(
          if (is.null(offset)) {
            glm(y ~ x, family = quasipoisson(link = "log"), data = df)
          } else {
            glm(y ~ x + offset(off), family = quasipoisson(link = "log"), data = df)
          },
          warning = function(w2) invokeRestart("muffleWarning")
        )
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(model)) return(NA_real_)
  coef(model)[["x"]]
}

# ── Pre-treatment trends — raw counts ────────────────────────────────────────

inj_trends <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI    = poisson_slope(Car_Van_KSI,       time),
    trend_car_slight = poisson_slope(Car_Van_Slight,    time),
    trend_cyc_KSI    = poisson_slope(Cyclist_KSI,       time),
    trend_cyc_slight = poisson_slope(Cyclist_Slight,    time),
    trend_ped_KSI    = poisson_slope(Pedestrian_KSI,    time),
    trend_ped_slight = poisson_slope(Pedestrian_Slight, time),
    trend_other_KSI    = poisson_slope(Other_KSI,       time),
    trend_other_slight = poisson_slope(Other_Slight,    time),
    trend_total      = poisson_slope(total_injuries,    time),
    .groups = "drop"
  )

# ── Baseline injury levels (mean quarterly, pre-period, zero-filled) ──────────

inj_baseline <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    mean_car_KSI    = mean(Car_Van_KSI,        na.rm = TRUE),
    mean_car_slight = mean(Car_Van_Slight,      na.rm = TRUE),
    mean_cyc_KSI    = mean(Cyclist_KSI,         na.rm = TRUE),
    mean_cyc_slight = mean(Cyclist_Slight,      na.rm = TRUE),
    mean_ped_KSI    = mean(Pedestrian_KSI,      na.rm = TRUE),
    mean_ped_slight = mean(Pedestrian_Slight,   na.rm = TRUE),
    mean_other_KSI    = mean(Other_KSI,        na.rm = TRUE),
    mean_other_slight = mean(Other_Slight,      na.rm = TRUE),
    mean_total      = mean(total_injuries,      na.rm = TRUE),
    .groups = "drop"
  )

# ── Per-km baselines ──────────────────────────────────────────────────────────

inj_per_km <- inj_baseline %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    mean_other_KSI_pkm    = mean_other_KSI    / road_length_km,
    mean_other_slight_pkm = mean_other_slight / road_length_km,
    mean_car_KSI_pkm    = mean_car_KSI    / road_length_km,
    mean_car_slight_pkm = mean_car_slight / road_length_km,
    mean_cyc_KSI_pkm    = mean_cyc_KSI    / road_length_km,
    mean_cyc_slight_pkm = mean_cyc_slight / road_length_km,
    mean_ped_KSI_pkm    = mean_ped_KSI    / road_length_km,
    mean_ped_slight_pkm = mean_ped_slight / road_length_km,
    mean_total_pkm      = mean_total      / road_length_km
  ) %>%
  dplyr::  select(OA, ends_with("_pkm"))

# ── Per-km trend slopes — quasi-Poisson GLM with offset ──────────────────────

inj_pre_offset <- inj_pre %>%
  left_join(road_lengths, by = "OA") %>%
  mutate(
    log_road_km = if_else(
      !is.na(road_length_km) & road_length_km > 0,
      log(road_length_km),
      NA_real_
    )
  )
inj_trends_pkm <- inj_pre_offset %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI_pkm     = poisson_slope(Car_Van_KSI,    time, log_road_km),
    trend_car_slight_pkm  = poisson_slope(Car_Van_Slight, time, log_road_km),
    trend_cyc_KSI_pkm     = poisson_slope(Cyclist_KSI,    time, log_road_km),
    trend_cyc_slight_pkm  = poisson_slope(Cyclist_Slight, time, log_road_km),
    trend_ped_KSI_pkm     = poisson_slope(Pedestrian_KSI,    time, log_road_km),
    trend_ped_slight_pkm  = poisson_slope(Pedestrian_Slight, time, log_road_km),
    trend_other_KSI_pkm   = poisson_slope(Other_KSI,    time, log_road_km),
    trend_other_slight_pkm = poisson_slope(Other_Slight, time, log_road_km),
    trend_total_pkm       = poisson_slope(total_injuries, time, log_road_km),
    .groups = "drop"
  )

# ── Road composition & network density ───────────────────────────────────────

oa_area <- oa_sub %>%
  st_make_valid() %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(OA, area_km2)

road_composition <- OA_roads_clean %>%
  select(OA, total_road_length, n_A, n_B, n_minor, n_roads) %>%
  left_join(oa_area, by = "OA") %>%
  mutate(
    total_road_length  = if_else(
      is.na(total_road_length) | total_road_length == 0,
      NA_real_, total_road_length),
    pct_A_road         = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_A     / n_roads, NA_real_),
    pct_B_road         = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_B     / n_roads, NA_real_),
    pct_minor_road     = if_else(
      !is.na(n_roads) & n_roads > 0, 100 * n_minor / n_roads, NA_real_),
    road_density_m_km2 = if_else(
      !is.na(area_km2) & area_km2 > 0,
      total_road_length / area_km2, NA_real_)
  ) %>%
  dplyr::  select(OA, pct_A_road, pct_B_road, pct_minor_road, road_density_m_km2)

# ── Assemble OA matching dataset ─────────────────────────────────────────────
# OA_analysis is the base — all other tables are one-row-per-OA so no
# row inflation can occur. OA_roads_clean has been stripped of any columns
# already in OA_analysis to prevent .x/.y suffix collisions.

OA_matching_dataset <- OA_analysis %>%
  arrange(OA, scheme) %>%
  distinct(OA, .keep_all = TRUE) %>%   # ← add this
  left_join(OA_roads_clean,   by = "OA") %>%
  left_join(inj_baseline,     by = "OA") %>%
  left_join(inj_trends,       by = "OA") %>%
  left_join(inj_per_km,       by = "OA") %>%
  left_join(inj_trends_pkm,   by = "OA") %>%
  left_join(road_composition, by = "OA") %>%
  mutate(road_length_km = total_road_length / 1000)

# Replace NAs with 0 for injury count variables only
# _pkm NAs kept — undefined where OA has no roads
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(across(starts_with("mean_"),  ~ replace_na(.x, 0))) %>%
  mutate(across(starts_with("trend_"), ~ replace_na(.x, 0)))

# ── Checks ────────────────────────────────────────────────────────────────────

# Unit counts by assignment group
OA_matching_dataset %>%
  count(treated_OA, control_group1_OA, control_group2_OA, buffer_OA) %>%
  print()

# No OAs missing from oa_sub
anti_join(
  oa_sub %>% st_drop_geometry() %>% select(OA),
  OA_matching_dataset, by = "OA"
) %>% nrow()

# Pre-period length by scheme —  using balanced panel ##

inj_pre %>%
  group_by(scheme, treated_OA) %>%
  summarise(
    n_OAs           = n_distinct(OA),
    min_quarters    = min(table(OA)),
    max_quarters    = max(table(OA)),
    median_quarters = median(as.numeric(table(OA))),
    .groups = "drop"
  ) %>%
  print()



# Short pre-period flag (should be 0 or very few after zero-filling)
short_pre <- inj_pre %>%
  filter(treated_OA == 1) %>%
  count(OA, name = "n_pre_quarters") %>%
  filter(n_pre_quarters < 4)

###### other

inj_pre %>%
  group_by(OA) %>%
  summarise(
    all_zero_other_KSI    = all(Other_KSI == 0),
    all_zero_other_slight = all(Other_Slight == 0),
    .groups = "drop"
  ) %>%
  summarise(
    pct_zero_other_KSI    = 100 * mean(all_zero_other_KSI),
    pct_zero_other_slight = 100 * mean(all_zero_other_slight)
  )

class(OA_injuries_balanced$quarter_year)
head(OA_injuries_balanced$quarter_year)
# Diagnostic version with calendar time — kept separate to avoid overwriting
# inj_pre (which uses row_number() time and was used for all trend estimation above)
inj_pre_cal <- OA_injuries_balanced %>%
  filter(period == "pre") %>%
  arrange(OA, quarter_year) %>%
  mutate(
    time_cal = as.numeric(quarter_year - min(quarter_year)) / 90
  )


cat("OAs with fewer than 4 pre-treatment quarters:", nrow(short_pre), "\n")

# Zero inflation in baseline injury variables — all three groups
inj_baseline %>%
  left_join(OA_analysis %>% select(OA, treated_OA, control_group1_OA, control_group2_OA),
            by = "OA") %>%
  group_by(treated_OA, control_group1_OA, control_group2_OA) %>%
  summarise(
    n                = n(),
    pct_zero_car_KSI = mean(mean_car_KSI == 0) * 100,
    pct_zero_cyc_KSI = mean(mean_cyc_KSI == 0) * 100,
    pct_zero_ped_KSI = mean(mean_ped_KSI == 0) * 100,
    pct_zero_total   = mean(mean_total   == 0) * 100,
    .groups = "drop"
  ) %>%
  print()

# NAs remaining in matching variables for treated OAs
OA_matching_dataset %>%
  filter(treated_OA == 1) %>%
  summarise(across(
    c(starts_with("trend_"), starts_with("mean_"), ends_with("_pkm")),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  filter(n_NA > 0) %>%
  arrange(desc(n_NA)) %>%
  print()

# Parallel trends plot — uses inj_pre directly (already has treated_OA etc.)
trend_plot <- inj_pre %>%
  mutate(group = case_when(
    treated_OA        == 1 ~ "Treated",
    control_group1_OA == 1 ~ "Control (Group 1)",
    control_group2_OA == 1 ~ "Control (Group 2)",
    TRUE                   ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(group, quarter_year) %>%
  summarise(mean_total = mean(total_injuries, na.rm = TRUE), .groups = "drop")

ggplot(trend_plot, aes(x = quarter_year, y = mean_total, colour = group)) +
  geom_line() +
  geom_point(size = 1) +
  labs(
    title = "Pre-treatment parallel trends check",
    x = "Quarter", y = "Mean total casualties per OA", colour = NULL
  ) +
  theme_minimal()

# 8. Baseline balance — treated vs control
OA_matching_dataset %>%
  filter(treated_OA == 1 | control_group2_OA == 1) %>%
  mutate(group = if_else(treated_OA == 1, "Treated", "Control")) %>%
  group_by(group) %>%
  summarise(
    n_OAs             = n(),
    mean_total_inj    = mean(mean_total,         na.rm = TRUE),
    mean_road_km      = mean(road_length_km,     na.rm = TRUE),
    mean_road_density = mean(road_density_m_km2, na.rm = TRUE),
    mean_pct_A        = mean(pct_A_road,         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# Duplicate OA check — must be zero
dupes <- OA_matching_dataset %>% count(OA) %>% filter(n > 1)
cat("Duplicate OAs:", nrow(dupes), "\n")
stopifnot(nrow(dupes) == 0)

#  Balanced panel verification
cat("Balanced panel rows:    ", nrow(inj_pre), "\n")
cat("Expected (n_OA × n_qtr):",
    n_distinct(inj_pre$OA) * n_distinct(inj_pre$quarter_year), "\n")



# Which OAs are in oa_scheme_lookup but never appeared in OA_injuries?
missing_oas <- anti_join(
  oa_scheme_lookup %>% select(OA),
  OA_injuries %>% select(OA),
  by = "OA"
)
nrow(missing_oas)
cat("OAs never in OA_injuries:", nrow(missing_oas), "\n")

# Are any of these treated?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  count(treated_OA, control_group2_OA)

# Where are these zero-injury treated OAs geographically?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  count(scheme) %>%
  arrange(desc(n))

# Do the OA codes look the same format?
missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  pull(OA) %>%
  head(10)

OA_injuries %>% pull(OA) %>% head(10)


missing_oas %>%
  left_join(oa_scheme_lookup, by = "OA") %>%
  filter(treated_OA == 1) %>%
  left_join(OA_roads_clean, by = "OA") %>%
  summarise(
    n_no_roads       = sum(is.na(n_roads) | n_roads == 0),
    n_has_roads      = sum(!is.na(n_roads) & n_roads > 0),
    mean_road_length = mean(total_road_length, na.rm = TRUE)
  )

## — these are CAZ OAs that happen to have no recorded STATS19 casualties across the whole study period.

#flag to excluded them from the matching pool — they cannot serve as useful donors or treated units 

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(zero_injury_OA = if_else(
    OA %in% (missing_oas %>% pull(OA)), 1L, 0L
  ))

# Summary of what gets flagged
OA_matching_dataset %>%
  count(treated_OA, control_group2_OA, zero_injury_OA) %>%
  arrange(treated_OA, control_group2_OA, zero_injury_OA)

# Effective analysis sample after flagging
OA_matching_dataset %>%
  filter(zero_injury_OA == 0) %>%
  count(treated_OA, control_group2_OA, buffer_OA)



miss_var_summary(OA_matching_dataset)


# replace the no roads OA with zeros in the 1752 NAs

OA_matching_dataset <- OA_matching_dataset    %>%   
  mutate(across(
    c(pct_A_road, pct_B_road, pct_minor_road, road_density_m_km2),
    ~ replace_na(.x, 0)
  ))



table(OA_matching_dataset$assignment)

saveRDS(
  OA_matching_dataset,
  here("data", "processed", "OA_matching_dataset.rds")
)
#OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_dataset.rds"))
names(OA_matching_dataset)


######################################################################################
#######################################################################################
####################################################################################
##################################################################################
##################################################################################
###########################################################################
# ── Variable descriptions ─────────────────────────────────────────────────────

var_description <- c(
  OA                   = "Output Area identifier (ONS 2011 OA code). Each row is one OA.",
  n_roads              = "Total number of road segments intersecting the OA.",
  total_road_length    = "Total length of roads within the OA (metres).",
  road_length_km       = "Total length of roads within the OA (kilometres).",
  n_A                  = "Number of A-road segments within the OA.",
  n_B                  = "Number of B-road segments within the OA.",
  n_motorway           = "Number of motorway segments within the OA.",
  n_minor              = "Number of minor/local road segments within the OA.",
  mean_car_KSI         = "Mean quarterly car/van KSI casualties, pre-treatment period (zero-filled).",
  mean_car_slight      = "Mean quarterly car/van slight casualties, pre-treatment period (zero-filled).",
  mean_cyc_KSI         = "Mean quarterly cyclist KSI casualties, pre-treatment period (zero-filled).",
  mean_cyc_slight      = "Mean quarterly cyclist slight casualties, pre-treatment period (zero-filled).",
  mean_ped_KSI         = "Mean quarterly pedestrian KSI casualties, pre-treatment period (zero-filled).",
  mean_ped_slight      = "Mean quarterly pedestrian slight casualties, pre-treatment period (zero-filled).",
  mean_total           = "Mean quarterly total casualties (all modes/severities), pre-treatment (zero-filled).",
  trend_car_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for car/van KSI, pre-treatment.",
  trend_car_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for car/van slight, pre-treatment.",
  trend_cyc_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for cyclist KSI, pre-treatment.",
  trend_cyc_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for cyclist slight, pre-treatment.",
  trend_ped_KSI        = "Quasi-Poisson GLM slope (log-rate/quarter) for pedestrian KSI, pre-treatment.",
  trend_ped_slight     = "Quasi-Poisson GLM slope (log-rate/quarter) for pedestrian slight, pre-treatment.",
  trend_total          = "Quasi-Poisson GLM slope (log-rate/quarter) for total casualties, pre-treatment.",
  mean_car_KSI_pkm     = "Mean quarterly car/van KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_car_slight_pkm  = "Mean quarterly car/van slight casualties per road-km, pre-treatment (zero-filled).",
  mean_cyc_KSI_pkm     = "Mean quarterly cyclist KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_cyc_slight_pkm  = "Mean quarterly cyclist slight casualties per road-km, pre-treatment (zero-filled).",
  mean_ped_KSI_pkm     = "Mean quarterly pedestrian KSI casualties per road-km, pre-treatment (zero-filled).",
  mean_ped_slight_pkm  = "Mean quarterly pedestrian slight casualties per road-km, pre-treatment (zero-filled).",
  mean_total_pkm       = "Mean quarterly total casualties per road-km, pre-treatment (zero-filled).",
  trend_car_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for car/van KSI rate.",
  trend_car_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for car/van slight rate.",
  trend_cyc_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for cyclist KSI rate.",
  trend_cyc_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for cyclist slight rate.",
  trend_ped_KSI_pkm    = "Quasi-Poisson GLM slope with log(road_km) offset for pedestrian KSI rate.",
  trend_ped_slight_pkm = "Quasi-Poisson GLM slope with log(road_km) offset for pedestrian slight rate.",
  trend_total_pkm      = "Quasi-Poisson GLM slope with log(road_km) offset for total casualty rate.",
  pct_A_road           = "Percentage of road segments in the OA that are A-roads.",
  pct_B_road           = "Percentage of road segments in the OA that are B-roads.",
  pct_minor_road       = "Percentage of road segments in the OA that are minor roads.",
  road_density_m_km2   = "Total road length per km² of OA area (metres per km²).",
  zero_injury_OA       = "OA with no injury"
  
)

tibble(
  variable    = names(var_description),
  description = unname(var_description)
) %>%
  print(n = Inf)



# ============================================================
# Merge OA Matching Dataset with Census Characteristics


library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(naniar)

### census OA data 
OA_char_raw <- read.csv(here("data","processed","outputArea_raw.csv"))
OA_char_percent <- read.csv(here("data","processed","outputArea_percent.csv"))

miss_var_summary(OA_char_raw)
miss_var_summary(OA_char_percent)


OA_char_percent <-  OA_char_percent %>% filter(!is.na(Bicycle))

### business  OA data 
OA__EW_businesses <- read.csv(here("data","processed","OA_EW_businesses.csv"))
OA__scot_businesses<- read.csv(here("data","processed","OA_Sco_businesses.csv"))


OA_EW_businesses_clean <- OA__EW_businesses %>%
  rename(OA = OA21CD)
OA_Sco_businesses_clean <- OA__scot_businesses %>%
  rename(OA = OA22) 

# Combine datasets
OA_businesses <- bind_rows(
  OA_EW_businesses_clean,
  OA_Sco_businesses_clean
)


OA_matching_dataset <- readRDS(
  here("data","processed","OA_matching_dataset.rds")
)

# 
# OA shapefile and compute area
# ------------------------------------------------------------

oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  mutate(
    area_m2 = as.numeric(st_area(.)),
    area_km2 = area_m2 / 1e6
  )

# --------------- 
# checks
# -------------- 

cat("OA_matching_dataset rows:", nrow(OA_matching_dataset), "\n")
cat("OA_char_raw rows:", nrow(OA_char_raw), "\n")
cat("OA_char_percent rows:", nrow(OA_char_percent), "\n")

cat("OA_char_raw dups:",
    OA_char_raw %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

cat("OA_char_percent dups:",
    OA_char_percent %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")
cat("OA__EW_businesses dups:",
    OA__EW_businesses %>% count(OA21CD) %>% filter(n > 1) %>% nrow(), "\n")

OA_businesses %>% count(OA) %>% filter(n > 1)%>% nrow()


OA_businesses %>% 
  count(OA, sort = TRUE)

OA_businesses %>% 
  filter(OA %in% (OA_businesses %>% count(OA) %>% filter(n > 1) %>% pull(OA))) %>%
  head(20)

OA_businesses <- OA_businesses %>% distinct()

OA_businesses %>% 
  count(OA) %>% 
  filter(n > 1)

### still 4 duplicates 
OA_businesses <- OA_businesses %>%
  group_by(OA) %>%
  summarise(
    business_retail_per_km2 = sum(business_retail_per_km2, na.rm = TRUE),
    business_accommodation_food_per_km2 = sum(business_accommodation_food_per_km2, na.rm = TRUE),
    .groups = "drop"
  )



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

# Merge census tables
# ------------------------

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_renamed, by="OA")

# ------------------------------------------------------------
# Merge census onto OA matching dataset
# ------------------------------------------------------------

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by="OA")


#'### ------- add busnesses 
OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by = "OA") %>%
  left_join(OA_businesses, by = "OA")


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



summary(OA_matching_census$pop_density)
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



naniar::miss_var_summary(OA_matching_census)

#### remoive the OA with nas 
OA_matching_census <- OA_matching_census %>%
  drop_na(
    White_pct, Mixed_pct, Asian_pct, Black_pct, Other_ethnicity_pct,
    X4under_pct, X5to9_pct, X10to14_pct, X15to19_pct
  )

OA_matching_census <- OA_matching_census %>%
  mutate(scheme = replace_na(scheme, "Control"))

# ============================================================
# PRE-MATCHING QA
# All checks must pass before saving. A hard stop is raised
# at the end if any fail.
# ============================================================

qa_pass <- TRUE
qa_fail <- function(msg) { cat("  FAIL \u2717", msg, "\n"); qa_pass <<- FALSE }
qa_ok   <- function(msg) { cat("  PASS \u2713", msg, "\n") }

cat("\n================================================\n")
cat("PRE-MATCHING QA\n")
cat("================================================\n")

# ── No duplicate OAs ───────────────────────────────────────────────────────
n_rows <- nrow(OA_matching_census)
n_oa   <- n_distinct(OA_matching_census$OA)
if (n_rows == n_oa) {
  qa_ok(sprintf("No duplicate OAs (%d rows)", n_rows))
} else {
  qa_fail(sprintf("%d rows but only %d distinct OA codes — duplicates present", n_rows, n_oa))
}

# ── Country  ─────────────────────────────────────────────────

country_tbl <- OA_matching_census |>
  mutate(.country = case_when(
    substr(LAD24CD, 1, 1) == "E" ~ "England",
    substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
    TRUE                         ~ "Unknown"
  )) |>
  count(.country)
print(country_tbl)
n_unknown <- country_tbl |> filter(.country == "Unknown") |> pull(n) |> sum()
if (n_unknown == 0) qa_ok("No unrecognised LAD24CD prefixes") else
  qa_fail(sprintf("%d OAs have unrecognised LAD24CD prefix (not E/S)", n_unknown))

# ── Assignment  counts by country ─────────────────────────────────────

OA_matching_census |>
  count(treated_OA, control_group1_OA, control_group2_OA, buffer_OA) |>
  arrange(treated_OA, control_group1_OA, control_group2_OA) |>
  print()

n_treated <- sum(OA_matching_census$treated_OA        == 1, na.rm = TRUE)
n_ctrl1   <- sum(OA_matching_census$control_group1_OA == 1, na.rm = TRUE)
n_ctrl2   <- sum(OA_matching_census$control_group2_OA == 1, na.rm = TRUE)
if (n_treated == 0) qa_fail("No treated OAs found")       else qa_ok(sprintf("Treated OAs:       %d", n_treated))
if (n_ctrl1   == 0) qa_fail("No control group 1 OAs")     else qa_ok(sprintf("Control group 1:   %d", n_ctrl1))
if (n_ctrl2   == 0) qa_fail("No control group 2 OAs")     else qa_ok(sprintf("Control group 2:   %d", n_ctrl2))

# injury profiles by group × country ─────────────────
# 

inj_profile <- OA_matching_census |>
  mutate(
    .country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Other"
    ),
    .group = case_when(
      treated_OA        == 1 ~ "Treated",
      control_group1_OA == 1 ~ "Control_1",
      control_group2_OA == 1 ~ "Control_2",
      buffer_OA         == 1 ~ "Buffer",
      TRUE                   ~ "Other"
    )
  ) |>
  group_by(.country, .group) |>
  summarise(
    n                    = n(),
    pct_zero_mean_total  = round(100 * mean(mean_total     == 0, na.rm = TRUE), 1),
    median_mean_total    = round(median(mean_total,              na.rm = TRUE), 4),
    pct_zero_trend_total = round(100 * mean(trend_total_pkm == 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(.country, .group)
print(inj_profile, n = Inf)

scot_c1 <- inj_profile |> filter(.country == "Scotland", .group == "Control_1")
if (nrow(scot_c1) == 0) {
  cat("  NOTE: No Scottish Control_1 OAs found.\n")
} else if (scot_c1$pct_zero_mean_total < 100) {
  qa_ok(sprintf(
    "Scottish Control_1 OAs: %.1f%% zero mean injuries (was 100%% before script-12 fix)",
    scot_c1$pct_zero_mean_total
  ))
} else {
  qa_fail("Scottish Control_1 OAs STILL ALL ZEROS — re-run script 12 with the fix applied")
}

eng_c1 <- inj_profile |> filter(.country == "England", .group == "Control_1")
if (nrow(eng_c1) > 0 && eng_c1$pct_zero_mean_total == 100)
  qa_fail("English Control_1 OAs are ALL ZEROS — script 12 fix may not be applied")

# Stage 1 variable completeness (treated OAs) ───────────────────────────
cat("\n[5] Stage 1 variable completeness (treated OAs)\n")

stage1_expected <- c(
  "road_density_m_km2", "road_length_km", "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)
missing_cols_s1 <- setdiff(stage1_expected, names(OA_matching_census))
if (length(missing_cols_s1) > 0)
  qa_fail(sprintf("Stage 1 columns absent from dataset: %s",
                  paste(missing_cols_s1, collapse = ", ")))

na_s1 <- OA_matching_census |>
  filter(treated_OA == 1) |>
  summarise(across(all_of(intersect(stage1_expected, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_s1) == 0) {
  qa_ok("No missing Stage 1 values in treated OAs")
} else {
  cat("  WARNING — Stage 1 NAs in treated OAs:\n"); print(na_s1)
  qa_fail("Missing Stage 1 values in treated OAs")
}

#  Stage 2 variable completeness (matching pool, roads > 0) ───────────────
cat("\n[6] Stage 2 variable completeness (matching pool, n_roads > 0)\n")

stage2_expected <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm", "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm", "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm", "trend_total_pkm",
  "mean_car_KSI_pkm", "mean_car_slight_pkm", "mean_cyc_KSI_pkm",
  "mean_cyc_slight_pkm", "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_other_KSI_pkm", "mean_other_slight_pkm", "mean_total_pkm"
)
na_s2 <- OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1,
         buffer_OA == 0, n_roads > 0) |>
  summarise(across(all_of(intersect(stage2_expected, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_s2) == 0) {
  qa_ok("No missing Stage 2 values in matching pool (OAs with roads)")
} else {
  cat("  INVESTIGATE — Stage 2 NAs in matching pool (pkm NAs unexpected for roaded OAs):\n")
  print(na_s2)
  qa_fail("Missing Stage 2 values for roaded OAs in matching pool")
}

# Road variable sanity ───────────────────────────────────────────────────


road_sanity <- OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1) |>
  summarise(
    n_negative_road_length = sum(road_length_km    < 0, na.rm = TRUE),
    n_negative_density     = sum(road_density_m_km2 < 0, na.rm = TRUE),
    n_pct_A_over_100       = sum(pct_A_road         > 100, na.rm = TRUE),
    n_zero_roads           = sum(n_roads == 0 | is.na(n_roads))
  )
print(road_sanity)
if (road_sanity$n_negative_road_length == 0 &&
    road_sanity$n_negative_density     == 0 &&
    road_sanity$n_pct_A_over_100       == 0) {
  qa_ok("Road variables have no impossible values")
} else {
  qa_fail("Impossible values in road variables — check OA_roads_clean pipeline")
}

# ── Per-km injury rate plausibility (treated, roads > 0) ──────────────────
cat("\n[8] Per-km injury rate plausibility (treated OAs with roads)\n")

pkm_sanity <- OA_matching_census |>
  filter(treated_OA == 1, n_roads > 0) |>
  summarise(
    n                    = n(),
    n_na_mean_total_pkm  = sum(is.na(mean_total_pkm)),
    n_negative_pkm       = sum(mean_total_pkm < 0,              na.rm = TRUE),
    median_mean_total    = round(median(mean_total_pkm,         na.rm = TRUE), 4),
    p99_mean_total       = round(quantile(mean_total_pkm, 0.99, na.rm = TRUE), 4),
    max_mean_total       = round(max(mean_total_pkm,           na.rm = TRUE), 4)
  )
print(pkm_sanity)
if (pkm_sanity$n_negative_pkm == 0)
  qa_ok("No negative per-km injury rates in treated OAs") else
    qa_fail("Negative per-km injury rates in treated OAs")
if (pkm_sanity$n_na_mean_total_pkm == 0)
  qa_ok("No NA mean_total_pkm for treated OAs with roads") else
    qa_fail(sprintf(
      "%d treated OAs with roads have NA mean_total_pkm", pkm_sanity$n_na_mean_total_pkm
    ))

# ──Zero-injury flag consistency ──────────────────────────────────────────
cat("\n[9] zero_injury_OA flag\n")

OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1) |>
  count(treated_OA, control_group1_OA, zero_injury_OA) |>
  print()

n_zero_treated <- sum(OA_matching_census$treated_OA   == 1 &
                        OA_matching_census$zero_injury_OA == 1, na.rm = TRUE)
cat(sprintf("  Zero-injury treated OAs (excluded from Analysis A): %d / %d (%.1f%%)\n",
            n_zero_treated, n_treated, 100 * n_zero_treated / max(n_treated, 1)))
qa_ok("zero_injury_OA flag present and consistent")

# ── Census variable completeness (treated OAs) ───────────────────────────
cat("\n[10] Census variable completeness (treated OAs)\n")

census_vars <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                 "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct",
                 "pop_density")
na_census <- OA_matching_census |>
  filter(treated_OA == 1) |>
  summarise(across(all_of(intersect(census_vars, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_census) == 0) qa_ok("No missing census variables in treated OAs") else {
  cat("  Missing census values in treated OAs:\n"); print(na_census)
}

# ── Final result ──────────────────────────────────────────────────────────────
cat("\n================================================\n")
if (qa_pass) {
  cat("ALL QA CHECKS PASSED — safe to save and proceed to matching\n")
} else {
  cat("ONE OR MORE CHECKS FAILED — fix issues above before matching\n")
}
cat("================================================\n\n")

stopifnot("Pre-matching QA failed — see output above" = qa_pass)

# ------------------------------------------------------------

saveRDS(
  OA_matching_census,
  here("data","processed","OA_matching_census.rds")
)

# OA_matching_census<- readRDS(here("data","processed","OA_matching_census.rds"))
names(OA_matching_census)


# ---------------
#  spatial version
# ---------------------------------

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by="OA") %>%
  st_as_sf()

st_write(
  OA_matching_census_sf,
  here("data","processed","shp_files","OA_matching_census.gpkg"),
  delete_dsn = TRUE
)







# ============================================================
# Panel Construction
# Road × Quarter Panel Dataset
# ============================================================

# This script constructs the analytical panel data 
# The treatment is defined at the road level:
#
#   A road is treated if:
#       • it lies inside a CAZ boundary
#       • AND the quarter occurs after CAZ implementation.
#
# 
# ------------------------------------------------------------
# CONTROL GROUP DEFINITIONS
# ------------------------------------------------------------
#
# Buffer Control (Spillover Zone)
# ------------------------------------------------------------
# Roads located within 1 km outside the CAZ boundary.
#
# Control Group 1: Same-City Controls
# ------------------------------------------------------------
# Roads located in the same CAZ cities but outside
# the 1 km CAZ buffer.
#
# These roads share:
#
#   • the same economic conditions
#   • the same road network
#   • the same city-specific policies
#
# This provides the **primary counterfactual**.
#
#
# Control Group 2:  Non-CAZ Cities
# ------------------------------------------------------------
# STEPS
# ------------------------------------------------------------
#
#   1. Identify treated roads inside CAZ
#   2. Construct a 1 km CAZ spillover buffer
#   3. Classify roads into treatment and control groups
#   4. Restrict dataset to relevant roads
#   5. Expand roads into a road × quarter panel
#   6. Aggregate injuries to road × quarter level
#   7. Merge treatment timing
## ------------------------------------------------------------
# OUTCOME VARIABLES
# ------------------------------------------------------------
#
# Road traffic injuries aggregated to road × quarter level.
#
# Separate outcomes are constructed for:
#
#   • KSI_adj       (Killed or Seriously Injured)
#   • Slight_adj    (Slight injuries)
#
# and by casualty type.
#
## ------------------------------------------------------------
# FINAL DATASET STRUCTURE
# ------------------------------------------------------------
#
# identifier        : Road identifier
# quarter_year      : Quarter of observation
# treated           : Post-treatment indicator
# treated_group     : Ever-treated road indicator
# buffer_control    : Roads within 1 km outside CAZ
# control_group1    : Same-city roads outside buffer
# control_group2    : Non-CAZ citys roads
#
# KSI_adj_*         : Severe injuries by road user type
# Slight_adj_*      : Slight injuries by road user type
# total_inj_adj_*   : Total injuries
#
# ============================================================


library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(zoo)
library(arrow)
library(units)

options(arrow.use_mmap = FALSE)


# ============================================================
#=================================================


# --data
# RTI matched data
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))      #from 7

# all roads with thier attribute  - @ road level 
road_attributes <- st_read(here("data","processed","road_attributes_OA.gpkg"))  # from 9


# roads inside the CAZs  == the treatment @ road level 
road_caz_prop<- readRDS(here("data", "processed", "roads_caz_props.rds"))  # from 10

glimpse(road_caz_prop)
#create treatment indicators 
road_caz_prop <- road_caz_prop %>%
  mutate(
    ever_treated_any   = 1,
    ever_treated_50pct = if_else(prop_in_caz >= 0.5, 1, 0)
  )


names(road_attributes)
names(road_caz_prop)


# -----------------------------

# - Treated Cities
# -----------------------------
treated_cities <- road_caz_prop %>%
  distinct(scheme) %>%
  pull(scheme)

# -----------------------------
# Control 1: Same city, outside CAZ
# -----------------------------
road_classification <- road_attributes %>%
  left_join(road_caz_prop %>% select(identifier, ever_treated_any, ever_treated_50pct), by="identifier") 

road_classification <- road_classification %>%
  mutate(
    ever_treated_any   = replace_na(ever_treated_any, 0),
    ever_treated_50pct = replace_na(ever_treated_50pct, 0),
    
    # Control groups for 50%+ definition (can also adapt for any)
    control_group1 = if_else(scheme %in% treated_cities & ever_treated_50pct == 0, 1, 0),
    control_group2 = if_else(!scheme %in% treated_cities & ever_treated_50pct == 0, 1, 0),
    
    treated_group_any   = ever_treated_any,
    treated_group_50pct = ever_treated_50pct,
    
    control_group3_mixed = if_else(control_group1 == 1 | control_group2 == 1, 1, 0)
  )


# Keep only relevant roads
analysis_roads <- road_classification %>%
  filter(treated_group_any == 1 | treated_group_50pct == 1 | control_group1 == 1 | control_group2 == 1) 

# -----------------------------
# Create Road × Quarter Panel
# -----------------------------
all_quarters <- unique(injuries$quarter_year)

road_panel <- expand_grid(
  identifier = analysis_roads$identifier,  # only analysis roads
  quarter_year = all_quarters
)

# -----------------------------
# Aggregate Injuries (KSI & Slight separate)
# -----------------------------
roadlevel_long <- injuries %>%
  group_by(identifier, quarter_year, casualty_type1) %>%
  summarise(
    KSI_adj = sum(KSI_adj, na.rm=TRUE),
    Slight_adj = sum(Slight_adj, na.rm=TRUE),
    KSI_unadj = sum(KSI_unadj, na.rm=TRUE),
    Slight_unadj = sum(Slight_unadj, na.rm=TRUE),
    total_inj_adj = sum(KSI_adj + Slight_adj, na.rm=TRUE),
    total_inj_unadj = sum(KSI_unadj + Slight_unadj, na.rm=TRUE),
    .groups = "drop"
  )

injury_wide <- roadlevel_long %>%
  pivot_wider(
    names_from = casualty_type1,
    values_from = c(KSI_adj, Slight_adj, KSI_unadj, Slight_unadj, total_inj_adj, total_inj_unadj),
    values_fill = 0
  )

# -----------------------------
# Merge panel + injuries + controls + treatment timing
# -----------------------------
road_panel_complete <- road_panel %>%
  left_join(injury_wide, by=c("identifier","quarter_year")) %>%
  mutate(across(starts_with(c("KSI","Slight","total_inj")), ~replace_na(.,0))) %>%
  left_join(analysis_roads, by="identifier") %>%
  left_join(road_caz_prop %>% select(identifier, caz_start_q), by="identifier") %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q),
    
    # Two post-treatment indicators
    treated_any   = if_else(treated_group_any == 1 & quarter_year >= caz_start_q, 1, 0),
    treated_50pct = if_else(treated_group_50pct == 1 & quarter_year >= caz_start_q, 1, 0)
  )

names(road_panel_complete)

road_panel_model <- road_panel_complete %>%
  st_drop_geometry() %>%
  mutate(
    # total across all modes — needed as primary outcome
    total_inj_adj_All   = rowSums(across(starts_with("total_inj_adj_")),   na.rm = TRUE),
    total_inj_unadj_All = rowSums(across(starts_with("total_inj_unadj_")), na.rm = TRUE),
    KSI_adj_All         = rowSums(across(starts_with("KSI_adj_")),         na.rm = TRUE),
    Slight_adj_All      = rowSums(across(starts_with("Slight_adj_")),      na.rm = TRUE)
  ) %>%
  select(
    identifier, quarter_year,
    treated_any, treated_50pct,
    treated_group_any, treated_group_50pct,
    scheme,
    caz_start_q,                          # ADD THIS
    control_group1, control_group2, control_group3_mixed,
    starts_with("KSI"), starts_with("Slight"), starts_with("total_inj")
  ) %>%
  rename_with(~ make.names(.x))







### add some road and city vars 

road_panel_model <- road_panel_model %>%
  rename_with(~make.names(.x))
names(road_panel_model)

write_dataset(
  road_panel_model,
  path = here("data","processed","road_panel_dataset"),
  format = "parquet"
)

saveRDS(analysis_roads, here("data", "processed", "analysis_roads.rds"))


# How many unique roads are in injuries vs panel?
length(unique(injuries$identifier))
length(unique(road_panel_model$identifier))

# How many injury rows are missing after join?
sum(!injuries$identifier %in% road_panel_model$identifier)






# road_panel_model<-arrow::open_dataset(here("data","processed","road_panel_dataset")) %>% collect()




road_panel_model <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

#  dimensions
cat("Roads:   ", n_distinct(road_panel_model$identifier), "\n")
cat("Quarters:", n_distinct(road_panel_model$quarter_year), "\n")
cat("Rows:    ", nrow(road_panel_model), "\n")

# Variable names
names(road_panel_model)

# Treatment group counts
road_panel_model %>%
  distinct(identifier, treated_group_any, treated_group_50pct,
           control_group1, control_group2, scheme) %>%
  count(treated_group_50pct, control_group1, control_group2) %>%
  print()

#Treatment timing — what does caz_start_q look like?
road_panel_model %>%
  filter(treated_group_50pct == 1) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q) %>%
  print()


class(road_panel_model$quarter_year)
head(road_panel_model$quarter_year, 5)

#  Outcome variable check
road_panel_model %>%
  summarise(
    mean_KSI   = mean(KSI_adj_All, na.rm = TRUE),
    mean_slight = mean(Slight_adj_All, na.rm = TRUE),
    mean_total  = mean(total_inj_adj_All, na.rm = TRUE),
    pct_zero    = mean(total_inj_adj_All == 0) * 100
  ) %>%
  print()

# Does the panel link back to OAs?
"OA" %in% names(road_panel_model)







# =============================================================================
# matching code  OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING
# =============================================================================
#
#
#   Construct matched comparison groups for a DiD evaluation of CAZ/LEZ
#   interventions using two-stage Mahalanobis Distance Matching.
#
# CONTROL GROUP STRATEGY (asymmetric by country):
#
#   ENGLAND — Other-city controls only (control_group2):
#     OAs from LADs that contain no treated zone. Same-city control OAs
#     (control_group1, from CAZ/LEZ LADs) are excluded to prevent indirect
#     treatment contamination from traffic displacement and CAZ spillover.
#
#   SCOTLAND — Other-city + same-city controls (control_group2 + control_group1):
#     All four major Scottish cities (Glasgow, Edinburgh, Aberdeen, Dundee)
#     implemented LEZs. Applying the same-city exclusion alone would restrict
#     Scotland to the two non-LEZ LADs (S12000014, S12000050). Retaining
#     same-city controls (outside the buffer zone) maximises pool size and
#     match quality for Scotland at the cost of potential downward bias:
#     same-city controls may have experienced indirect LEZ effects (e.g.,
#     suppressed feeder-road traffic), which would attenuate Scottish
#     estimates toward zero.
#     LIMITATION: Scottish results should be interpreted as conservative
#     lower-bound estimates and reported separately from England.
#
#   Compare 16b (other-city only for both countries) against this script to
#   assess sensitivity of Scottish estimates to same-city control inclusion.
#
# OUTPUTS:
#   OA_matched_treated_mixed.rds        — treated OA IDs + weights + stratum
#   OA_matched_donors_mixed.rds         — control OA IDs + weights
#   OA_matched_full_mixed.rds           — combined matched dataset (Eng + Scot)
#   OA_matched_full_mixed_England.rds   — England only
#   OA_matched_full_mixed_Scotland.rds  — Scotland only
#   OA_common_support_flags_mixed.rds   — structurally isolated OA flags
#   OA_ratio_selection_mixed.rds        — ratio selection diagnostics by country
#   OA_balance_tests_mixed.rds          — balance improvement test results
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(here)
library(MASS)
library(purrr)
library(sf)
library(ggrepel)
library(tidyverse)
library(patchwork)

set.seed(222)   #  governs all MatchIt nearest-neighbour calls


select <- dplyr::select
filter <- dplyr::filter

dir.create(here("output", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics")

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))


# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

stage1_road     <- c("road_density_m_km2", "road_length_km",
                     "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban    <- c("dist_citycentre", "pop_density", "area_km2")
stage1_business <- c("business_retail_per_km2")
stage1_socdem   <- c(
  "IMD",
  "cars_one_pct", "cars_twoPlus_pct",
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
  "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
  "age_45to64_pct", "age_65to84_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_business, stage1_socdem)

stage2_trends <- c(
  "trend_car_KSI_pkm",   "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",   "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",   "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm",   "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",   "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",   "mean_ped_slight_pkm",
  "mean_other_KSI_pkm", "mean_other_slight_pkm",
  "mean_total_pkm"
)
stage2_vars <- c(stage2_trends, stage2_levels)

log_transform_s1        <- c("road_length_km", "pop_density", "dist_citycentre",
                             "road_density_m_km2", "business_retail_per_km2")
log_nozero_s1           <- c("area_km2")
log_transform_s2_levels <- stage2_levels

log_names_s1        <- paste0("log1p_", log_transform_s1)
log_nozero_names_s1 <- paste0("log_",   log_nozero_s1)
log_names_s2        <- paste0("log1p_", log_transform_s2_levels)

stage1_vars_log <- c(
  log_names_s1, log_nozero_names_s1,
  setdiff(stage1_vars, c(log_transform_s1, log_nozero_s1))
)
stage2_vars_log <- c(stage2_trends, log_names_s2)

# =============================================================================
# BUILD DATASETS — ENGLAND AND SCOTLAND SEPARATELY
# =============================================================================

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                          ~ "Unknown"
    )
  )

treated_lads_england <- OA_matching_dataset %>%
  filter(treated_OA == 1, country == "England") %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

treated_lads_scotland <- OA_matching_dataset %>%
  filter(treated_OA == 1, country == "Scotland") %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

cat("English CAZ LADs (same-city controls excluded):\n")
print(treated_lads_england)
cat("\nScottish LEZ LADs (same-city controls RETAINED — see header note):\n")
print(treated_lads_scotland)

# --- ENGLAND: other-city controls only (control_group2) ----------------------
# control_group1 OAs (from treated LADs) are fully excluded to prevent
# contamination from traffic displacement and CAZ spillover effects.
# control_group1_OA == 0 guards against any OA being double-classified.

data_england <- OA_matching_dataset %>%
  filter(
    country           == "England",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA         == 0,
    n_roads            > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("\n=== ENGLAND dataset (other-city controls only) ===\n")
cat("Treated:", sum(data_england$treat_indicator == 1),
    "| Controls:", sum(data_england$treat_indicator == 0), "\n")
cat("Control-to-treated ratio:",
    round(sum(data_england$treat_indicator == 0) /
            sum(data_england$treat_indicator == 1), 1), "\n\n")

eng_samecity_leak <- data_england %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_england) %>%
  nrow()
cat("Same-city control leak (should be 0):", eng_samecity_leak, "\n\n")

# --- SCOTLAND: other-city + same-city controls --------------------------------
# Both control_group1 (same-city, outside buffer) and control_group2
# (other-city, non-LEZ LADs) are retained to maximise pool size and match
# quality. This introduces potential downward bias if same-city controls
# were indirectly affected by the LEZ; results should be read as conservative.
# Compare with 16b (other-city only) to assess sensitivity.

data_scotland <- OA_matching_dataset %>%
  filter(
    country == "Scotland",
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1)
    # same-city exclusion deliberately omitted — see header
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cat("=== SCOTLAND dataset (other-city + same-city controls) ===\n")
cat("Treated:", sum(data_scotland$treat_indicator == 1),
    "| Controls:", sum(data_scotland$treat_indicator == 0), "\n")
cat("Control-to-treated ratio:",
    round(sum(data_scotland$treat_indicator == 0) /
            sum(data_scotland$treat_indicator == 1), 1), "\n")
cat("LIMITATION: same-city controls retained — potential attenuation bias.\n\n")

scotland_samecity_controls <- data_scotland %>%
  filter(treat_indicator == 0, LAD24CD %in% treated_lads_scotland) %>%
  nrow()
scotland_othercity_controls <- data_scotland %>%
  filter(treat_indicator == 0, !(LAD24CD %in% treated_lads_scotland)) %>%
  nrow()
cat("Scottish same-city controls:", scotland_samecity_controls, "\n")
cat("Scottish other-city controls:", scotland_othercity_controls, "\n\n")

# =============================================================================
# AGE BAND AGGREGATION + WINSORISE + LOG-TRANSFORM
# =============================================================================

prep_dataset <- function(data) {
  data %>%
    mutate(
      age_under15_pct = X4under_pct  + X5to9_pct   + X10to14_pct,
      age_15to24_pct  = X15to19_pct  + X20to24_pct,
      age_25to44_pct  = X25to29_pct  + X30to34_pct + X35to39_pct + X40to44_pct,
      age_45to64_pct  = X45to49_pct  + X50to54_pct + X55to59_pct + X60to64_pct,
      age_65to84_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct
    )
}

data_england  <- prep_dataset(data_england)
data_scotland <- prep_dataset(data_scotland)

winsorise_and_log_s1 <- function(data, raw_vars, log_vars,
                                 log_nozero_vars = character(0)) {
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  for (v in intersect(log_nozero_vars, names(data))) {
    data[[paste0("log_", v)]] <- log(data[[v]])
  }
  data
}

winsorise_and_log_s2 <- function(data, raw_vars, log_vars) {
  for (v in intersect(raw_vars, names(data))) {
    q <- quantile(data[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
    data[[v]] <- pmin(pmax(data[[v]], q[1]), q[2])
  }
  for (v in intersect(log_vars, names(data))) {
    data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  data
}

# Transform each country independently so winsorisation thresholds
# reflect each country's own distribution, not the pooled distribution.
data_england_clean  <- winsorise_and_log_s1(data_england,  stage1_vars,
                                            log_transform_s1, log_nozero_s1)
data_england_clean  <- winsorise_and_log_s2(data_england_clean,  stage2_levels,
                                            log_transform_s2_levels)

data_scotland_clean <- winsorise_and_log_s1(data_scotland, stage1_vars,
                                            log_transform_s1, log_nozero_s1)
data_scotland_clean <- winsorise_and_log_s2(data_scotland_clean, stage2_levels,
                                            log_transform_s2_levels)

check_vars <- function(data, vars, label) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0)
    cat("WARNING —", label, "missing:", paste(missing, collapse = ", "), "\n")
  vcheck <- data %>%
    summarise(across(all_of(intersect(vars, names(data))),
                     ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")
  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)
  if (length(low) > 0)
    cat("Dropping near-zero variance (", label, "):",
        paste(low, collapse = ", "), "\n")
  setdiff(intersect(vars, names(data)), low)
}

s1_vars_england  <- check_vars(data_england_clean,  stage1_vars_log, "S1 England")
s2_vars_england  <- check_vars(data_england_clean,  stage2_vars_log, "S2 England")
s1_vars_scotland <- check_vars(data_scotland_clean, stage1_vars_log, "S1 Scotland")
s2_vars_scotland <- check_vars(data_scotland_clean, stage2_vars_log, "S2 Scotland")


flow_counts <- tibble(
  stage = character(),
  treated_OAs = integer(),
  control_OAs = integer()
)




# =============================================================================
# BALANCE TEST FUNCTION
# =============================================================================

balance_test_log <- list()

run_balance_tests <- function(matchit_obj, trend_vars, label) {
  bt  <- bal.tab(matchit_obj, un = TRUE,
                 stats = c("mean.diffs", "variance.ratios"))
  bal <- bt$Balance
  
  mean_un  <- mean(abs(bal$Diff.Un),  na.rm = TRUE)
  mean_adj <- mean(abs(bal$Diff.Adj), na.rm = TRUE)
  test_a   <- mean_adj < mean_un
  cat(sprintf("  [TEST a] Mean |SMD|: %.3f → %.3f  %s\n",
              mean_un, mean_adj, if (test_a) "PASS ✓" else "FAIL ✗"))
  
  trend_in_bal  <- intersect(trend_vars, rownames(bal))
  max_trend_smd <- if (length(trend_in_bal) > 0)
    max(abs(bal[trend_in_bal, "Diff.Adj"]), na.rm = TRUE) else NA_real_
  test_b <- !is.na(max_trend_smd) && max_trend_smd < 0.1
  cat(sprintf("  [TEST b] Max trend |SMD|: %.4f  %s\n",
              max_trend_smd, if (test_b) "PASS ✓" else "FAIL ✗"))
  
  vr_col  <- if ("Var.Ratio.Adj" %in% names(bal)) "Var.Ratio.Adj" else NULL
  vr_fail <- character(0)
  if (!is.null(vr_col)) {
    vr      <- bal[[vr_col]]
    vr_fail <- rownames(bal)[!is.na(vr) & (vr < 0.5 | vr > 2.0)]
    test_c  <- length(vr_fail) == 0
    cat(sprintf("  [TEST c] Variance ratio [0.5, 2.0]: %d/%d pass  %s\n",
                sum(is.na(vr) | (vr >= 0.5 & vr <= 2.0), na.rm = TRUE),
                sum(!is.na(vr)),
                if (test_c) "PASS ✓" else "FAIL ✗"))
  } else { test_c <- NA }
  
  balance_test_log[[label]] <<- tibble(
    label         = label,
    mean_smd_un   = round(mean_un,       4),
    mean_smd_adj  = round(mean_adj,      4),
    max_trend_smd = round(max_trend_smd, 4),
    test_a_pass   = test_a,
    test_b_pass   = test_b,
    test_c_pass   = if (is.na(test_c)) NA else test_c,
    vr_fail_vars  = paste(vr_fail, collapse = "; ")
  )
  invisible(balance_test_log[[label]])
}

# =============================================================================
# RATIO SELECTION FUNCTION
# =============================================================================

prepare_s2_for_ratio <- function(data_clean, s1_vars, s2_vars, label) {
  s1v        <- check_vars(data_clean, s1_vars, paste("S1 ratio prep", label))
  formula_s1 <- reformulate(s1v, response = "treat_indicator")
  m_s1 <- matchit(formula_s1, data = data_clean, method = "nearest",
                  distance = "mahalanobis", ratio = 10, replace = TRUE)
  mm_s1       <- m_s1$match.matrix
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  s2_raw      <- bind_rows(treated_s1, controls_s1)
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  for (v in intersect(stage2_levels, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <-
      log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }
  for (v in intersect(stage2_trends, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  s2v <- check_vars(s2_raw, s2_vars, paste("S2 ratio prep", label))
  list(data = s2_raw, s2_vars = s2v)
}

select_ratio <- function(data, s2_vars, trend_vars, label,
                         ratios_to_test = 1:10) {
  n_t <- sum(data$treat_indicator == 1)
  n_c <- sum(data$treat_indicator == 0)
  ratios_to_test <- ratios_to_test[ratios_to_test <= floor(n_c / n_t)]
  
  cat("\n--- Ratio selection:", label, "---\n")
  cat("  Treated:", n_t, "| Controls:", n_c, "\n")
  
  formula_s2 <- reformulate(s2_vars, response = "treat_indicator")
  
  ratio_results <- map_df(ratios_to_test, function(r) {
    m <- tryCatch(
      matchit(formula_s2, data = data, method = "nearest",
              distance = "mahalanobis", ratio = r, replace = TRUE),
      error = function(e) NULL
    )
    if (is.null(m)) return(tibble(ratio = r, mean_smd = NA,
                                  max_trend_smd = NA, max_level_smd = NA))
    bt          <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
    trend_rows  <- rownames(bt)[rownames(bt) %in% trend_vars]
    level_rows  <- rownames(bt)[!rownames(bt) %in% trend_vars]
    tibble(
      ratio         = r,
      mean_smd      = round(mean(abs(bt$Diff.Adj), na.rm = TRUE), 4),
      max_trend_smd = round(max(abs(bt[trend_rows, "Diff.Adj"]),
                                na.rm = TRUE), 4),
      max_level_smd = round(max(abs(bt[level_rows, "Diff.Adj"]),
                                na.rm = TRUE), 4)
    )
  })
  
  print(ratio_results)
  
  best <- ratio_results %>%
    filter(!is.na(max_trend_smd)) %>%
    arrange(max_trend_smd, mean_smd) %>%
    slice(1)
  
  cat(sprintf("  Selected: 1:%d (max trend |SMD| = %.4f)\n",
              best$ratio, best$max_trend_smd))
  
  list(optimal_ratio = best$ratio, ratio_results = ratio_results, label = label)
}

# =============================================================================
# RATIO SELECTION — ENGLAND AND SCOTLAND SEPARATELY
# =============================================================================

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("RATIO SELECTION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

s2_prep_eng <- prepare_s2_for_ratio(data_england_clean,  s1_vars_england,
                                    s2_vars_england,  "England")
s2_prep_sco <- prepare_s2_for_ratio(data_scotland_clean, s1_vars_scotland,
                                    s2_vars_scotland, "Scotland")

trend_vars_eng <- intersect(s2_vars_england,  stage2_trends)
trend_vars_sco <- intersect(s2_vars_scotland, stage2_trends)

ratio_eng <- select_ratio(s2_prep_eng$data, s2_prep_eng$s2_vars,
                          trend_vars_eng, "England",  1:10)
# Scotland tested up to 1:10; same-city controls expand the pool so higher
# ratios are feasible, though match quality may plateau earlier.
ratio_sco <- select_ratio(s2_prep_sco$data, s2_prep_sco$s2_vars,
                          trend_vars_sco, "Scotland", 1:10)

optimal_ratio_england  <- ratio_eng$optimal_ratio
optimal_ratio_scotland <- ratio_sco$optimal_ratio

cat("\nOptimal ratio — England:  1:", optimal_ratio_england, "\n")
cat("Optimal ratio — Scotland: 1:", optimal_ratio_scotland, "\n\n")

# Save and plot ratio selection
ratio_combined <- bind_rows(
  ratio_eng$ratio_results %>% mutate(
    country = "England (other-city controls only)"),
  ratio_sco$ratio_results %>% mutate(
    country = "Scotland (other-city + same-city controls)")
)
saveRDS(ratio_combined,
        here("data", "processed", "OA_ratio_selection_mixed.rds"))

p_ratio <- ratio_combined %>%
  filter(!is.na(max_trend_smd)) %>%
  ggplot(aes(x = ratio, y = max_trend_smd,
             colour = country, group = country)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#888888", linewidth = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             colour = "#555555", linewidth = 0.5) +
  scale_colour_manual(
    values = c(
      "England (other-city controls only)"            = "#2E6FAB",
      "Scotland (other-city + same-city controls)"    = "#6B3FA0"
    )
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title    = "Ratio Selection — Maximum Trend |SMD| by Country",
    subtitle = paste0(
      "England: optimal = 1:", optimal_ratio_england,
      " (other-city only)",
      " | Scotland: optimal = 1:", optimal_ratio_scotland,
      " (other-city + same-city)\n",
      "Scottish estimates may be attenuated if same-city controls were",
      " indirectly affected by the LEZ"
    ),
    x       = "Matching ratio (1:k)",
    y       = "Maximum trend |SMD| after matching",
    colour  = NULL,
    caption = paste0(
      "England: control_group2 only (same-city exclusion applied).\n",
      "Scotland: control_group1 + control_group2 (same-city exclusion not applied).\n",
      "Compare with 16b (other-city only for both) to assess sensitivity."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1A2E5A"),
    plot.subtitle   = element_text(colour = "#555555", size = 10),
    plot.caption    = element_text(colour = "#888888", hjust = 0, size = 9),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "fig08_ratio_selection_by_country.png"),
       p_ratio, width = 13, height = 7, dpi = 300, bg = "white")
cat("Saved: fig08_ratio_selection_by_country.png\n\n")

# =============================================================================
# MATCHING FUNCTION
# =============================================================================

run_matching <- function(data_clean, s1_vars, s2_vars, ratio,
                         label, trend_vars) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("MATCHING —", label, "| ratio 1:", ratio, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # ---- STAGE 1 ---------------------------------------------------------------
  cat("--- Stage 1:", label, "---\n")
  s1v        <- check_vars(data_clean, s1_vars, paste("S1", label))
  formula_s1 <- reformulate(s1v, response = "treat_indicator")
  
  m_s1 <- tryCatch(
    matchit(formula_s1, data = data_clean, method = "nearest",
            distance = "mahalanobis", ratio = 10, replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s1)) return(NULL)
  
  mm_s1       <- m_s1$match.matrix
  treated_s1  <- data_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
  controls_s1 <- data_clean[ctrl_idx, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated:", nrow(treated_s1),
      "| Unique controls in pool:", nrow(controls_s1), "\n")
  
  # Common support
  pool_idx <- c(as.integer(rownames(mm_s1)), ctrl_idx)
  S_s1 <- cov(data_clean[pool_idx, s1v], use = "pairwise.complete.obs")
  dist_s1 <- map_df(seq_len(nrow(mm_s1)), function(i) {
    t_idx     <- as.integer(rownames(mm_s1)[i])
    trow      <- data_clean[t_idx, , drop = FALSE]
    c_indices <- mm_s1[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    map_df(seq_along(c_indices), function(j) {
      crow <- data_clean[as.integer(c_indices[j]), , drop = FALSE]
      tibble(treated_OA = trow[["OA"]], control_OA = crow[["OA"]],
             mdist = mahalanobis(as.numeric(crow[s1v]),
                                 as.numeric(trow[s1v]), S_s1))
    })
  })
  
  # Flag isolated OAs by their minimum Stage 1 Mahalanobis distance to any
  # control. No arbitrary percentile threshold: every treated OA gets its
  # min_dist_s1 recorded and structurally_isolated = TRUE if it exceeds the
  # median absolute deviation (MAD) rule  (> median + 3*MAD), a
  # distribution-free outlier criterion that adapts to each country's pool.
  min_dist <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  
  mad_threshold  <- median(min_dist$min_dist_s1) +
    3 * mad(min_dist$min_dist_s1, constant = 1)
  
  isolated_OAs <- min_dist %>%
    mutate(
      structurally_isolated = min_dist_s1 > mad_threshold,
      flag_threshold        = round(mad_threshold, 4)
    )
  
  cat(sprintf("  Isolated OAs (median + 3*MAD threshold = %.2f): %d / %d\n",
              mad_threshold,
              sum(isolated_OAs$structurally_isolated),
              nrow(isolated_OAs)))
  
  
  min_dist     <- dist_s1 %>% group_by(treated_OA) %>%
    summarise(min_dist_s1 = min(mdist), .groups = "drop")
  threshold_95 <- quantile(min_dist$min_dist_s1, 0.95)
  isolated_OAs <- min_dist %>% filter(min_dist_s1 > threshold_95) %>%
    mutate(structurally_isolated = TRUE)
  
  cat("  Isolated OAs (95th pct):", nrow(isolated_OAs), "/",
      nrow(treated_s1), "\n")
  
  cat("\n  Stage 1 balance:\n")
  run_balance_tests(m_s1, trend_vars = character(0),
                    label = paste0("S1_", label))
  
  
  
  ######
  
  ####################
  
  
  # ---- STAGE 2 ---------------------------------------------------------------
  cat("\n--- Stage 2:", label, "---\n")
  s2_raw      <- bind_rows(treated_s1, controls_s1) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  for (v in intersect(stage2_trends, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
  }
  for (v in intersect(stage2_levels, names(s2_raw))) {
    q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
    q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
    s2_raw[[paste0("log1p_", v)]] <-
      log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
  }
  
  s2v        <- check_vars(s2_raw, s2_vars, paste("S2", label))
  formula_s2 <- reformulate(s2v, response = "treat_indicator")
  
  m_s2 <- tryCatch(
    matchit(formula_s2, data = s2_raw, method = "nearest",
            distance = "mahalanobis", ratio = ratio, replace = TRUE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m_s2)) return(NULL)
  
  matched_data <- match.data(m_s2)
  
  cat("  Treated:", sum(matched_data$treat_indicator == 1),
      "| Controls:", sum(matched_data$treat_indicator == 0), "\n")
  
  # Distances
  mm_s2 <- m_s2$match.matrix
  S_s2  <- cov(s2_raw[as.integer(rownames(mm_s2)), s2v],
               use = "pairwise.complete.obs")
  dist_s2 <- map_df(seq_len(nrow(mm_s2)), function(i) {
    t_idx     <- as.integer(rownames(mm_s2)[i])
    trow      <- s2_raw[t_idx, , drop = FALSE]
    c_indices <- mm_s2[i, ]; c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- s2_raw[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(as.numeric(crow[s2v]), as.numeric(trow[s2v]), S_s2)
    })
    tibble(OA = trow[["OA"]], mdist = mean(dists))
  })
  matched_data <- matched_data %>% left_join(dist_s2, by = "OA")
  
  cat("\n  Stage 2 balance:\n")
  run_balance_tests(m_s2, trend_vars = trend_vars,
                    label = paste0("S2_", label))
  
  # Extract treated→control OA pairs (used by scripts 17 and 18)
  mm_pairs    <- m_s2$match.matrix
  treated_oas <- s2_raw$OA[as.integer(rownames(mm_pairs))]
  pairs <- map_df(seq_len(nrow(mm_pairs)), function(i) {
    t_oa  <- treated_oas[i]
    c_idx <- as.integer(mm_pairs[i, ])
    c_idx <- c_idx[!is.na(c_idx)]
    if (length(c_idx) == 0) return(tibble())
    tibble(treated_OA = t_oa, control_OA = s2_raw$OA[c_idx])
  })
  
  list(matched_data = matched_data, isolated_OAs = isolated_OAs,
       matchit_s2 = m_s2, dist_s2 = dist_s2, country = label,
       pairs = pairs)
}

# =============================================================================
# RUN MATCHING
# =============================================================================

result_england  <- run_matching(
  data_clean  = data_england_clean,
  s1_vars     = s1_vars_england,
  s2_vars     = s2_vars_england,
  ratio       = optimal_ratio_england,
  label       = "England",
  trend_vars  = trend_vars_eng
)

result_scotland <- run_matching(
  data_clean  = data_scotland_clean,
  s1_vars     = s1_vars_scotland,
  s2_vars     = s2_vars_scotland,
  ratio       = optimal_ratio_scotland,
  label       = "Scotland",
  trend_vars  = trend_vars_sco
)

# =============================================================================
# COMBINE
# =============================================================================

matched_england  <- result_england$matched_data  %>% mutate(country = "England")
matched_scotland <- result_scotland$matched_data %>% mutate(country = "Scotland")
matched_full_mixed   <- bind_rows(matched_england, matched_scotland)

isolated_combined <- bind_rows(
  result_england$isolated_OAs  %>% mutate(country = "England"),
  result_scotland$isolated_OAs %>% mutate(
    country = "Scotland",
    note    = "Same-city controls included — potential downward bias"
  )
)

cat("\n=== FINAL COMBINED DATASET ===\n")
cat("England  — Treated:", sum(matched_england$treat_indicator  == 1),
    "| Controls:", sum(matched_england$treat_indicator  == 0),
    "(other-city only)\n")
cat("Scotland — Treated:", sum(matched_scotland$treat_indicator == 1),
    "| Controls:", sum(matched_scotland$treat_indicator == 0),
    "(other-city + same-city)\n")
cat("Total    — Treated:", sum(matched_full_mixed$treat_indicator   == 1),
    "| Controls:", sum(matched_full_mixed$treat_indicator   == 0), "\n\n")
cat("NOTE: Scottish controls include same-city OAs (outside buffer).\n")
cat("Compare Scotland results against 16b (other-city only) for sensitivity.\n")

# =============================================================================
# BASELINE INJURY STRATIFICATION
# =============================================================================

q_breaks <- quantile(
  matched_full_mixed %>%
    filter(treat_indicator == 1) %>%
    pull(log1p_mean_total_pkm),
  probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE
)

matched_full_mixed <- matched_full_mixed %>%
  mutate(baseline_injury_stratum = case_when(
    treat_indicator == 0                ~ NA_integer_,
    log1p_mean_total_pkm <= q_breaks[2] ~ 1L,
    log1p_mean_total_pkm <= q_breaks[3] ~ 2L,
    log1p_mean_total_pkm <= q_breaks[4] ~ 3L,
    TRUE                                ~ 4L
  ))

# =============================================================================
# EXTRACT + INTEGRITY CHECKS + SAVE
# =============================================================================

matched_mixed_treated  <- matched_full_mixed %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights, baseline_injury_stratum, country)

matched_mixed_controls <- matched_full_mixed %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights, country)

stopifnot(
  "treated weights == 1"     = all(matched_mixed_treated$weights == 1),
  "no NA control weights"    = !anyNA(matched_mixed_controls$weights),
  "no duplicate treated OAs" = !anyDuplicated(matched_mixed_treated$OA)
)
cat("\nAll integrity checks passed.\n")

saveRDS(matched_mixed_treated,
        here("data", "processed", "OA_matched_treated_mixed.rds"))
saveRDS(matched_mixed_controls,
        here("data", "processed", "OA_matched_donors_mixed.rds"))
saveRDS(matched_full_mixed,
        here("data", "processed", "OA_matched_full_mixed.rds"))
saveRDS(matched_england,
        here("data", "processed", "OA_matched_full_mixed_England.rds"))
saveRDS(matched_scotland,
        here("data", "processed", "OA_matched_full_mixed_Scotland.rds"))
saveRDS(isolated_combined,
        here("data", "processed", "OA_common_support_flags_mixed.rds"))
saveRDS(bind_rows(balance_test_log),
        here("data", "processed", "OA_balance_tests_mixed.rds"))

# --- Overall balance summary (Stage 1 & Stage 2) ----------------------------
balance_summary <- bind_rows(balance_test_log) |>
  mutate(stage = if_else(str_starts(label, "S1"), "Stage 1", "Stage 2")) |>
  group_by(stage) |>
  summarise(
    mean_abs_smd_unmatched = mean(mean_smd_un,  na.rm = TRUE),
    mean_abs_smd_matched   = mean(mean_smd_adj, na.rm = TRUE),
    max_trend_smd          = max(max_trend_smd, na.rm = TRUE),
    all_test_a_pass        = all(test_a_pass,   na.rm = TRUE),
    all_test_b_pass        = all(test_b_pass,   na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== OVERALL BALANCE SUMMARY ===\n")
print(balance_summary)

pairs_mixed <- bind_rows(
  result_england$pairs  %>% mutate(country = "England"),
  result_scotland$pairs %>% mutate(country = "Scotland")
)
saveRDS(pairs_mixed,
        here("data", "processed", "OA_matching_pairs_mixed.rds"))

cat("\n=== OUTPUTS SAVED ===\n")
cat("  OA_matched_full_mixed.rds          — combined England + Scotland\n")
cat("  OA_matched_full_mixed_England.rds  — England (other-city controls)\n")
cat("  OA_matched_full_mixed_Scotland.rds — Scotland (other-city + same-city)\n")
cat("  OA_matched_treated_mixed.rds\n")
cat("  OA_matched_donors_mixed.rds\n")
cat("  OA_common_support_flags_mixed.rds\n")
cat("  OA_ratio_selection_mixed.rds\n")
cat("  OA_balance_tests_mixed.rds\n")
cat("  OA_matching_pairs_mixed.rds\n")
cat("  fig08_ratio_selection_by_country.png\n")







