# ==========================================================
# STAT19 – OS Open Roads Linkage Framework
# Script: 01_prepare_os_roads.R
# Purpose: Prepare OS Open Roads data within selected LADs
#          defined in Script 00 (large city subset)
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
cities_path  <- here("../stat19-os-road-linkage-data", "big_cities_with_LADs.rds")

# ----------------------------------------------------------
# Load LAD Boundaries
# ----------------------------------------------------------
## LADs with CAZ 

CAZ_LADs <- c(
  # CAZ
  "E06000022","E08000025","E08000032","E06000023",
  "E06000044","E08000019","E08000021","E08000018", 
  "S12000049","S12000036","S12000033", "S12000042")

LADs <- st_read(lads_path, quiet = TRUE)

View(caz_LADS<-LADs %>% filter(LAD24CD %in% CAZ_LADs))

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
  
  ### save the LADs sub 

saveRDS(LADs_sub, here("data", "processed", "LADs_sub.rds"))

lads_union <- st_union(LADs_sub)

# ----------------------------------------------------------
# Load OS Open Roads
# ----------------------------------------------------------

roads <- st_read(roads_path, quiet = TRUE)

# for speed ## dealing with .gpkg is faster than .shp 

roads <- roads %>% st_zm(roads, drop = TRUE, what = "ZM") %>% select(-fid)
st_write(roads, "data/processed/roads.gpkg", delete_dsn = TRUE)
rm(roads)

roads <- st_read(
  "data/processed/roads.gpkg",
  wkt_filter = st_as_text(lads_union),
  quiet = TRUE
)
# ----------------------------------------------------------
# Recode Road Classification
# ----------------------------------------------------------

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

# ----------------------------------------------------------
# Save 
# ----------------------------------------------------------
output_path <- here("data", "processed", "roads_filtered.rds")



saveRDS(roads, output_path)

