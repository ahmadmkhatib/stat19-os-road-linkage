# ==========================================================
# STAT19 – OS Open Roads Linkage Framework
# Script: 01_prepare_os_roads.R
# Purpose: Prepare OS Open Roads data within selected LADs
# ==========================================================

library(sf)
library(tidyverse)
library(stringr)
library(here)

# ----------------------------------------------------------
# USER INPUT: Set Data Paths
# ----------------------------------------------------------

# Users must download:
# - OS Open Roads
# - LAD boundaries (ONS)

roads_path <- here("data", "raw", "OS highways all.shp")
lads_path  <- here("data", "raw", "LAD_DEC_24_UK_BGC.shp")

# ----------------------------------------------------------
# Load LAD Boundaries
# ----------------------------------------------------------

LADs <- st_read(lads_path, quiet = TRUE)

# ----------------------------------------------------------
# Define LAD Selection (modify this section as needed)
# ----------------------------------------------------------

selected_lads <- c(
  # England CAZ
  "E06000022","E08000025","E08000032","E06000023",
  "E06000044","E08000019","E08000021","E08000018",
  # Controls
  "E06000043","E08000026","E06000015","E06000010",
  "E08000035","E06000016","E08000012","E06000032",
  "E08000003","E06000062","E07000148","E06000018",
  "E06000026","E06000038","E06000045","E06000021",
  "E08000031",
  # Scotland
  "S12000049","S12000036","S12000033",
  "S12000042","S12000038","S12000029",
  # Wales
  "W06000015","W06000022","W06000011"
)

LADs_sub <- LADs %>%
  filter(LAD24CD %in% selected_lads)

lads_union <- st_union(LADs_sub)

# ----------------------------------------------------------
# Load OS Open Roads
# ----------------------------------------------------------

roads <- st_read(roads_path, quiet = TRUE)

# for speed 
dir.create(here("data", "processed"), showWarnings = FALSE)
roads <- roads %>% st_zm(roads, drop = TRUE, what = "ZM") %>% select(-fid)
st_write(roads, "data/processed/roads.gpkg", delete_dsn = TRUE)
## dealing with .gpkg is faster than .shp 
rm(roads)
gc()
## n = 824790

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

dir.create(here("data", "processed"), showWarnings = FALSE)

saveRDS(roads, output_path)

