
# ============================================================
# Panel Construction
# a road × quarter panel dataset

# The treatment is defined at the road level:
#   A road is treated if it lies inside a CAZ boundary
#   AND the quarter is after the CAZ implementation date.
#
# The key identifying variation comes from comparing:
##   • Roads inside CAZ boundaries (treated)
#   • Roads outside CAZ boundaries (controls)
#
# The counterfactual question is:
#   What would have happened to treated roads in the absence of CAZ?
#
#  three control groups.
#
# ------------------------------------------------------------
# CONTROL GROUP DEFINITIONS
# ------------------------------------------------------------
#
# Control Group 1: Same-city controls # -----------------------------------------------------
# Roads located in the same CAZ cities but outside CAZ boundaries.
#
##   - Same macroeconomic conditions
#   - Same contex  and infrastructure

# Control Group 2: City-centre roads in non-CAZ cities
# -----------------------------------------------------
# proxy "city centre" 
#
#   1. Compute the centroid of each non-CAZ city
#   2. Create a 1 km buffer around the centroid
#   3. Classify roads within this buffer as pseudo-city-centre roads
#
#  - Ensures comparable traffic density and urban structure
#   - Avoids comparing treated urban cores to suburban/rural roads
#
#
# Control Group 3: Mixed controls
# --------------------------------
# Combination of Control Group 1 and Control Group 2.
#
## ------------------------------------------------------------
# # Steps:
#   1. Classify roads as treated / control at road level
#   2. Keep only roads relevant for analysis
#   3. Expand selected roads to quarter panel
#   4. Aggregate injuries to road × quarter level
#   5. Merge treatment timing and construct post indicator
#   6. Export as Parquet dataset using Arrow
#
# ------------------------------------------------------------
# OUTCOME VARIABLES
# ------------------------------------------------------------
#
# Injuries  aggregated to road × quarter level.
#
# Outcomes separate for:
#   - KSI_adj     (Killed or Seriously Injured, adjusted)
#   - Slight_adj  (Slight injuries, adjusted)
#
# This allows estimating heterogeneous effects by severity.
#
#
# ------------------------------------------------------------
# FINAL DATASET STRUCTURE
# ------------------------------------------------------------
#
#   identifier              : Road ID
#   quarter_year            : Time period
#   treated                 : Post-treatment indicator
#   treated_group           : Ever-treated road indicator
#   control_group1          : Same-city outer roads
#   control_group2          : Non-CAZ city-centre roads
#   control_group3_mixed    : Combined controls
#   KSI_adj                 : Severe injuries
#   Slight_adj              : Slight injuries
#
#
# 
# ============================================================


library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(zoo)
library(arrow)
library(units)

# --data
# RTI matched data
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))
# all roads with thier attribute  - @ road level 
road_attributes<- readRDS(here("data", "processed", "road_attributes.rds"))
# roads inside the CAZs  == the treatment @ road level 
road_caz_prop<- readRDS(here("data", "processed", "roads_caz_props.rds"))
#create treatment indicator 
road_caz_prop<- road_caz_prop %>%
  mutate(
    ever_treated = 1
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
  left_join(road_caz_prop %>% select(identifier, scheme, ever_treated),
            by = "identifier") %>%
  mutate(
    ever_treated = replace_na(ever_treated, 0),
    scheme = scheme, # assign scheme for non-CAZ roads
    control_group1 = if_else(scheme %in% treated_cities & ever_treated == 0, 1, 0)
  )

# -----------------------------
#  Control 2: City centre roads in non-CAZ cities
# -----------------------------
# Non-CAZ roads
non_caz_roads <- road_classification %>%
  filter(!scheme %in% treated_cities) %>%
  st_as_sf()

# Compute city centroids
city_centroids <- non_caz_roads %>%
  group_by(scheme) %>%
  summarise(geometry = st_union(geom)) %>%
  st_centroid()

# Create 3km buffer around centroid
city_buffers <- st_buffer(city_centroids, dist = 1000)

# Roads within buffer → pseudo city centres
roads_in_centres <- st_join(non_caz_roads, city_buffers, join = st_within, left = FALSE) %>%
  pull(identifier)

# Add indicator
road_classification <- road_classification %>%
  mutate(
    control_group2 = if_else(identifier %in% roads_in_centres, 1, 0)
  )

# -----------------------------
# Control 3: Mixed (1 + 2)
# -----------------------------
road_classification <- road_classification %>%
  mutate(
    treated_group = ever_treated,
    control_group3_mixed = if_else(control_group1 == 1 | control_group2 == 1, 1, 0)
  )

# -----------------------------
# Restrict to Relevant Roads for Analysis
# -----------------------------
analysis_roads <- road_classification %>%
  filter(treated_group == 1 | control_group1 == 1 | control_group2 == 1) %>%
  select(identifier, treated_group, control_group1, control_group2, control_group3_mixed)

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
injuries <- injuries %>%
  st_drop_geometry() 


roadlevel_long <- injuries  %>%
  group_by(identifier, quarter_year, casualty_type1) %>%
  summarise(
    KSI_adj      = sum(KSI_adj, na.rm = TRUE),
    Slight_adj   = sum(Slight_adj, na.rm = TRUE),
    KSI_unadj    = sum(KSI_unadj, na.rm = TRUE),
    Slight_unadj = sum(Slight_unadj, na.rm = TRUE),
    total_inj_adj   = sum(KSI_adj + Slight_adj, na.rm = TRUE),
    total_inj_unadj = sum(KSI_unadj + Slight_unadj, na.rm = TRUE),
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
  left_join(injury_wide, by = c("identifier","quarter_year")) %>%
  mutate(across(starts_with(c("KSI","Slight","total_inj")), ~replace_na(.,0))) %>%
  left_join(analysis_roads, by = "identifier") %>%
  left_join(road_caz_prop %>% select(identifier, caz_start_q), by = "identifier") %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q),
    treated = if_else(treated_group == 1 & quarter_year >= caz_start_q, 1, 0)
  )

# -----------------------------
# Arrow Dataset
# -----------------------------
road_panel_model <- road_panel_complete %>%
  st_drop_geometry() %>%
  select(
    identifier,
    quarter_year,
    treated,
    treated_group,
    control_group1,
    control_group2,
    control_group3_mixed,
    starts_with("KSI"),
    starts_with("Slight"),
    starts_with("total_inj")
  )

write_dataset(
  road_panel_model,
  path = here("data","processed","road_panel_dataset"),
  format = "parquet"
)



table(road_panel_model$`KSI_adj_Car/Van`)
sum(road_panel_model$`KSI_adj_Car/Van`)

sum (injuries$KSI_adj)


# How many unique roads are in injuries vs panel?
length(unique(injuries$identifier))
length(unique(road_panel_model$identifier))

# How many injury rows are missing after join?
sum(!injuries$identifier %in% road_panel_model$identifier)



sum(road_panel_model$`KSI_adj_Car/Van`)
sum(road_panel_model$KSI_adj_Pedestrian)
sum(road_panel_model$KSI_adj_Cyclist)
sum(road_panel_model$KSI_adj)
 