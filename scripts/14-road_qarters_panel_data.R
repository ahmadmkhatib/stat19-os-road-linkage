
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
    ever_treated_50pct = if_else(prop_inside >= 0.5, 1, 0)
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
  left_join(road_caz_prop %>% select(identifier, scheme, ever_treated_any, ever_treated_50pct), by="identifier") 

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


road_panel_model <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

#1#1#1#1#1#1#1#1#1#1#1#1#11#
# These three outputs are what I need to see tomorrow
class(road_panel_model$quarter_year)
head(road_panel_model$quarter_year, 5)
road_panel_model %>%
  filter(treated_group_50pct == 1) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q)
# -----------------------------
# Arrow Dataset
# -----------------------------


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

# 1. Basic dimensions
cat("Roads:   ", n_distinct(road_panel_model$identifier), "\n")
cat("Quarters:", n_distinct(road_panel_model$quarter_year), "\n")
cat("Rows:    ", nrow(road_panel_model), "\n")

# 2. Variable names
names(road_panel_model)

# 3. Treatment group counts
road_panel_model %>%
  distinct(identifier, treated_group_any, treated_group_50pct,
           control_group1, control_group2, scheme) %>%
  count(treated_group_50pct, control_group1, control_group2) %>%
  print()

# 4. Treatment timing — what does caz_start_q look like?
road_panel_model %>%
  filter(treated_group_50pct == 1) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q) %>%
  print()

# 5. How is quarter_year formatted?
class(road_panel_model$quarter_year)
head(road_panel_model$quarter_year, 5)

# 6. Outcome variable check
road_panel_model %>%
  summarise(
    mean_KSI   = mean(KSI_adj_All, na.rm = TRUE),
    mean_slight = mean(Slight_adj_All, na.rm = TRUE),
    mean_total  = mean(total_inj_adj_All, na.rm = TRUE),
    pct_zero    = mean(total_inj_adj_All == 0) * 100
  ) %>%
  print()

# 7. Does the panel link back to OAs?
"OA" %in% names(road_panel_model)







 