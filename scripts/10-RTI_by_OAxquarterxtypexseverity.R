# ============================================================
# OA-Level Injury Dataset Construction
# ============================================================
#
# Purpose
# ------------------------------------------------------------
#  aggregates road traffic injuries from the
# matched STATS19–road dataset to the Output Area (OA) level.
#
# Injuries are aggregated by:
#   • OA
#   • Quarter  
#   • Casualty type
#   • Injury severity
#
# The script produces an OA × quarter panel dataset of injury
# counts that is later used to construct:
#
#   • Pre-treatment injury baselines
#   • Pre-treatment injury trends
#   • Matching variables for synthetic control construction
#
#  Steps
# ------------------------------------------------------------
# 
# Aggregate injuries by OA × quarter × casualty type × severity.
# Pivot the dataset to wide format 
# Get total injuries per OA × quarter across all
# casualty types and severities.
#
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



injuries <- readRDS( here("data", "processed","injuries_matched_final.rds"))


table(injuries$casualty_type1)

## remove the / from the Car/Van 
injuries <- injuries %>%
  mutate(
    casualty_type1 = str_replace_all(casualty_type1, "/", "_")
  )


table(injuries$casualty_severity1)

# ------------------------------------------------------------
# Aggregate injuries by OA × quarter × type × severity
# ------------------------------------------------------------
OA_injuries_q <- injuries %>%
  group_by(OA, quarter_year, casualty_type1, casualty_severity1) %>%
  summarise(n_injuries = n(), .groups = "drop")

glimpse(OA_injuries_q)

# ------------------------------------------------------------
# Pivot to wide format: type + severity columns
# ------------------------------------------------------------
OA_injuries_wide <- OA_injuries_q %>%
  unite(type_severity, casualty_type1, casualty_severity1) %>%
  pivot_wider(
    names_from = type_severity,
    values_from = n_injuries,
    values_fill = 0
  )

# ------------------------------------------------------------
# Compute total injuries per OA × quarter
# ------------------------------------------------------------
# List of columns to sum for total injuries
inj_cols <- c(
  "Car_Van_KSI", "Car_Van_Slight",
  "Cyclist_KSI", "Cyclist_Slight",
  "Pedestrian_KSI", "Pedestrian_Slight",
  "Other_KSI", "Other_Slight"
)

# Compute total injuries safely
OA_injuries_wide <- OA_injuries_wide %>%
  mutate(
    total_injuries = rowSums(across(any_of(inj_cols)), na.rm = TRUE)
  )




glimpse(OA_injuries_wide)
# ------------------------------------------------------------
# OA × quarter injury dataset
# ------------------------------------------------------------
saveRDS(
  OA_injuries_wide,
  here("data","processed","OA_injuries_quarterly.rds"))