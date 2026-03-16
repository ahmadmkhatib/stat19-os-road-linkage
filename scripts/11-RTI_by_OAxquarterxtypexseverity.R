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

# ------------------------------------------------------------
# Load matched injuries dataset
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Save OA × quarter injury dataset
# ------------------------------------------------------------
saveRDS(
  OA_injuries_wide,
  here("data","processed","OA_injuries_quarterly.rds")
)