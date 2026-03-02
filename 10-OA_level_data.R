#Construct a cross-sectional Output Area (OA) dataset


library(tidyverse)
library(sf)
library(here)


# Road attributes with OA linkage       ###       (from script 8)
road_attributes <- readRDS(
  here("data", "processed", "road_attributes.rds")
)

# Road-level treatment classification       ###   (from script 9)
road_classification <- readRDS(
  here("data", "processed", "road_classification.rds")
)

names(road_classification)
# ------------------------------
# Attach road treatment classification to OA
# -------------------------------------------------

road_OA_classification <- road_attributes %>% 
  st_drop_geometry() %>%
  select(identifier, OA_CODE) %>%
  left_join(
    road_classification %>%
      st_drop_geometry() %>%
      select(identifier,
             treated_group,
             control_group1,
             control_group2,
             control_group3_mixed,
             scheme),
    by = "identifier"
  )  %>%
  mutate(
    across(
      c(treated_group,
        control_group1,
        control_group2,
        control_group3_mixed),
      ~replace_na(., 0)
    )
  )

# -------------------------------------------------
#  Collapse to OA level
# -------------------------------------------------

OA_level <- road_OA_classification %>%
  group_by(OA_CODE) %>%
  summarise(
    treated_OA  = max(treated_group, na.rm = TRUE),
    control_group1_OA = max(control_group1, na.rm = TRUE),
    control_group2_OA = max(control_group2, na.rm = TRUE),
    control_group3_mixed_OA = max(control_group3_mixed, na.rm = TRUE),
    LAD24CD = first(LAD24CD),
    scheme  = first(scheme),
    .groups = "drop"
  )


#
####  Define Spillover Buffer OAs (Same LAD as Treated)
# ------------------------------------------

treated_LADs <- OA_level %>%
  filter(treated_OA == 1) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

OA_level <- OA_level %>%
  mutate(
    buffer_OA = if_else(
      LAD24CD %in% treated_LADs & treated_OA == 0,
      1,
      0
    )
  )

# Restrict to Relevant OAs
# ---------------------------

OA_analysis <- OA_level %>%
  filter(
    treated_OA == 1 |
      buffer_OA == 1 |
      control_group2_OA == 1
  )



write_rds(
  OA_analysis,
  here("data", "processed", "OA_cross_sectios.rds")
)


