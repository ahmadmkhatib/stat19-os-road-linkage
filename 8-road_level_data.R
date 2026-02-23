library(tidyverse)
library(lubridate)
library(here)

# ------------------------------------------------------------
# Load matched RTI data
# -------------------------------
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))
injuries <- injuries %>%
  st_drop_geometry()  %>%
  rename(identifier = matched_roadID)

n_distinct(injuries$matched_roadID)

table(injuries$casualty_type1)

roads_filtered<- read_rds(here("data", "processed", "roads_filtered.rds"))
names(roads_filtered)

n_distinct(roads_filtered$identifier)
table(roads_filtered$road_class)


road_attributes <- roads_filtered %>%
  select(identifier, road_class, geom) %>%
  distinct(identifier, .keep_all = TRUE)


# Aggregate at road × quarter x type
#------------------------------------------------------------
#-----

roadlevel_long <- injuries %>%
    group_by(identifier, quarter_year, casualty_type1) %>%
    summarise(
      KSI_adj     = sum(KSI_adj, na.rm = TRUE),
      Slight_adj  = sum(Slight_adj, na.rm = TRUE),
      KSI_unadj   = sum(KSI_unadj, na.rm = TRUE),
      Slight_unadj = sum(Slight_unadj, na.rm = TRUE),
      total_inj_adj = sum(KSI_adj + Slight_adj, na.rm = TRUE),
      total_inj_unadj = sum(KSI_unadj + Slight_unadj, na.rm = TRUE),
      .groups = "drop"
    )
  
# Pivot to wide (mode-specific columns)
# ------------------------------------------------------------
injury_wide <- roadlevel_long %>%
  pivot_wider(
    names_from = casualty_type1,
    values_from = c(KSI_adj, Slight_adj,
                    KSI_unadj, Slight_unadj,
                    total_inj_adj, total_inj_unadj),
    values_fill = 0
  )

  # Create balanced road × quarter panel ---all roads subset
# ------------------------------------------------------------
all_roads <- unique(road_attributes$identifier)
all_quarters <- unique(injuries$quarter_year)

road_panel <- expand_grid(
  identifier = all_roads,
  quarter_year = all_quarters
)

road_panel_complete <- road_panel %>%
  left_join(injury_wide, by = c("identifier", "quarter_year")) %>%
  mutate(
    across(starts_with(c("KSI", "Slight", "total_inj")), 
           ~replace_na(., 0))
  )

# ------------------------------------------------------------

# Attach road attributes 

road_panel_complete <- road_panel_complete %>%
  left_join(road_attributes, by = "identifier")

# ------------------------------------------------------------
#  checks

stopifnot(
  all(
    road_panel_complete %>%
      st_drop_geometry() %>%
      select(starts_with("KSI_adj")) %>%
      as.matrix() >= 0
  )
)

 saveRDS(road_panel_complete, here("data", "processed", "road_panel_complete.rds"))
 