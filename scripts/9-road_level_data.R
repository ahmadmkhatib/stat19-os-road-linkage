library(tidyverse)
library(lubridate)
library(here)
library(sf)

# ------------------------------------------------------------
# Load matched RTI data
# -------------------------------
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))

injuries <- injuries %>%
  st_drop_geometry()  %>%
  rename(identifier = matched_roadID)

n_distinct(injuries$identifier)

table(injuries$casualty_type1)
injuries <- injuries %>%  
  mutate(casualty_type1=if_else(casualty_type1 =="Car or van driver or occupant","Car/Van", casualty_type1))




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
  
# Pivot to wide
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
      select(starts_with("KSI_adj")) %>%
      as.matrix() >= 0
  )
)




------------------------------------------------------------
  # Prepare panel quarter variable
  # ------------------------------------------------------------

road_panel_complete <- road_panel_complete %>%
  mutate(
    quarter_yq = as.yearqtr(quarter_year, format = "%Y Q%q")
  )

# ------------------------------------------------------------
# 8. Join CAZ exposure to panel
# ------------------------------------------------------------

road_panel_treated <- road_panel_complete %>%
  left_join(road_caz_first, by = "identifier") %>%
  mutate(
    treated = if_else(
      !is.na(caz_start_q) & quarter_yq >= caz_start_q,
      1, 0
    ),
    ever_treated = if_else(!is.na(caz_start_q), 1, 0)
  )

# ------------------------------------------------------------
# Create optional DiD event-time variable
# ------------------------------------------------------------

road_panel_treated <- road_panel_treated %>%
  mutate(
    event_time = if_else(
      ever_treated == 1,
      as.numeric(quarter_yq - caz_start_q),
      NA_real_
    )
  )

# ------------------------------------------------------------
# 10. Final checks
# ------------------------------------------------------------

summary(road_panel_treated$treated)
table(road_panel_treated$ever_treated)

# ------------------------------------------------------------
# 11. Save output
# ------------------------------------------------------------

saveRDS(
  road_panel_treated,
  here("data", "processed", "road_panel_with_CAZ_treatment.rds")
  
  
  
  
  

 saveRDS(road_panel_complete %>% 
           st_drop_geometry(), 
         here("data", "processed", "road_panel_complete.rds"))
 
saveRDS(road_attributes, here("data", "processed", "road_attributes.rds"))
 