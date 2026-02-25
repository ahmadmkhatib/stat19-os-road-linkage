


library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(zoo)

# ------------------------------------------------------------
# RTI matched data
injuries <- read_rds(here("data", "processed", "injuries_matched_final.rds"))
# all roads with thier attribute  - roads as rows 
road_attributes<- readRDS(here("data", "processed", "road_attributes.rds"))
# roads inside the CAZs  == the treatment @ road level 
road_caz_prop<- readRDS(here("data", "processed", "roads_caz_props.rds"))
#create treatment indicator 
road_caz_prop<- road_caz_prop %>%
  mutate(
    ever_treated = 1
  )
    
## recode the injuries for the joins     
injuries <- injuries %>%
  st_drop_geometry()  %>%
  rename(identifier = matched_roadID)

n_distinct(injuries$identifier)

table(injuries$casualty_type1)
injuries <- injuries %>%  
  mutate(casualty_type1=if_else(casualty_type1 =="Car or van driver or occupant","Car/Van", casualty_type1))

# Aggregate injuries to road × quarter x type
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
  
# Pivot - road * quarter to wide
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

# ------------------------------------------------------------
#  Join CAZ ttt to the panel
# ------------------------------------------------------------
### 
road_panel_complete <- road_panel_complete %>%
  left_join(road_caz_prop, by = "identifier") 


road_panel_complete <- road_panel_complete %>%
  mutate(
    ever_treated = replace_na(ever_treated, 0),
    treated = if_else(
      ever_treated == 1 & quarter_yq >= caz_start_q,
      1, 0
    )
  )



sum(road_panel_complete$treated)
sum(road_panel_complete$ever_treated)



saveRDS(road_panel_complete %>% 
          st_drop_geometry(), 
        here("data", "processed", "road_panel_complete.rds"))


# road_panel_complete <- readRDS(here("data", "processed", "road_panel_complete.rds"))

 