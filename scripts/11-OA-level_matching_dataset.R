



# ============================================================
# Create OA-Level Matching Dataset for Synthetic Control Analysis
# ============================================================

# This script prepares a dataset of Output Areas (OAs) for use in
# constructing synthetic control units in the CAZ evaluation.
# . Restrict to pre-treatment period (up to 2020-12-31).
# . Calculate baseline injury levels per OA (mean injuries by type/severity).
# . Estimate pre-treatment trends in injuries for each OA (slope per quarter).
# . Combine OA characteristics, road characteristics, baseline injuries, and trends
#    into a single OA-level matching dataset for synthetic control construction.
#=========================================================



library(tidyverse)
library(lubridate)
library(here)

# OA characteristics
OA_analysis <- readRDS(
  here("data","processed","OA_level_from_polygons.rds")
)

# road characteristics aggregated to OA
OA_roads <- readRDS(
  here("data","processed","OA_roads.rds")
)

# injuries aggregated by OA + quarter
OA_injuries <- readRDS(
  here("data","processed","OA_injuries_quarterly.rds")
)



#pre-treatment period## 
pre_period_end <- as.Date("2020-12-31")
inj_pre <- OA_injuries %>%
  filter(quarter_year <= pre_period_end)
summary(OA_injuries$quarter_year)
#


####Create baseline injury levels 

inj_baseline <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    
    mean_car_KSI = mean(Car_Van_KSI, na.rm = TRUE),
    mean_car_slight = mean(Car_Van_Slight, na.rm = TRUE),
    mean_cyc_KSI = mean(Cyclist_KSI, na.rm = TRUE),
    mean_cyc_slight = mean(Cyclist_Slight, na.rm = TRUE),
    mean_ped_KSI = mean(Pedestrian_KSI, na.rm = TRUE),
    mean_ped_slight = mean(Pedestrian_Slight, na.rm = TRUE),
    mean_total = mean(total_injuries, na.rm = TRUE),
    .groups = "drop"
  )

#############################################
#####Create pre-treatment trends
#############################################
#Convert quarter to numeric time
inj_pre <- inj_pre %>%
  arrange(OA, quarter_year) %>%
  group_by(OA) %>%
  mutate(
    time = row_number()
  ) %>%
  ungroup()


#Function to estimate slopes
trend_fun <- function(df, var){
  model <- lm(df[[var]] ~ df$time)
  coef(model)[2]
}

#### slopes per OA 



inj_trends <- inj_pre %>%
  group_by(OA) %>%
  summarise(
    trend_car_KSI    = lm(.data[["Car_Van_KSI"]] ~ time)$coefficients[2],
    trend_car_slight = lm(.data[["Car_Van_Slight"]] ~ time)$coefficients[2],
    trend_cyc_KSI    = lm(.data[["Cyclist_KSI"]] ~ time)$coefficients[2],
    trend_cyc_slight = lm(.data[["Cyclist_Slight"]] ~ time)$coefficients[2],
    trend_ped_KSI    = lm(.data[["Pedestrian_KSI"]] ~ time)$coefficients[2],
    trend_ped_slight = lm(.data[["Pedestrian_Slight"]] ~ time)$coefficients[2],
    trend_total      = lm(.data[["total_injuries"]] ~ time)$coefficients[2],
    .groups = "drop"
  )

#### # # ## #     OA matching Dataset

OA_matching_dataset <- OA_analysis %>%
  left_join(OA_roads, by = "OA") %>%
  left_join(inj_baseline, by = "OA") %>%
  left_join(inj_trends, by = "OA")

all(OA_matching_dataset$LAD24CD.x == OA_matching_dataset$LAD24CD.y)
OA_matching_dataset %>%
  summarise(
    same_LAD = sum(LAD24CD.x == LAD24CD.y, na.rm = TRUE),
    total = n()
  )

OA_matching_dataset %>%
  summarise(
    LAD_y_missing = sum(is.na(LAD24CD.y)),
    LAD_x_missing = sum(is.na(LAD24CD.x))
  )
#### #684 OAs do not have roads 
### OA_analysis is the master dataset 
### some clean up 

OA_matching_dataset <- OA_matching_dataset %>%
  select(
    OA,
    
    LAD24CD = LAD24CD.x,
    LAD24NM = LAD24NM.x,
    
    treated_OA = treated_OA.x,
    buffer_OA = buffer_OA.x,
    control_group1_OA = control_group1_OA.x,
    control_group2_OA = control_group2_OA.x,
    
    dist_citycentre = dist_citycentre.x,
    assignment = assignment.x,
    
    everything(),
    
    -ends_with(".x"),
    -ends_with(".y")
  )


names(OA_matching_dataset)

nrow(OA_matching_dataset)
length(unique(OA_matching_dataset$OA))


OA_matching_dataset %>%
  summarise(across(starts_with("mean"), ~sum(is.na(.))))

# #  the NAs because the was no injuries in these OAs 
#### replace with 0s 

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(across(starts_with("mean"), ~replace_na(.x, 0))) %>% 
  mutate(across(starts_with("trend"), ~replace_na(.x, 0)))


OA_matching_dataset %>% summarise(across(starts_with("mean"), ~sum(is.na(.))))
OA_matching_dataset %>% summarise(across(starts_with("trend"), ~sum(is.na(.))))



summary(OA_matching_dataset$mean_total)
summary(OA_matching_dataset$trend_total)   

saveRDS(
  OA_matching_dataset,
  here("data","processed","OA_matching_dataset.rds")
)
#









# ============================================================
# Variable descriptions for OA_matching_dataset
# ============================================================

var_description <- c(
  
  OA = "Output Area identifier (ONS 2011 OA code). Each row represents one OA.",
  
  n_roads = "Total number of road segments intersecting the OA.",
  
  total_road_length = "Total length of roads within the OA (meters).",
  
  n_A = "Number of A-road segments within the OA.",
  
  n_B = "Number of B-road segments within the OA.",
  
  n_motorway = "Number of motorway segments within the OA.",
  
  n_minor = "Number of minor/local road segments within the OA.",
  
  mean_car_KSI = "Mean quarterly number of car/van KSI (Killed or Seriously Injured) casualties in the OA during the pre-treatment period.",
  
  mean_car_slight = "Mean quarterly number of car/van slight injury casualties in the OA during the pre-treatment period.",
  
  mean_cyc_KSI = "Mean quarterly number of cyclist KSI casualties in the OA during the pre-treatment period.",
  
  mean_cyc_slight = "Mean quarterly number of cyclist slight injury casualties in the OA during the pre-treatment period.",
  
  mean_ped_KSI = "Mean quarterly number of pedestrian KSI casualties in the OA during the pre-treatment period.",
  
  mean_ped_slight = "Mean quarterly number of pedestrian slight injury casualties in the OA during the pre-treatment period.",
  
  mean_total = "Mean quarterly total number of casualties (all modes and severities) in the OA during the pre-treatment period.",
  
  trend_car_KSI = "Linear pre-treatment trend (slope per quarter) in car/van KSI casualties in the OA.",
  
  trend_car_slight = "Linear pre-treatment trend (slope per quarter) in car/van slight casualties in the OA.",
  
  trend_cyc_KSI = "Linear pre-treatment trend (slope per quarter) in cyclist KSI casualties in the OA.",
  
  trend_cyc_slight = "Linear pre-treatment trend (slope per quarter) in cyclist slight casualties in the OA.",
  
  trend_ped_KSI = "Linear pre-treatment trend (slope per quarter) in pedestrian KSI casualties in the OA.",
  
  trend_ped_slight = "Linear pre-treatment trend (slope per quarter) in pedestrian slight casualties in the OA.",
  
  trend_total = "Linear pre-treatment trend (slope per quarter) in total casualties in the OA."
)

# table 
var_description_table <- tibble(
  variable = names(var_description),
  description = unname(var_description)
)

print(var_description_table)



