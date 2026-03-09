# ============================================================
# Exploratory & Descriptive Summaries
# Road-Level Panel + OA-Level Data
# ============================================================

library(tidyverse)
library(here)
library(arrow)
library(sf)
library(zoo)
library(knitr)

# -----------------------------
# Load Relevant Data
# -----------------------------

# Road-level panel and classification
road_panel_model<-arrow::open_dataset(here("data","processed","road_panel_dataset")) %>% collect()
analysis_roads   <- readRDS(here("data", "processed", "analysis_roads.rds"))

# OA-level data
OA_analysis <- readRDS(here("data", "processed", "OA_level_from_polygons.rds"))

# Injuries data
injuries <- readRDS(here("data", "processed", "injuries_matched_final.rds"))


# -----------------------------
# Road-Level Counts
# -----------------------------
road_summary <- analysis_roads %>%
  summarise(
    n_roads_total = n(),
    n_treated_any = sum(treated_group_any),
    n_treated50 = sum(treated_group_50pct),
    n_control1 = sum(control_group1),
    n_control2 = sum(control_group2)
  )

# -----------------------------
#  Injury Summaries by Treatment/Control
inj_summary <- road_panel_model %>%
  group_by(treated_group_any,treated_group_50pct, control_group1, control_group2) %>%
  summarise(
    total_KSI = sum(across(starts_with("KSI_adj")), na.rm = TRUE),
    total_Slight = sum(across(starts_with("Slight_adj")), na.rm = TRUE),
    total_injuries = sum(across(starts_with("total_inj_adj")), na.rm = TRUE),
    n_roads = n_distinct(identifier),
    .groups = "drop"
  )
# -----------------------------
# Quarterly Injury Distribution
# -----------------------------
quarterly_summary <- road_panel_model %>%
  group_by(quarter_year) %>%
  summarise(
    total_KSI = sum(across(starts_with("KSI_adj")), na.rm = TRUE),
    total_Slight = sum(across(starts_with("Slight_adj")), na.rm = TRUE),
    total_injuries = sum(across(starts_with("total_inj_adj")), na.rm = TRUE)
  )

# -----------------------------
#  OA-Level Descriptive Statistics
# -----------------------------
OA_summary <- OA_analysis %>%
  summarise(
    n_OAs_total = n(),
    n_treated_any = sum(treated_OA_any),
    n_treated_50 =sum(treated_OA_50pct),
    n_buffer = sum(buffer_OA),
    n_control1 = sum(control_group1_OA),
    n_control2 = sum(control_group2_OA),
    median_dist_to_centre = median(dist_citycentre),
    mean_dist_to_centre = mean(dist_citycentre),
    min_dist_to_centre = min(dist_citycentre),
    max_dist_to_centre = max(dist_citycentre)
  )






road_summary %>% st_drop_geometry() %>%kable(
  caption = "Road-Level Counts by Treatment and Control Groups"
)


inj_summary %>% st_drop_geometry() %>%kable(
  caption = "Total Injuries by Road Group"
)

quarterly_summary %>% st_drop_geometry() %>% kable(
  caption = "Quarterly Injury Distribution"
)

OA_summary %>% st_drop_geometry() %>%kable(
  caption = "OA-Level Summary Statistics"
)





library(flextable)
library(officer)

doc <- read_docx()

doc <- body_add_flextable(doc,
                          flextable(st_drop_geometry(road_summary)))

doc <- body_add_flextable(doc,
                          flextable(st_drop_geometry(inj_summary)))

print(doc, target = "tables.docx")
















glimpse(road_panel_model)

road_panel_model %>% select(treated_any) %>% table(useNA = "always")

# check the dimensions
n_distinct(road_panel_model$identifier) * n_distinct(road_panel_model$quarter_year)


road_panel_model %>% 
  summarise(
    n_roads = n_distinct(identifier),
    n_quarters = n_distinct(quarter_year),
    total_rows = n()
  )



## balance - panel per road 

road_panel_model %>%
  count(identifier) %>%
  summarise(
    min_obs = min(n),
    max_obs = max(n)
  )


## ever treated roads 
road_panel_model %>%
  distinct(identifier, treated_group_any) %>%
  count(treated_group_any
        )





#ttt timmimg 
road_panel_model %>%
  filter(treated_group_any == 1) %>%
  distinct(identifier, caz_start_q) %>%
  count(caz_start_q)

#Visualise treatment rollout 
road_panel %>%
  group_by(quarter_year) %>%
  summarise(
    n_treated = sum(treated)
  ) %>%
  ggplot(aes(quarter_year, n_treated)) +
  geom_line() +
  labs(title = "Number of Treated Roads Over Time")


# distibution 
road_panel %>%
  summarise(
    mean_inj = mean(total_inj_adj),
    var_inj  = var(total_inj_adj),
    zero_prop = mean(total_inj_adj == 0)
  )




road_panel %>%
  ggplot(aes(total_inj_adj)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0, 10)) +
  labs(title = "Distribution of Quarterly Injuries per Road")

#road avarages 

road_panel %>%
  group_by(identifier) %>%
  summarise(avg_inj = mean(total_inj_adj)) %>%
  ggplot(aes(avg_inj)) +
  geom_histogram(bins = 40) +
  labs(title = "Average Injuries per Road")



road_panel %>%
  group_by(quarter_year, ever_treated) %>%
  summarise(mean_inj = mean(total_inj_adj), .groups = "drop") %>%
  ggplot(aes(quarter_year, mean_inj, colour = factor(ever_treated))) +
  geom_line(size = 1) +
  labs(
    title = "Pre-Treatment Trends: Treated vs Control Roads",
    colour = "Ever Treated"
  )


min_caz <- min(road_panel$caz_start_q, na.rm = TRUE)

pre_data <- road_panel %>%
  filter(quarter_year < min_caz)

pre_data %>%
  group_by(ever_treated) %>%
  summarise(
    mean_inj = mean(total_inj_adj),
    mean_length = mean(length, na.rm = TRUE),
    prop_A_road = mean(road_class == "A", na.rm = TRUE),
    n_roads = n_distinct(identifier)
  )

## hetrogenity byu calss 


road_panel %>%
  group_by(road_class, ever_treated) %>%
  summarise(mean_inj = mean(total_inj_adj), .groups = "drop") %>%
  ggplot(aes(road_class, mean_inj, fill = factor(ever_treated))) +
  geom_col(position = "dodge")



## serial corelation check 

top_road <- road_panel %>%
  group_by(identifier) %>%
  summarise(avg = mean(total_inj_adj)) %>%
  arrange(desc(avg)) %>%
  slice(1) %>%
  pull(identifier)

road_panel %>%
  filter(identifier == top_road) %>%
  ggplot(aes(quarter_year, total_inj_adj)) +
  geom_line()


#   event structure 


road_panel <- road_panel %>%
  mutate(
    event_time = as.numeric(quarter_year - caz_start_q)
  )


road_panel %>%
  filter(ever_treated == 1) %>%
  ggplot(aes(event_time)) +
  geom_histogram(bins = 30)


### zero inf 

road_panel %>%
  group_by(ever_treated) %>%
  summarise(
    zero_rate = mean(total_inj_adj == 0)
  )



## agregation 

road_panel %>%
  mutate(year = year(as.Date(quarter_year))) %>%
  group_by(identifier, year) %>%
  summarise(yearly_inj = sum(total_inj_adj), .groups = "drop") %>%
  summarise(
    mean_yearly = mean(yearly_inj),
    zero_prop_yearly = mean(yearly_inj == 0)
  )




