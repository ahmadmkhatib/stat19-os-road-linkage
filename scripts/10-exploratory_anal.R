library(tidyverse)
library(lubridate)
library(here)

road_panel_complete<- readRDS(here("data", "processed", "road_panel_complete.rds"))


glimpse(road_panel_complete)

road_panel_complete %>% select(treated) %>% table(useNA = "always")

# check the dimensions  in the matrix
n_distinct(road_panel_complete$identifier) * n_distinct(road_panel_complete$quarter_year)


road_panel %>% 
  summarise(
    n_roads = n_distinct(identifier),
    n_quarters = n_distinct(quarter_year),
    total_rows = n()
  )



## balance - panel per road 


road_panel %>%
  count(identifier) %>%
  summarise(
    min_obs = min(n),
    max_obs = max(n)
  )


## ttt roads 
road_panel %>%
  distinct(identifier, ever_treated) %>%
  count(ever_treated
        )

#ttt timmimg 
road_panel %>%
  filter(ever_treated == 1) %>%
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




