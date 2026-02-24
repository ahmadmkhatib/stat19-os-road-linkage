library(tidyverse)
library(lubridate)
library(here)

road_panel_complete<- readRDS(here("data", "processed", "road_panel_complete.rds"))


glimpse(road_panel_complete)

n_distinct(road_panel_complete$identifier) * n_distinct(road_panel_complete$quarter_year)
