# ============================================================
# Merge OA Matching Dataset with Census Characteristics
# ============================================================

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(naniar)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
### census OA data 
OA_char_raw <- read.csv(here("data","processed","outputArea_raw.csv"))
OA_char_percent <- read.csv(here("data","processed","outputArea_percent.csv"))

miss_var_summary(OA_char_raw)
miss_var_summary(OA_char_percent)


OA_char_percent <-  OA_char_percent %>% filter(!is.na(Bicycle))

### business  OA data 
OA__EW_businesses <- read.csv(here("data","processed","OA_EW_businesses.csv"))
OA__scot_businesses<- read.csv(here("data","processed","OA_Sco_businesses.csv"))


OA_EW_businesses_clean <- OA__EW_businesses %>%
  rename(OA = OA21CD)
OA_Sco_businesses_clean <- OA__scot_businesses %>%
  rename(OA = OA22) 

# Combine datasets
OA_businesses <- bind_rows(
  OA_EW_businesses_clean,
  OA_Sco_businesses_clean
)


OA_matching_dataset <- readRDS(
  here("data","processed","OA_matching_dataset.rds")
)

# ------------------------------------------------------------
# Load OA shapefile and compute area
# ------------------------------------------------------------

oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  mutate(
    area_m2 = as.numeric(st_area(.)),
    area_km2 = area_m2 / 1e6
  )

# ------------------------------------------------------------
# Pre-merge checks
# ------------------------------------------------------------

cat("OA_matching_dataset rows:", nrow(OA_matching_dataset), "\n")
cat("OA_char_raw rows:", nrow(OA_char_raw), "\n")
cat("OA_char_percent rows:", nrow(OA_char_percent), "\n")

cat("OA_char_raw dups:",
    OA_char_raw %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

cat("OA_char_percent dups:",
    OA_char_percent %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")
cat("OA__EW_businesses dups:",
    OA__EW_businesses %>% count(OA21CD) %>% filter(n > 1) %>% nrow(), "\n")

OA_businesses %>% count(OA) %>% filter(n > 1)%>% nrow()


OA_businesses %>% 
  count(OA, sort = TRUE)

OA_businesses %>% 
  filter(OA %in% (OA_businesses %>% count(OA) %>% filter(n > 1) %>% pull(OA))) %>%
  head(20)

OA_businesses <- OA_businesses %>% distinct()

OA_businesses %>% 
  count(OA) %>% 
  filter(n > 1)

### still 4 duplicates 
OA_businesses <- OA_businesses %>%
  group_by(OA) %>%
  summarise(
    business_retail_per_km2 = sum(business_retail_per_km2, na.rm = TRUE),
    business_accommodation_food_per_km2 = sum(business_accommodation_food_per_km2, na.rm = TRUE),
    .groups = "drop"
  )



# ------------------------------------------------------------
# Rename census variables
# ------------------------------------------------------------

vars_to_rename <- setdiff(
  names(OA_char_raw),
  c("OA","country","Total","IMD")
)

OA_char_raw_renamed <- OA_char_raw %>%
  rename_with(~ paste0(.x,"_n"), all_of(vars_to_rename))

names(OA_char_percent)
OA_char_pct_renamed <- OA_char_percent %>%
  rename_with(~ paste0(.x,"_pct"), all_of(vars_to_rename)) %>%
  select(-country,-Total,-IMD)

# ------------------------------------------------------------
# Merge census tables
# ------------------------------------------------------------

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_renamed, by="OA")

# ------------------------------------------------------------
# Merge census onto OA matching dataset
# ------------------------------------------------------------

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by="OA")


#'### ------- add busnesses 
OA_matching_census <- OA_matching_dataset %>%
left_join(OA_census, by = "OA") %>%
  left_join(OA_businesses, by = "OA")


# ------------------------------------------------------------
# Add OA area + population density
# ------------------------------------------------------------

OA_matching_census <- OA_matching_census %>%
  left_join(
    oa_sub %>%
      st_drop_geometry() %>%
      select(OA, area_km2),
    by = "OA"
  ) %>%
  mutate(
    pop_density = Total / area_km2,
    log_pop_density = log1p(pop_density)
  )



table(OA_matching_census$pop_density)
# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

stopifnot(nrow(OA_matching_census) == nrow(OA_matching_dataset))

stopifnot(
  OA_matching_census %>%
    count(OA) %>%
    filter(n > 1) %>%
    nrow() == 0
)



naniar::miss_var_summary(OA_matching_census)

#### remoive the OA with nas 
OA_matching_census <- OA_matching_census %>%
  drop_na(
    White_pct, Mixed_pct, Asian_pct, Black_pct, Other_ethnicity_pct,
    X4under_pct, X5to9_pct, X10to14_pct, X15to19_pct
     )

OA_matching_census <- OA_matching_census %>%
  mutate(scheme = replace_na(scheme, "Control"))
# ------------------------------------------------------------
# ------------------------------------------------------------

saveRDS(
  OA_matching_census,
  here("data","processed","OA_matching_census.rds")
)

# OA_matching_census<- readRDS(here("data","processed","OA_matching_census.rds"))
names(OA_matching_census)


# ------------------------------------------------------------
# Create spatial version
# ------------------------------------------------------------

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by="OA") %>%
  st_as_sf()

# ------------------------------------------------------------
# Save spatial dataset
# ------------------------------------------------------------

st_write(
  OA_matching_census_sf,
  here("data","processed","shp_files","OA_matching_census.gpkg"),
  delete_dsn = TRUE
)


