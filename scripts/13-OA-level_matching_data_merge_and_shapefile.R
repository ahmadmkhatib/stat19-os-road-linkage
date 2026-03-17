# ============================================================
# Merge OA Matching Dataset with Census Characteristics
# ============================================================
# This script merges the OA-level matching dataset with 2021 Census
# characteristics for use in synthetic DiD matching.
#
# Two census files are merged:
#   - outputArea_raw.csv    : absolute counts (population, households, etc.)
#   - outputArea_percent.csv: percentage equivalents of the same variables
#
# Variables cover:
#   - Ethnicity (White, Mixed, Asian, Black, Other)
#   - Age groups (under 4 through 85+, plus broad bands)
#   - Car ownership (none, one, two, two+, three+)
#   - Travel to work mode (drive, walk, cycle, bus, train, etc.)
#   - Employment (all in work, work at home)
#   - IMD score
#
# Output: data/processed/OA_matching_census.rds
#         data/processed/shp_files/OA_matching_census.gpkg
# ============================================================

#=========================================================
library(tidyverse)
library(lubridate)
library(here)
library(sf)

#### Jamies data numbers and percent in each OA 
OA_char_raw<-read.csv(here("data","processed","outputArea_raw.csv"))
OA_char_percent<-read.csv(here("data","processed","outputArea_percent.csv"))
names(OA_char_raw)
names(OA_char_percent)


OA_matching_dataset <- readRDS(
  here("data", "processed", "OA_matching_dataset.rds")
)



# ── Pre-merge checks ──────────────────────────────────────────────────────────

cat("OA_matching_dataset rows:", nrow(OA_matching_dataset), "\n")
cat("OA_char_raw rows:        ", nrow(OA_char_raw), "\n")
cat("OA_char_percent rows:    ", nrow(OA_char_percent), "\n")

# Duplicates in census files
cat("OA_char_raw dups:    ", OA_char_raw     %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")
cat("OA_char_percent dups:", OA_char_percent %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

# OA code format check
cat("Matching dataset OA sample:", head(OA_matching_dataset$OA, 3), "\n")
cat("Census raw OA sample:      ", head(OA_char_raw$OA, 3), "\n")
cat("Census percent OA sample:  ", head(OA_char_percent$OA, 3), "\n")

# How many matching dataset OAs are in each census file?
cat("OAs in matching + raw:    ",
    sum(OA_matching_dataset$OA %in% OA_char_raw$OA), "/",
    nrow(OA_matching_dataset), "\n")
cat("OAs in matching + percent:",
    sum(OA_matching_dataset$OA %in% OA_char_percent$OA), "/",
    nrow(OA_matching_dataset), "\n")

# OAs in matching dataset but missing from census
missing_from_raw <- anti_join(
  OA_matching_dataset %>% select(OA),
  OA_char_raw, by = "OA"
)
missing_from_pct <- anti_join(
  OA_matching_dataset %>% select(OA),
  OA_char_percent, by = "OA"
)
cat("OAs missing from raw:    ", nrow(missing_from_raw), "\n")
cat("OAs missing from percent:", nrow(missing_from_pct), "\n")


# ── Rename census columns before merging ─────────────────────────────────────
# Raw counts get _n suffix, percentages get _pct suffix
# OA, country, Total, IMD are kept without suffix (same in both or non-numeric)

vars_to_rename <- setdiff(
  names(OA_char_raw),
  c("OA", "country", "Total", "IMD")
)

OA_char_raw_renamed <- OA_char_raw %>%
  rename_with(~ paste0(.x, "_n"), all_of(vars_to_rename))

OA_char_pct_renamed <- OA_char_percent %>%
  rename_with(~ paste0(.x, "_pct"), all_of(vars_to_rename)) %>%
  select(-country, -Total, -IMD)   # already in raw, don't duplicate

# ── Merge census files together ───────────────────────────────────────────────
# Join raw and percent census into one wide census table first,
# then join onto matching dataset

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_renamed, by = "OA")

cat("OA_census rows after joining raw + pct:", nrow(OA_census), "\n")
cat("OA_census dups:", OA_census %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

# ── Merge onto OA matching dataset ───────────────────────────────────────────

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by = "OA")

cat("OA_matching_census rows:", nrow(OA_matching_census), "\n")
cat("Duplicate OAs:          ", OA_matching_census %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

# ── Post-merge checks ─────────────────────────────────────────────────────────

# 1. Row count unchanged
stopifnot(nrow(OA_matching_census) == nrow(OA_matching_dataset))

# 2. No new duplicate OAs introduced
stopifnot(OA_matching_census %>% count(OA) %>% filter(n > 1) %>% nrow() == 0)

# 3. NA rates in census variables for treated OAs
OA_matching_census %>%
  filter(treated_OA == 1) %>%
  summarise(across(
    c(ends_with("_n"), ends_with("_pct")),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  filter(n_NA > 0) %>%
  arrange(desc(n_NA)) %>%
  print()

# 4. Sense check key census variables for treated vs control
OA_matching_census %>%
  filter(treated_OA == 1 | control_group2_OA == 1) %>%
  mutate(group = if_else(treated_OA == 1, "Treated", "Control")) %>%
  group_by(group) %>%
  summarise(
    n_OAs            = n(),
    mean_population  = mean(Total,           na.rm = TRUE),
    mean_pct_drive   = mean(Drive_Car_pct,   na.rm = TRUE),
    mean_pct_walk    = mean(Walk_pct,        na.rm = TRUE),
    mean_pct_cycle   = mean(Bicycle_pct,     na.rm = TRUE),
    mean_pct_no_car  = mean(cars_none_pct,   na.rm = TRUE),
    mean_IMD         = mean(IMD,             na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()




# ── Save ──────────────────────────────────────────────────────────────────────

saveRDS(
  OA_matching_census,
  here("data", "processed", "OA_matching_census.rds")
)





# Spatial version
oa_sub <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by = "OA") %>%
  st_as_sf()

st_write(
  OA_matching_census_sf,
  here("data", "processed", "shp_files", "OA_matching_census.gpkg"),
  delete_dsn = TRUE
)

cat("Saved: OA_matching_census.rds and OA_matching_census.gpkg\n")






OA_matching_dataset<- readRDS(here("data","processed","OA_matching_dataset.rds"))
glimpse(OA_matching_dataset)

####### as a shapefile   # # # # 



# OA shpfile 
oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid()



# Join geometry to matching dataset
OA_matching_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_dataset, by = "OA") %>%
  st_as_sf()

# Quick check
print(OA_matching_sf)


table(st_is_valid(OA_matching_sf))
st_geometry_type(OA_matching_sf) %>% table()


# Save as GeoPackage
st_write(
  OA_matching_sf,
  here("data","processed","shp_files","OA_matching_dataset.gpkg"),
  delete_dsn = TRUE
)


OA_matching_sf<- st_read(here("data","processed","shp_files","OA_matching_dataset.gpkg"))

OA_matching_sf %>%
  ggplot() +
  geom_sf(aes(fill = treated_OA), colour = NA) +
  theme_minimal()
## check the CAZ shapfile 


