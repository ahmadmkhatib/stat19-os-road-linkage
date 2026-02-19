library(sf)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(here)



injuries_sf <- readRDS(here("data","processed","injuries_final.rds"))

#structure check
print(injuries_sf)
glimpse(injuries_sf)
st_crs(injuries_sf)

# Check invalid geometries
invalid_geom <- sum(!st_is_valid(injuries_sf))
cat("Invalid geometries:", invalid_geom, "\n")

# Check missing geometry
missing_geom <- sum(is.na(st_geometry(injuries_sf)))
cat("Missing geometries:", missing_geom, "\n")

# ==========================================================
# MISSINGN SUMMARY
# ==========================================================

missing_summary <- injuries_sf %>%
  st_drop_geometry() %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to = "variable",
               values_to = "prop_missing") %>%
  arrange(desc(prop_missing))

print(missing_summary, n=113)

# ==========================================================
#  DUPLICATE CHECKS
# ==========================================================

# Duplicate collision IDs
if ("collision_ref_no" %in% names(injuries_sf)) {
  dup_collisions <- injuries_sf %>%
    st_drop_geometry() %>%
    count(collision_ref_no) %>%
    filter(n > 1)
  
  cat("Duplicate collision IDs:", nrow(dup_collisions), "\n")
}

# Duplicate collision + casualty combination
if (all(c("collision_ref_no", "casualty_reference") %in% names(injuries_sf))) {
  dup_casualty <- injuries_sf %>%
    st_drop_geometry() %>%
    count(collision_ref_no, casualty_reference) %>%
    filter(n > 1)
  
  cat("Duplicate casualty rows:", nrow(dup_casualty), "\n")
}



injuries_sf %>%
  st_drop_geometry() %>%
  count(injury_id) %>%
  filter(n > 1)

# ==========================================================
# 6. COORDINATE RANGE CHECK (IF OSGB)
# ==========================================================

coords <- st_coordinates(injuries_sf)

coord_summary <- tibble(
  min_easting  = min(coords[,1], na.rm = TRUE),
  max_easting  = max(coords[,1], na.rm = TRUE),
  min_northing = min(coords[,2], na.rm = TRUE),
  max_northing = max(coords[,2], na.rm = TRUE)
)

print(coord_summary)

# Check zero coordinates
zero_coords <- sum(coords[,1] == 0 | coords[,2] == 0, na.rm = TRUE)
cat("Zero coordinate rows:", zero_coords, "\n")

# ==========================================================
# TEMPORAL CONSISTENCY
# ==========================================================

if ("date" %in% names(injuries_sf)) {
  
  injuries_sf <- injuries_sf %>%
    mutate(date = as.Date(date))
  
  cat("Date range:",
      min(injuries_sf$date, na.rm = TRUE), "to",
      max(injuries_sf$date, na.rm = TRUE), "\n")
  
  print(table(year(injuries_sf$date)))
}

# ==========================================================
# ROAD CLASS DISTRIBUTIONS
# ==========================================================

if ("first_road_class1" %in% names(injuries_sf)) {
  print(injuries_sf %>% st_drop_geometry() %>% count(first_road_class1))
}

if ("second_road_class1" %in% names(injuries_sf)) {
  print(injuries_sf %>% st_drop_geometry() %>% count(second_road_class1))
}

if ("atJunction" %in% names(injuries_sf)) {
  print(injuries_sf %>% st_drop_geometry() %>% count(atJunction))
}

# ==========================================================
#  CASUALTY SEVERITY CHECK
# ==========================================================

if (all(c("casualty_severity", "casualty_severity1") %in% names(injuries_sf))) {
  print(
    injuries_sf %>%
      st_drop_geometry() %>%
      count(casualty_severity, casualty_severity1)
  )
}

# ==========================================================
#  MATCHING KEY COMPLETENESS
# ==========================================================

matching_keys <- c("USRN", "TOID")

existing_keys <- matching_keys[matching_keys %in% names(injuries_sf)]

if (length(existing_keys) > 0) {
  key_missing <- injuries_sf %>%
    st_drop_geometry() %>%
    summarise(across(all_of(existing_keys),
                     ~ mean(is.na(.))))
  
  print(key_missing)
}

# ==========================================================
# IDENTICAL LOCATION + DATE CHECK
# ==========================================================

if ("date" %in% names(injuries_sf)) {
  
  duplicate_spacetime <- injuries_sf %>%
    mutate(x = coords[,1],
           y = coords[,2]) %>%
    st_drop_geometry() %>%
    count(x, y, date) %>%
    filter(n > 1)
  
  cat("Identical space-time duplicates:",
      nrow(duplicate_spacetime), "\n")
}
 
# ==========================================================
# SUMMARY QA REPORT
# ==========================================================

qa_report <- injuries_sf %>%
  mutate(x = coords[,1],
         y = coords[,2]) %>%
  st_drop_geometry() %>%
  summarise(
    n_rows = n(),
    prop_missing_geometry = missing_geom / nrow(injuries_sf),
    prop_zero_coords = zero_coords / nrow(injuries_sf),
    min_year = if ("date" %in% names(.)) min(year(date), na.rm = TRUE) else NA,
    max_year = if ("date" %in% names(.)) max(year(date), na.rm = TRUE) else NA
  )

print(qa_report)


