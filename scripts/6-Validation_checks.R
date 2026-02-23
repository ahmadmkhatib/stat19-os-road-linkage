# ============================================================
# Script 05: Validation Checks for STAT19–OS Linkage
# ============================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)

# ------------------------------------------------------------
# Load matched injuries_with_oa data
# ------------------------------------------------------------

injuries_with_oa <- read_rds(here("data", "processed", "injuries_matched_OA.rds"))

# ------------------------------------------------------------
# Basic Structure Checks
# ------------------------------------------------------------

cat("Total injuries_with_oa:", nrow(injuries_with_oa), "\n")

# Check duplicates
dup_count <- sum(duplicated(injuries_with_oa$injury_id))
cat("Duplicate injury_id count:", dup_count, "\n")

# ------------------------------------------------------------
# Coordinate / Geometry Checks
# ------------------------------------------------------------

# Check if any geometries are empty
empty_geom <- sum(st_is_empty(injuries_with_oa))
cat("Empty geometries:", empty_geom, "\n")

# Check invalid geometries
invalid_geom <- sum(!st_is_valid(injuries_with_oa))
cat("Invalid geometries:", invalid_geom, "\n")

# CRS check
cat("CRS (EPSG):", st_crs(injuries_with_oa)$epsg, "\n")

# ------------------------------------------------------------
# Output Area (OA) Assignment Validation
# ------------------------------------------------------------

oa_summary <- injuries_with_oa %>%
  summarise(
    total = n(),
    missing_oa = sum(is.na(OA_CODE)),
    pct_missing = 100 * missing_oa / total
  )

print(oa_summary)

# ------------------------------------------------------------
# Road Linkage Distance Check (if matched road geometry available)
# ------------------------------------------------------------

if("road_geom" %in% names(injuries_with_oa)) {
  cat("Computing distances to matched roads...\n")
  
  distances <- st_distance(
    injuries_with_oa,
    injuries_with_oa$road_geom,
    by_element = TRUE
  )
  
  injuries_with_oa$road_distance_m <- as.numeric(distances)
  
  print(summary(injuries_with_oa$road_distance_m))
  
  extreme_matches <- injuries_with_oa %>%
    filter(road_distance_m > 100)
  
  cat("Injuries >100m from matched road:", nrow(extreme_matches), "\n")
}

# ------------------------------------------------------------
# Temporal Validation
# ------------------------------------------------------------

date_check <- injuries_with_oa %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE)
  )

print(date_check)

# Check for future dates
future_dates <- injuries_with_oa %>%
  filter(date > Sys.Date())

cat("Future-dated injuries:", nrow(future_dates), "\n")

# ------------------------------------------------------------
# 7. Optional: OA fallback / distance check
# ------------------------------------------------------------

# If needed, you can check nearest OA fallback distances:
# (This is just a rough spatial sanity check)
if(!all(is.na(injuries_with_oa$OA_CODE))) {
  cat("Checking OA centroid distances...\n")
  # centroids of OA polygons (approx)
  # only meaningful if you have OA polygons stored separately
  # injuries_with_oa %>% st_join(oa_2011) %>% st_distance(...)
}

# ------------------------------------------------------------
# Save Summary
# ------------------------------------------------------------

validation_report <- list(
  total_injuries = nrow(injuries_with_oa),
  duplicate_injuries = dup_count,
  missing_oa = oa_summary$missing_oa,
  pct_missing_oa = oa_summary$pct_missing,
  empty_geometries = empty_geom,
  invalid_geometries = invalid_geom,
  min_date = date_check$min_date,
  max_date = date_check$max_date,
  future_dates = nrow(future_dates)
)

view(validation_report)

write_rds(
  validation_report,
  here("data", "processed", "validation_summary.rds")
)

cat("Validation complete. Summary saved to 'validation_summary.rds'.\n")
