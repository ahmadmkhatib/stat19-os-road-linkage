# ==========================================================
# STAT19 – OS Open Roads Linkage Framework
# Script: 03_match_injuries_to_roads_by_type.R
# Purpose: Match STATS19 injuries to nearest plausible road link
# CRS: British National Grid (EPSG:27700)
# ==========================================================

library(sf)
library(dplyr)
library(here)

# ----------------------------------------------------------
# Load Processed Data
# ----------------------------------------------------------

injuries_path <- here("data", "processed", "injuries.rds")
roads_path    <- here("data", "processed", "roads_filtered.rds")

injuries <- readRDS(injuries_path)
roads    <- readRDS(roads_path)

# Ensure CRS consistency
stopifnot(st_crs(injuries) == st_crs(roads))

## checks
table(injuries$injury_class, useNA = "ifany")
any(duplicated(injuries$injury_id))
# ----------------------------------------------------------
# Split Roads by Class
# --------------------------------------------------------

roads_A        <- roads %>% filter(road_class == "A")
roads_B        <- roads %>% filter(road_class == "B")
roads_motorway <- roads %>% filter(road_class == "Motorway")
roads_minor    <- roads %>% filter(road_class == "minor")

# ----------------------------------------------------------
# Matching Function
# ----------------------------------------------------------

match_one_class <- function(injuries, roads_same, roads_any,
                            max_same = 50, max_any = 100) {
  
  # nearest ANY class
  idx_any <- st_nearest_feature(injuries, roads_any)
  dist_any <- as.numeric(
    st_distance(injuries, roads_any[idx_any, ], by_element = TRUE)
  )
  
  # nearest SAME class
  idx_same <- st_nearest_feature(injuries, roads_same)
  dist_same <- as.numeric(
    st_distance(injuries, roads_same[idx_same, ], by_element = TRUE)
  )
  
  injuries %>%
    mutate(
      matched_roadID = case_when(
        dist_same <= max_same ~ roads_same$identifier[idx_same],
        dist_any  <= max_any  ~ roads_any$identifier[idx_any],
        TRUE ~ NA_character_
      ),
      match_type = case_when(
        dist_same <= max_same ~ "same_class",
        dist_any  <= max_any  ~ "fallback_any",
        TRUE ~ "dropped"
      ),
      dist_same = dist_same,
      dist_any  = dist_any
    ) %>%
    filter(match_type != "dropped")
}
# ----------------------------------------------------------
# Perform Matching by Injury Class
# ----------------------------------------------------------

matched_A <- match_one_class(
  injuries %>% filter(injury_class == "A"),
  roads_A,
  roads
)

matched_B <- match_one_class(
  injuries %>% filter(injury_class == "B"),
  roads_B,
  roads
)

matched_motorway <- match_one_class(
  injuries %>% filter(injury_class == "Motorway"),
  roads_motorway,
  roads
)

matched_minor <- match_one_class(
  injuries %>% filter(injury_class == "minor"),
  roads_minor,
  roads
)

matched <- bind_rows(
  matched_A,
  matched_B,
  matched_motorway,
  matched_minor
)

# ----------------------------------------------------------
# Matching Diagnostics (injury-level)
# ----------------------------------------------------------

# total unique injuries entering matching
n_total <- injuries %>%
  distinct(injury_id) %>%
  nrow()

# unique injuries successfully matched
n_matched <- matched %>%
  distinct(injury_id) %>%
  nrow()

n_unmatched <- n_total - n_matched

summary_table <- tibble(
  injuries_total    = n_total,
  injuries_matched  = n_matched,
  injuries_unmatched = n_unmatched,
  pct_unmatched     = 100 * n_unmatched / n_total
)

print(summary_table)


####  -------
# Check duplicates
any(duplicated(matched$injury_id))

sum(duplicated(matched$injury_id))




# ----------------------------------------------------------
# Save Output
# ----------------------------------------------------------

output_path <- here("data", "processed", "injuries_matched.rds")

saveRDS(matched, output_path)

cat("Matching complete. Output saved to:", output_path)
