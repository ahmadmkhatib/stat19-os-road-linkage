# ==============================================================================
# STATS19 - OS Open Roads linkage
# Script: 4_match_injuries_to_roads_by_type.R
# Purpose: Link each STATS19 injury to the nearest plausible OS Open Roads link
# CRS: OSGB36 / British National Grid (EPSG:27700)
#
# Linkage hierarchy (adapted from Furlong et al., 2025):
#   1. Match to the nearest road link of the corresponding STATS19 road class
#      when that link is no more than 50 m away.
#   2. Otherwise, match to the nearest road link of any class when it is no more
#      than 100 m away.
#   3. Exclude records that do not meet either rule.
#
# C and Unclassified STATS19 roads are combined as "minor" to correspond with
# the minor-road category in the processed OS Open Roads data.
#
# Junctions: STATS19 first road class identifies the road on which the collision
# occurred, or the priority road when the collision occurred at a junction.
# Therefore, first_road_class_label1 is also the appropriate class for the
# initial class-restricted search at junctions. The second road class is not used.
# ==============================================================================

library(sf)
library(tidyverse)
library(here)

injuries <- readRDS(here("data", "processed", "injuries_final.rds"))
roads    <- readRDS(here("data", "processed", "roads_filtered.rds"))

required_injury_fields <- c(
  "injury_id", "first_road_class_label1", "LAD24CD", "LAD24NM"
)
required_road_fields <- c("identifier", "road_class")

missing_injury_fields <- setdiff(required_injury_fields, names(injuries))
missing_road_fields   <- setdiff(required_road_fields, names(roads))

if (length(missing_injury_fields) > 0) {
  stop(
    "Missing required injury field(s): ",
    paste(missing_injury_fields, collapse = ", ")
  )
}

if (length(missing_road_fields) > 0) {
  stop(
    "Missing required road field(s): ",
    paste(missing_road_fields, collapse = ", ")
  )
}

if (anyDuplicated(injuries$injury_id) > 0) {
  stop("injury_id must uniquely identify rows before spatial linkage.")
}

if (is.na(st_crs(injuries)) || is.na(st_crs(roads))) {
  stop("Both spatial datasets must have a defined coordinate reference system.")
}

if (st_crs(injuries)$epsg != 27700 || st_crs(roads)$epsg != 27700) {
  stop("Both datasets must use British National Grid (EPSG:27700).")
}

if (any(st_is_empty(injuries)) || any(st_is_empty(roads))) {
  stop("Empty geometries must be resolved before spatial linkage.")
}

if (!all(st_geometry_type(injuries) %in% c("POINT", "MULTIPOINT"))) {
  stop("The injury dataset must contain point geometries.")
}

if (!all(st_geometry_type(roads) %in% c("LINESTRING", "MULTILINESTRING"))) {
  stop("The road dataset must contain line geometries.")
}

# ------------------------------------------------------------------------------
# 2. Harmonise road classifications
# ------------------------------------------------------------------------------

recode_road_class <- function(x) {
  x_clean <- str_to_lower(str_squish(as.character(x)))
  
  case_when(
    x_clean %in% c("motorway", "m")             ~ "Motorway",
    x_clean %in% c("a", "a road")              ~ "A",
    x_clean %in% c("b", "b road")              ~ "B",
    x_clean %in% c(
      "c", "c road", "unclassified", "minor", "minor road"
    )                                             ~ "minor",
    is.na(x_clean) | x_clean == ""               ~ NA_character_,
    TRUE                                           ~ NA_character_
  )
}

injuries <- injuries %>%
  mutate(road_class = recode_road_class(first_road_class_label1))

roads <- roads %>%
  mutate(road_class = recode_road_class(road_class))

# Fail visibly if a populated source value was not recognised by the recode.
unknown_injury_classes <- injuries %>%
  st_drop_geometry() %>%
  filter(
    !is.na(first_road_class_label1),
    str_squish(as.character(first_road_class_label1)) != "",
    is.na(road_class)
  ) %>%
  distinct(first_road_class_label1) %>%
  pull(first_road_class_label1)

if (length(unknown_injury_classes) > 0) {
  stop(
    "Unrecognised STATS19 road class value(s): ",
    paste(unknown_injury_classes, collapse = ", ")
  )
}

if (any(is.na(roads$road_class))) {
  stop("Some OS road classes are missing or unrecognised after recoding.")
}

expected_classes <- c("Motorway", "A", "B", "minor")
missing_road_classes <- setdiff(expected_classes, unique(roads$road_class))

if (length(missing_road_classes) > 0) {
  stop(
    "No OS road links are available for class(es): ",
    paste(missing_road_classes, collapse = ", ")
  )
}

roads_by_class <- split(roads, roads$road_class)

# ------------------------------------------------------------------------------
# 3. Matching function
# ------------------------------------------------------------------------------

match_one_class <- function(injury_subset, roads_same, roads_any,
                            max_same_m = 50, max_any_m = 100) {
  if (nrow(injury_subset) == 0) {
    return(injury_subset %>%
             mutate(
               matched_roadID = character(),
               matched_road_class = character(),
               match_type = character(),
               matched_distance_m = double(),
               nearest_same_distance_m = double(),
               nearest_any_distance_m = double()
             ))
  }
  
  if (nrow(roads_same) == 0 || nrow(roads_any) == 0) {
    stop("A road candidate set supplied to match_one_class() is empty.")
  }
  
  nearest_same_index <- st_nearest_feature(injury_subset, roads_same)
  nearest_any_index  <- st_nearest_feature(injury_subset, roads_any)
  
  nearest_same_distance_m <- as.numeric(
    st_distance(
      injury_subset,
      roads_same[nearest_same_index, ],
      by_element = TRUE
    )
  )
  
  nearest_any_distance_m <- as.numeric(
    st_distance(
      injury_subset,
      roads_any[nearest_any_index, ],
      by_element = TRUE
    )
  )
  
  use_same     <- nearest_same_distance_m <= max_same_m
  use_fallback <- !use_same & nearest_any_distance_m <= max_any_m
  
  injury_subset %>%
    mutate(
      matched_roadID = case_when(
        use_same ~ as.character(roads_same$identifier[nearest_same_index]),
        use_fallback ~ as.character(roads_any$identifier[nearest_any_index]),
        TRUE ~ NA_character_
      ),
      matched_road_class = case_when(
        use_same ~ as.character(roads_same$road_class[nearest_same_index]),
        use_fallback ~ as.character(roads_any$road_class[nearest_any_index]),
        TRUE ~ NA_character_
      ),
      match_type = case_when(
        use_same ~ "same_class",
        use_fallback ~ "fallback_any",
        TRUE ~ "excluded_over_100m"
      ),
      # This is the distance to the road actually selected. It must be used for
      # linkage diagnostics instead of nearest_any_distance_m for every record.
      matched_distance_m = case_when(
        use_same ~ nearest_same_distance_m,
        use_fallback ~ nearest_any_distance_m,
        TRUE ~ NA_real_
      ),
      nearest_same_distance_m = nearest_same_distance_m,
      nearest_any_distance_m = nearest_any_distance_m
    )
}

# ------------------------------------------------------------------------------
# 4. Perform matching separately by STATS19 first road class
# ------------------------------------------------------------------------------

classified_injuries <- injuries %>%
  filter(!is.na(road_class))

unclassified_injuries <- injuries %>%
  filter(is.na(road_class)) %>%
  mutate(
    matched_roadID = NA_character_,
    matched_road_class = NA_character_,
    match_type = "excluded_missing_road_class",
    matched_distance_m = NA_real_,
    nearest_same_distance_m = NA_real_,
    nearest_any_distance_m = NA_real_
  )

match_results <- map_dfr(
  expected_classes,
  function(class_name) {
    match_one_class(
      injury_subset = classified_injuries %>%
        filter(road_class == class_name),
      roads_same = roads_by_class[[class_name]],
      roads_any = roads,
      max_same_m = 50,
      max_any_m = 100
    )
  }
)

# Retain all outcomes until diagnostics are complete.
linkage_audit <- bind_rows(match_results, unclassified_injuries)

# Final analytical linkage contains successful matches only.
matched <- linkage_audit %>%
  filter(match_type %in% c("same_class", "fallback_any"))

excluded <- linkage_audit %>%
  filter(!match_type %in% c("same_class", "fallback_any"))

# ------------------------------------------------------------------------------
# 5. Diagnostics printed to the console
# ------------------------------------------------------------------------------

linkage_summary <- linkage_audit %>%
  st_drop_geometry() %>%
  summarise(
    injuries_total = n_distinct(injury_id),
    injuries_with_valid_class = n_distinct(injury_id[!is.na(road_class)]),
    injuries_matched = n_distinct(
      injury_id[match_type %in% c("same_class", "fallback_any")]
    ),
    excluded_missing_class = n_distinct(
      injury_id[match_type == "excluded_missing_road_class"]
    ),
    excluded_over_100m = n_distinct(
      injury_id[match_type == "excluded_over_100m"]
    )
  ) %>%
  mutate(
    injuries_unmatched = injuries_total - injuries_matched,
    pct_matched = 100 * injuries_matched / injuries_total,
    pct_unmatched = 100 * injuries_unmatched / injuries_total
  )

match_type_summary <- linkage_audit %>%
  st_drop_geometry() %>%
  count(match_type, name = "n") %>%
  mutate(pct_of_all = 100 * n / sum(n))

distance_summary <- matched %>%
  st_drop_geometry() %>%
  summarise(
    mean_m = mean(matched_distance_m, na.rm = TRUE),
    median_m = median(matched_distance_m, na.rm = TRUE),
    p90_m = as.numeric(quantile(matched_distance_m, 0.90, na.rm = TRUE)),
    max_m = max(matched_distance_m, na.rm = TRUE)
  )

class_summary <- matched %>%
  st_drop_geometry() %>%
  count(road_class, match_type, name = "n") %>%
  group_by(road_class) %>%
  mutate(pct_within_class = 100 * n / sum(n)) %>%
  ungroup()

cat("\n--- Linkage summary ---\n")
print(linkage_summary)

cat("\n--- Match types ---\n")
print(match_type_summary)

cat("\n--- Distance to the selected road link ---\n")
print(distance_summary)

cat("\n--- Match type by STATS19 road class ---\n")
print(class_summary, n = Inf)

cat("\n--- Matched injuries by local authority ---\n")
matched %>%
  st_drop_geometry() %>%
  count(LAD24CD, LAD24NM, name = "n") %>%
  arrange(LAD24NM) %>%
  print(n = Inf)

cat("\n--- Scottish local authorities present after linkage ---\n")
matched %>%
  st_drop_geometry() %>%
  filter(str_starts(LAD24CD, "S")) %>%
  distinct(LAD24CD, LAD24NM) %>%
  arrange(LAD24NM) %>%
  print(n = Inf)

# These checks ensure that every input injury has exactly one audited outcome.
stopifnot(nrow(linkage_audit) == nrow(injuries))
stopifnot(n_distinct(linkage_audit$injury_id) == nrow(injuries))
stopifnot(all(!is.na(matched$matched_roadID)))
stopifnot(all(matched$matched_distance_m <= 100))
stopifnot(all(
  matched$matched_distance_m[matched$match_type == "same_class"] <= 50
))

# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------

saveRDS(matched, here("data", "processed", "injuries_matched.rds"))

# The audit file documents both successful and unsuccessful linkage attempts.
saveRDS(linkage_audit, here("data", "processed", "injury_road_linkage_audit.rds"))

