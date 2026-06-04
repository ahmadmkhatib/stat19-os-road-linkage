# ============================================================
# Post-Matching Road × Quarter Panel Construction
# ============================================================
#
# Purpose:
#   Restricts the road-link panel to roads within matched OAs
#   (treated + matched controls) and attaches:
#     - OA-level matching weights
#     - treatment indicators and timing
#     - road and area characteristics
#
# Inputs:
#   road_panel_dataset/           (parquet, from panel construction script)
#   road_attributes_OA.gpkg       (road-level attributes with OA assignment)
#   OA_matched_full_mixed.rds     (matched OAs with weights, filtered to England)
#   OA_matching_pairs_mixed.rds   (treated → control OA pairs)
#   roads_caz_props.rds           (road-level CAZ proportions and timing)
#
# Output:
#   road_panel_matched.parquet    — road × quarter panel, matched sample only
#
# Key columns in output:
#   identifier         : Road link TOID
#   quarter_year       : Quarter of observation
#   OA                 : Output Area the road is assigned to
#   treat_group        : 1 = ever-treated road (≥50% in CAZ/LEZ)
#   post              : 1 if quarter >= caz_start_q and treat_group == 1
#   scheme             : CAZ/LEZ scheme name
#   caz_start_q        : Quarter CAZ/LEZ became operational
#   oa_weight          : Matching weight from design-stage MDM (treated = 1)
#   road_class         : A / B / Motorway / minor
#   length             : Road link length (metres)
#   KSI_adj_All        : Total adjusted KSI injuries (all modes)
#   Slight_adj_All     : Total adjusted slight injuries (all modes)
#   total_inj_adj_All  : Total adjusted injuries (all modes)
#   [mode-specific columns carried through from road_panel_dataset]
#
# ============================================================

library(tidyverse)
library(arrow)
library(sf)
library(zoo)
library(here)

# ============================================================
# Load inputs
# ============================================================

# Full unmatched road panel (parquet)
road_panel <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

cat("Road panel loaded:", nrow(road_panel), "rows |",
    n_distinct(road_panel$identifier), "roads |",
    n_distinct(road_panel$quarter_year), "quarters\n")

# Road-level attributes (OA assignment, road class, length)
road_attrs <- st_read(
  here("data", "processed", "road_attributes_OA.gpkg"),
  quiet = TRUE
) %>%
  st_drop_geometry() %>%
  select(identifier, OA, road_class, length)

# Matched OA dataset (England only, per-scheme matching from script 16)
# Controls already carry their scheme from per-scheme matching.
matched_oas <- readRDS(
  here("data", "processed", "OA_matched_full_mixed.rds")
) %>%
  select(OA, weights, treat_indicator, scheme,
         baseline_injury_stratum)

matching_pairs <- readRDS(
  here("data", "processed", "OA_matching_pairs_mixed.rds")
)

cat("Matched OAs:", nrow(matched_oas), "\n")
cat("  Treated:", sum(matched_oas$treat_indicator == 1),
    "| Controls:", sum(matched_oas$treat_indicator == 0), "\n")

# CAZ proportions and timing at road level (England only)
scottish_schemes <- c("Aberdeen", "Dundee", "Edinburgh", "Glasgow")

road_caz_props <- readRDS(
  here("data", "processed", "roads_caz_props.rds")
) %>%
  filter(!scheme %in% scottish_schemes) %>%
  select(identifier, scheme, prop_in_caz, caz_start_q) %>%
  mutate(caz_start_q = as.yearqtr(caz_start_q))

# ============================================================
# Attach OA to road panel
# ============================================================

road_panel <- road_panel %>%
  left_join(
    road_attrs %>% select(identifier, OA, road_class, length),
    by = "identifier"
  )

n_missing_oa <- sum(is.na(road_panel$OA))
cat("Roads missing OA assignment:", n_missing_oa,
    sprintf("(%.2f%%)\n", 100 * n_missing_oa / nrow(road_panel)))

# ============================================================
# Restrict to matched OAs only
# ============================================================

matched_oa_ids <- matched_oas %>% pull(OA)

road_panel_matched <- road_panel %>%
  filter(OA %in% matched_oa_ids)

cat("\nAfter restricting to matched OAs:\n")
cat("  Rows:  ", nrow(road_panel_matched), "\n")
cat("  Roads: ", n_distinct(road_panel_matched$identifier), "\n")

# ============================================================
# Attach matching weights and OA-level metadata
# ============================================================

# Per-scheme matching: control OAs may appear in multiple schemes
# (matched to different treated schemes). Include scheme in the join
# so the panel has one entry per road × quarter × scheme.
# Scheme-specific DiD subsets by scheme; pooled DiD should deduplicate.

road_panel_matched <- road_panel_matched %>%
  left_join(
    matched_oas %>%
      select(OA, oa_weight = weights, treat_indicator,
             baseline_injury_stratum, oa_scheme = scheme),
    by = "OA",
    relationship = "many-to-many"
  )

# For treated roads, keep their actual scheme (from road_caz_props).
# For control roads, use the scheme from matching.
road_panel_matched <- road_panel_matched %>%
  mutate(
    scheme = if_else(treat_indicator == 1, scheme, oa_scheme)
  ) %>%
  select(-oa_scheme)

# Check duplication from multi-scheme controls
n_panel_rows <- nrow(road_panel_matched)
n_unique_roadqtr <- road_panel_matched %>%
  distinct(identifier, quarter_year) %>%
  nrow()
cat("Panel rows:", n_panel_rows,
    "| Unique road × quarter:", n_unique_roadqtr, "\n")
if (n_panel_rows > n_unique_roadqtr) {
  cat("  Note:", n_panel_rows - n_unique_roadqtr,
      "duplicated rows from multi-scheme control OAs\n")
}

# Sanity: all rows should have a weight
stopifnot(
  "Missing OA weights after join" = !anyNA(road_panel_matched$oa_weight)
)

# ============================================================
# Rebuild clean treatment indicators
# ============================================================
# treat_group: road is inside a CAZ/LEZ (>=50% of length)
# post:        treated road AND quarter >= scheme start

road_panel_matched <- road_panel_matched %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q),
    
    # Treated group indicator (road-level, ever-treated)
    treat_group = as.integer(treated_group_50pct == 1),
    
    # Post-treatment indicator (road × quarter)
    post = case_when(
      treat_group == 1 & !is.na(caz_start_q) &
        quarter_year >= caz_start_q ~ 1L,
      TRUE ~ 0L
    ),
    
    # Relative time to treatment (quarters; NA for controls)
    rel_time = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(quarter_year - caz_start_q),
      TRUE ~ NA_real_
    )
  )

# ============================================================
# Attach CAZ proportions
# ============================================================
# Controls already carry scheme from per-scheme matching (script 16).

road_panel_matched <- road_panel_matched %>%
  left_join(
    road_caz_props %>%
      select(identifier, prop_in_caz),
    by = "identifier"
  ) %>%
  mutate(
    prop_in_caz = replace_na(prop_in_caz, 0)
  )

# ============================================================
# Final column selection and ordering
# ============================================================

road_panel_matched <- road_panel_matched %>%
  select(
    # Identifiers
    identifier,
    OA,
    quarter_year,
    
    # Treatment
    treat_group,
    post,
    rel_time,
    scheme,
    caz_start_q,
    prop_in_caz,
    
    # Matching metadata
    oa_weight,
    treat_indicator,
    baseline_injury_stratum,
    
    # Road characteristics
    road_class,
    length,
    
    # Outcomes — adjusted (primary)
    KSI_adj_All,
    Slight_adj_All,
    total_inj_adj_All,
    
    # Outcomes — by mode (adjusted)
    starts_with("KSI_adj_"),
    starts_with("Slight_adj_"),
    starts_with("total_inj_adj_")
  ) %>%
  rename_with(~ make.names(.x))


##############
road_panel_matched %>%
  distinct(identifier, treat_group, scheme) %>%
  count(treat_group)



# ============================================================
# Diagnostics
# ============================================================

cat("\n=== POST-MATCHING PANEL DIAGNOSTICS ===\n")

cat("\nTreatment group breakdown (road level):\n")
road_panel_matched %>%
  distinct(identifier, treat_group, scheme) %>%
  count(treat_group) %>%
  print()

cat("\nScheme × treatment timing:\n")
road_panel_matched %>%
  filter(treat_group == 1) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q) %>%
  print()

cat("\nWeight distribution by treatment status:\n")
road_panel_matched %>%
  distinct(OA, oa_weight, treat_indicator) %>%
  group_by(treat_indicator) %>%
  summarise(
    n_OAs      = n(),
    mean_wt    = round(mean(oa_weight), 3),
    median_wt  = round(median(oa_weight), 3),
    max_wt     = round(max(oa_weight), 3),
    .groups    = "drop"
  ) %>%
  print()

cat("\nOutcome summary (treated roads, post period):\n")
road_panel_matched %>%
  filter(treat_group == 1, post == 1) %>%
  summarise(
    n_obs           = n(),
    mean_KSI        = round(mean(KSI_adj_All,       na.rm = TRUE), 4),
    mean_slight     = round(mean(Slight_adj_All,     na.rm = TRUE), 4),
    mean_total      = round(mean(total_inj_adj_All,  na.rm = TRUE), 4),
    pct_zero_total  = round(mean(total_inj_adj_All == 0) * 100, 1)
  ) %>%
  print()

cat("\nQuarter coverage:\n")
cat("  Min quarter:", as.character(min(road_panel_matched$quarter_year)), "\n")
cat("  Max quarter:", as.character(max(road_panel_matched$quarter_year)), "\n")
cat("  N quarters: ", n_distinct(road_panel_matched$quarter_year), "\n")

# Panel balance check
n_roads <- n_distinct(road_panel_matched$identifier)
n_qtrs  <- n_distinct(road_panel_matched$quarter_year)
cat("\nPanel balance:\n")
cat("  Roads:    ", n_roads, "\n")
cat("  Quarters: ", n_qtrs, "\n")
cat("  Rows:     ", nrow(road_panel_matched), "\n")
cat("  Expected: ", n_roads * n_qtrs, "\n")
cat("  Balanced: ",
    nrow(road_panel_matched) == n_roads * n_qtrs, "\n")

# ============================================================
# Save
# ============================================================

arrow::write_parquet(
  road_panel_matched,
  here("data", "processed", "road_panel_matched.parquet")
)


road_panel_matched<- arrow::read_parquet(
  here("data", "processed", "road_panel_matched.parquet")
)
cat("Rows:", nrow(road_panel_matched),
    "| Roads:", n_distinct(road_panel_matched$identifier),
    "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n")



glimpse(road_panel_matched)
        
