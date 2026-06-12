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
matched_oas <- readRDS(
  here("data", "processed", "OA_matched_full_mixed.rds")
)

matching_pairs <- readRDS(
  here("data", "processed", "OA_matching_pairs_mixed.rds")
)

# Treated OA metadata
treated_oas <- matched_oas %>%
  filter(treat_indicator == 1) %>%
  select(OA, treat_indicator, scheme, baseline_injury_stratum)

# Control OA → scheme mapping with match counts (for row expansion).
# A control matched N times within a scheme gets N rows in the panel.
# This replaces frequency weights with explicit row duplication.
ctrl_match_counts <- matching_pairs %>%
  count(control_OA, scheme, name = "match_count")

cat("Treated OAs:", nrow(treated_oas), "\n")
cat("Matching pairs:", nrow(matching_pairs), "\n")
cat("Unique control OAs:", n_distinct(ctrl_match_counts$control_OA), "\n")

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
# Build panel: treated roads + expanded control roads
# ============================================================
# Instead of frequency weights, control roads are DUPLICATED
# proportional to how many times their OA was matched (with
# replacement). Each copy gets a unique panel_id so att_gt()
# treats them as separate units. OA-level clustering absorbs
# the within-OA correlation from duplication.

# --- Treated roads: one entry per road × quarter ---
all_matched_oas <- unique(c(treated_oas$OA, ctrl_match_counts$control_OA))

road_panel_matched <- road_panel %>%
  filter(OA %in% all_matched_oas)

cat("\nRoads in matched OAs:", n_distinct(road_panel_matched$identifier), "\n")

panel_treated <- road_panel_matched %>%
  semi_join(treated_oas, by = "OA") %>%
  mutate(treat_indicator = 1L, match_copy = 1L)

# --- Control roads: expand by match_count per scheme ---
panel_controls <- road_panel_matched %>%
  inner_join(
    ctrl_match_counts %>% rename(OA = control_OA),
    by = "OA",
    relationship = "many-to-many"
  ) %>%
  mutate(treat_indicator = 0L) %>%
  # Expand: repeat each road × quarter × scheme row match_count times
  uncount(match_count) %>%
  group_by(identifier, quarter_year, scheme) %>%
  mutate(match_copy = row_number()) %>%
  ungroup()

# --- Combine ---
road_panel_matched <- bind_rows(panel_treated, panel_controls)

# Unique panel identifier (road × copy, for att_gt idname)
road_panel_matched <- road_panel_matched %>%
  mutate(panel_id = paste0(identifier, "_", match_copy))

# Attach treated OA metadata (scheme, baseline_injury_stratum)
road_panel_matched <- road_panel_matched %>%
  left_join(
    treated_oas %>% select(OA, baseline_injury_stratum),
    by = "OA"
  )

cat("\nPanel constructed:\n")
cat("  Treated rows:", sum(road_panel_matched$treat_indicator == 1), "\n")
cat("  Control rows:", sum(road_panel_matched$treat_indicator == 0), "\n")
cat("  Unique roads:", n_distinct(road_panel_matched$identifier), "\n")
cat("  Unique panel_ids:", n_distinct(road_panel_matched$panel_id), "\n")

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
    panel_id,
    identifier,
    OA,
    quarter_year,
    match_copy,
    
    # Treatment
    treat_group,
    post,
    rel_time,
    scheme,
    caz_start_q,
    prop_in_caz,
    
    # Matching metadata
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



# ============================================================
# Diagnostics
# ============================================================

cat("\n=== POST-MATCHING PANEL DIAGNOSTICS ===\n")

cat("\nTreatment group breakdown (panel_id level):\n")
road_panel_matched %>%
  distinct(panel_id, treat_group, scheme) %>%
  count(treat_group) %>%
  print()

cat("\nPer-scheme breakdown:\n")
road_panel_matched %>%
  distinct(panel_id, treat_group, scheme) %>%
  group_by(scheme) %>%
  summarise(
    treated  = sum(treat_group == 1),
    controls = sum(treat_group == 0),
    .groups  = "drop"
  ) %>%
  print()

cat("\nScheme × treatment timing:\n")
road_panel_matched %>%
  filter(treat_group == 1) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q) %>%
  print()

cat("\nRow expansion from replacement matching:\n")
cat("  Unique roads:", n_distinct(road_panel_matched$identifier), "\n")
cat("  Unique panel_ids:", n_distinct(road_panel_matched$panel_id), "\n")
cat("  Expansion factor:",
    round(n_distinct(road_panel_matched$panel_id) /
            n_distinct(road_panel_matched$identifier), 2), "\n")

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

# Panel balance check (on panel_id, not identifier)
n_units <- n_distinct(road_panel_matched$panel_id)
n_qtrs  <- n_distinct(road_panel_matched$quarter_year)
cat("\nPanel balance (panel_id × quarter):\n")
cat("  Panel units:", n_units, "\n")
cat("  Quarters:   ", n_qtrs, "\n")
cat("  Rows:       ", nrow(road_panel_matched), "\n")
cat("  Expected:   ", n_units * n_qtrs, "\n")
cat("  Balanced:   ",
    nrow(road_panel_matched) == n_units * n_qtrs, "\n")

# ============================================================
# Save
# ============================================================

arrow::write_parquet(
  road_panel_matched,
  here("data", "processed", "road_panel_matched.parquet")
)

cat("\nSaved: road_panel_matched.parquet\n")
cat("  Rows:", nrow(road_panel_matched),
    "| Panel units:", n_units,
    "| Quarters:", n_qtrs, "\n")
