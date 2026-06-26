# ============================================================
# Post-Matching Road x Quarter Panel - POOLED MATCHING
# ============================================================
#
# Control roads are row-duplicated according to how many times
# their OA was matched with replacement. Each copy gets a unique
# panel_id so panel models treat them as separate matched units.
# OA-level clustering absorbs within-OA correlation.
#
# IMPORTANT FIX:
#   The treatment role in the matched panel is defined by the
#   matching output, not by the original road-level CAZ flag.
#
#   treat_indicator == 1 means treated matched row
#   treat_indicator == 0 means matched control row
#
#   Therefore treat_group is set from treat_indicator. This prevents
#   matched control rows whose source road was originally CAZ-treated
#   from being reclassified as treated in the final model panel.
#
# Inputs:
#   road_panel_dataset/           (parquet)
#   road_attributes_OA.gpkg
#   OA_matched_full_pooled.rds
#   OA_matching_pairs_pooled.rds
#   roads_caz_props.rds
#
# Output:
#   road_panel_matched_pooled.parquet
#
# Notes:
#   - The trajectory-shape matching variables are used upstream during OA
#     matching. They do not need to be carried into this road x quarter panel.
#   - CAZ road proportions are joined by both identifier and scheme to avoid
#     accidental row duplication for roads appearing in multiple scheme records.
#
# ============================================================

library(tidyverse)
library(arrow)
library(sf)
library(zoo)
library(here)


road_panel <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>%
  dplyr::select(
    identifier, quarter_year,
    treated_any, treated_50pct,
    treated_group_any, treated_group_50pct,
    scheme, caz_start_q,
    control_group1, control_group2, control_group3_mixed,
    starts_with("total_inj_adj_")
  ) %>%
  collect()

cat("Road panel loaded:", nrow(road_panel), "rows |",
    n_distinct(road_panel$identifier), "roads |",
    n_distinct(road_panel$quarter_year), "quarters\n")

road_attrs <- st_read(
  here("data", "processed", "road_attributes_OA.gpkg"),
  quiet = TRUE
) %>%
  st_drop_geometry() %>%
  dplyr::select(identifier, OA, road_class, length)

matched_oas <- readRDS(
  here("data", "processed", "OA_matched_full_pooled.rds")
)

matching_pairs <- readRDS(
  here("data", "processed", "OA_matching_pairs_pooled.rds")
)

# Treated OA metadata from matching output.
treated_oas <- matched_oas %>%
  filter(treat_indicator == 1) %>%
  distinct(OA, scheme, .keep_all = TRUE) %>%
  dplyr::select(OA, treat_indicator, scheme, baseline_injury_stratum)

# Control OA -> scheme mapping with match counts for row expansion.
ctrl_match_counts <- matching_pairs %>%
  count(control_OA, scheme, name = "match_count")

# Extra guard: remove any same-scheme treated OA from the control expansion.
treated_scheme_oas <- treated_oas %>%
  distinct(scheme, OA)

ctrl_match_counts <- ctrl_match_counts %>%
  rename(OA = control_OA) %>%
  anti_join(treated_scheme_oas, by = c("scheme", "OA"))

cat("Treated OAs:", nrow(treated_oas), "\n")
cat("Matching pairs after same-scheme treated-OA exclusion:",
    sum(ctrl_match_counts$match_count), "\n")
cat("Unique control OAs:", n_distinct(ctrl_match_counts$OA), "\n")

# CAZ proportions (England only).
scottish_schemes <- c("Aberdeen", "Dundee", "Edinburgh", "Glasgow")
road_caz_props <- readRDS(
  here("data", "processed", "roads_caz_props.rds")
) %>%
  filter(!scheme %in% scottish_schemes) %>%
  dplyr::select(identifier, scheme, prop_in_caz, caz_start_q) %>%
  mutate(caz_start_q = as.yearqtr(caz_start_q)) %>%
  distinct(identifier, scheme, .keep_all = TRUE)

# Scheme timing is needed for treated rows after the matched role is assigned.
scheme_timing <- road_caz_props %>%
  filter(!is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

# ============================================================
# Attach OA to road panel
# ============================================================

road_panel <- road_panel %>%
  left_join(
    road_attrs %>% dplyr::select(identifier, OA, road_class, length),
    by = "identifier"
  )

n_missing_oa <- sum(is.na(road_panel$OA))
cat("Roads missing OA assignment:", n_missing_oa,
    sprintf("(%.2f%%)\n", 100 * n_missing_oa / nrow(road_panel)))

# ============================================================
# Build panel: treated roads + expanded control roads
# ============================================================

all_matched_oas <- unique(c(treated_oas$OA, ctrl_match_counts$OA))

road_panel_matched_base <- road_panel %>%
  filter(OA %in% all_matched_oas)

cat("\nRoads in matched OAs:",
    n_distinct(road_panel_matched_base$identifier), "\n")

# --- Treated roads: one entry per road x quarter x treated scheme ---
# Use the scheme from the matching output, not any raw road-panel scheme value.
panel_treated <- road_panel_matched_base %>%
  dplyr::select(-scheme, -caz_start_q) %>%
  inner_join(
    treated_oas %>%
      dplyr::select(OA, scheme, baseline_injury_stratum),
    by = "OA",
    relationship = "many-to-many"
  ) %>%
  left_join(scheme_timing, by = "scheme") %>%
  mutate(
    treat_indicator = 1L,
    match_copy = 1L
  )

# --- Control roads: expand by match_count per matched scheme ---
# Use the scheme from ctrl_match_counts, not the raw road-panel scheme.
panel_controls <- road_panel_matched_base %>%
  dplyr::select(-scheme, -caz_start_q) %>%
  inner_join(
    ctrl_match_counts,
    by = "OA",
    relationship = "many-to-many"
  ) %>%
  left_join(scheme_timing, by = "scheme") %>%
  mutate(
    treat_indicator = 0L,
    baseline_injury_stratum = factor(NA, levels = levels(treated_oas$baseline_injury_stratum))
  ) %>%
  uncount(match_count) %>%
  group_by(identifier, quarter_year, scheme) %>%
  mutate(match_copy = row_number()) %>%
  ungroup()

# --- Combine ---
road_panel_matched <- bind_rows(panel_treated, panel_controls)

# Guard against same-scheme OA contamination after road expansion.
treated_scheme_oas_panel <- road_panel_matched %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme, OA)

road_panel_matched <- bind_rows(
  road_panel_matched %>% filter(treat_indicator == 1),
  road_panel_matched %>%
    filter(treat_indicator == 0) %>%
    anti_join(treated_scheme_oas_panel, by = c("scheme", "OA"))
)

# Unique panel identifier. Include scheme so a road reused across schemes remains
# separate, and include match_copy so replacement copies remain separate.
road_panel_matched <- road_panel_matched %>%
  mutate(panel_id = paste0(identifier, "_", scheme, "_", match_copy))

cat("\nPanel constructed:\n")
cat("  Treated rows:", sum(road_panel_matched$treat_indicator == 1), "\n")
cat("  Control rows:", sum(road_panel_matched$treat_indicator == 0), "\n")
cat("  Unique roads:", n_distinct(road_panel_matched$identifier), "\n")
cat("  Unique panel_ids:", n_distinct(road_panel_matched$panel_id), "\n")

# ============================================================
# Treatment indicators
# ============================================================

road_panel_matched <- road_panel_matched %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q = as.yearqtr(caz_start_q),
    
    # Critical: matched role, not raw road CAZ flag.
    treat_group = as.integer(treat_indicator == 1),
    
    post = case_when(
      treat_group == 1 & !is.na(caz_start_q) &
        quarter_year >= caz_start_q ~ 1L,
      TRUE ~ 0L
    ),
    rel_time = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(quarter_year - caz_start_q),
      TRUE ~ NA_real_
    )
  )

# ============================================================
# Attach CAZ proportions
# ============================================================

n_rows_before_caz_join <- nrow(road_panel_matched)

road_panel_matched <- road_panel_matched %>%
  left_join(
    road_caz_props %>% dplyr::select(identifier, scheme, prop_in_caz),
    by = c("identifier", "scheme")
  ) %>%
  mutate(prop_in_caz = replace_na(prop_in_caz, 0))

stopifnot(
  "CAZ proportion join changed row count" =
    nrow(road_panel_matched) == n_rows_before_caz_join
)

# ============================================================
# Final column selection
# ============================================================

road_panel_matched <- road_panel_matched %>%
  dplyr::select(
    panel_id,
    identifier,
    OA,
    quarter_year,
    match_copy,
    
    treat_group,
    post,
    rel_time,
    scheme,
    caz_start_q,
    prop_in_caz,
    
    treat_indicator,
    baseline_injury_stratum,
    
    road_class,
    length,
    
    starts_with("total_inj_adj_")
  ) %>%
  mutate(
    total_inj_adj_Vehicle = total_inj_adj_Car.Van + total_inj_adj_Other,
    total_inj_adj_ActiveTravel = total_inj_adj_Cyclist + total_inj_adj_Pedestrian
  )

# ============================================================
# Diagnostics
# ============================================================

cat("\n=== POST-MATCHING PANEL DIAGNOSTICS (POOLED) ===\n")

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
    treated = sum(treat_group == 1),
    controls = sum(treat_group == 0),
    .groups = "drop"
  ) %>%
  print()

cat("\nSame-scheme OA contamination check:\n")
contamination_check <- road_panel_matched %>%
  group_by(scheme, OA) %>%
  summarise(
    appears_treated = any(treat_group == 1),
    appears_control = any(treat_group == 0),
    .groups = "drop"
  ) %>%
  filter(appears_treated & appears_control)

print(contamination_check)
cat("  Contaminated scheme-OA pairs:", nrow(contamination_check), "\n")

cat("\nScheme x treatment timing:\n")
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
    n_obs = n(),
    mean_total = round(mean(total_inj_adj_All, na.rm = TRUE), 4),
    pct_zero_total = round(mean(total_inj_adj_All == 0) * 100, 1)
  ) %>%
  print()

cat("\nQuarter coverage:\n")
cat("  Min quarter:", as.character(min(road_panel_matched$quarter_year)), "\n")
cat("  Max quarter:", as.character(max(road_panel_matched$quarter_year)), "\n")
cat("  N quarters: ", n_distinct(road_panel_matched$quarter_year), "\n")

n_units <- n_distinct(road_panel_matched$panel_id)
n_qtrs <- n_distinct(road_panel_matched$quarter_year)
cat("\nPanel balance (panel_id x quarter):\n")
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
  here("data", "processed", "road_panel_matched_pooled.parquet")
)

cat("\nSaved: road_panel_matched_pooled.parquet\n")
cat("  Rows:", nrow(road_panel_matched),
    "| Panel units:", n_units,
    "| Quarters:", n_qtrs, "\n")



