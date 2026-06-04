# ============================================================
# Post-Matching Road × Quarter Panel — POOLED MATCHING
# ============================================================
#
# Same structure as main script 20, but reads from pooled
# matching outputs (total injuries matching only).
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
# ============================================================

library(tidyverse)
library(arrow)
library(sf)
library(zoo)
library(here)

# ============================================================
# Load inputs
# ============================================================

road_panel <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

cat("Road panel loaded:", nrow(road_panel), "rows |",
    n_distinct(road_panel$identifier), "roads |",
    n_distinct(road_panel$quarter_year), "quarters\n")

road_attrs <- st_read(
  here("data", "processed", "road_attributes_OA.gpkg"),
  quiet = TRUE
) %>%
  st_drop_geometry() %>%
  select(identifier, OA, road_class, length)

# Pooled matching outputs (England only, no Scotland filtering needed)
matched_oas <- readRDS(
  here("data", "processed", "OA_matched_full_pooled.rds")
) %>%
  select(OA, weights, treat_indicator, scheme,
         baseline_injury_stratum)

matching_pairs <- readRDS(
  here("data", "processed", "OA_matching_pairs_pooled.rds")
)

cat("Matched OAs:", nrow(matched_oas), "\n")
cat("  Treated:", sum(matched_oas$treat_indicator == 1),
    "| Controls:", sum(matched_oas$treat_indicator == 0), "\n")

# Map each control OA to its paired treated OA's scheme
treated_scheme <- matched_oas %>%
  filter(treat_indicator == 1) %>%
  select(OA, scheme)

ctrl_scheme_lookup <- matching_pairs %>%
  left_join(treated_scheme, by = c("treated_OA" = "OA")) %>%
  select(OA = control_OA, scheme_from_treated = scheme) %>%
  distinct()

cat("Control OA → scheme lookup:", nrow(ctrl_scheme_lookup), "entries\n")

# CAZ proportions (England only)
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
# Restrict to matched OAs
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

road_panel_matched <- road_panel_matched %>%
  left_join(
    matched_oas %>%
      select(OA, oa_weight = weights, treat_indicator,
             baseline_injury_stratum),
    by = "OA"
  )

stopifnot(
  "Missing OA weights after join" = !anyNA(road_panel_matched$oa_weight)
)

# ============================================================
# Treatment indicators
# ============================================================

road_panel_matched <- road_panel_matched %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q),
    treat_group = as.integer(treated_group_50pct == 1),
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
# Attach CAZ proportions and scheme for control roads
# ============================================================

road_panel_matched <- road_panel_matched %>%
  left_join(
    road_caz_props %>% select(identifier, prop_in_caz),
    by = "identifier"
  ) %>%
  mutate(prop_in_caz = replace_na(prop_in_caz, 0)) %>%
  left_join(ctrl_scheme_lookup, by = "OA") %>%
  mutate(
    scheme = if_else(
      treat_group == 0 & is.na(scheme),
      scheme_from_treated,
      scheme
    )
  ) %>%
  select(-scheme_from_treated)

# ============================================================
# Final column selection
# ============================================================

road_panel_matched <- road_panel_matched %>%
  select(
    identifier, OA, quarter_year,
    treat_group, post, rel_time,
    scheme, caz_start_q, prop_in_caz,
    oa_weight, treat_indicator, baseline_injury_stratum,
    road_class, length,
    KSI_adj_All, Slight_adj_All, total_inj_adj_All,
    starts_with("KSI_adj_"), starts_with("Slight_adj_"),
    starts_with("total_inj_adj_")
  ) %>%
  rename_with(~ make.names(.x))

# ============================================================
# Diagnostics
# ============================================================

cat("\n=== POST-MATCHING PANEL DIAGNOSTICS (POOLED) ===\n")

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
    n_OAs     = n(),
    mean_wt   = round(mean(oa_weight), 3),
    median_wt = round(median(oa_weight), 3),
    max_wt    = round(max(oa_weight), 3),
    .groups   = "drop"
  ) %>%
  print()

n_roads <- n_distinct(road_panel_matched$identifier)
n_qtrs  <- n_distinct(road_panel_matched$quarter_year)
cat("\nPanel balance:\n")
cat("  Roads:    ", n_roads, "\n")
cat("  Quarters: ", n_qtrs, "\n")
cat("  Rows:     ", nrow(road_panel_matched), "\n")
cat("  Expected: ", n_roads * n_qtrs, "\n")
cat("  Balanced: ", nrow(road_panel_matched) == n_roads * n_qtrs, "\n")

# ============================================================
# Save
# ============================================================

arrow::write_parquet(
  road_panel_matched,
  here("data", "processed", "road_panel_matched_pooled.parquet")
)

cat("\nSaved: road_panel_matched_pooled.parquet\n")
cat("Rows:", nrow(road_panel_matched),
    "| Roads:", n_roads, "| Quarters:", n_qtrs, "\n")
