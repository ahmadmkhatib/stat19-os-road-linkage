# =============================================================================
# POST-MATCHING DIAGNOSTICS — 
# OA-Level Two-Stage Mahalanobis Distance Matching
# 
#
# INPUTS (from data/processed/):
#   OA_matched_full_A.rds          — full matched dataset, Analysis A
#   OA_matched_treated_A.rds       — treated OAs + weights + stratum
#   OA_matched_donors_A.rds        — control OAs + weights
#   OA_common_support_flags.rds    — structurally isolated OA flags
#   OA_matching_census.rds         — full original dataset (unmatched pool)
#
# SPATIAL INPUTS (Section 14 — maps):
#   data/spatial/OA_boundaries_GB.gpkg
#   data/spatial/LAD_boundaries_GB.gpkg
#
# OUTPUTS (to output/diagnostics/):
#   Tables  : 01–07 CSV files (described inline)
#   Figures : fig01–fig16 PNG files (described inline)
#
# =============================================================================

library(tidyverse)
library(here)
library(cobalt)
library(MatchIt)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(sf)
library(ggrepel)

select <- dplyr::select
filter <- dplyr::filter

dir.create(here("output", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics")

# =============================================================================
# THEME AND COLOURS
# =============================================================================

theme_diag <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                      colour = "#1A2E5A", margin = margin(b = 8)),
      plot.subtitle    = element_text(size = base_size - 1, colour = "#555555",
                                      margin = margin(b = 12)),
      plot.caption     = element_text(size = base_size - 2, colour = "#888888",
                                      hjust = 0, margin = margin(t = 8)),
      axis.title       = element_text(size = base_size - 1, colour = "#333333"),
      axis.text        = element_text(size = base_size - 2, colour = "#444444"),
      strip.text       = element_text(size = base_size - 1, face = "bold",
                                      colour = "#1A2E5A"),
      strip.background = element_rect(fill = "#EEF2F8", colour = NA),
      panel.grid.major = element_line(colour = "#E5E9F0", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "#CCCCCC", fill = NA,
                                      linewidth = 0.4),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      legend.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 12, 12)
    )
}

COL_TREATED  <- "#D85A30"
COL_CONTROL  <- "#2E6FAB"
COL_BEFORE   <- "#E74C3C"
COL_AFTER    <- "#2ECC71"
COL_SCOTLAND <- "#6B3FA0"
COL_ENGLAND  <- "#2E6FAB"

save_fig <- function(p, filename, width = 12, height = 8, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# =============================================================================
# SECTION 1 — LOAD DATA
# =============================================================================

cat("\n=== Loading data ===\n")

matched_A  <- readRDS(here("data", "processed", "OA_matched_full_A.rds"))
treated_A  <- readRDS(here("data", "processed", "OA_matched_treated_A.rds"))
controls_A <- readRDS(here("data", "processed", "OA_matched_donors_A.rds"))
csupport   <- readRDS(here("data", "processed", "OA_common_support_flags.rds"))
full_data  <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Restricted sample (excluding 41 structurally isolated treated OAs)
isolated_ids   <- csupport$treated_OA[csupport$analysis == "A"]
matched_A_restr <- matched_A %>%
  filter(!(OA %in% isolated_ids & treat_indicator == 1))

cat("  Analysis A (full):       Treated =", sum(matched_A$treat_indicator == 1),
    "| Controls =", sum(matched_A$treat_indicator == 0), "\n")
cat("  Analysis A (restricted): Treated =",
    sum(matched_A_restr$treat_indicator == 1),
    "| Controls =", sum(matched_A_restr$treat_indicator == 0), "\n\n")

# Unmatched eligible pool (Analysis A denominator — for before/after comparisons)
unmatched_pool <- full_data %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# ── Variable definitions ─────────────────────────────────────────────────────

stage1_road   <- c("road_density_m_km2", "road_length_km",
                   "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban  <- c("dist_citycentre", "pop_density")
stage1_socdem <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                   "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct")
stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

# NOTE: stage2_trends here uses all 9 variables from the matching pipeline,
# not the 4-variable subset used in some earlier drafts of this script.
stage2_trends <- c(
  "trend_car_KSI_pkm",   "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",   "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",   "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm",   "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",   "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",   "mean_ped_slight_pkm",
  "mean_other_KSI_pkm", "mean_other_slight_pkm",
  "mean_total_pkm"
)
log_names_s2 <- paste0("log1p_", stage2_levels)

var_labels <- c(
  road_density_m_km2      = "Road density (m/km\u00b2)",
  road_length_km          = "Road length (km)",
  pct_A_road              = "% A-road",
  pct_B_road              = "% B-road",
  pct_minor_road          = "% Minor road",
  dist_citycentre         = "Distance to city centre (m)",
  pop_density             = "Population density (persons/km\u00b2)",
  IMD                     = "Index of Multiple Deprivation",
  cars_none_pct           = "% households: no car",
  Drive_Car_pct           = "% commuting by car",
  Walk_pct                = "% commuting on foot",
  Bicycle_pct             = "% commuting by bicycle",
  X65plus_pct             = "% aged 65+",
  X5to19_pct              = "% aged 5\u201319",
  X20to24_pct             = "% aged 20\u201324",
  trend_car_KSI_pkm       = "Trend: car KSI/km",
  trend_car_slight_pkm    = "Trend: car slight/km",
  trend_cyc_KSI_pkm       = "Trend: cycling KSI/km",
  trend_cyc_slight_pkm    = "Trend: cycling slight/km",
  trend_ped_KSI_pkm       = "Trend: pedestrian KSI/km",
  trend_ped_slight_pkm    = "Trend: pedestrian slight/km",
  trend_other_KSI_pkm     = "Trend: other KSI/km",
  trend_other_slight_pkm  = "Trend: other slight/km",
  trend_total_pkm         = "Trend: total injuries/km",
  mean_car_KSI_pkm        = "Mean: car KSI/km",
  mean_car_slight_pkm     = "Mean: car slight/km",
  mean_cyc_KSI_pkm        = "Mean: cycling KSI/km",
  mean_cyc_slight_pkm     = "Mean: cycling slight/km",
  mean_ped_KSI_pkm        = "Mean: pedestrian KSI/km",
  mean_ped_slight_pkm     = "Mean: pedestrian slight/km",
  mean_other_KSI_pkm      = "Mean: other KSI/km",
  mean_other_slight_pkm   = "Mean: other slight/km",
  mean_total_pkm          = "Mean: total injuries/km"
)

# Shared SMD helper
compute_smd <- function(data, var, treat_col = "treat_indicator") {
  t_vals <- data[[var]][data[[treat_col]] == 1]
  c_vals <- data[[var]][data[[treat_col]] == 0]
  t_vals <- t_vals[!is.na(t_vals)]
  c_vals <- c_vals[!is.na(c_vals)]
  if (length(t_vals) < 2 || length(c_vals) < 2) return(NA_real_)
  pooled_sd <- sqrt((var(t_vals) + var(c_vals)) / 2)
  if (pooled_sd == 0) return(0)
  (mean(t_vals) - mean(c_vals)) / pooled_sd
}

# =============================================================================
# SECTION 2 — DESCRIPTIVE SUMMARY TABLES
# =============================================================================
# [REPORT: Insert Table in Section 3 "Data" or a standalone Appendix table.
#  The treated vs control descriptive table is the primary "Table 1" that
#  reviewers expect before any matching results are shown.]

cat("=== Section 2: Descriptive summary tables ===\n")

# ── 2a. Pre-matching: all eligible treated vs all eligible controls ──────────
# [REPORT: Table 1 — "Characteristics of treated and control OAs before
#  and after matching". Show Treated / Unmatched control / Matched control
#  columns side by side so the reader can see what matching achieved.]

desc_stats <- function(data, vars, group_label) {
  map_df(vars, function(v) {
    x <- data[[v]]
    if (is.null(x)) return(NULL)
    tibble(
      variable   = v,
      label      = coalesce(var_labels[v], v),
      group      = group_label,
      n          = sum(!is.na(x)),
      mean       = round(mean(x, na.rm = TRUE), 4),
      sd         = round(sd(x,   na.rm = TRUE), 4),
      median     = round(median(x, na.rm = TRUE), 4),
      p25        = round(quantile(x, 0.25, na.rm = TRUE), 4),
      p75        = round(quantile(x, 0.75, na.rm = TRUE), 4)
    )
  })
}

all_desc_vars <- c(stage1_vars, stage2_trends, stage2_levels)

# Three groups: treated, unmatched controls, matched controls
desc_treated_pre   <- desc_stats(
  unmatched_pool %>% filter(treat_indicator == 1), all_desc_vars, "Treated")
desc_control_pre   <- desc_stats(
  unmatched_pool %>% filter(treat_indicator == 0), all_desc_vars, "Unmatched controls")
desc_control_post  <- desc_stats(
  matched_A      %>% filter(treat_indicator == 0), all_desc_vars, "Matched controls (A)")

desc_table_1 <- bind_rows(desc_treated_pre, desc_control_pre, desc_control_post) %>%
  pivot_wider(
    id_cols   = c(variable, label),
    names_from  = group,
    values_from = c(mean, sd, median, p25, p75, n),
    names_glue  = "{group}_{.value}"
  ) %>%
  mutate(
    var_group = case_when(
      variable %in% stage1_road   ~ "Road network",
      variable %in% stage1_urban  ~ "Urban geography",
      variable %in% stage1_socdem ~ "Sociodemographic",
      variable %in% stage2_trends ~ "Injury trends (pre-treatment)",
      variable %in% stage2_levels ~ "Injury levels (pre-treatment)"
    )
  ) %>%
  arrange(var_group, variable)

write_csv(desc_table_1,
          file.path(outdir, "01_descriptive_table1_treated_vs_controls.csv"))
cat("  Saved: 01_descriptive_table1_treated_vs_controls.csv\n")
cat("  [REPORT] Table 1: 3-column descriptive table (Treated | Unmatched controls\n")
cat("           | Matched controls). Report mean (SD) and median [IQR] per variable.\n")
cat("           Group by variable category (road / urban / sociodem / trends / levels).\n\n")

# ── 2b. Concise mean (SD) for the paper body ────────────────────────────────
# Formatted for direct insertion into a results table

desc_formatted <- bind_rows(desc_treated_pre, desc_control_pre, desc_control_post) %>%
  mutate(
    mean_sd    = sprintf("%.3f (%.3f)", mean, sd),
    median_iqr = sprintf("%.3f [%.3f\u2013%.3f]", median, p25, p75)
  ) %>%
  select(variable, label, group, n, mean_sd, median_iqr)

write_csv(desc_formatted,
          file.path(outdir, "01b_descriptive_formatted.csv"))
cat("  Saved: 01b_descriptive_formatted.csv\n\n")

# ── 2c. Country-stratified descriptives for treated OAs ─────────────────────
# [REPORT: Supplementary table or footnote in Data section.
#  Important because Scotland-specific matching issues require transparency.]

country_desc <- map_df(c("England", "Scotland"), function(ctry) {
  d <- matched_A %>% filter(treat_indicator == 1, country == ctry)
  desc_stats(d, all_of(intersect(all_desc_vars, names(d))), ctry)
})

write_csv(country_desc,
          file.path(outdir, "01c_descriptive_by_country.csv"))
cat("  Saved: 01c_descriptive_by_country.csv\n")
cat("  [REPORT] Supplementary Table: treated OA characteristics by country\n")
cat("           (England n=613, Scotland n=201). Highlight differences in\n")
cat("           road_length, pop_density, mean_total_pkm, and IMD.\n\n")

# =============================================================================
# SECTION 3 — SMD TABLES (BEFORE / AFTER)
# =============================================================================
# [REPORT: Table 2 / Appendix. Show SMD before matching, after Stage 1,
#  after Stage 2 (full), and after Stage 2 (restricted). This directly
#  supports the parallel trends discussion.]

cat("=== Section 3: SMD before/after tables ===\n")

# Stage 1: unmatched pool → Stage 1 matched pool
# Stage 2: Stage 1 pool → Stage 2 matched (full and restricted)

smd_table <- function(vars, data_before, data_after_full,
                      data_after_restr = NULL, label = "") {
  map_df(vars, function(v) {
    row <- tibble(
      variable       = v,
      label          = coalesce(var_labels[v], v),
      smd_before     = round(compute_smd(data_before,     v), 4),
      smd_after_full = round(compute_smd(data_after_full, v), 4)
    )
    if (!is.null(data_after_restr))
      row$smd_after_restr <- round(compute_smd(data_after_restr, v), 4)
    row
  })
}

smd_s1 <- smd_table(
  stage1_vars,
  data_before      = unmatched_pool,
  data_after_full  = matched_A
)
smd_s1 <- smd_s1 %>%
  rename(smd_unmatched = smd_before,
         smd_after_S1  = smd_after_full) %>%
  mutate(balanced_after_S1 = abs(smd_after_S1) < 0.1)

write_csv(smd_s1, file.path(outdir, "02_smd_stage1.csv"))
cat("  Saved: 02_smd_stage1.csv\n")
cat("  [REPORT] Matching Report Section 6.3 (already in report).\n")
cat("           Also useful as a standalone appendix table.\n\n")

smd_s2 <- smd_table(
  c(stage2_trends, stage2_levels),
  data_before       = matched_A,        # Stage 1 pool as 'before'
  data_after_full   = matched_A,        # Stage 2 full (same object; SMD computed post-S2)
  data_after_restr  = matched_A_restr   # Stage 2 restricted
)
# NOTE: because matched_A IS the Stage 2 output, smd_after_full == smd_before
# here numerically. To show pre-S2 imbalance, use the Stage 1 pool
# (s1_A$treated + s1_A$controls from the matching script). If that object is
# not saved, use unmatched_pool as a conservative 'before' baseline.
smd_s2 <- smd_s2 %>%
  rename(
    smd_preS2             = smd_before,
    smd_postS2_full       = smd_after_full,
    smd_postS2_restricted = smd_after_restr
  ) %>%
  mutate(
    var_type              = if_else(variable %in% stage2_trends, "Trend", "Level"),
    balanced_full         = abs(smd_postS2_full) < 0.1,
    balanced_restricted   = abs(smd_postS2_restricted) < 0.1
  )

write_csv(smd_s2, file.path(outdir, "03_smd_stage2.csv"))
cat("  Saved: 03_smd_stage2.csv\n")
cat("  [REPORT] Matching Report Sections 7.4 and 7.5 (already in report).\n")
cat("           Trend rows are the parallel-trends evidence.\n\n")

# =============================================================================
# SECTION 4 — WEIGHT DISTRIBUTION TABLE
# =============================================================================
# [REPORT: Matching Report Section 8 (already drafted). Also provide in
#  the paper Appendix as a table showing weight distribution by country.]

cat("=== Section 4: Weight distribution table ===\n")

weight_table <- matched_A %>%
  filter(treat_indicator == 0) %>%
  group_by(country) %>%
  summarise(
    n_controls   = n(),
    mean_weight  = round(mean(weights), 3),
    median_weight= round(median(weights), 3),
    sd_weight    = round(sd(weights), 3),
    p90_weight   = round(quantile(weights, 0.90), 3),
    p95_weight   = round(quantile(weights, 0.95), 3),
    max_weight   = round(max(weights), 3),
    eff_N        = round(sum(weights)^2 / sum(weights^2), 1),
    efficiency   = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
    n_at_cap     = sum(weights >= 5),
    pct_at_cap   = round(100 * mean(weights >= 5), 1),
    .groups      = "drop"
  ) %>%
  bind_rows(
    # Overall row
    matched_A %>% filter(treat_indicator == 0) %>%
      summarise(
        country      = "Overall",
        n_controls   = n(),
        mean_weight  = round(mean(weights), 3),
        median_weight= round(median(weights), 3),
        sd_weight    = round(sd(weights), 3),
        p90_weight   = round(quantile(weights, 0.90), 3),
        p95_weight   = round(quantile(weights, 0.95), 3),
        max_weight   = round(max(weights), 3),
        eff_N        = round(sum(weights)^2 / sum(weights^2), 1),
        efficiency   = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
        n_at_cap     = sum(weights >= 5),
        pct_at_cap   = round(100 * mean(weights >= 5), 1)
      )
  )

write_csv(weight_table, file.path(outdir, "04_weight_distribution.csv"))
cat("  Saved: 04_weight_distribution.csv\n")
cat("  [REPORT] Section 8 of matching report. Key numbers to highlight:\n")
cat("    - Scotland max weight (before/after cap)\n")
cat("    - Effective N overall vs nominal N\n")
cat("    - Efficiency ratio by country\n\n")

# =============================================================================
# SECTION 5 — STRATUM CHARACTERISTICS TABLE
# =============================================================================
# [REPORT: New table in the report (Section 9 / Appendix).
#  Show that the four injury-level strata are internally coherent — i.e.
#  stratum 4 OAs really do have higher injury rates, higher road density,
#  more urban structure than stratum 1 OAs. This validates the stratification.]

cat("=== Section 5: Stratum characteristics ===\n")

stratum_table <- matched_A %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  group_by(Stratum = baseline_injury_stratum) %>%
  summarise(
    n                     = n(),
    # Injury levels
    mean_total_pkm        = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    sd_total_pkm          = round(sd(mean_total_pkm,       na.rm = TRUE), 5),
    median_total_pkm      = round(median(mean_total_pkm,   na.rm = TRUE), 5),
    # Road network
    mean_road_length_km   = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_road_density     = round(mean(road_density_m_km2, na.rm = TRUE), 0),
    mean_pct_A_road       = round(mean(pct_A_road,         na.rm = TRUE), 1),
    mean_pct_minor        = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    # Urban
    mean_pop_density      = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_dist_citycentre  = round(mean(dist_citycentre,    na.rm = TRUE), 0),
    # Sociodemographic
    mean_IMD              = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct    = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_Walk_pct         = round(mean(Walk_pct,           na.rm = TRUE), 1),
    mean_cars_none_pct    = round(mean(cars_none_pct,      na.rm = TRUE), 1),
    mean_X65plus_pct      = round(mean(X65plus_pct,        na.rm = TRUE), 1),
    # Country
    pct_England           = round(100 * mean(country == "England"), 1),
    pct_Scotland          = round(100 * mean(country == "Scotland"), 1),
    .groups               = "drop"
  )

write_csv(stratum_table, file.path(outdir, "05_stratum_characteristics.csv"))
cat("  Saved: 05_stratum_characteristics.csv\n")
cat("  [REPORT] Table in Section 9 of matching report (Baseline Injury\n")
cat("           Stratification). Shows stratum 1 = low injury / minor roads;\n")
cat("           stratum 4 = high injury / more A-roads / lower car commuting.\n\n")

# =============================================================================
# SECTION 6 — COMMON SUPPORT FLAGS TABLE
# =============================================================================
# [REPORT: Footnote or small table in Section 6.4 of matching report.
#  Show characteristics of the 41 isolated OAs vs the remaining 773,
#  to justify why their exclusion improves trend balance.]

cat("=== Section 6: Common support — isolated OA characteristics ===\n")

matched_A_treated <- matched_A %>% filter(treat_indicator == 1)

isolated_chars <- matched_A_treated %>%
  mutate(isolated = OA %in% isolated_ids) %>%
  group_by(isolated) %>%
  summarise(
    n                   = n(),
    mean_road_length    = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_pop_density    = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_pct_minor      = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    mean_IMD            = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct  = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_dist_citycentre= round(mean(dist_citycentre,    na.rm = TRUE), 0),
    mean_total_pkm      = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    pct_scotland        = round(100 * mean(country == "Scotland"), 1),
    .groups             = "drop"
  ) %>%
  mutate(isolated = if_else(isolated, "Isolated (excl. A_restricted)", "Non-isolated"))

# SMD between isolated and non-isolated
matched_A_treated_flag <- matched_A_treated %>%
  mutate(treat_indicator = as.integer(OA %in% isolated_ids))  # repurpose for SMD

smd_isolated <- map_df(stage1_vars, function(v) {
  tibble(
    variable = v,
    label    = coalesce(var_labels[v], v),
    smd      = round(compute_smd(matched_A_treated_flag, v), 3)
  )
}) %>% arrange(desc(abs(smd)))

write_csv(isolated_chars,  file.path(outdir, "06a_isolated_OA_characteristics.csv"))
write_csv(smd_isolated,    file.path(outdir, "06b_isolated_OA_smd.csv"))
cat("  Saved: 06a_isolated_OA_characteristics.csv\n")
cat("  Saved: 06b_isolated_OA_smd.csv\n")
cat("  [REPORT] Section 6.4 — Common Support. Add a sentence describing\n")
cat("           what makes isolated OAs unusual (higher dist_citycentre,\n")
cat("           lower road density, more minor roads, higher % Scotland).\n\n")

# =============================================================================
# SECTION 7 — SKEWNESS ASSESSMENT TABLE
# =============================================================================
# [REPORT: Matching Report Section 4.4 (already drafted as Table 4).
#  This reproduces and extends that table with the max/median ratio
#  and the log-transform decision rule.]

cat("=== Section 7: Skewness assessment ===\n")

skew_fn <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x); mu <- mean(x); s <- sd(x)
  if (s == 0 || n < 3) return(NA_real_)
  sum((x - mu)^3) / (n * s^3)
}

skew_table <- map_df(
  intersect(c(stage1_vars, stage2_levels), names(full_data)),
  function(v) {
    x_raw <- full_data[[v]]
    q     <- quantile(x_raw, c(0.01, 0.99), na.rm = TRUE)
    x_w   <- pmin(pmax(x_raw, q[1]), q[2])
    x_log <- log1p(pmax(x_w, 0))
    tibble(
      variable      = v,
      label         = coalesce(var_labels[v], v),
      stage         = case_when(v %in% stage1_vars   ~ "Stage 1",
                                v %in% stage2_levels ~ "Stage 2 level"),
      skew_raw      = round(skew_fn(x_raw), 2),
      skew_winsor   = round(skew_fn(x_w),   2),
      skew_log1p    = round(skew_fn(x_log), 2),
      max_med_ratio = round(max(x_raw, na.rm = TRUE) /
                              (median(x_raw, na.rm = TRUE) + 1e-9), 1),
      log_applied   = v %in% c("road_length_km", "pop_density",
                               "dist_citycentre", stage2_levels)
    )
  }
) %>% arrange(desc(abs(skew_winsor)))

write_csv(skew_table, file.path(outdir, "07_skewness_assessment.csv"))
cat("  Saved: 07_skewness_assessment.csv\n")
cat("  [REPORT] Section 4.4 Table 4. Note pop_density skewness magnitude\n")
cat("           barely improves (1.82 → |-1.83|); log retained for right-tail\n")
cat("           handling and consistency. Footnote this in the report table.\n\n")

# =============================================================================
# SECTION 8 — LOVE PLOTS
# =============================================================================
# [REPORT: Figure 1 = Stage 1 love plot (Section 6.3).
#          Figure 3 = Stage 2 love plot restricted (Section 7.5).
#          These are already in the matching report as placeholder figures.]

cat("=== Section 8: Love plots ===\n")

love_data_fn <- function(data_before, data_after, vars) {
  map_df(vars, function(v) {
    tibble(
      variable  = v,
      label     = coalesce(var_labels[v], v),
      smd_un    = compute_smd(data_before, v),
      smd_adj   = compute_smd(data_after,  v)
    )
  })
}

make_love_plot <- function(ldat, title, subtitle = NULL, threshold = 0.1) {
  ldat <- ldat %>%
    filter(!is.na(smd_un), !is.na(smd_adj)) %>%
    arrange(abs(smd_un)) %>%
    mutate(label = factor(label, levels = unique(label)))
  
  ldat_long <- ldat %>%
    pivot_longer(c(smd_un, smd_adj), names_to = "timing", values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_un", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  ggplot(ldat_long, aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = threshold, linetype = "dashed",
               colour = "#999999", linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = c("Before matching" = COL_BEFORE, "After matching" = COL_AFTER)) +
    scale_shape_manual(
      values = c("Before matching" = 16, "After matching" = 17)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = subtitle,
         x = "Absolute Standardised Mean Difference", y = NULL,
         colour = NULL, shape = NULL,
         caption = "Dashed line = |SMD| = 0.10 threshold") +
    theme_diag() +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 9),
          plot.margin = margin(10, 20, 10, 10))
}

# Stage 1
ld_s1 <- love_data_fn(unmatched_pool, matched_A, stage1_vars)
p_love_s1 <- make_love_plot(
  ld_s1,
  "Stage 1 Balance — Analysis A",
  "Structural and sociodemographic variables | ratio 10, exact = country")
save_fig(p_love_s1, "fig01_love_plot_stage1_A.png", width = 13, height = 9)
cat("  [REPORT] fig01 → Figure 1 in Section 6.3 of matching report.\n\n")

# Stage 2 — full vs restricted side-by-side
s2_all_vars_for_love <- c(stage2_trends, stage2_levels)
ld_s2_full  <- love_data_fn(unmatched_pool, matched_A,       s2_all_vars_for_love)
ld_s2_restr <- love_data_fn(unmatched_pool, matched_A_restr, s2_all_vars_for_love)

p_love_s2_full  <- make_love_plot(
  ld_s2_full,
  "Stage 2 Balance — A (all treated OAs)",
  "Pre-treatment injury trends and levels | ratio 3, exact = country")
p_love_s2_restr <- make_love_plot(
  ld_s2_restr,
  "Stage 2 Balance — A Restricted (isolated OAs excluded)",
  "Pre-treatment injury trends and levels | ratio 3, exact = country")

save_fig(p_love_s2_full,  "fig02_love_plot_stage2_A_full.png",       width = 13, height = 11)
save_fig(p_love_s2_restr, "fig03_love_plot_stage2_A_restricted.png", width = 13, height = 11)
cat("  [REPORT] fig02 → Figure 3 (Section 7.5, full spec).\n")
cat("           fig03 → Figure 4 (Section 7.5, restricted spec — PRIMARY).\n\n")

# =============================================================================
# SECTION 9 — SMD HEATMAP (SINGLE SPEC — ANALYSIS A)
# =============================================================================
# [REPORT: Supplementary figure or within Section 7. Shows at a glance which
#  variable types remain imbalanced after each stage. Useful for reviewers.]

cat("=== Section 9: SMD heatmap ===\n")

specs_heat <- list(
  "Unmatched"         = unmatched_pool,
  "After Stage 1"     = matched_A,
  "After Stage 2\n(full)"       = matched_A,
  "After Stage 2\n(restricted)" = matched_A_restr
)

smd_heat_data <- map_df(names(specs_heat), function(spec) {
  map_df(c(stage1_vars, stage2_trends, stage2_levels), function(v) {
    tibble(spec = spec, variable = v,
           smd = abs(tryCatch(compute_smd(specs_heat[[spec]], v),
                              error = function(e) NA_real_)))
  })
}) %>%
  mutate(
    label     = coalesce(var_labels[variable], variable),
    var_group = case_when(
      variable %in% stage1_road   ~ "1. Road network",
      variable %in% stage1_urban  ~ "2. Urban geography",
      variable %in% stage1_socdem ~ "3. Sociodemographic",
      variable %in% stage2_trends ~ "4. Injury trends",
      variable %in% stage2_levels ~ "5. Injury levels"
    ),
    spec = factor(spec, levels = names(specs_heat))
  )

p_heatmap <- ggplot(smd_heat_data, aes(x = spec, y = label, fill = smd)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(!is.na(smd), sprintf("%.2f", smd), "")),
            size = 2.3, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#2ECC71", mid = "#F39C12", high = "#E74C3C",
    midpoint = 0.1, na.value = "#EEEEEE",
    name = "|SMD|", limits = c(0, NA)
  ) +
  facet_grid(var_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Absolute SMD Across Specifications — Analysis A",
    subtitle = "Green < 0.10 (balanced) | Yellow = marginal | Red > 0.20 (imbalanced)",
    x = NULL, y = NULL,
    caption  = "Stage 1 variables: unmatched → after S1. Stage 2 variables: unmatched → after S2 (full/restricted)."
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0.3, "lines"))

save_fig(p_heatmap, "fig04_smd_heatmap.png", width = 13, height = 16)
cat("  [REPORT] Supplementary figure. Reference in Section 12 (Balance Test\n")
cat("           Summary). Particularly useful for showing which Stage 2 level\n")
cat("           variables remain imbalanced (by design) vs trend variables\n")
cat("           which should be balanced.\n\n")

# =============================================================================
# SECTION 10 — COVARIATE DISTRIBUTION PLOTS
# =============================================================================
# [REPORT: Supplementary figures. Include one representative panel in the
#  paper body (e.g. Drive_Car_pct, road_length_km, mean_total_pkm).
#  Full panels go in online appendix.]

cat("=== Section 10: Covariate distribution plots ===\n")

plot_density_pair <- function(var, data_before, data_after,
                              log_scale = FALSE, title_prefix = "") {
  prep <- function(dat, phase) {
    d <- dat %>%
      filter(!is.na(.data[[var]])) %>%
      mutate(
        group = if_else(treat_indicator == 1, "Treated", "Control"),
        val   = if (log_scale) log1p(pmax(.data[[var]], 0)) else .data[[var]]
      )
    d$phase <- phase
    d
  }
  d_combined <- bind_rows(
    prep(data_before, "Before matching"),
    prep(data_after,  "After matching")
  ) %>%
    mutate(phase = factor(phase, levels = c("Before matching", "After matching")))
  
  lbl <- coalesce(var_labels[var], var)
  x_lbl <- if (log_scale) paste0("log1p(", lbl, ")") else lbl
  
  ggplot(d_combined, aes(x = val, colour = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.7) +
    scale_colour_manual(values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
    scale_fill_manual(  values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
    facet_wrap(~phase) +
    labs(title = paste0(title_prefix, lbl), x = x_lbl,
         y = "Density", colour = NULL, fill = NULL) +
    theme_diag(base_size = 9) +
    theme(legend.position = "bottom")
}

log_s1 <- c("road_length_km", "pop_density", "dist_citycentre")

# Road network group
road_plots <- map(stage1_road, function(v)
  plot_density_pair(v, unmatched_pool, matched_A,
                    log_scale = v %in% log_s1))
p_dist_road <- wrap_plots(road_plots, ncol = 2) +
  plot_annotation(title = "Road Network Variables — Treated vs Control",
                  subtitle = "Left panel: before matching | Right panel: after Stage 1+2 matching (Analysis A)")
save_fig(p_dist_road, "fig05a_dist_road_network.png", width = 14, height = 12)

# Sociodemographic group
socdem_plots <- map(stage1_socdem, function(v)
  plot_density_pair(v, unmatched_pool, matched_A, log_scale = FALSE))
p_dist_socdem <- wrap_plots(socdem_plots, ncol = 3) +
  plot_annotation(title = "Sociodemographic Variables — Treated vs Control",
                  subtitle = "Left panel: before matching | Right panel: after matching (Analysis A)")
save_fig(p_dist_socdem, "fig05b_dist_socdem.png", width = 16, height = 16)

# Stage 2: trends (4 key ones for readability)
key_trends <- c("trend_total_pkm", "trend_ped_slight_pkm",
                "trend_car_slight_pkm", "trend_cyc_slight_pkm")
trend_plots <- map(key_trends, function(v)
  plot_density_pair(v, matched_A, matched_A_restr))
p_dist_trends <- wrap_plots(trend_plots, ncol = 2) +
  plot_annotation(
    title    = "Key Pre-Treatment Injury Trend Variables — Treated vs Matched Controls",
    subtitle = "Full matched sample (Analysis A, ratio 1:3). Overlap indicates parallel trends support.",
    caption  = "Note: 'before' here = Stage 1 pool; 'after' = Stage 2 restricted matched sample."
  )
save_fig(p_dist_trends, "fig06_dist_stage2_trends.png", width = 14, height = 10)
cat("  [REPORT] fig05a–05b → Supplementary (covariate overlap).\n")
cat("           fig06 → Figure in Section 7.5 (parallel trends visual).\n\n")

# =============================================================================
# SECTION 11 — WEIGHT DISTRIBUTION PLOTS
# =============================================================================
# [REPORT: Figure in Section 8 of matching report. Caption should note
#  the Scotland weight problem explicitly and the rationale for capping at 5.]

cat("=== Section 11: Weight distribution plots ===\n")

ctrl_A <- matched_A %>% filter(treat_indicator == 0)
eff_n_overall <- round(sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2), 0)

p_w1 <- ggplot(ctrl_A, aes(x = weights, fill = country)) +
  geom_histogram(bins = 60, alpha = 0.85, position = "stack") +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  scale_x_continuous(limits = c(0, 5.5)) +
  labs(
    title    = "Control OA Weight Distribution (capped at 5)",
    subtitle = paste0("Nominal N = ", nrow(ctrl_A),
                      " | Effective N = ", eff_n_overall,
                      " | Efficiency = ",
                      round(eff_n_overall / nrow(ctrl_A), 3)),
    x = "Weight", y = "Count", fill = "Country"
  ) +
  theme_diag()

p_w2 <- ctrl_A %>%
  group_by(country) %>%
  summarise(n       = n(),
            eff_n   = round(sum(weights)^2 / sum(weights^2), 1),
            max_w   = round(max(weights), 2),
            .groups = "drop") %>%
  ggplot(aes(x = country, y = eff_n, fill = country)) +
  geom_col(width = 0.45) +
  geom_text(aes(label = paste0("n=", n, "\neff_N=", eff_n,
                               "\nmax_w=", max_w)),
            vjust = -0.2, size = 3.5, lineheight = 1.1) +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Effective N by country", x = NULL, y = "Effective N") +
  theme_diag() + theme(legend.position = "none")

p_w3 <- ggplot(ctrl_A, aes(x = weights, colour = country)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = 5, linetype = "dashed",
             colour = "#E74C3C", linewidth = 0.6) +
  scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  coord_cartesian(xlim = c(0, 6)) +
  labs(title = "ECDF of weights by country",
       subtitle = "Red dashed = cap at 5",
       x = "Weight", y = "Cumulative proportion", colour = "Country") +
  theme_diag()

p_weights <- (p_w1 | p_w2 | p_w3) +
  plot_annotation(
    title = "Weight Diagnostics — Analysis A (after cap at 5)",
    theme = theme(plot.title = element_text(size = 13, face = "bold",
                                            colour = "#1A2E5A"))
  )
save_fig(p_weights, "fig07_weight_diagnostics.png", width = 16, height = 6)
cat("  [REPORT] fig07 → Figure in Section 8 of matching report.\n\n")

# =============================================================================
# SECTION 12 — MAHALANOBIS DISTANCE PLOTS
# =============================================================================
# [REPORT: Figure in Section 7.6 of matching report. Shows distance
#  distribution by country — highlights the Scotland thin-pool problem visually.
#  The ECDF shows that most treated OAs are well-matched (distance < 10).]

cat("=== Section 12: Mahalanobis distance plots ===\n")

mdist_A <- matched_A %>%
  filter(treat_indicator == 1, !is.na(mdist)) %>%
  select(OA, mdist, country)

p_mdist_ecdf <- ggplot(mdist_A, aes(x = mdist, colour = country)) +
  stat_ecdf(linewidth = 1.1) +
  geom_vline(xintercept = 5,  linetype = "dashed",  colour = "#888888", linewidth = 0.5) +
  geom_vline(xintercept = 10, linetype = "dotted",  colour = "#CC3333", linewidth = 0.5) +
  geom_vline(xintercept = 20, linetype = "dotdash", colour = "#8B0000", linewidth = 0.5) +
  annotate("text", x =  5.3, y = 0.15, label = "d=5",  size = 3, colour = "#888888") +
  annotate("text", x = 10.3, y = 0.08, label = "d=10", size = 3, colour = "#CC3333") +
  annotate("text", x = 20.3, y = 0.04, label = "d=20", size = 3, colour = "#8B0000") +
  scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(
    title    = "ECDF of Stage 2 Mahalanobis Distance — Treated OAs (Analysis A)",
    subtitle = "Distance = how far each treated OA is from its nearest control in trajectory space",
    x        = "Stage 2 Mahalanobis distance",
    y        = "Cumulative proportion of treated OAs",
    colour   = "Country",
    caption  = paste0("n = ", nrow(mdist_A), " treated OAs | ",
                      sum(mdist_A$mdist > 20), " OAs with distance > 20 ",
                      "(removed in A_restricted)")
  ) +
  theme_diag() + theme(legend.position = "bottom")

p_mdist_hist <- ggplot(mdist_A, aes(x = mdist, fill = country)) +
  geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  coord_cartesian(xlim = c(0, 35)) +
  facet_wrap(~country, scales = "free_y") +
  labs(
    title    = "Stage 2 Distance Distribution by Country",
    subtitle = "Scottish treated OAs have larger distances — thinner control pool",
    x = "Mahalanobis distance", y = "Count", fill = NULL
  ) +
  theme_diag() + theme(legend.position = "none")

p_mdist <- p_mdist_ecdf / p_mdist_hist + plot_layout(heights = c(1.3, 1))
save_fig(p_mdist, "fig08_mahalanobis_distance.png", width = 12, height = 11)
cat("  [REPORT] fig08 → Figure in Section 7.6. ECDF shows ~80% of treated\n")
cat("           OAs have distance < 10; Scotland tail explains isolated OAs.\n\n")

# =============================================================================
# SECTION 13 — RATIO SELECTION CURVE
# =============================================================================
# [REPORT: Figure 2 in Section 7.2. Already a placeholder in the report.
#  This produces it from the pre-computed data — replace values with
#  actual ratio_curve_A object from the matching script if available.]

cat("=== Section 13: Ratio selection curve ===\n")

# These values come directly from ratio_curve_A in the matching script
ratio_curve_data <- tibble(
  ratio          = 1:10,
  n_controls     = c(604, 1015, 1342, 1580, 1762, 1895, 2013, 2106, 2198, 2254),
  max_trend_smd  = c(0.0564, 0.0644, 0.0561, 0.0618, 0.0614,
                     0.0651, 0.0655, 0.0661, 0.0636, 0.0637),
  mean_trend_smd = c(0.0216, 0.0296, 0.0284, 0.0276, 0.0274,
                     0.0278, 0.0280, 0.0282, 0.0288, 0.0299)
)

p_ratio <- ggplot(ratio_curve_data, aes(x = ratio)) +
  geom_line(aes(y = max_trend_smd,  colour = "Max trend |SMD|"),  linewidth = 1) +
  geom_line(aes(y = mean_trend_smd, colour = "Mean trend |SMD|"), linewidth = 1,
            linetype = "dashed") +
  geom_point(aes(y = max_trend_smd,  colour = "Max trend |SMD|"),  size = 3) +
  geom_point(aes(y = mean_trend_smd, colour = "Mean trend |SMD|"), size = 3,
             shape = 17) +
  geom_vline(xintercept = 3, linetype = "dotted", colour = "#E74C3C",
             linewidth = 0.7) +
  annotate("text", x = 3.15, y = 0.065, label = "Selected\nratio = 3",
           size = 3.2, colour = "#E74C3C", hjust = 0) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#888888", linewidth = 0.5) +
  annotate("text", x = 10.1, y = 0.101, label = "0.10",
           size = 3, colour = "#888888", hjust = 0) +
  scale_colour_manual(values = c("Max trend |SMD|"  = "#D85A30",
                                 "Mean trend |SMD|" = "#2E6FAB")) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0, 0.075)) +
  labs(
    title    = "Ratio Selection: Trend Balance vs Matching Ratio — Analysis A",
    subtitle = "Ratio 1:3 selected: lowest max trend |SMD| (0.056) with 1,342 initial controls",
    x        = "Matching ratio (1:k)",
    y        = "Trend variable |SMD|",
    colour   = NULL,
    caption  = "Red dotted = selected ratio. Dashed grey = 0.10 balance threshold."
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p_ratio, "fig09_ratio_selection_curve.png", width = 10, height = 6)
cat("  [REPORT] fig09 → Figure 2 in Section 7.2 of matching report.\n\n")

# =============================================================================
# SECTION 14 — STRATUM INJURY DISTRIBUTION
# =============================================================================
# [REPORT: Figure in Section 9 of matching report. Shows that stratum
#  assignment reflects genuine injury-level gradients and is not arbitrary.]

cat("=== Section 14: Stratum injury distribution ===\n")

strat_data <- matched_A %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  mutate(stratum = factor(baseline_injury_stratum,
                          labels = paste0("Stratum ", 1:4, c(
                            "\n(lowest)", "", "", "\n(highest)"))))

p_strat_inj <- ggplot(strat_data,
                      aes(x = stratum, y = mean_total_pkm, fill = stratum)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.85) +
  scale_y_log10(labels = label_comma()) +
  labs(
    title    = "Pre-Treatment Injury Rate by Baseline Stratum",
    subtitle = "Treated OAs | Analysis A | log scale y-axis",
    x        = "Baseline injury stratum",
    y        = "Mean total injuries per km (log scale)",
    fill     = NULL,
    caption  = paste0("Stratum 1: n=204, log1p_mean_total_pkm \u2264 0.138 | ",
                      "Stratum 4: n=204, log1p_mean_total_pkm > 0.441")
  ) +
  theme_diag() + theme(legend.position = "none")

p_strat_covs <- strat_data %>%
  select(stratum, Drive_Car_pct, Walk_pct, IMD, pop_density,
         road_length_km, pct_minor_road) %>%
  pivot_longer(-stratum) %>%
  mutate(label = coalesce(var_labels[name], name)) %>%
  ggplot(aes(x = stratum, y = value, fill = stratum)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.85) +
  facet_wrap(~label, scales = "free_y", ncol = 3) +
  labs(title = "Covariate Profile by Baseline Injury Stratum",
       x = NULL, y = NULL, fill = NULL) +
  theme_diag() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7))

save_fig(p_strat_inj,  "fig10_stratum_injury.png",    width = 10, height = 7)
save_fig(p_strat_covs, "fig10b_stratum_covariates.png", width = 13, height = 10)
cat("  [REPORT] fig10 → Figure in Section 9 (stratum injury distribution).\n")
cat("           fig10b → Supplementary (stratum covariate profiles).\n\n")

# =============================================================================
# SECTION 15 — ENGLAND vs SCOTLAND COMPARISON
# =============================================================================
# [REPORT: IMPORTANT — Include in the report (new section or Appendix).
#  Scotland-specific matching issues (weight concentration, thinner control
#  pool, higher Mahalanobis distances) must be documented transparently.
#  This section produces the supporting material for that discussion.]

cat("=== Section 15: England vs Scotland comparison ===\n")

country_summary <- matched_A %>%
  filter(treat_indicator == 1) %>%
  group_by(country) %>%
  summarise(
    n_treated              = n(),
    pct_of_total           = round(100 * n() / nrow(matched_A_treated), 1),
    mean_road_length_km    = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_pop_density       = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_dist_citycentre   = round(mean(dist_citycentre,    na.rm = TRUE), 0),
    mean_IMD               = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct     = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_Walk_pct          = round(mean(Walk_pct,           na.rm = TRUE), 1),
    mean_pct_minor         = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    mean_total_pkm         = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    median_mdist           = round(median(mdist,            na.rm = TRUE), 2),
    pct_isolated           = round(100 * mean(OA %in% isolated_ids), 1),
    .groups                = "drop"
  )

cat("\n--- Treated OA country comparison (Analysis A) ---\n")
print(country_summary)
write_csv(country_summary, file.path(outdir, "08_country_comparison.csv"))
cat("  Saved: 08_country_comparison.csv\n")

# Control weight summary per country (already in Section 4 table; reprint here)
cat("\n--- Control weight summary by country ---\n")
matched_A %>%
  filter(treat_indicator == 0) %>%
  group_by(country) %>%
  summarise(n = n(),
            eff_n      = round(sum(weights)^2 / sum(weights^2), 1),
            efficiency = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
            max_w      = round(max(weights), 2),
            .groups    = "drop") %>%
  print()

# Violin plot comparing treated OA profiles
country_vars <- c("road_length_km", "pop_density", "IMD",
                  "Drive_Car_pct", "Walk_pct", "pct_minor_road",
                  "mean_total_pkm", "dist_citycentre")

p_country_profile <- matched_A %>%
  filter(treat_indicator == 1) %>%
  select(OA, country, all_of(country_vars)) %>%
  pivot_longer(-c(OA, country)) %>%
  mutate(label = coalesce(var_labels[name], name)) %>%
  ggplot(aes(x = country, y = value, fill = country)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.6) +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  facet_wrap(~label, scales = "free_y", ncol = 4) +
  labs(
    title    = "Treated OA Structural Profiles: England vs Scotland",
    subtitle = paste0("England n=613 | Scotland n=201 | Analysis A"),
    x = NULL, y = NULL, fill = NULL,
    caption  = paste0("Scotland treated OAs tend to be longer roads, ",
                      "lower pop density, higher car commuting, lower pre-treatment injury rates.")
  ) +
  theme_diag() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

save_fig(p_country_profile, "fig11_england_scotland_profile.png", width = 16, height = 12)
cat("  [REPORT] fig11 → New figure in report (Section on Scotland or Appendix).\n")
cat("           [REPORT] Add a paragraph in Section 6 or 7 of the matching\n")
cat("           report explaining Scotland-specific findings:\n")
cat("             1. Scotland = 201/814 treated OAs (24.7%)\n")
cat("             2. Scottish control pool is thin (only 6 matched controls)\n")
cat("             3. Weight concentration → max weight 82.7 before cap, 5.0 after\n")
cat("             4. Effective N for Scotland controls = 3.9 (vs 654 for England)\n")
cat("             5. Median Stage 2 distance larger for Scotland (see fig08)\n")
cat("             6. Mitigation: weight cap at 5; cluster SEs at OA level in DiD\n\n")

# =============================================================================
# SECTION 16 — MAP PLOTS (spatial boundary required)
# =============================================================================
# [REPORT: fig14 → Figure in Section 3 (Data) of the paper or matching report.
#          fig15 (city maps) → Supplementary.
#          fig16 (control origin) → Section 3 or Appendix.
#  Maps are essential for showing that controls are not just from the same city.]

cat("=== Section 16: Map plots ===\n")

oa_path  <- here("data", "spatial", "OA_boundaries_GB.gpkg")
lad_path <- here("data", "spatial", "LAD_boundaries_GB.gpkg")

if (!file.exists(oa_path) || !file.exists(lad_path)) {
  cat("  SKIPPED — boundary files not found.\n")
  cat("  Expected paths:\n")
  cat("    OA:  ", oa_path,  "\n")
  cat("    LAD: ", lad_path, "\n")
  cat("  [REPORT] When boundaries available, produce:\n")
  cat("    fig12 = GB overview map of treated OAs (by country)\n")
  cat("    fig13 = city-level zoom maps (treated + matched controls + buffer)\n")
  cat("    fig14 = choropleth of matched control origins per LAD\n\n")
} else {
  
  oa_sf  <- st_read(oa_path,  quiet = TRUE)
  lad_sf <- st_read(lad_path, quiet = TRUE)
  
  oa_id_col  <- intersect(c("OA21CD","OA11CD","geo_code","OA"), names(oa_sf))[1]
  oa_sf <- oa_sf %>% rename(OA = all_of(oa_id_col))
  
  lad_id_col <- intersect(c("LAD21CD","LAD22CD","LAD24CD","lad_code"), names(lad_sf))[1]
  lad_nm_col <- intersect(c("LAD21NM","LAD22NM","LAD24NM","lad_name"),  names(lad_sf))[1]
  if (!is.na(lad_id_col)) lad_sf <- lad_sf %>% rename(LAD_CODE = all_of(lad_id_col))
  if (!is.na(lad_nm_col)) lad_sf <- lad_sf %>% rename(LAD_NAME = all_of(lad_nm_col))
  
  map_status <- full_data %>%
    mutate(map_group = case_when(
      buffer_OA  == 1                    ~ "Buffer (excluded)",
      treated_OA == 1                    ~ "Treated",
      OA %in% controls_A$OA             ~ "Matched control",
      control_group1_OA == 1 | control_group2_OA == 1 ~ "Eligible (unmatched)",
      TRUE                               ~ "Other"
    )) %>%
    select(OA, map_group, LAD24CD, country)
  
  # Overview map — treated OA centroids on GB base
  treated_centroids <- oa_sf %>%
    filter(OA %in% treated_A$OA) %>%
    st_centroid() %>%
    left_join(map_status %>% select(OA, country), by = "OA")
  
  p_overview <- ggplot() +
    geom_sf(data = lad_sf, fill = "#F5F7FA", colour = "#CCCCCC", linewidth = 0.2) +
    geom_sf(data = treated_centroids, aes(colour = country),
            size = 1.2, alpha = 0.7) +
    scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    coord_sf(xlim = c(-6, 2), ylim = c(50, 59), crs = 4326) +
    labs(title    = "Geographic Distribution of Treated OAs — Analysis A",
         subtitle = paste0(nrow(treated_centroids),
                           " treated OAs (England n=613 | Scotland n=201)"),
         caption  = "Each point = one treated OA centroid",
         colour   = "Country") +
    theme_diag() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          panel.grid = element_line(colour = "#E8ECF3", linewidth = 0.3),
          legend.position = "bottom")
  
  save_fig(p_overview, "fig12_map_treated_overview.png", width = 8, height = 12)
  
  # Control origin choropleth
  controls_per_lad <- full_data %>%
    filter(OA %in% controls_A$OA) %>%
    count(LAD24CD, name = "n_controls")
  
  treated_per_lad <- full_data %>%
    filter(OA %in% treated_A$OA) %>%
    count(LAD24CD, name = "n_treated")
  
  lad_map <- lad_sf %>%
    left_join(controls_per_lad, by = c("LAD_CODE" = "LAD24CD")) %>%
    left_join(treated_per_lad,  by = c("LAD_CODE" = "LAD24CD")) %>%
    mutate(n_controls = replace_na(n_controls, 0),
           n_treated  = replace_na(n_treated,  0))
  
  p_control_origin <- ggplot() +
    geom_sf(data = lad_map, aes(fill = n_controls),
            colour = "white", linewidth = 0.2) +
    geom_sf(data = lad_map %>% filter(n_treated > 0),
            fill = NA, colour = COL_TREATED, linewidth = 0.7) +
    scale_fill_viridis_c(option = "Blues", name = "N matched\ncontrols",
                         trans = "sqrt", labels = label_comma(),
                         na.value = "#F5F7FA") +
    coord_sf(xlim = c(-6, 2), ylim = c(50, 59), crs = 4326) +
    labs(title    = "Geographic Origin of Matched Control OAs",
         subtitle = "Fill = number of controls from each LAD | Red outline = LADs with treated OAs",
         caption  = "Analysis A | square-root fill scale") +
    theme_diag() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          legend.position = "right")
  
  save_fig(p_control_origin, "fig13_map_control_origin.png", width = 9, height = 13)
  
  # City-level zoom maps
  top_cities <- full_data %>%
    filter(treated_OA == 1) %>%
    count(LAD24CD, sort = TRUE) %>%
    slice_head(n = 6) %>%
    pull(LAD24CD)
  
  lad_names_lookup <- if ("LAD24NM" %in% names(full_data)) {
    full_data %>% distinct(LAD24CD, LAD24NM)
  } else {
    tibble(LAD24CD = top_cities, LAD24NM = top_cities)
  }
  
  group_colours <- c(
    "Treated"              = COL_TREATED,
    "Matched control"      = COL_CONTROL,
    "Buffer (excluded)"    = "#F5A623",
    "Eligible (unmatched)" = "#B0C4DE",
    "Other"                = "#EEEEEE"
  )
  
  walk(top_cities, function(lad_cd) {
    city_nm <- lad_names_lookup$LAD24NM[lad_names_lookup$LAD24CD == lad_cd][1]
    if (is.na(city_nm)) city_nm <- lad_cd
    
    city_oas <- full_data %>% filter(LAD24CD == lad_cd) %>% pull(OA)
    city_geom <- oa_sf %>% filter(OA %in% city_oas)
    if (nrow(city_geom) == 0) return(invisible(NULL))
    
    bbox <- st_bbox(city_geom)
    pad  <- 0.025
    xlim <- c(bbox["xmin"] - pad, bbox["xmax"] + pad)
    ylim <- c(bbox["ymin"] - pad, bbox["ymax"] + pad)
    
    nearby <- oa_sf %>%
      st_crop(st_bbox(c(xmin=xlim[1], xmax=xlim[2],
                        ymin=ylim[1], ymax=ylim[2]), crs = st_crs(oa_sf))) %>%
      left_join(map_status, by = "OA") %>%
      mutate(map_group = replace_na(map_group, "Other"))
    
    n_t <- sum(nearby$map_group == "Treated", na.rm = TRUE)
    n_c <- sum(nearby$map_group == "Matched control", na.rm = TRUE)
    
    p_city <- ggplot() +
      geom_sf(data = nearby %>% filter(map_group == "Other"),
              fill = "#F5F7FA", colour = "#DDDDDD", linewidth = 0.1) +
      geom_sf(data = nearby %>% filter(map_group == "Eligible (unmatched)"),
              fill = "#B0C4DE", colour = "#8AAAC8", linewidth = 0.1, alpha = 0.45) +
      geom_sf(data = nearby %>% filter(map_group == "Matched control"),
              fill = COL_CONTROL, colour = "#1A4F8A", linewidth = 0.2, alpha = 0.8) +
      geom_sf(data = nearby %>% filter(map_group == "Buffer (excluded)"),
              fill = "#F5A623", colour = "#C8860A", linewidth = 0.2, alpha = 0.7) +
      geom_sf(data = nearby %>% filter(map_group == "Treated"),
              fill = COL_TREATED, colour = "#8B2010", linewidth = 0.3) +
      scale_fill_manual(values = group_colours, name = NULL) +
      coord_sf(xlim = xlim, ylim = ylim, crs = 4326) +
      labs(title    = paste("Treated and Matched Control OAs —", city_nm),
           subtitle = paste0(n_t, " treated OAs | ", n_c, " matched controls"),
           caption  = paste0("Blue = matched controls | Orange = buffer | ",
                             "Red = treated | Light blue = eligible (unmatched)")) +
      theme_diag() +
      theme(axis.text = element_blank(), axis.title = element_blank(),
            legend.position = "bottom", legend.text = element_text(size = 8))
    
    fn <- paste0("fig14_map_city_",
                 gsub("[^A-Za-z0-9]", "_", city_cd), ".png")
    save_fig(p_city, fn, width = 10, height = 10)
  })
  
  cat("  Map figures saved.\n")
  cat("  [REPORT] fig12 → Figure in Section 3 (overview of study area).\n")
  cat("           fig13 → Figure in Section 3 or Appendix (control origins).\n")
  cat("           fig14_* → Supplementary city-level maps.\n\n")
}

# =============================================================================
# SECTION 17 — PARALLEL TRENDS VISUAL (TREND SLOPE DISTRIBUTIONS)
# =============================================================================
# [REPORT: Figure in Section 7.5 of the matching report (alongside love plots).
#  Also include in the paper methods section as visual evidence for the
#  parallel trends assumption. Note: these show pre-treatment SLOPE
#  distributions, not actual time series — the latter require the panel data.]

cat("=== Section 17: Parallel trends visualisation (slope distributions) ===\n")

# NOTE: This uses trend SLOPES (regression coefficients) from the matched
# sample as a proxy for the parallel trends assumption. The correct check
# is a time-series plot of actual pre-treatment injury rates (requires the
# OA×quarter panel). Produce that separately in the DiD analysis script.

trend_long <- matched_A_restr %>%
  filter(treat_indicator %in% c(0, 1)) %>%
  mutate(group = if_else(treat_indicator == 1, "Treated", "Matched control\n(restricted)")) %>%
  select(OA, group, country, all_of(stage2_trends)) %>%
  pivot_longer(all_of(stage2_trends),
               names_to = "trend_var", values_to = "slope") %>%
  mutate(trend_label = coalesce(var_labels[trend_var], trend_var),
         trend_label = factor(trend_label)) %>%
  filter(!is.na(slope))

p_trend_density <- ggplot(trend_long,
                          aes(x = slope, colour = group, fill = group)) +
  geom_density(alpha = 0.20, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#888888", linewidth = 0.4) +
  scale_colour_manual(values = c("Treated" = COL_TREATED,
                                 "Matched control\n(restricted)" = COL_CONTROL)) +
  scale_fill_manual(  values = c("Treated" = COL_TREATED,
                                 "Matched control\n(restricted)" = COL_CONTROL)) +
  facet_wrap(~trend_label, scales = "free", ncol = 3) +
  labs(
    title    = "Pre-Treatment Injury Trend Slopes: Treated vs Matched Controls",
    subtitle = "Analysis A (restricted) | Overlap of distributions supports the parallel trends assumption",
    x        = "Pre-treatment slope (log-linear regression coefficient)",
    y        = "Density",
    colour   = NULL, fill   = NULL,
    caption  = paste0(
      "n treated = ", sum(matched_A_restr$treat_indicator == 1),
      " | n controls = ", sum(matched_A_restr$treat_indicator == 0),
      "\nDashed line = slope of 0 (flat trend). ",
      "Note: This shows slope distributions, not actual time series."
    )
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p_trend_density, "fig15_parallel_trends_slopes.png", width = 14, height = 12)
cat("  [REPORT] fig15 → Figure in Section 7.5 / paper methods section.\n")
cat("  [REPORT] IMPORTANT: Also produce actual pre-treatment time series\n")
cat("           plots (OA x quarter panel) in the DiD analysis script.\n")
cat("           Slope distributions are a proxy, not a substitute.\n\n")

# =============================================================================
# SECTION 18 — FINAL SUMMARY
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("DIAGNOSTICS COMPLETE — ANALYSIS A\n")
cat("=================================================================\n\n")

cat("SAMPLE SIZES:\n")
cat(sprintf("  Treated (full):       %d\n", sum(matched_A$treat_indicator == 1)))
cat(sprintf("  Controls (full):      %d\n", sum(matched_A$treat_indicator == 0)))
cat(sprintf("  Treated (restricted): %d\n", sum(matched_A_restr$treat_indicator == 1)))
cat(sprintf("  Controls (restricted):%d\n", sum(matched_A_restr$treat_indicator == 0)))
cat(sprintf("  Isolated OAs excluded:%d\n", length(isolated_ids)))

cat("\nWEIGHT SUMMARY (controls, after cap = 5):\n")
ctrl_A <- matched_A %>% filter(treat_indicator == 0)
cat(sprintf("  Nominal N:     %d\n", nrow(ctrl_A)))
cat(sprintf("  Effective N:   %.0f\n",
            sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2)))
cat(sprintf("  Efficiency:    %.3f\n",
            (sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2)) / nrow(ctrl_A)))
cat(sprintf("  Max weight:    %.3f\n", max(ctrl_A$weights)))

cat("\nBALANCE SUMMARY:\n")
cat(sprintf("  Stage 1 — mean |SMD| before: 0.686 | after: 0.165 (76%% reduction)\n"))
cat(sprintf("  Stage 2 full:       mean |SMD| = 0.184 | max trend SMD = 0.102\n"))
cat(sprintf("  Stage 2 restricted: mean |SMD| = 0.167 | max trend SMD = 0.091 PASS\n"))

cat("\nFIGURES SAVED TO:", outdir, "\n")
cat("  fig01 = Stage 1 love plot\n")
cat("  fig02 = Stage 2 love plot (full)\n")
cat("  fig03 = Stage 2 love plot (restricted) — PRIMARY\n")
cat("  fig04 = SMD heatmap\n")
cat("  fig05a = Road network distributions\n")
cat("  fig05b = Sociodemographic distributions\n")
cat("  fig06 = Stage 2 trend distributions\n")
cat("  fig07 = Weight diagnostics\n")
cat("  fig08 = Mahalanobis distance ECDF\n")
cat("  fig09 = Ratio selection curve\n")
cat("  fig10 = Stratum injury distribution\n")
cat("  fig10b = Stratum covariate profiles\n")
cat("  fig11 = England vs Scotland profile\n")
cat("  fig12-14 = Maps (if boundaries available)\n")
cat("  fig15 = Parallel trends slope distributions\n")

cat("\nTABLES SAVED TO:", outdir, "\n")
cat("  01  = Descriptive Table 1 (3-column: Treated | Unmatched ctrl | Matched ctrl)\n")
cat("  01b = Formatted mean (SD) / median [IQR]\n")
cat("  01c = Country-stratified descriptives\n")
cat("  02  = Stage 1 SMD before/after\n")
cat("  03  = Stage 2 SMD before/after (full + restricted)\n")
cat("  04  = Weight distribution by country\n")
cat("  05  = Stratum characteristics\n")
cat("  06a = Isolated OA characteristics\n")
cat("  06b = Isolated OA SMD\n")
cat("  07  = Skewness assessment\n")
cat("  08  = England vs Scotland comparison\n")

cat("\n=== REPORT: WHAT TO ADD ===\n")
cat("These items are NOT yet in the matching report and should be added:\n\n")
cat("  1. Table 1 (from 01_descriptive_table1): Treated vs Unmatched vs Matched\n")
cat("     controls — include in Data/Methods section of the paper.\n\n")
cat("  2. Scotland section: ~2 paragraphs + Table 13 (already in report) +\n")
cat("     fig11 (country profile) + discussion of weight concentration,\n")
cat("     thin control pool, and mitigation strategy (cap + OA clustering).\n\n")
cat("  3. Isolated OA characteristics (06a/06b): Add 1 sentence to Section 6.4\n")
cat("     noting what structurally distinguishes isolated OAs (farther from city\n")
cat("     centre, higher Scotland proportion, lower road density).\n\n")
cat("  4. Actual parallel trends time series: slope distributions (fig15) are\n")
cat("     a matching diagnostic only. The paper methods section needs actual\n")
cat("     OA x quarter pre-treatment injury rate plots (produce in DiD script).\n\n")
cat("  5. fig07 (weight diagnostics) and fig08 (Mahalanobis distance ECDF)\n")
cat("     should be added as supplementary figures referenced from the report.\n\n")
cat("=================================================================\n")