# =============================================================================
# POST-MATCHING DIAGNOSTICS
# OA-Level Two-Stage Mahalanobis Distance Matching 
#
# INPUTS (from data/processed/):
#   OA_matched_full_mixed.rds          — full matched dataset (England + Scotland)
#   OA_matched_treated_mixed.rds       — treated OAs + weights + stratum
#   OA_matched_donors_mixed.rds        — control OAs + weights
#   OA_common_support_flags_mixed.rds  — structurally isolated OA flags
#   OA_matching_pairs_mixed.rds        — treated→control OA pairs
#   OA_matching_census.rds             — full original dataset (unmatched pool)
#
# SOURCE: 16_Matching_England_othercityControlsScotland_mix.R
#   England: other-city controls only
#   Scotland: other-city + same-city controls (see script 16 header)
#
# SPATIAL INPUTS (Section 16 — maps):
#   data/spatial/OA_boundaries_GB.gpkg
#   data/spatial/LAD_boundaries_GB.gpkg
#
# OUTPUTS (to output/diagnostics/):
#   Tables  : 01–08 CSV files
#   Figures : fig01–fig15 PNG files
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

theme_diag <- function(base_size = 14) {
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

save_fig <- function(p, filename, width = 14, height = 10, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# =============================================================================
# LOAD DATA
# =============================================================================

matched_A  <- readRDS(here("data", "processed", "OA_matched_full_mixed.rds"))
treated_A  <- readRDS(here("data", "processed", "OA_matched_treated_mixed.rds"))
controls_A <- readRDS(here("data", "processed", "OA_matched_donors_mixed.rds"))
csupport   <- readRDS(here("data", "processed", "OA_common_support_flags_mixed.rds"))
pairs_mixed <- readRDS(here("data", "processed", "OA_matching_pairs_mixed.rds"))
full_data  <- readRDS(here("data", "processed", "OA_matching_census.rds"))

cat("  Treated =", sum(matched_A$treat_indicator == 1),
    "| Controls =", sum(matched_A$treat_indicator == 0), "\n\n")

# Unmatched eligible pool (denominator for before/after comparisons)
unmatched_pool <- full_data %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# =========
# VARIABLES


stage1_road   <- c("log1p_road_length_km", "log1p_road_density_m_km2",
                   "log_area_km2",
                   "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban  <- c("log1p_dist_citycentre", "log1p_pop_density",
                   "log1p_business_retail_per_km2")
stage1_socdem <- c("IMD",
                   "cars_one_pct", "cars_twoPlus_pct",
                   "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
                   "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
                   "Taxi_pct", "workAthome_pct", "Other_pct",
                   "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
                   "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
                   "age_45to64_pct", "age_65plus_pct")
stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

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
  # Road infrastructure (log-transformed names used in matching)
  log1p_road_length_km                 = "Road length (km)",
  log1p_road_density_m_km2             = "Road density (m/km\u00b2)",
  log_area_km2                         = "Area (km\u00b2)",
  pct_A_road                           = "% A-road",
  pct_B_road                           = "% B-road",
  pct_minor_road                       = "% Minor road",
  # Urban / area (log-transformed)
  log1p_dist_citycentre                = "Distance to city centre (m)",
  log1p_pop_density                    = "Population density (persons/km\u00b2)",
  log1p_business_retail_per_km2        = "Retail businesses (per km\u00b2)",
  # Raw names (for density plots and other uses)
  road_density_m_km2                   = "Road density (m/km\u00b2)",
  road_length_km                       = "Road length (km)",
  dist_citycentre                      = "Distance to city centre (m)",
  pop_density                          = "Population density (persons/km\u00b2)",
  area_km2                             = "Area (km\u00b2)",
  business_retail_per_km2              = "Retail businesses (per km\u00b2)",
  # Socioeconomic
  IMD                                  = "Index of Multiple Deprivation",
  cars_one_pct                         = "% households: 1 car",
  cars_twoPlus_pct                     = "% households: 2+ cars",
  Drive_Car_pct                        = "% commuting: drive car",
  Passenger_Car_pct                    = "% commuting: car passenger",
  Walk_pct                             = "% commuting: walk",
  Bicycle_pct                          = "% commuting: bicycle",
  bus_Coach_pct                        = "% commuting: bus/coach",
  Train_pct                            = "% commuting: train",
  Underground_train_tram_pct           = "% commuting: underground/tram",
  Taxi_pct                             = "% commuting: taxi",
  workAthome_pct                       = "% working from home",
  Other_pct                            = "% commuting: other mode",
  # Ethnicity
  White_pct                            = "% White",
  Mixed_pct                            = "% Mixed ethnicity",
  Asian_pct                            = "% Asian",
  Black_pct                            = "% Black",
  # Age structure (derived groups used in matching)
  age_under15_pct                      = "% aged under 15",
  age_15to24_pct                       = "% aged 15\u201324",
  age_25to44_pct                       = "% aged 25\u201344",
  age_45to64_pct                       = "% aged 45\u201364",
  age_65plus_pct                       = "% aged 65+",
  # Age structure (original columns, used for distribution plots)
  X4under_pct                          = "% aged under 5",
  X5to19_pct                           = "% aged 5\u201319",
  X20to24_pct                          = "% aged 20\u201324",
  X65plus_pct                          = "% aged 65+",
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
#DESCRIPTIVE SUMMARY TABLES

cat("=== Section 2: Descriptive summary tables ===\n")

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

desc_treated_pre  <- desc_stats(
  unmatched_pool %>% filter(treat_indicator == 1), all_desc_vars, "Treated")
desc_control_pre  <- desc_stats(
  unmatched_pool %>% filter(treat_indicator == 0), all_desc_vars, "Unmatched controls")
desc_control_post <- desc_stats(
  matched_A      %>% filter(treat_indicator == 0), all_desc_vars, "Matched controls")

desc_table_1 <- bind_rows(desc_treated_pre, desc_control_pre, desc_control_post) %>%
  pivot_wider(
    id_cols     = c(variable, label),
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

desc_formatted <- bind_rows(desc_treated_pre, desc_control_pre, desc_control_post) %>%
  mutate(
    mean_sd    = sprintf("%.3f (%.3f)", mean, sd),
    median_iqr = sprintf("%.3f [%.3f\u2013%.3f]", median, p25, p75)
  ) %>%
  select(variable, label, group, n, mean_sd, median_iqr)

write_csv(desc_formatted,
          file.path(outdir, "01b_descriptive_formatted.csv"))
cat("  Saved: 01b_descriptive_formatted.csv\n")

country_desc <- map_df(c("England", "Scotland"), function(ctry) {
  d <- matched_A %>% filter(treat_indicator == 1, country == ctry)
  desc_stats(d, all_of(intersect(all_desc_vars, names(d))), ctry)
})

write_csv(country_desc,
          file.path(outdir, "01c_descriptive_by_country.csv"))
cat("  Saved: 01c_descriptive_by_country.csv\n\n")

# ==================================
# — SMD TABLES (BEFORE / AFTER)
# =============================
cat("=== Section 3: SMD before/after tables ===\n")

smd_table <- function(vars, data_before, data_after) {
  map_df(vars, function(v) {
    tibble(
      variable       = v,
      label          = coalesce(var_labels[v], v),
      smd_before     = round(compute_smd(data_before, v), 4),
      smd_after      = round(compute_smd(data_after,  v), 4)
    )
  })
}

smd_s1 <- smd_table(stage1_vars, unmatched_pool, matched_A) %>%
  rename(smd_unmatched = smd_before,
         smd_after_S1  = smd_after) %>%
  mutate(balanced_after_S1 = abs(smd_after_S1) < 0.1)

write_csv(smd_s1, file.path(outdir, "02_smd_stage1.csv"))
cat("  Saved: 02_smd_stage1.csv\n")

smd_s2 <- smd_table(c(stage2_trends, stage2_levels), unmatched_pool, matched_A) %>%
  rename(
    smd_preS2  = smd_before,
    smd_postS2 = smd_after
  ) %>%
  mutate(
    var_type  = if_else(variable %in% stage2_trends, "Trend", "Level"),
    balanced  = abs(smd_postS2) < 0.1
  )

write_csv(smd_s2, file.path(outdir, "03_smd_stage2.csv"))
cat("  Saved: 03_smd_stage2.csv\n\n")

# =====================================
#  — WEIGHT DISTRIBUTION 
# ============================

cat("=== Section 4: Weight distribution table ===\n")

weight_table <- matched_A %>%
  filter(treat_indicator == 0) %>%
  group_by(country) %>%
  summarise(
    n_controls    = n(),
    mean_weight   = round(mean(weights), 3),
    median_weight = round(median(weights), 3),
    sd_weight     = round(sd(weights), 3),
    p90_weight    = round(quantile(weights, 0.90), 3),
    p95_weight    = round(quantile(weights, 0.95), 3),
    max_weight    = round(max(weights), 3),
    eff_N         = round(sum(weights)^2 / sum(weights^2), 1),
    efficiency    = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
    n_at_cap      = sum(weights >= 5),
    pct_at_cap    = round(100 * mean(weights >= 5), 1),
    .groups       = "drop"
  ) %>%
  bind_rows(
    matched_A %>% filter(treat_indicator == 0) %>%
      summarise(
        country       = "Overall",
        n_controls    = n(),
        mean_weight   = round(mean(weights), 3),
        median_weight = round(median(weights), 3),
        sd_weight     = round(sd(weights), 3),
        p90_weight    = round(quantile(weights, 0.90), 3),
        p95_weight    = round(quantile(weights, 0.95), 3),
        max_weight    = round(max(weights), 3),
        eff_N         = round(sum(weights)^2 / sum(weights^2), 1),
        efficiency    = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
        n_at_cap      = sum(weights >= 5),
        pct_at_cap    = round(100 * mean(weights >= 5), 1)
      )
  )

write_csv(weight_table, file.path(outdir, "04_weight_distribution.csv"))
cat("  Saved: 04_weight_distribution.csv\n\n")

# =============================================================================
# SECTION 5 — STRATUM CHARACTERISTICS TABLE
# =============================================================================

cat("=== Section 5: Stratum characteristics ===\n")

stratum_table <- matched_A %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  group_by(Stratum = baseline_injury_stratum) %>%
  summarise(
    n                     = n(),
    mean_total_pkm        = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    sd_total_pkm          = round(sd(mean_total_pkm,       na.rm = TRUE), 5),
    median_total_pkm      = round(median(mean_total_pkm,   na.rm = TRUE), 5),
    mean_road_length_km   = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_road_density     = round(mean(road_density_m_km2, na.rm = TRUE), 0),
    mean_pct_A_road       = round(mean(pct_A_road,         na.rm = TRUE), 1),
    mean_pct_minor        = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    mean_pop_density      = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_dist_citycentre  = round(mean(dist_citycentre,    na.rm = TRUE), 0),
    mean_IMD              = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct    = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_Walk_pct         = round(mean(Walk_pct,           na.rm = TRUE), 1),
    mean_cars_none_pct    = round(mean(cars_none_pct,      na.rm = TRUE), 1),
    mean_X65plus_pct      = round(mean(X65plus_pct,        na.rm = TRUE), 1),
    pct_England           = round(100 * mean(country == "England"), 1),
    pct_Scotland          = round(100 * mean(country == "Scotland"), 1),
    .groups               = "drop"
  )

write_csv(stratum_table, file.path(outdir, "05_stratum_characteristics.csv"))
cat("  Saved: 05_stratum_characteristics.csv\n\n")

# ========================================
# SECTION 6 — COMMON SUPPORT FLAGS TABLE
# =======================================

cat("=== Section 6: Common support — isolated OA characteristics ===\n")

isolated_ids <- csupport$treated_OA

matched_A_treated <- matched_A %>% filter(treat_indicator == 1)

isolated_chars <- matched_A_treated %>%
  mutate(isolated = OA %in% isolated_ids) %>%
  group_by(isolated) %>%
  summarise(
    n                    = n(),
    mean_road_length     = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_pop_density     = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_pct_minor       = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    mean_IMD             = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct   = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_dist_citycentre = round(mean(dist_citycentre,    na.rm = TRUE), 0),
    mean_total_pkm       = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    pct_scotland         = round(100 * mean(country == "Scotland"), 1),
    .groups              = "drop"
  ) %>%
  mutate(isolated = if_else(isolated, "Isolated", "Non-isolated"))

matched_A_treated_flag <- matched_A_treated %>%
  mutate(treat_indicator = as.integer(OA %in% isolated_ids))

smd_isolated <- map_df(stage1_vars, function(v) {
  tibble(
    variable = v,
    label    = coalesce(var_labels[v], v),
    smd      = round(compute_smd(matched_A_treated_flag, v), 3)
  )
}) %>% arrange(desc(abs(smd)))

write_csv(isolated_chars, file.path(outdir, "06a_isolated_OA_characteristics.csv"))
write_csv(smd_isolated,   file.path(outdir, "06b_isolated_OA_smd.csv"))
cat("  Saved: 06a_isolated_OA_characteristics.csv\n")
cat("  Saved: 06b_isolated_OA_smd.csv\n\n")

# =============================================================================
# SECTION 7 — SKEWNESS ASSESSMENT TABLE
# =============================================================================

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
cat("  Saved: 07_skewness_assessment.csv\n\n")

# =============================================================================
# SECTION 8 — LOVE PLOTS
# =============================================================================

cat("=== Section 8: Love plots ===\n")

love_data_fn <- function(data_before, data_after, vars) {
  map_df(vars, function(v) {
    tibble(
      variable = v,
      label    = coalesce(var_labels[v], v),
      smd_un   = compute_smd(data_before, v),
      smd_adj  = compute_smd(data_after,  v)
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
    geom_point(size = 4) +
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
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 13),
          axis.text.y     = element_text(size = 13),
          axis.text.x     = element_text(size = 12),
          axis.title.x    = element_text(size = 13),
          plot.title      = element_text(size = 16),
          plot.subtitle   = element_text(size = 13),
          plot.margin     = margin(10, 25, 10, 10))
}

# Add log-transformed columns to unmatched_pool so love plot SMDs
# are computed on the same scale as the matching (log-space)
unmatched_pool_log <- unmatched_pool |>
  mutate(
    log1p_road_length_km          = log1p(pmax(road_length_km,          0)),
    log1p_road_density_m_km2      = log1p(pmax(road_density_m_km2,      0)),
    log_area_km2                  = log(area_km2),
    log1p_dist_citycentre         = log1p(pmax(dist_citycentre,         0)),
    log1p_pop_density             = log1p(pmax(pop_density,             0)),
    log1p_business_retail_per_km2 = log1p(pmax(business_retail_per_km2, 0)),
    age_under15_pct = X4under_pct  + X5to9_pct   + X10to14_pct,
    age_15to24_pct  = X15to19_pct  + X20to24_pct,
    age_25to44_pct  = X25to29_pct  + X30to34_pct + X35to39_pct + X40to44_pct,
    age_45to64_pct  = X45to49_pct  + X50to54_pct + X55to59_pct + X60to64_pct,
    age_65plus_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct
  )

ld_s1 <- love_data_fn(unmatched_pool_log, matched_A, stage1_vars)
p_love_s1 <- make_love_plot(
  ld_s1,
  "Stage 1 Balance ",
  "Structural and sociodemographic variables | ratio 10, exact = country")
save_fig(p_love_s1, "fig01_love_plot_stage1_A.png", width = 16, height = 18)

ld_s2 <- love_data_fn(unmatched_pool, matched_A, c(stage2_trends, stage2_levels))
p_love_s2 <- make_love_plot(
  ld_s2,
  "Stage 2 Balance ",
  "Pre-treatment injury trends and levels | ratio 3, exact = country")
save_fig(p_love_s2, "fig02_love_plot_stage2_A.png", width = 16, height = 16)

cat("  Saved: fig01_love_plot_stage1_A.png\n")
cat("  Saved: fig02_love_plot_stage2_A.png\n\n")

# =============================================================================
# SECTION 9 — SMD HEATMAP
# =============================================================================

cat("=== Section 9: SMD heatmap ===\n")

specs_heat <- list(
  "Unmatched"     = unmatched_pool,
  "After Stage 1" = matched_A,
  "After Stage 2" = matched_A
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
    title    = "Absolute SMD Across SpecificationsA",
    subtitle = "Green < 0.10 (balanced) | Yellow = marginal | Red > 0.20 (imbalanced)",
    x = NULL, y = NULL,
    caption  = "Stage 1 variables: unmatched → after S1. Stage 2 variables: unmatched → after S2."
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0.3, "lines"))

save_fig(p_heatmap, "fig03_smd_heatmap.png", width = 14, height = 20)
cat("  Saved: fig03_smd_heatmap.png\n\n")

# =============================================================================
# SECTION 10 — COVARIATE DISTRIBUTION PLOTS
# =============================================================================

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
  
  lbl   <- coalesce(var_labels[var], var)
  x_lbl <- if (log_scale) paste0("log1p(", lbl, ")") else lbl
  
  ggplot(d_combined, aes(x = val, colour = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.7) +
    scale_colour_manual(values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
    scale_fill_manual(  values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
    facet_wrap(~phase) +
    labs(title = paste0(title_prefix, lbl), x = x_lbl,
         y = "Density", colour = NULL, fill = NULL) +
    theme_diag(base_size = 12) +
    theme(legend.position = "bottom")
}

# Raw variable names for distribution plots (use original columns, not log-transformed)
dist_road_vars   <- c("road_length_km", "road_density_m_km2", "area_km2",
                      "pct_A_road", "pct_B_road", "pct_minor_road")
dist_urban_vars  <- c("dist_citycentre", "pop_density", "business_retail_per_km2")
dist_socdem_vars <- c("IMD",
                      "cars_one_pct", "cars_twoPlus_pct",
                      "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
                      "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
                      "Taxi_pct", "workAthome_pct", "Other_pct",
                      "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
                      "X4under_pct", "X5to19_pct", "X20to24_pct", "X65plus_pct")
log_s1 <- c("road_length_km", "road_density_m_km2", "area_km2",
            "dist_citycentre", "pop_density", "business_retail_per_km2")

road_plots <- map(c(dist_road_vars, dist_urban_vars), function(v)
  plot_density_pair(v, unmatched_pool, matched_A, log_scale = v %in% log_s1))
p_dist_road <- wrap_plots(road_plots, ncol = 2) +
  plot_annotation(
    title    = "Road Network & Urban Variables — Treated vs Control",
    subtitle = "Left panel: before matching | Right panel: after matching)")
save_fig(p_dist_road, "fig04a_dist_road_network.png", width = 16, height = 22)

socdem_plots <- map(dist_socdem_vars, function(v)
  plot_density_pair(v, unmatched_pool, matched_A, log_scale = FALSE))
p_dist_socdem <- wrap_plots(socdem_plots, ncol = 3) +
  plot_annotation(
    title    = "Sociodemographic Variables — Treated vs Control",
    subtitle = "Left panel: before matching | Right panel: after matchin")
save_fig(p_dist_socdem, "fig04b_dist_socdem.png", width = 16, height = 24)

key_trends <- c("trend_total_pkm", "trend_ped_slight_pkm",
                "trend_car_slight_pkm", "trend_cyc_slight_pkm")
trend_plots <- map(key_trends, function(v)
  plot_density_pair(v, unmatched_pool, matched_A))
p_dist_trends <- wrap_plots(trend_plots, ncol = 2) +
  plot_annotation(
    title    = "Key Pre-Treatment Injury Trend Variables — Treated vs Matched Controls",
    subtitle = " Overlap of distributions supports the parallel trends assumption."
  )
save_fig(p_dist_trends, "fig05_dist_stage2_trends.png", width = 14, height = 10)

cat("  Saved: fig04a_dist_road_network.png\n")
cat("  Saved: fig04b_dist_socdem.png\n")
cat("  Saved: fig05_dist_stage2_trends.png\n\n")

# =============================================================================
# SECTION 11 — WEIGHT DISTRIBUTION PLOTS
# =============================================================================

cat("=== Section 11: Weight distribution plots ===\n")

ctrl_A     <- matched_A %>% filter(treat_indicator == 0)
eff_n_overall <- round(sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2), 0)

p_w1 <- ggplot(ctrl_A, aes(x = weights, fill = country)) +
  geom_histogram(bins = 60, alpha = 0.85, position = "stack") +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  scale_x_continuous(limits = c(0, 5.5)) +
  labs(
    title    = "Control OA Weight Distribution",
    subtitle = paste0("Nominal N = ", nrow(ctrl_A),
                      " | Effective N = ", eff_n_overall,
                      " | Efficiency = ",
                      round(eff_n_overall / nrow(ctrl_A), 3)),
    x = "Weight", y = "Count", fill = "Country"
  ) +
  theme_diag()

p_w2 <- ctrl_A %>%
  group_by(country) %>%
  summarise(n     = n(),
            eff_n = round(sum(weights)^2 / sum(weights^2), 1),
            max_w = round(max(weights), 2),
            .groups = "drop") %>%
  ggplot(aes(x = country, y = eff_n, fill = country)) +
  geom_col(width = 0.45) +
  geom_text(aes(label = paste0("n=", n, "\neff_N=", eff_n, "\nmax_w=", max_w)),
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
  labs(title    = "ECDF of weights by country",
       subtitle = "Red dashed = cap at 5",
       x = "Weight", y = "Cumulative proportion", colour = "Country") +
  theme_diag()

p_weights <- (p_w1 | p_w2 | p_w3) +
  plot_annotation(
    title = "Weight Diagnostics (after cap at 5)",
    theme = theme(plot.title = element_text(size = 13, face = "bold",
                                            colour = "#1A2E5A"))
  )
save_fig(p_weights, "fig06_weight_diagnostics.png", width = 16, height = 6)
cat("  Saved: fig06_weight_diagnostics.png\n\n")




control_reuse <- matched_A %>%
  filter(treat_indicator == 0) %>%
  count(OA, name = "times_used") %>%
  count(times_used, name = "n_controls") %>%
  mutate(pct = round(100 * n_controls / sum(n_controls), 1))




# =============================================================================
# SECTION 12 — MAHALANOBIS DISTANCE PLOTS
# =============================================================================

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
    title    = "ECDF of Stage 2 Mahalanobis Distance — Treated OAs",
    subtitle = "Distance = how far each treated OA is from its nearest control in trajectory space",
    x        = "Stage 2 Mahalanobis distance",
    y        = "Cumulative proportion of treated OAs",
    colour   = "Country",
    caption  = paste0("n = ", nrow(mdist_A), " treated OAs | ",
                      sum(mdist_A$mdist > 20), " OAs with distance > 20")
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
save_fig(p_mdist, "fig07_mahalanobis_distance.png", width = 12, height = 11)
cat("  Saved: fig07_mahalanobis_distance.png\n\n")

# =============================================================================
# SECTION 12b — SMD BY MAHALANOBIS DISTANCE QUARTILE
# Does residual imbalance concentrate in poorly-matched (high-distance) OAs?
# =============================================================================

cat("=== Section 12b: SMD by Mahalanobis distance quartile ===\n")

# Assign each treated OA to a distance quartile based on Stage 2 Mahalanobis
# distance (stored in matched_A from script 16)
quartile_labels <- c("Q1\n(best matched)", "Q2", "Q3", "Q4\n(worst matched)")

treated_quartiles <- matched_A |>
  filter(treat_indicator == 1, !is.na(mdist)) |>
  mutate(dist_quartile       = ntile(mdist, 4),
         dist_quartile_label = factor(quartile_labels[dist_quartile],
                                      levels = quartile_labels))

# Build treated→control linkage using saved OA-code pairs
# (pairs_mixed loaded at top of script from OA_matching_pairs_mixed.rds)
pair_list <- pairs_mixed |>
  left_join(
    treated_quartiles |> select(OA, dist_quartile, dist_quartile_label),
    by = c("treated_OA" = "OA")
  ) |>
  filter(!is.na(dist_quartile)) |>
  pivot_longer(c(treated_OA, control_OA),
               names_to = "role", values_to = "OA") |>
  select(OA, dist_quartile, dist_quartile_label) |>
  distinct()

matched_with_q <- matched_A |>
  inner_join(pair_list, by = "OA", relationship = "many-to-many") |>
  filter(!is.na(dist_quartile))

# Variables to assess — Stage 2 trends + key Stage 1 structural vars
smd_q_vars <- c(
  stage2_trends,
  "Drive_Car_pct", "Walk_pct", "log1p_pop_density",
  "log1p_road_density_m_km2", "log1p_business_retail_per_km2", "IMD"
)
smd_q_vars <- intersect(smd_q_vars, names(matched_A))

# Compute |SMD| within each quartile group (treated vs their matched controls)
smd_by_quartile <- map_df(1:4, function(q) {
  d <- matched_with_q |> filter(dist_quartile == q)
  map_df(smd_q_vars, function(v) {
    tibble(
      dist_quartile       = q,
      dist_quartile_label = factor(quartile_labels[q], levels = quartile_labels),
      variable            = v,
      label               = coalesce(var_labels[v], v),
      smd                 = abs(compute_smd(d, v))
    )
  })
}) |>
  mutate(var_type = if_else(variable %in% stage2_trends,
                            "Stage 2 trend variable", "Stage 1 structural variable"))

# Quartile n sizes for subtitle
q_sizes <- treated_quartiles |>
  count(dist_quartile_label) |>
  mutate(lbl = paste0(gsub("\n", " ", as.character(dist_quartile_label)),
                      " n=", n)) |>
  pull(lbl) |>
  paste(collapse = " | ")

p_smd_quartile <- ggplot(
  smd_by_quartile,
  aes(x = dist_quartile_label, y = smd, group = label, colour = label)
) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "turbo", name = "Variable") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ var_type, ncol = 1, scales = "free_y") +
  labs(
    title    = "Residual Imbalance by Mahalanobis Distance Quartile — Analysis A",
    subtitle = paste0("Each line = one variable | Dashed = |SMD| = 0.10 threshold\n",
                      q_sizes),
    x        = "Distance quartile (Stage 2 Mahalanobis distance)",
    y        = "Absolute SMD",
    caption  = paste0("Rising lines toward Q4 indicate that residual imbalance ",
                      "is concentrated in poorly-matched OAs.")
  ) +
  theme_diag() +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 10),
        legend.key.height = unit(0.6, "cm"),
        strip.text       = element_text(size = 13, face = "bold"))

save_fig(p_smd_quartile, "fig07b_smd_by_distance_quartile.png", width = 16, height = 14)
cat("  Saved: fig07b_smd_by_distance_quartile.png\n\n")

# =============================================================================
#  — STRATUM INJURY DISTRIBUTION
# =============================================================================

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
    subtitle = "Treated OAs  | log scale y-axis",
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

save_fig(p_strat_inj,  "fig09_stratum_injury.png",      width = 10, height = 7)
save_fig(p_strat_covs, "fig09b_stratum_covariates.png", width = 13, height = 10)
cat("  Saved: fig09_stratum_injury.png\n")
cat("  Saved: fig09b_stratum_covariates.png\n\n")

# =============================================================================
# — ENGLAND vs SCOTLAND COMPARISON
# =============================================================================

cat("=== Section 15: England vs Scotland comparison ===\n")

country_summary <- matched_A %>%
  filter(treat_indicator == 1) %>%
  group_by(country) %>%
  summarise(
    n_treated            = n(),
    pct_of_total         = round(100 * n() / nrow(matched_A_treated), 1),
    mean_road_length_km  = round(mean(road_length_km,     na.rm = TRUE), 3),
    mean_pop_density     = round(mean(pop_density,        na.rm = TRUE), 0),
    mean_dist_citycentre = round(mean(dist_citycentre,    na.rm = TRUE), 0),
    mean_IMD             = round(mean(IMD,                na.rm = TRUE), 1),
    mean_Drive_Car_pct   = round(mean(Drive_Car_pct,      na.rm = TRUE), 1),
    mean_Walk_pct        = round(mean(Walk_pct,           na.rm = TRUE), 1),
    mean_pct_minor       = round(mean(pct_minor_road,     na.rm = TRUE), 1),
    mean_total_pkm       = round(mean(mean_total_pkm,     na.rm = TRUE), 5),
    median_mdist         = round(median(mdist,            na.rm = TRUE), 2),
    pct_isolated         = round(100 * mean(OA %in% isolated_ids), 1),
    .groups              = "drop"
  )

write_csv(country_summary, file.path(outdir, "08_country_comparison.csv"))
cat("  Saved: 08_country_comparison.csv\n")

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
    subtitle = paste0("England n=613 | Scotland n=201"),
    x = NULL, y = NULL, fill = NULL,
    caption  = paste0("Scotland treated OAs tend to be longer roads, ",
                      "lower pop density, higher car commuting, lower pre-treatment injury rates.")
  ) +
  theme_diag() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

save_fig(p_country_profile, "fig10_england_scotland_profile.png", width = 16, height = 12)
cat("  Saved: fig10_england_scotland_profile.png\n\n")

# =============================================================================
# SECTION 16 — MAP PLOTS (spatial boundary required)
# =============================================================================

cat("=== Section 16: Map plots ===\n")

oa_path  <- here("data", "spatial", "OA_boundaries_GB.gpkg")
lad_path <- here("data", "spatial", "LAD_boundaries_GB.gpkg")

if (!file.exists(oa_path) || !file.exists(lad_path)) {
  cat("  SKIPPED — boundary files not found.\n")
  cat("  Expected paths:\n")
  cat("    OA: ", oa_path,  "\n")
  cat("    LAD:", lad_path, "\n\n")
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
      buffer_OA  == 1                                  ~ "Buffer (excluded)",
      treated_OA == 1                                  ~ "Treated",
      OA %in% controls_A$OA                           ~ "Matched control",
      control_group1_OA == 1 | control_group2_OA == 1 ~ "Eligible (unmatched)",
      TRUE                                             ~ "Other"
    )) %>%
    select(OA, map_group, LAD24CD, country)
  
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
    labs(title    = "Geographic Distribution of Treated OAs",
         subtitle = paste0(nrow(treated_centroids),
                           " treated OAs (England n=613 | Scotland n=201)"),
         caption  = "Each point = one treated OA centroid",
         colour   = "Country") +
    theme_diag() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          panel.grid = element_line(colour = "#E8ECF3", linewidth = 0.3),
          legend.position = "bottom")
  
  save_fig(p_overview, "fig11_map_treated_overview.png", width = 8, height = 12)
  
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
         caption  = " square-root fill scale") +
    theme_diag() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          legend.position = "right")
  
  save_fig(p_control_origin, "fig12_map_control_origin.png", width = 9, height = 13)
  
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
    
    city_oas  <- full_data %>% filter(LAD24CD == lad_cd) %>% pull(OA)
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
    
    n_t <- sum(nearby$map_group == "Treated",         na.rm = TRUE)
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
    
    fn <- paste0("fig13_map_city_", gsub("[^A-Za-z0-9]", "_", lad_cd), ".png")
    save_fig(p_city, fn, width = 10, height = 10)
  })
  
  cat("  Map figures saved.\n\n")
}

# =============================================================================
# SECTION 17 — PARALLEL TRENDS (TREND SLOPE DISTRIBUTIONS)
# =============================================================================

cat("=== Section 17: Parallel trends visualisation (slope distributions) ===\n")

trend_long <- matched_A %>%
  filter(treat_indicator %in% c(0, 1)) %>%
  mutate(group = if_else(treat_indicator == 1, "Treated", "Matched control")) %>%
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
                                 "Matched control" = COL_CONTROL)) +
  scale_fill_manual(  values = c("Treated" = COL_TREATED,
                                 "Matched control" = COL_CONTROL)) +
  facet_wrap(~trend_label, scales = "free", ncol = 3) +
  labs(
    title    = "Pre-Treatment Injury Trend Slopes: Treated vs Matched Controls",
    subtitle = " Overlap of distributions supports the parallel trends assumption",
    x        = "Pre-treatment slope (log-linear regression coefficient)",
    y        = "Density",
    colour   = NULL, fill = NULL,
    caption  = paste0(
      "n treated = ", sum(matched_A$treat_indicator == 1),
      " | n controls = ", sum(matched_A$treat_indicator == 0),
      "\nDashed line = slope of 0 (flat trend). ",
      "Note: slope distributions are a proxy — also produce actual time series in the DiD script."
    )
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p_trend_density, "fig14_parallel_trends_slopes.png", width = 14, height = 12)
cat("  Saved: fig14_parallel_trends_slopes.png\n\n")



cat("SAMPLE SIZES:\n")
cat(sprintf("  Treated:  %d\n", sum(matched_A$treat_indicator == 1)))
cat(sprintf("  Controls: %d\n", sum(matched_A$treat_indicator == 0)))
cat(sprintf("  Isolated OAs (flagged, included in this analysis): %d\n",
            length(isolated_ids)))

cat("\nWEIGHT SUMMARY (controls, after cap = 5):\n")
ctrl_A <- matched_A %>% filter(treat_indicator == 0)
cat(sprintf("  Nominal N:   %d\n", nrow(ctrl_A)))
cat(sprintf("  Effective N: %.0f\n",
            sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2)))
cat(sprintf("  Efficiency:  %.3f\n",
            (sum(ctrl_A$weights)^2 / sum(ctrl_A$weights^2)) / nrow(ctrl_A)))
cat(sprintf("  Max weight:  %.3f\n", max(ctrl_A$weights)))

cat("\nBALANCE SUMMARY:\n")
cat(sprintf("  Stage 1 — mean |SMD| before: %.3f | after: %.3f\n",
            mean(abs(smd_s1$smd_unmatched), na.rm = TRUE),
            mean(abs(smd_s1$smd_after_S1),  na.rm = TRUE)))
cat(sprintf("  Stage 2 — mean |SMD|: %.3f | max trend SMD: %.3f\n",
            mean(abs(smd_s2$smd_postS2), na.rm = TRUE),
            max(abs(smd_s2$smd_postS2[smd_s2$var_type == "Trend"]), na.rm = TRUE)))

cat("\nFIGURES SAVED TO:", outdir, "\n")
cat("  fig01  = Stage 1 love plot\n")
cat("  fig02  = Stage 2 love plot\n")
cat("  fig03  = SMD heatmap\n")
cat("  fig04a = Road network distributions\n")
cat("  fig04b = Sociodemographic distributions\n")
cat("  fig05  = Stage 2 trend distributions\n")
cat("  fig06  = Weight diagnostics\n")
cat("  fig07  = Mahalanobis distance ECDF\n")
cat("  fig08  = Stratum injury distribution\n")
cat("  fig09 = Stratum covariate profiles\n")
cat("  fig10  = England vs Scotland profile\n")
cat("  fig11-13 = Maps (if boundaries available)\n")
cat("  fig14  = Parallel trends slope distributions\n")

cat("\nTABLES SAVED TO:", outdir, "\n")
cat("  01  = Descriptive Table 1 (Treated | Unmatched ctrl | Matched ctrl)\n")
cat("  01b = Formatted mean (SD) / median [IQR]\n")
cat("  01c = Country-stratified descriptives\n")
cat("  02  = Stage 1 SMD before/after\n")
cat("  03  = Stage 2 SMD before/after\n")
cat("  04  = Weight distribution by country\n")
cat("  05  = Stratum characteristics\n")
cat("  06a = Isolated OA characteristics\n")
cat("  06b = Isolated OA SMD\n")
cat("  07  = Skewness assessment\n")
cat("  08  = England vs Scotland comparison\n")
cat("=================================================================\n")