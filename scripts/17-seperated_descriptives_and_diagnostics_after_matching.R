# =============================================================================
# POST-MATCHING DIAGNOSTICS — ENGLAND / SCOTLAND SEPARATED
# OA-Level Two-Stage Mahalanobis Distance Matching
#
# INPUTS (from data/processed/):
#   OA_matched_full_mixed.rds               — full matched dataset (England + Scotland)
#   OA_matched_full_mixed_England.rds       — England matched dataset
#   OA_matched_full_mixed_Scotland.rds      — Scotland matched dataset
#   OA_matched_treated_mixed.rds            — treated OAs + weights + stratum
#   OA_matched_donors_mixed.rds             — control OAs + weights
#   OA_common_support_flags_mixed.rds       — structurally isolated OA flags
#   OA_matching_pairs_mixed.rds             — treated->control OA pairs
#   OA_matching_census.rds                  — full original dataset (unmatched pool)
#
# SOURCE: 16_Matching_England_othercityControlsScotland_mix.R
#   England: other-city controls only (ratio 1:2)
#   Scotland: other-city + same-city controls (ratio 1:1)
#
# CHANGE LOG:
#   Stage 2 love plots and SMD tables now computed on the MATCHING scale:
#     - Injury trends: raw scale (as matched)
#     - Injury levels: log1p scale (as matched)
#   Previously raw-scale levels were shown, producing inflated post-match SMDs
#   due to right-skew. Log-scale SMDs reflect what the matcher actually achieved.
#
# OUTPUTS (to output/diagnostics/):
#   Combined  : pooled results
#   England   : suffix _england
#   Scotland  : suffix _scotland
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
      plot.title        = element_text(size = base_size + 2, face = "bold",
                                       colour = "#1A2E5A", margin = margin(b = 8)),
      plot.subtitle     = element_text(size = base_size - 1, colour = "#555555",
                                       margin = margin(b = 12)),
      plot.caption      = element_text(size = base_size - 2, colour = "#888888",
                                       hjust = 0, margin = margin(t = 8)),
      axis.title        = element_text(size = base_size - 1, colour = "#333333"),
      axis.text         = element_text(size = base_size - 2, colour = "#444444"),
      strip.text        = element_text(size = base_size - 1, face = "bold",
                                       colour = "#1A2E5A"),
      strip.background  = element_rect(fill = "#EEF2F8", colour = NA),
      panel.grid.major  = element_line(colour = "#E5E9F0", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "#CCCCCC", fill = NA,
                                       linewidth = 0.4),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      legend.text       = element_text(size = base_size - 2),
      legend.background = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      plot.margin       = margin(12, 16, 12, 12)
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



matched_A        <- readRDS(here("data", "processed", "OA_matched_full_mixed.rds"))
treated_A        <- readRDS(here("data", "processed", "OA_matched_treated_mixed.rds"))
controls_A       <- readRDS(here("data", "processed", "OA_matched_donors_mixed.rds"))
csupport         <- readRDS(here("data", "processed", "OA_common_support_flags_mixed.rds"))
pairs_mixed      <- readRDS(here("data", "processed", "OA_matching_pairs_mixed.rds"))
full_data        <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Per-scheme matching: England only. matched_eng = matched_A (no Scotland).
matched_eng <- matched_A
matched_sco <- matched_A %>% filter(FALSE)   # empty — Scotland removed
has_scotland <- nrow(matched_sco) > 0

cat("=== Sample sizes ===\n")
cat("  Treated:", sum(matched_A$treat_indicator == 1),
    "| Controls:", sum(matched_A$treat_indicator == 0), "\n\n")

# Controls
ctrl_eng <- matched_eng %>% filter(treat_indicator == 0)
ctrl_sco <- matched_sco %>% filter(treat_indicator == 0)

# Control → scheme lookup (controls have scheme == "Control" in matched_full;
# per-scheme matching pairs carry the actual scheme)
ctrl_scheme_lookup <- pairs_mixed %>%
  select(OA = control_OA, scheme) %>%
  distinct()

# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

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
                   "age_45to64_pct", "age_65to84_pct")



stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

# Stage 2: 4 trends + 4 levels (KSI + slight combined per mode; Other in total only)
stage2_trends <- c(
  "trend_car_pkm", "trend_cyc_pkm", "trend_ped_pkm", "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_pkm", "mean_cyc_pkm", "mean_ped_pkm", "mean_total_pkm"
)


# Raw Stage 1 names — for descriptive tables only (human-readable scale)
stage1_road_raw   <- c("road_length_km", "road_density_m_km2", "area_km2",
                       "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban_raw  <- c("dist_citycentre", "pop_density", "business_retail_per_km2")
stage1_socdem_raw <- c("IMD",
                       "cars_one_pct", "cars_twoPlus_pct",
                       "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
                       "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
                       "Taxi_pct", "workAthome_pct", "Other_pct",
                       "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
                       "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
                       "age_45to64_pct", "age_65to84_pct")
stage1_vars_raw <- c(stage1_road_raw, stage1_urban_raw, stage1_socdem_raw)

all_desc_vars <- c(stage1_vars_raw, stage2_trends, stage2_levels)



# Log-scale versions of level variables — these match what the matcher used
stage2_levels_log <- paste0("log1p_", stage2_levels)

# =============================================================================
# VAR LABELS
# =============================================================================

var_labels <- c(
  # Stage 1 — road
  log1p_road_length_km                 = "Road length (km)",
  log1p_road_density_m_km2             = "Road density (m/km\u00b2)",
  log_area_km2                         = "Area (km\u00b2)",
  pct_A_road                           = "% A-road",
  pct_B_road                           = "% B-road",
  pct_minor_road                       = "% Minor road",
  # Stage 1 — urban
  log1p_dist_citycentre                = "Distance to city centre (m)",
  log1p_pop_density                    = "Population density (persons/km\u00b2)",
  log1p_business_retail_per_km2        = "Retail businesses (per km\u00b2)",
  # Stage 1 — raw (for desc tables)
  road_density_m_km2                   = "Road density (m/km\u00b2)",
  road_length_km                       = "Road length (km)",
  dist_citycentre                      = "Distance to city centre (m)",
  pop_density                          = "Population density (persons/km\u00b2)",
  area_km2                             = "Area (km\u00b2)",
  business_retail_per_km2              = "Retail businesses (per km\u00b2)",
  # Stage 1 — sociodemographic
  IMD                                  = "Index of Multiple Deprivation",
  cars_one_pct                         = "% households: 1 car",
  cars_twoPlus_pct                     = "% households: 2+ cars",
  Drive_Car_pct                        = "% commuting: drive car",
  Passenger_Car_pct                    = "% commuting: car passenger",
  Walk_pct                             = "% commuting: walk",
  Bicycle_pct                          = "% commuting: bicycle",
  Bus_Coach_pct                        = "% commuting: bus/coach",
  Train_pct                            = "% commuting: train",
  Underground_train_tram_pct           = "% commuting: underground/tram",
  Taxi_pct                             = "% commuting: taxi",
  workAthome_pct                       = "% working from home",
  Other_pct                            = "% commuting: other mode",
  White_pct                            = "% White",
  Mixed_pct                            = "% Mixed ethnicity",
  Asian_pct                            = "% Asian",
  Black_pct                            = "% Black",
  age_under15_pct                      = "% aged under 15",
  age_15to24_pct                       = "% aged 15\u201324",
  age_25to44_pct                       = "% aged 25\u201344",
  age_45to64_pct                       = "% aged 45\u201364",
  age_65to84_pct                       = "% aged 65\u201384",
  X4under_pct                          = "% aged under 5",
  X5to19_pct                           = "% aged 5\u201319",
  X20to24_pct                          = "% aged 20\u201324",
  X65to84_pct                          = "% aged 65\u201324",
  # Stage 2 — trends (raw scale)
  # Stage 2 — trends (KSI + slight combined per mode)
  trend_car_pkm                        = "Trend: car/km",
  trend_cyc_pkm                        = "Trend: cycling/km",
  trend_ped_pkm                        = "Trend: pedestrian/km",
  trend_total_pkm                      = "Trend: total injuries/km",
  # Stage 2 — levels (raw scale, for descriptive tables only)
  mean_car_pkm                         = "Mean: car/km",
  mean_cyc_pkm                         = "Mean: cycling/km",
  mean_ped_pkm                         = "Mean: pedestrian/km",
  mean_total_pkm                       = "Mean: total injuries/km",
  # Stage 2 — levels (log1p scale, for balance diagnostics — matching scale)
  log1p_mean_car_pkm                   = "Mean: car/km (log)",
  log1p_mean_cyc_pkm                   = "Mean: cycling/km (log)",
  log1p_mean_ped_pkm                   = "Mean: pedestrian/km (log)",
  log1p_mean_total_pkm                 = "Mean: total injuries/km (log)"
)

# Descriptor vars use raw levels (for summary tables)
all_desc_vars <- c(stage1_vars, stage2_trends, stage2_levels)

# =============================================================================
# UNMATCHED POOL (with derived variables)
# =============================================================================

# add_log_vars: produces log-transformed Stage 1 vars AND log1p level vars
# so that unmatched_*_log can be used as the "before" reference for Stage 2
# balance diagnostics on the matching scale.
add_log_vars <- function(df) {
  df %>% mutate(
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
    age_65to84_pct  = X65to69_pct  + X70to74_pct + X75to79_pct + X80to84_pct,
    # Log1p-transformed injury levels — matching scale
    # Combined _pkm vars created upstream in script 12 (Other in total only)
    log1p_mean_car_pkm            = log1p(pmax(mean_car_pkm,            0)),
    log1p_mean_cyc_pkm            = log1p(pmax(mean_cyc_pkm,            0)),
    log1p_mean_ped_pkm            = log1p(pmax(mean_ped_pkm,            0)),
    log1p_mean_total_pkm          = log1p(pmax(mean_total_pkm,          0))
  )
}

unmatched_pool <- full_data %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(
    treat_indicator = as.integer(treated_OA == 1),
    country = if_else(country == "EnglandWales", "England", country)
  )

unmatched_pool_log <- add_log_vars(unmatched_pool)
unmatched_eng_log  <- unmatched_pool_log %>% filter(country == "England")
unmatched_sco_log  <- unmatched_pool_log %>% filter(country == "Scotland")

# Raw unmatched (for descriptive tables — raw scale only)
unmatched_eng <- unmatched_pool %>% filter(country == "England")
unmatched_sco <- unmatched_pool %>% filter(country == "Scotland")

# =============================================================================
# HELPERS
# =============================================================================

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

desc_stats <- function(data, vars, group_label) {
  map_df(vars, function(v) {
    x <- data[[v]]
    if (is.null(x)) return(NULL)
    tibble(
      variable = v,
      label    = coalesce(var_labels[v], v),
      group    = group_label,
      n        = sum(!is.na(x)),
      mean     = round(mean(x,           na.rm = TRUE), 4),
      sd       = round(sd(x,             na.rm = TRUE), 4),
      median   = round(median(x,         na.rm = TRUE), 4),
      p25      = round(quantile(x, 0.25, na.rm = TRUE), 4),
      p75      = round(quantile(x, 0.75, na.rm = TRUE), 4)
    )
  })
}

smd_table <- function(vars, data_before, data_after) {
  map_df(vars, function(v) {
    tibble(
      variable   = v,
      label      = coalesce(var_labels[v], v),
      smd_before = round(compute_smd(data_before, v), 4),
      smd_after  = round(compute_smd(data_after,  v), 4)
    )
  })
}

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
  
  ggplot(ldat_long,
         aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = threshold, linetype = "dashed",
               colour = "#999999", linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.5) +
    geom_point(size = 4) +
    scale_colour_manual(
      values = c("Before matching" = COL_BEFORE, "After matching" = COL_AFTER)) +
    scale_shape_manual(
      values = c("Before matching" = 16, "After matching" = 17)) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = subtitle,
         x = "Absolute Standardised Mean Difference", y = NULL,
         colour = NULL, shape = NULL,
         caption = "Dashed line = |SMD| = 0.10 threshold") +
    theme_diag() +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 13),
          axis.text.y      = element_text(size = 13),
          axis.text.x      = element_text(size = 12),
          axis.title.x     = element_text(size = 13),
          plot.title       = element_text(size = 16),
          plot.subtitle    = element_text(size = 13),
          plot.margin      = margin(10, 25, 10, 10))
}

# =============================================================================
# SECTION 2 — DESCRIPTIVE SUMMARY TABLES (raw scale; combined + by country)
# =============================================================================

cat("=== Section 2: Descriptive summary tables ===\n")

# --- 2a: Combined treated vs matched controls --------------------------------
desc_table_combined <- bind_rows(
  desc_stats(unmatched_pool %>% filter(treat_indicator == 1), all_desc_vars, "Treated"),
  desc_stats(matched_A      %>% filter(treat_indicator == 0), all_desc_vars, "Matched")
) %>%
  pivot_wider(id_cols = c(variable, label), names_from = group,
              values_from = c(mean, sd, median), names_glue = "{group}_{.value}") %>%
  mutate(var_group = case_when(
    variable %in% stage2_trends ~ "Injury trends",
    variable %in% stage2_levels ~ "Injury levels",
    TRUE                        ~ "Demographic & road characteristics"
  )) %>%
  select(var_group, label, Treated_mean, Matched_mean,
         Treated_sd, Matched_sd, Treated_median, Matched_median) %>%
  arrange(var_group, label)

write_csv(desc_table_combined,
          file.path(outdir, "01_descriptive_table_combined.csv"))
cat("  Saved: 01_descriptive_table_combined.csv\n")

# --- 2b: England -----------------------------------------------------------
desc_table_eng <- bind_rows(
  desc_stats(unmatched_eng_log %>% filter(treat_indicator == 1),
             all_desc_vars, "Treated"),
  desc_stats(matched_eng %>% filter(treat_indicator == 0),
             all_desc_vars, "Matched")
) %>%
  pivot_wider(id_cols = c(variable, label), names_from = group,
              values_from = c(mean, sd, median), names_glue = "{group}_{.value}") %>%
  mutate(var_group = case_when(
    variable %in% stage2_trends ~ "Injury trends",
    variable %in% stage2_levels ~ "Injury levels",
    TRUE                        ~ "Demographic & road characteristics"
  )) %>%
  select(var_group, label, Treated_mean, Matched_mean,
         Treated_sd, Matched_sd, Treated_median, Matched_median) %>%
  arrange(var_group, label)

write_csv(desc_table_eng,
          file.path(outdir, "01_descriptive_table_england.csv"))
cat("  Saved: 01_descriptive_table_england.csv\n")

# --- 2c: Scotland (skipped if no Scotland data) ----------------------------
if (has_scotland) {
  desc_table_sco <- bind_rows(
    desc_stats(unmatched_sco_log %>% filter(treat_indicator == 1),
               all_desc_vars, "Treated"),
    desc_stats(matched_sco %>% filter(treat_indicator == 0),
               all_desc_vars, "Matched")
  ) %>%
    pivot_wider(id_cols = c(variable, label), names_from = group,
                values_from = c(mean, sd, median), names_glue = "{group}_{.value}") %>%
    mutate(var_group = case_when(
      variable %in% stage2_trends ~ "Injury trends",
      variable %in% stage2_levels ~ "Injury levels",
      TRUE                        ~ "Demographic & road characteristics"
    )) %>%
    select(var_group, label, Treated_mean, Matched_mean,
           Treated_sd, Matched_sd, Treated_median, Matched_median) %>%
    arrange(var_group, label)
  write_csv(desc_table_sco, file.path(outdir, "01_descriptive_table_scotland.csv"))
  cat("  Saved: 01_descriptive_table_scotland.csv\n")
}
cat("\n")

# =============================================================================
# SECTION 3 — SMD TABLES
#
# Stage 1: log-scale structural vars (as matched)
# Stage 2: raw trends + log1p levels (as matched)
# NOTE: descriptive tables above remain on raw scale for interpretability.
# =============================================================================

cat("=== Section 3: SMD before/after tables ===\n")

# --- Stage 1: combined -------------------------------------------------------
smd_s1 <- smd_table(stage1_vars, unmatched_pool_log, matched_A) %>%
  rename(smd_unmatched = smd_before, smd_after_S1 = smd_after) %>%
  mutate(balanced_after_S1 = abs(smd_after_S1) < 0.1)
write_csv(smd_s1, file.path(outdir, "02_smd_stage1_combined.csv"))
cat("  Saved: 02_smd_stage1_combined.csv\n")

# --- Stage 1: England --------------------------------------------------------
smd_s1_eng <- smd_table(stage1_vars, unmatched_eng_log, matched_eng) %>%
  rename(smd_unmatched = smd_before, smd_after_S1 = smd_after) %>%
  mutate(balanced_after_S1 = abs(smd_after_S1) < 0.1)
write_csv(smd_s1_eng, file.path(outdir, "02_smd_stage1_england.csv"))
cat("  Saved: 02_smd_stage1_england.csv\n")

if (has_scotland) {
  smd_s1_sco <- smd_table(stage1_vars, unmatched_sco_log, matched_sco) %>%
    rename(smd_unmatched = smd_before, smd_after_S1 = smd_after) %>%
    mutate(balanced_after_S1 = abs(smd_after_S1) < 0.1)
  write_csv(smd_s1_sco, file.path(outdir, "02_smd_stage1_scotland.csv"))
  cat("  Saved: 02_smd_stage1_scotland.csv\n")
} else { smd_s1_sco <- tibble() }

# --- Stage 2: combined (log-scale levels + raw trends) -----------------------
smd_s2 <- smd_table(
  c(stage2_trends, stage2_levels_log),
  unmatched_pool_log, matched_A
) %>%
  rename(smd_preS2 = smd_before, smd_postS2 = smd_after) %>%
  mutate(
    var_type = if_else(variable %in% stage2_trends, "Trend", "Level (log)"),
    balanced = abs(smd_postS2) < 0.1,
    note     = if_else(
      variable %in% stage2_levels_log,
      "log1p scale — matches optimisation scale",
      "raw scale"
    )
  )
write_csv(smd_s2, file.path(outdir, "03_smd_stage2_combined.csv"))
cat("  Saved: 03_smd_stage2_combined.csv\n")

# --- Stage 2: England --------------------------------------------------------
smd_s2_eng <- smd_table(
  c(stage2_trends, stage2_levels_log),
  unmatched_eng_log, matched_eng
) %>%
  rename(smd_preS2 = smd_before, smd_postS2 = smd_after) %>%
  mutate(
    var_type = if_else(variable %in% stage2_trends, "Trend", "Level (log)"),
    balanced = abs(smd_postS2) < 0.1,
    note     = if_else(
      variable %in% stage2_levels_log,
      "log1p scale — matches optimisation scale",
      "raw scale"
    )
  )
write_csv(smd_s2_eng, file.path(outdir, "03_smd_stage2_england.csv"))
cat("  Saved: 03_smd_stage2_england.csv\n")

if (has_scotland) {
  smd_s2_sco <- smd_table(
    c(stage2_trends, stage2_levels_log),
    unmatched_sco_log, matched_sco
  ) %>%
    rename(smd_preS2 = smd_before, smd_postS2 = smd_after) %>%
    mutate(
      var_type = if_else(variable %in% stage2_trends, "Trend", "Level (log)"),
      balanced = abs(smd_postS2) < 0.1,
      note     = if_else(
        variable %in% stage2_levels_log,
        "log1p scale — matches optimisation scale",
        "raw scale"
      )
    )
  write_csv(smd_s2_sco, file.path(outdir, "03_smd_stage2_scotland.csv"))
  cat("  Saved: 03_smd_stage2_scotland.csv\n")
} else { smd_s2_sco <- tibble() }
cat("\n")

# =============================================================================
# SECTION 4 — WEIGHT DISTRIBUTION + CONCENTRATION DIAGNOSTICS (by country)
# =============================================================================

cat("=== Section 4: Weight diagnostics ===\n")

ctrl_A            <- matched_A %>% filter(treat_indicator == 0)
total_ctrl_weight <- sum(ctrl_A$weights)

weight_summary_fn <- function(ctrl_df, label) {
  ctrl_df %>%
    summarise(
      country       = label,
      n_controls    = n(),
      mean_weight   = round(mean(weights),           3),
      median_weight = round(median(weights),         3),
      sd_weight     = round(sd(weights),             3),
      p90_weight    = round(quantile(weights, 0.90), 3),
      p95_weight    = round(quantile(weights, 0.95), 3),
      max_weight    = round(max(weights),            3),
      eff_N         = round(sum(weights)^2 / sum(weights^2), 1),
      efficiency    = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
      n_at_cap      = sum(weights >= 5),
      pct_at_cap    = round(100 * mean(weights >= 5), 1)
    )
}

weight_table <- bind_rows(
  weight_summary_fn(ctrl_eng, "England"),
  if (has_scotland) weight_summary_fn(ctrl_sco, "Scotland") else NULL,
  weight_summary_fn(ctrl_A,   "Overall")
)
write_csv(weight_table, file.path(outdir, "04_weight_distribution.csv"))
cat("  Saved: 04_weight_distribution.csv\n")

reuse_fn <- function(ctrl_df) {
  ctrl_df %>%
    count(OA, name = "times_used") %>%
    count(times_used, name = "n_controls") %>%
    mutate(pct = round(100 * n_controls / sum(n_controls), 1))
}
write_csv(reuse_fn(ctrl_eng), file.path(outdir, "04b_control_reuse_england.csv"))
if (has_scotland) write_csv(reuse_fn(ctrl_sco), file.path(outdir, "04b_control_reuse_scotland.csv"))
cat("  Saved: 04b_control_reuse_england.csv\n")

top_ctrl_fn <- function(ctrl_df, total_w) {
  ctrl_df %>%
    select(OA, weights, country, LAD24CD) %>%
    arrange(desc(weights)) %>%
    slice_head(n = 20) %>%
    mutate(rank = row_number(),
           cumulative_weight_share_pct = round(100 * cumsum(weights) / total_w, 1))
}
write_csv(top_ctrl_fn(ctrl_eng, sum(ctrl_eng$weights)),
          file.path(outdir, "04c_top_weight_controls_england.csv"))
if (has_scotland) {
  write_csv(top_ctrl_fn(ctrl_sco, sum(ctrl_sco$weights)),
            file.path(outdir, "04c_top_weight_controls_scotland.csv"))
}
cat("  Saved: 04c_top_weight_controls_england.csv\n\n")

# =============================================================================
# SECTION 5 — STRATUM CHARACTERISTICS
# =============================================================================

cat("=== Section 5: Stratum characteristics ===\n")

matched_A_treated <- matched_A %>% filter(treat_indicator == 1)

stratum_table <- matched_A_treated %>%
  filter(!is.na(baseline_injury_stratum)) %>%
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
    mean_X65to84_pct      = round(mean(age_65to84_pct,        na.rm = TRUE), 1),
    pct_England           = round(100 * mean(country == "England"),        1),
    pct_Scotland          = round(100 * mean(country == "Scotland"),       1),
    .groups               = "drop"
  )

write_csv(stratum_table, file.path(outdir, "05_stratum_characteristics.csv"))
cat("  Saved: 05_stratum_characteristics.csv\n\n")

# =============================================================================
# SECTION 6 — COMMON SUPPORT FLAGS
# =============================================================================

cat("=== Section 6: Common support — isolated OA characteristics ===\n")

isolated_ids <- csupport$treated_OA

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
    pct_scotland         = round(100 * mean(country == "Scotland"),       1),
    .groups              = "drop"
  ) %>%
  mutate(isolated = if_else(isolated, "Isolated", "Non-isolated"))

matched_A_treated_flag <- matched_A_treated %>%
  mutate(treat_indicator = as.integer(OA %in% isolated_ids))

smd_isolated <- map_df(stage1_vars, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd = round(compute_smd(matched_A_treated_flag, v), 3))
}) %>% arrange(desc(abs(smd)))

write_csv(isolated_chars, file.path(outdir, "06a_isolated_OA_characteristics.csv"))
write_csv(smd_isolated,   file.path(outdir, "06b_isolated_OA_smd.csv"))
cat("  Saved: 06a/06b isolated OA tables\n\n")

# =============================================================================
# SECTION 8 — LOVE PLOTS
#
# Stage 1: log-scale structural vars (as matched)
# Stage 2: raw trends + log1p levels (as matched)
#   Subtitle explicitly states the scale so readers are not misled.
# =============================================================================

cat("=== Section 8: Love plots ===\n")

# --- Stage 1 love plots ------------------------------------------------------
ld_s1_combined <- love_data_fn(unmatched_pool_log, matched_A,   stage1_vars)
ld_s1_eng      <- love_data_fn(unmatched_eng_log,  matched_eng, stage1_vars)
if (has_scotland) ld_s1_sco <- love_data_fn(unmatched_sco_log, matched_sco, stage1_vars)

p_love_s1_combined <- make_love_plot(
  ld_s1_combined,
  "Stage 1 Balance \u2014 Combined",
  "Structural & sociodemographic variables (log-scale where matched) | England 1:2, Scotland 1:1")
p_love_s1_eng <- make_love_plot(
  ld_s1_eng,
  "Stage 1 Balance \u2014 England",
  "Structural & sociodemographic variables (log-scale where matched) | ratio 1:2, other-city controls only")
save_fig(p_love_s1_combined, "fig01_love_plot_stage1_combined.png", width = 16, height = 18)
save_fig(p_love_s1_eng,      "fig01_love_plot_stage1_england.png",  width = 16, height = 18)
if (has_scotland) {
  p_love_s1_sco <- make_love_plot(
    ld_s1_sco,
    "Stage 1 Balance \u2014 Scotland",
    "Structural & sociodemographic variables (log-scale where matched) | ratio 1:1, other-city + same-city controls")
  save_fig(p_love_s1_sco, "fig01_love_plot_stage1_scotland.png", width = 16, height = 18)
}
cat("  Saved: fig01 love plots\n")

# --- Stage 2 love plots — ON MATCHING SCALE ----------------------------------
# Trends: raw (as matched)
# Levels: log1p (as matched)
# Using unmatched_*_log as the "before" reference so both timing points
# are on the same scale and the comparison is valid.

ld_s2_combined <- love_data_fn(unmatched_pool_log, matched_A,
                               c(stage2_trends, stage2_levels_log))
ld_s2_eng      <- love_data_fn(unmatched_eng_log,  matched_eng,
                               c(stage2_trends, stage2_levels_log))
if (has_scotland) ld_s2_sco <- love_data_fn(unmatched_sco_log, matched_sco,
                                            c(stage2_trends, stage2_levels_log))

p_love_s2_combined <- make_love_plot(
  ld_s2_combined,
  "Stage 2 Balance \u2014 Combined",
  "Injury trends (raw) + mean levels (log1p) \u2014 both on matching scale")
p_love_s2_eng <- make_love_plot(
  ld_s2_eng,
  "Stage 2 Balance \u2014 England",
  "Injury trends (raw) + mean levels (log1p) \u2014 both on matching scale | ratio 1:2")
save_fig(p_love_s2_combined, "fig02_love_plot_stage2_combined.png", width = 16, height = 16)
save_fig(p_love_s2_eng,      "fig02_love_plot_stage2_england.png",  width = 16, height = 16)
if (has_scotland) {
  p_love_s2_sco <- make_love_plot(
    ld_s2_sco,
    "Stage 2 Balance \u2014 Scotland",
    "Injury trends (raw) + mean levels (log1p) \u2014 both on matching scale | ratio 1:1")
  save_fig(p_love_s2_sco, "fig02_love_plot_stage2_scotland.png", width = 16, height = 16)
}
cat("  Saved: fig02 love plots\n\n")

# =============================================================================
# SECTION 9 — SMD HEATMAP (by country)
#
# Stage 1 vars: log-scale (matching scale)
# Stage 2 trends: raw (matching scale)
# Stage 2 levels: log1p (matching scale)
# =============================================================================

cat("=== Section 9: SMD heatmap ===\n")

make_heatmap_data <- function(matched_df, unmatched_log_df, label) {
  map_df(c(stage1_vars, stage2_trends, stage2_levels_log), function(v) {
    tibble(
      country  = label,
      variable = v,
      smd_un   = abs(tryCatch(compute_smd(unmatched_log_df, v),
                              error = function(e) NA_real_)),
      smd_adj  = abs(tryCatch(compute_smd(matched_df, v),
                              error = function(e) NA_real_))
    )
  })
}

heatmap_long <- bind_rows(
  make_heatmap_data(matched_eng, unmatched_eng_log, "England"),
  if (has_scotland) make_heatmap_data(matched_sco, unmatched_sco_log, "Scotland") else NULL
) %>%
  pivot_longer(c(smd_un, smd_adj), names_to = "timing", values_to = "smd") %>%
  mutate(
    label     = coalesce(var_labels[variable], variable),
    var_group = case_when(
      variable %in% stage1_road        ~ "1. Road network",
      variable %in% stage1_urban       ~ "2. Urban geography",
      variable %in% stage1_socdem      ~ "3. Sociodemographic",
      variable %in% stage2_trends      ~ "4. Injury trends (raw)",
      variable %in% stage2_levels_log  ~ "5. Injury levels (log1p)"
    ),
    timing  = if_else(timing == "smd_un", "Unmatched", "After matching"),
    timing  = factor(timing, levels = c("Unmatched", "After matching")),
    country = factor(country, levels = c("England", "Scotland"))
  )

p_heatmap <- ggplot(heatmap_long, aes(x = timing, y = label, fill = smd)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(!is.na(smd), sprintf("%.2f", smd), "")),
            size = 2.3, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#2ECC71", mid = "#F39C12", high = "#E74C3C",
    midpoint = 0.1, na.value = "#EEEEEE", name = "|SMD|", limits = c(0, NA)) +
  facet_grid(var_group ~ country, scales = "free_y", space = "free_y") +
  labs(
    title    = "Absolute SMD: England vs Scotland \u2014 Matching Scale",
    subtitle = "Green < 0.10 (balanced) | Yellow = marginal | Red > 0.20 (imbalanced)\nLevels shown on log1p scale; trends on raw scale \u2014 both as matched",
    x = NULL, y = NULL,
    caption  = "Left = unmatched pool; Right = after matching"
  ) +
  theme_diag() +
  theme(axis.text.x   = element_text(angle = 20, hjust = 1, size = 9),
        axis.text.y   = element_text(size = 8),
        panel.spacing = unit(0.3, "lines"))

save_fig(p_heatmap, "fig03_smd_heatmap_by_country.png", width = 16, height = 20)
cat("  Saved: fig03_smd_heatmap_by_country.png\n\n")

# =============================================================================
# SECTION 11 — WEIGHT DISTRIBUTION PLOTS
# =============================================================================

cat("=== Section 11: Weight distribution plots ===\n")

make_weight_panels <- function(ctrl_df, ctry_label, ctry_col) {
  eff_n   <- round(sum(ctrl_df$weights)^2 / sum(ctrl_df$weights^2), 0)
  total_w <- sum(ctrl_df$weights)
  
  p1 <- ggplot(ctrl_df, aes(x = weights)) +
    geom_histogram(bins = 50, fill = ctry_col, alpha = 0.85) +
    scale_x_continuous(limits = c(0, 5.5)) +
    labs(title    = paste0(ctry_label, ": Weight Distribution"),
         subtitle = paste0("Nominal N = ", nrow(ctrl_df),
                           " | Effective N = ", eff_n,
                           " | Efficiency = ",
                           round(eff_n / nrow(ctrl_df), 3)),
         x = "Weight", y = "Count") +
    theme_diag()
  
  p2 <- ggplot(ctrl_df, aes(x = weights)) +
    stat_ecdf(linewidth = 1, colour = ctry_col) +
    geom_vline(xintercept = 5, linetype = "dashed",
               colour = "#E74C3C", linewidth = 0.6) +
    coord_cartesian(xlim = c(0, 6)) +
    labs(title    = paste0(ctry_label, ": ECDF of Weights"),
         subtitle = "Red dashed = cap at 5",
         x = "Weight", y = "Cumulative proportion") +
    theme_diag()
  
  list(p1 = p1, p2 = p2)
}

eng_wp <- make_weight_panels(ctrl_eng, "England", COL_ENGLAND)

if (has_scotland) {
  sco_wp <- make_weight_panels(ctrl_sco, "Scotland", COL_SCOTLAND)
  p_weights_country <- (eng_wp$p1 | sco_wp$p1) / (eng_wp$p2 | sco_wp$p2) +
    plot_annotation(
      title = "Weight Diagnostics by Country (after cap at 5)",
      theme = theme(plot.title = element_text(size = 13, face = "bold",
                                              colour = "#1A2E5A"))
    )
} else {
  p_weights_country <- eng_wp$p1 / eng_wp$p2 +
    plot_annotation(
      title = "Weight Diagnostics (after cap at 5)",
      theme = theme(plot.title = element_text(size = 13, face = "bold",
                                              colour = "#1A2E5A"))
    )
}
save_fig(p_weights_country, "fig06_weight_diagnostics_by_country.png",
         width = 16, height = 12)
cat("  Saved: fig06_weight_diagnostics_by_country.png\n\n")

# =============================================================================
# SECTION 12 — MAHALANOBIS DISTANCE PLOTS (by country)
#
# mdist reflects Stage 2 matching on 18 variables (9 injury trends + 9 injury
# levels). Diagnostics here focus on the TREND component of match quality —
# specifically whether the distance distribution is consistent with good
# pre-treatment trend balance, which is the key assumption for DiD estimation.
#
# Reference distribution: with 9 trend variables, distances from trend-only
# matching would follow chi-squared(9) (95th pct ≈ 17). The composite 18-var
# distance is larger, so d = 20 is used as a conservative flagging threshold.
#
# If mdist absent (older .rds files), section is skipped with a warning.
# =============================================================================

cat("=== Section 12: Mahalanobis distance plots (trend-focus) ===\n")

if ("mdist" %in% names(matched_eng) && (has_scotland && "mdist" %in% names(matched_sco) || !has_scotland)) {
  
  mdist_eng <- matched_eng %>%
    filter(treat_indicator == 1, !is.na(mdist)) %>%
    select(OA, mdist, country)
  if (has_scotland) {
    mdist_sco <- matched_sco %>%
      filter(treat_indicator == 1, !is.na(mdist)) %>%
      select(OA, mdist, country)
  } else {
    mdist_sco <- mdist_eng %>% filter(FALSE)
  }
  mdist_all <- bind_rows(mdist_eng, mdist_sco)
  
  # --- summary stats for caption annotation ---
  n_total   <- nrow(mdist_all)
  n_gt20    <- sum(mdist_all$mdist > 20)
  pct_gt20  <- round(100 * n_gt20 / n_total, 1)
  pct_lt5   <- round(100 * mean(mdist_all$mdist <= 5),  1)
  pct_lt10  <- round(100 * mean(mdist_all$mdist <= 10), 1)
  pct_lt20  <- round(100 * mean(mdist_all$mdist <= 20), 1)
  
  n_gt20_eng <- sum(mdist_eng$mdist > 20)
  n_gt20_sco <- if (nrow(mdist_sco) > 0) sum(mdist_sco$mdist > 20) else 0L
  pct_gt20_eng <- round(100 * n_gt20_eng / nrow(mdist_eng), 1)
  pct_gt20_sco <- if (nrow(mdist_sco) > 0) round(100 * n_gt20_sco / nrow(mdist_sco), 1) else 0
  
  cat(sprintf("  Distance summary (trend-focused interpretation):\n"))
  cat(sprintf("    Median: %.1f  |  Mean: %.1f\n",
              median(mdist_all$mdist), mean(mdist_all$mdist)))
  cat(sprintf("    %% <= 5:  %.1f%%  |  %% <= 10: %.1f%%  |  %% <= 20: %.1f%%\n",
              pct_lt5, pct_lt10, pct_lt20))
  cat(sprintf("    OAs with d > 20: %d (%.1f%%)  —  England: %d (%.1f%%)  Scotland: %d (%.1f%%)\n\n",
              n_gt20, pct_gt20, n_gt20_eng, pct_gt20_eng, n_gt20_sco, pct_gt20_sco))
  
  # Subtitle note explaining what the distance covers and why trends are the focus
  ecdf_subtitle <- paste0(
    "Distance reflects matching on 9 injury trends + 9 injury levels  |  ",
    "trend balance is the key DiD assumption\n",
    "Scotland's smaller control pool produces systematically larger distances"
  )
  
  p_mdist_ecdf <- ggplot(mdist_all, aes(x = mdist, colour = country)) +
    stat_ecdf(linewidth = 1.1) +
    # Reference lines: d=5 (good), d=10 (acceptable), d=20 (flagging threshold)
    geom_vline(xintercept = c(5, 10, 20),
               linetype   = c("dashed", "dotted", "dotdash"),
               colour     = c("#888888", "#CC3333", "#8B0000"),
               linewidth  = 0.5) +
    annotate("text", x = c(5.3, 10.3, 20.3), y = c(0.12, 0.12, 0.12),
             label  = c("d=5", "d=10", "d=20 flag"),
             colour = c("#888888", "#CC3333", "#8B0000"),
             size   = 3.2, hjust = 0) +
    scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    coord_cartesian(xlim = c(0, 35)) +
    labs(
      title    = "ECDF of Stage 2 Mahalanobis Distance \u2014 Treated OAs by Country",
      subtitle = ecdf_subtitle,
      x        = "Stage 2 Mahalanobis distance (9 trends + 9 levels)",
      y        = "Cumulative proportion of treated OAs",
      colour   = "Country",
      caption  = paste0(
        "England n=", nrow(mdist_eng), " | Scotland n=", nrow(mdist_sco),
        "  |  OAs with d > 20: England ", n_gt20_eng, " (", pct_gt20_eng, "%)",
        ", Scotland ", n_gt20_sco, " (", pct_gt20_sco, "%)"
      )
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  p_mdist_hist <- ggplot(mdist_all, aes(x = mdist, fill = country)) +
    geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
    geom_vline(xintercept = 20, linetype = "dotdash",
               colour = "#8B0000", linewidth = 0.5) +
    scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    coord_cartesian(xlim = c(0, 35)) +
    facet_wrap(~country, scales = "free_y") +
    labs(
      title   = "Stage 2 Distance Distribution by Country",
      subtitle = "England shows clean right-skew; Scotland has a secondary tail reflecting thin control pool",
      x       = "Mahalanobis distance (9 trends + 9 levels)",
      y       = "Count",
      fill    = NULL,
      caption = "Dark red line = d=20 flagging threshold"
    ) +
    theme_diag() +
    theme(legend.position = "none")
  
  p_mdist <- p_mdist_ecdf / p_mdist_hist +
    plot_layout(heights = c(1.3, 1)) +
    plot_annotation(
      caption = paste0(
        "Note: Mahalanobis distance computed on 18 Stage 2 variables ",
        "(9 pre-treatment injury trends + 9 mean injury levels). ",
        "Trend variables — covering car, cycling, pedestrian and other KSI and slight injuries ",
        "per km — are the primary focus as they underpin the parallel trends assumption for DiD estimation."
      ),
      theme = theme(
        plot.caption = element_text(size = 9, colour = "#666666",
                                    hjust = 0, margin = margin(t = 8))
      )
    )
  
  save_fig(p_mdist, "fig07_mahalanobis_distance_by_country.png",
           width = 12, height = 11)
  cat("  Saved: fig07_mahalanobis_distance_by_country.png\n\n")
  
} else {
  cat("  WARNING: mdist column not found in matched data — Section 12 skipped.\n")
  cat("  Re-run matching script to regenerate .rds files with mdist included.\n\n")
}
# =============================================================================
# SECTION 15 — ENGLAND vs SCOTLAND COMPARISON TABLE
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
    median_mdist         = if ("mdist" %in% names(matched_A))
      round(median(mdist, na.rm = TRUE), 2) else NA_real_,
    pct_isolated         = round(100 * mean(OA %in% isolated_ids), 1),
    .groups              = "drop"
  )
write_csv(country_summary, file.path(outdir, "08_country_comparison.csv"))
cat("  Saved: 08_country_comparison.csv\n\n")

# =============================================================================
# SECTION 17 — PARALLEL TRENDS (by country)
# =============================================================================

cat("=== Section 17: Parallel trends slope distributions ===\n")

make_trend_plot <- function(matched_df, ctry_label) {
  trend_long <- matched_df %>%
    filter(treat_indicator %in% c(0, 1)) %>%
    mutate(group = if_else(treat_indicator == 1, "Treated", "Matched control")) %>%
    select(OA, group, all_of(stage2_trends)) %>%
    pivot_longer(all_of(stage2_trends), names_to = "trend_var", values_to = "slope") %>%
    mutate(trend_label = coalesce(var_labels[trend_var], trend_var),
           trend_label = factor(trend_label)) %>%
    filter(!is.na(slope))
  
  ggplot(trend_long, aes(x = slope, colour = group, fill = group)) +
    geom_density(alpha = 0.20, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "#888888", linewidth = 0.4) +
    scale_colour_manual(
      values = c("Treated" = COL_TREATED, "Matched control" = COL_CONTROL)) +
    scale_fill_manual(
      values = c("Treated" = COL_TREATED, "Matched control" = COL_CONTROL)) +
    facet_wrap(~trend_label, scales = "free", ncol = 3) +
    labs(
      title    = paste0("Pre-Treatment Injury Trend Slopes \u2014 ", ctry_label),
      subtitle = "Overlap of distributions supports the parallel trends assumption",
      x        = "Pre-treatment slope (log-linear regression coefficient)",
      y        = "Density", colour = NULL, fill = NULL,
      caption  = paste0("n treated = ",   sum(matched_df$treat_indicator == 1),
                        " | n controls = ", sum(matched_df$treat_indicator == 0))
    ) +
    theme_diag() + theme(legend.position = "bottom")
}

p_trends_eng <- make_trend_plot(matched_eng, "England")
save_fig(p_trends_eng, "fig14_parallel_trends_england.png", width = 14, height = 12)
cat("  Saved: fig14_parallel_trends_england.png\n")

if (nrow(matched_sco) > 0) {
  p_trends_sco <- make_trend_plot(matched_sco, "Scotland")
  save_fig(p_trends_sco, "fig14_parallel_trends_scotland.png", width = 14, height = 12)
  cat("  Saved: fig14_parallel_trends_scotland.png\n")
} else {
  cat("  Skipped: Scotland (no data)\n")
}
cat("\n")

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("================================================================\n")
cat("FINAL DIAGNOSTIC SUMMARY\n")
cat("================================================================\n\n")

cat("SAMPLE SIZES:\n")
cat(sprintf("  Combined — Treated: %d | Controls: %d\n",
            sum(matched_A$treat_indicator   == 1),
            sum(matched_A$treat_indicator   == 0)))
cat(sprintf("  England  — Treated: %d | Controls: %d\n",
            sum(matched_eng$treat_indicator == 1),
            sum(matched_eng$treat_indicator == 0)))
if (has_scotland) {
  cat(sprintf("  Scotland — Treated: %d | Controls: %d\n",
              sum(matched_sco$treat_indicator == 1),
              sum(matched_sco$treat_indicator == 0)))
}
cat(sprintf("  Isolated OAs (flagged, included): %d\n", length(isolated_ids)))

cat("\nWEIGHT SUMMARY (after cap = 5):\n")
for (nm in c("England", if (has_scotland) "Scotland", "Overall")) {
  df  <- switch(nm, England = ctrl_eng, Scotland = ctrl_sco, Overall = ctrl_A)
  eff <- round(sum(df$weights)^2 / sum(df$weights^2), 0)
  cat(sprintf("  %-8s — N: %d | Eff N: %d | Efficiency: %.3f | Max weight: %.3f\n",
              nm, nrow(df), eff, eff / nrow(df), max(df$weights)))
}

cat("\nBALANCE SUMMARY (matching scale):\n")
cat(sprintf("  England  S1 — mean |SMD| after: %.3f\n",
            mean(abs(smd_s1_eng$smd_after_S1), na.rm = TRUE)))
if (has_scotland) {
  cat(sprintf("  Scotland S1 — mean |SMD| after: %.3f\n",
              mean(abs(smd_s1_sco$smd_after_S1), na.rm = TRUE)))
}
cat(sprintf("  England  S2 — max trend |SMD|: %.3f | max level (log) |SMD|: %.3f\n",
            max(abs(smd_s2_eng$smd_postS2[smd_s2_eng$var_type == "Trend"]),    na.rm = TRUE),
            max(abs(smd_s2_eng$smd_postS2[smd_s2_eng$var_type == "Level (log)"]), na.rm = TRUE)))
if (has_scotland) {
  cat(sprintf("  Scotland S2 — max trend |SMD|: %.3f | max level (log) |SMD|: %.3f\n",
              max(abs(smd_s2_sco$smd_postS2[smd_s2_sco$var_type == "Trend"]),    na.rm = TRUE),
              max(abs(smd_s2_sco$smd_postS2[smd_s2_sco$var_type == "Level (log)"]), na.rm = TRUE)))
}

cat("\nOUTPUTS SAVED TO:", outdir, "\n")
cat("  Descriptive tables  : 01_*_combined/england/scotland.csv  [raw scale]\n")
cat("  SMD Stage 1         : 02_smd_stage1_*.csv                 [log scale, as matched]\n")
cat("  SMD Stage 2         : 03_smd_stage2_*.csv                 [trends raw, levels log1p]\n")
cat("  Weight diagnostics  : 04_*.csv\n")
cat("  Stratum table       : 05_stratum_characteristics.csv\n")
cat("  Isolated OAs        : 06a/06b_*.csv\n")
cat("  Country comparison  : 08_country_comparison.csv\n")
cat("  Love plots S1       : fig01_* [log scale, as matched]\n")
cat("  Love plots S2       : fig02_* [trends raw, levels log1p — matching scale]\n")
cat("  SMD heatmap         : fig03_smd_heatmap_by_country.png    [matching scale]\n")
cat("  Weight diagnostics  : fig06_weight_diagnostics_by_country.png\n")
cat("  Mahalanobis dist    : fig07_mahalanobis_distance_by_country.png\n")
cat("  Parallel trends     : fig14_parallel_trends_england/scotland.png\n")
cat("================================================================\n")


# Verify all matching-scale vars exist in unmatched_eng_log before running
required_s1  <- stage1_vars          # already log-transformed names
required_s2  <- c(stage2_trends, stage2_levels_log)

missing_s1 <- setdiff(required_s1, names(unmatched_eng_log))
missing_s2 <- setdiff(required_s2, names(unmatched_eng_log))

cat("Missing S1 vars:", if (length(missing_s1) == 0) "none" else paste(missing_s1, collapse=", "), "\n")
cat("Missing S2 vars:", if (length(missing_s2) == 0) "none" else paste(missing_s2, collapse=", "), "\n")

# ╔═══════════════════════════════════════════════════════════════════════╗
# ║          PART B — PER-SCHEME DIAGNOSTICS (replaces 18b)             ║
# ╚═══════════════════════════════════════════════════════════════════════╝

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART B — PER-SCHEME DIAGNOSTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

schemes <- matched_A %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("Schemes:", paste(schemes, collapse = ", "), "\n\n")

# =============================================================================
# 18b-1. PER-SCHEME SMD (trends only — the key DiD assumption)
# =============================================================================

cat("--- Per-scheme SMD (trends) ---\n")

smd_for_scheme <- function(scheme_name, vars) {
  treated_oas <- matched_A %>%
    filter(scheme == scheme_name, treat_indicator == 1) %>%
    pull(OA)

  control_oas <- ctrl_scheme_lookup %>%
    filter(scheme == scheme_name) %>%
    pull(OA)

  # Before: this scheme's treated vs full unmatched control pool
  before_df <- unmatched_eng_log %>%
    filter(treated_OA == 1 & OA %in% treated_oas | control_group2_OA == 1) %>%
    mutate(treat_indicator = as.integer(OA %in% treated_oas))

  # After: this scheme's treated + their matched controls
  after_df <- matched_A %>%
    filter(OA %in% c(treated_oas, control_oas)) %>%
    mutate(treat_indicator = as.integer(OA %in% treated_oas))

  vars_avail <- intersect(vars, intersect(names(before_df), names(after_df)))

  map_df(vars_avail, function(v) {
    tibble(
      scheme     = scheme_name,
      variable   = v,
      label      = coalesce(var_labels[v], v),
      smd_before = round(compute_smd(before_df, v), 4),
      smd_after  = round(compute_smd(after_df,  v), 4)
    )
  })
}

smd_all_schemes <- map_df(schemes, function(s) {
  cat("  Computing SMD:", s, "\n")
  smd_for_scheme(s, c(stage1_vars, stage2_trends, stage2_levels_log))
})

write_csv(smd_all_schemes, file.path(outdir, "09_smd_per_scheme.csv"))

# =============================================================================
# 18b-2. PER-SCHEME TREND BALANCE HEATMAP (trends only — same as original 18b)
# =============================================================================

cat("\n--- Per-scheme trend balance heatmap ---\n")

s2_trend_heatmap <- smd_all_schemes %>%
  filter(variable %in% stage2_trends)

p_heatmap_scheme <- ggplot(s2_trend_heatmap,
                           aes(x = scheme, y = label, fill = abs(smd_after))) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(!is.na(smd_after),
                                sprintf("%.3f", smd_after), "")),
            size = 2.8, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#2ECC71", mid = "#F39C12", high = "#E74C3C",
    midpoint = 0.1, na.value = "#EEEEEE", name = "|SMD|", limits = c(0, NA)) +
  labs(
    title    = "Pre-Treatment Trend Balance by Scheme \u2014 After Matching",
    subtitle = "Green < 0.10 (balanced) | Orange = marginal | Red > 0.20 (imbalanced)\nTrends on raw scale (as matched)",
    x = NULL, y = NULL,
    caption  = "Per-scheme matching — each scheme matched independently"
  ) +
  theme_diag() +
  theme(axis.text.x   = element_text(angle = 30, hjust = 1, size = 11),
        axis.text.y   = element_text(size = 10),
        panel.grid    = element_blank(),
        panel.spacing = unit(0.3, "lines"))

save_fig(p_heatmap_scheme, "fig18b_trend_heatmap_per_scheme.png",
         width = 14, height = 10)

# =============================================================================
# 18b-3. PER-SCHEME LOVE PLOTS (faceted — trends + key S1 vars)
# =============================================================================

cat("--- Per-scheme love plots ---\n")

key_vars_scheme <- c(stage2_trends,
                     "log1p_road_density_m_km2", "log1p_pop_density", "IMD")

love_scheme_data <- smd_all_schemes %>%
  filter(variable %in% key_vars_scheme) %>%
  pivot_longer(c(smd_before, smd_after),
               names_to = "timing", values_to = "smd") %>%
  mutate(
    timing = if_else(timing == "smd_before", "Before", "After"),
    timing = factor(timing, levels = c("Before", "After"))
  )

p_love_scheme <- ggplot(
  love_scheme_data,
  aes(x = abs(smd), y = label, colour = timing, shape = timing)
) +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "#999999") +
  geom_vline(xintercept = 0, colour = "#DDDDDD") +
  geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.3) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Before" = COL_BEFORE, "After" = COL_AFTER)) +
  scale_shape_manual(values = c("Before" = 16, "After" = 17)) +
  facet_wrap(~scheme, ncol = 4) +
  labs(
    title    = "Per-Scheme Balance: Injury Trends + Key Covariates",
    subtitle = "Before vs after per-scheme matching | dashed = 0.10 threshold",
    x = "|SMD|", y = NULL, colour = NULL, shape = NULL,
    caption = "Each scheme matched independently against the shared other-city control pool"
  ) +
  theme_diag(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8))

save_fig(p_love_scheme, "fig18b_love_plots_per_scheme.png",
         width = 16, height = 12)

# =============================================================================
# 18b-4. PER-SCHEME BALANCE SUMMARY TABLE
# =============================================================================

cat("\n--- Per-scheme balance summary ---\n\n")

scheme_balance <- smd_all_schemes %>%
  filter(variable %in% stage2_trends) %>%
  group_by(scheme) %>%
  summarise(
    mean_trend_smd_before = round(mean(abs(smd_before), na.rm = TRUE), 4),
    mean_trend_smd_after  = round(mean(abs(smd_after),  na.rm = TRUE), 4),
    max_trend_smd_after   = round(max(abs(smd_after),   na.rm = TRUE), 4),
    n_trends_imbalanced   = sum(abs(smd_after) >= 0.10, na.rm = TRUE),
    all_trends_balanced   = all(abs(smd_after) < 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    matched_A %>%
      filter(treat_indicator == 1) %>%
      count(scheme, name = "n_treated"),
    by = "scheme"
  ) %>%
  left_join(
    ctrl_scheme_lookup %>% count(scheme, name = "n_controls"),
    by = "scheme"
  )

write_csv(scheme_balance, file.path(outdir, "18b_scheme_balance_summary.csv"))

cat("=== PER-SCHEME TREND BALANCE SUMMARY ===\n")
print(scheme_balance, n = Inf)

cat("\n================================================================\n")
cat("PART B COMPLETE — per-scheme diagnostics saved\n")
cat("================================================================\n")
cat("  09_smd_per_scheme.csv\n")
cat("  18b_scheme_balance_summary.csv\n")
cat("  fig18b_trend_heatmap_per_scheme.png\n")
cat("  fig18b_love_plots_per_scheme.png\n")

