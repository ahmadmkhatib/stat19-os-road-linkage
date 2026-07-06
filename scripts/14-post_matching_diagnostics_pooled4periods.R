# =============================================================================
# POST-MATCHING DIAGNOSTICS & DESCRIPTIVES â€” POOLED MATCHING
# =============================================================================
#
## matching pipeline (total injuries only, England only).
#
# PART A â€” Overall diagnostics (script 17 equivalent):
#   1. Descriptive summary tables (treated vs control)
#   2. SMD tables (Stage 1 + Stage 2)
#   3. Love plots (Stage 1 + Stage 2)
#   4. Weight distribution diagnostics + plots
#   5. Stratum characteristics
#   6. Common support / isolated OAs
#   7. Mahalanobis distance diagnostics
#   8. Parallel trends density plots
#
# PART B â€” Per-scheme diagnostics (script 18b equivalent):
#   9. Per-scheme SMD computation
#  10. Per-scheme love plots (faceted)
#  11. SMD heatmap across schemes
#  12. Per-scheme balance summary table
#
# INPUTS (from data/processed/):
#   OA_matched_full_pooled.rds
#   OA_matched_treated_pooled.rds
#   OA_matched_donors_pooled.rds
#   OA_common_support_flags_pooled.rds
#   OA_matching_pairs_pooled.rds
#   OA_matching_census.rds
#
# OUTPUTS:
#   output/diagnostics/pooled/  (CSV tables + PNG figures)
#
# =============================================================================

library(tidyverse)
library(here)
library(patchwork)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "pooled")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, filename, width = 14, height = 10, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# =============================================================================
# THEME
# =============================================================================

COL_TREATED <- "#D85A30"
COL_CONTROL <- "#2E6FAB"
COL_BEFORE  <- "#E74C3C"
COL_AFTER   <- "#2ECC71"

# Central text-size constants â€” change these to rescale everything uniformly
BASE_SIZE       <- 16   # base for theme_diag (axis text, strip text, etc.)
TITLE_SIZE      <- 19   # plot titles
SUBTITLE_SIZE   <- 15   # plot subtitles
CAPTION_SIZE    <- 13   # captions
AXIS_TITLE_SIZE <- 15   # axis titles
AXIS_TEXT_SIZE  <- 14   # axis tick labels
STRIP_TEXT_SIZE <- 15   # facet strip labels
LEGEND_TEXT_SIZE <- 14  # legend text
CELL_LABEL_SIZE <- 5    # geom_text inside heatmap tiles (in ggplot units)
POINT_SIZE      <- 4    # geom_point size for love plots / distance plots

theme_diag <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Titles
      plot.title       = element_text(size = TITLE_SIZE,    face = "bold",
                                      colour = "#1A2E5A",   margin = margin(b = 8)),
      plot.subtitle    = element_text(size = SUBTITLE_SIZE, colour = "#555555",
                                      margin = margin(b = 12)),
      plot.caption     = element_text(size = CAPTION_SIZE,  colour = "#888888",
                                      hjust = 0,            margin = margin(t = 8)),
      # Axes
      axis.title       = element_text(size = AXIS_TITLE_SIZE),
      axis.text        = element_text(size = AXIS_TEXT_SIZE),
      # Strips
      strip.text       = element_text(size = STRIP_TEXT_SIZE, face = "bold",
                                      colour = "#1A2E5A"),
      strip.background = element_rect(fill = "#EEF2F8", colour = NA),
      # Legend
      legend.text      = element_text(size = LEGEND_TEXT_SIZE),
      legend.title     = element_text(size = LEGEND_TEXT_SIZE, face = "bold"),
      legend.key.size  = unit(1.1, "lines"),
      # Grid
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

# Stage 1 (same as script 16)
stage1_vars_raw <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_BUA_centroid", "pop_density", "area_km2",
  "business_retail_per_km2", "IMD",
  "ethnic_minority_pct",
  "no_car_households_pct", "two_plus_car_households_pct",
  "car_commute_pct", "public_transport_to_work_pct",
  "active_travel_to_work_pct", "work_from_home_pct",
  "children_0_19_pct", "young_people_5_19_pct",
  "young_adults_20_34_pct", "working_age_20_64_pct",
  "older_65plus_pct"
)

log_transform_s1 <- c("road_length_km", "pop_density", "dist_BUA_centroid",
                      "road_density_m_km2", "business_retail_per_km2")
log_nozero_s1    <- c("area_km2")

stage1_vars_log <- c(
  paste0("log1p_", log_transform_s1),
  paste0("log_", log_nozero_s1),
  setdiff(stage1_vars_raw, c(log_transform_s1, log_nozero_s1))
)

# Stage 2 - pooled outcome-history variables
stage2_trends <- c(
  "trend_total_pkm",
  "recent_minus_mid_total_pkm",
  "mid_minus_early_total_pkm",
  "covid_minus_precovid_total_pkm",
  "recovery_minus_covid_total_pkm"
)
stage2_levels <- c(
  "mean_total_pkm",
  "mean_precovid_total_pkm",
  "mean_lockdown_total_pkm",
  "mean_recovery_total_pkm",
  "zero_quarter_share_pre"
)
stage2_vars_log <- c(stage2_trends, paste0("log1p_", stage2_levels))

all_match_vars <- c(stage1_vars_log, stage2_vars_log)

# Variable labels for display
var_labels <- c(
  "log1p_road_length_km"       = "Road length (log)",
  "log1p_road_density_m_km2"   = "Road density (log)",
  "log_area_km2"               = "Area (log)",
  "pct_A_road"                 = "% A road",
  "pct_B_road"                 = "% B road",
  "pct_minor_road"             = "% Minor road",
  "log1p_dist_BUA_centroid"      = "Dist. to BUA centroid (log)",
  "log1p_pop_density"          = "Pop. density (log)",
  "log1p_business_retail_per_km2" = "Retail density (log)",
  "IMD"                        = "IMD",
  "ethnic_minority_pct"         = "% Ethnic minority",
  "no_car_households_pct"       = "% No car household",
  "two_plus_car_households_pct" = "% 2+ car household",
  "car_commute_pct"             = "% Car commute",
  "public_transport_to_work_pct" = "% Public/paid transport to work",
  "active_travel_to_work_pct"   = "% Active travel to work",
  "work_from_home_pct"          = "% Work from home",
  "children_0_19_pct"           = "% Aged 0-19",
  "young_people_5_19_pct"       = "% Aged 5-19",
  "young_adults_20_34_pct"      = "% Aged 20-34",
  "working_age_20_64_pct"       = "% Aged 20-64",
  "older_65plus_pct"            = "% Aged 65+",
  "cars_none_pct"              = "% No car",
  "no_cars_pct"                = "% No car",
  "cars_zero_pct"              = "% No car",
  "cars_one_pct"               = "% 1 car",
  "cars_twoPlus_pct"           = "% 2+ cars",
  "Drive_Car_pct"              = "% Drive to work",
  "Passenger_Car_pct"          = "% Passenger car",
  "Walk_pct"                   = "% Walk to work",
  "Bicycle_pct"                = "% Cycle to work",
  "bus_Coach_pct"              = "% Bus/coach",
  "Train_pct"                  = "% Train",
  "Underground_train_tram_pct" = "% Underground/tram",
  "Taxi_pct"                   = "% Taxi",
  "Motorcycle_pct"             = "% Motorcycle",
  "workAthome_pct"             = "% Work from home",
  "Other_pct"                  = "% Other commute",
  "White_pct"                  = "% White",
  "Mixed_pct"                  = "% Mixed ethnicity",
  "Asian_pct"                  = "% Asian",
  "Black_pct"                  = "% Black",
  "Other_ethnic_pct"           = "% Other ethnicity",
  "Other_ethnicity_pct"        = "% Other ethnicity",
  "Other_ethnic_group_pct"     = "% Other ethnicity",
  "age_under15_pct"            = "% Under 15",
  "age_15to24_pct"             = "% 15-24",
  "age_25to44_pct"             = "% 25-44",
  "age_45to64_pct"             = "% 45-64",
  "age_65to84_pct"             = "% 65-84",
  "X85plus_pct" = "% 85+",
  "trend_total_pkm"            = "Trend: total injuries/km",
  "recent_minus_mid_total_pkm" = "Recent - mid injuries/km",
  "mid_minus_early_total_pkm"  = "Mid - early injuries/km",
  "covid_minus_precovid_total_pkm" = "COVID - pre-COVID injuries/km",
  "recovery_minus_covid_total_pkm" = "Recovery - COVID injuries/km",
  "log1p_mean_total_pkm"       = "Level: mean total injuries/km (log)",
  "log1p_mean_precovid_total_pkm" = "Level: pre-COVID injuries/km (log)",
  "log1p_mean_lockdown_total_pkm" = "Level: lockdown injuries/km (log)",
  "log1p_mean_recovery_total_pkm" = "Level: recovery injuries/km (log)",
  "log1p_zero_quarter_share_pre"  = "Zero-injury quarter share (log)"
)

# I report balance in compact groups for the main text, while keeping the full
# variable-by-variable SMDs for the appendix.
balance_groups <- list(
  "Road/network covariates" = c(
    "log1p_road_length_km", "log1p_road_density_m_km2", "log_area_km2",
    "pct_A_road", "pct_B_road", "pct_minor_road",
    "log1p_dist_BUA_centroid", "log1p_pop_density",
    "log1p_business_retail_per_km2"
  ),
  "Socio-demographic covariates" = c(
    "IMD", "ethnic_minority_pct", "no_car_households_pct",
    "two_plus_car_households_pct", "car_commute_pct",
    "public_transport_to_work_pct", "active_travel_to_work_pct",
    "work_from_home_pct", "children_0_19_pct", "young_people_5_19_pct",
    "young_adults_20_34_pct", "working_age_20_64_pct", "older_65plus_pct"
  ),
  "Overall pre-treatment trend" = "trend_total_pkm",
  "Within-pre-period trend shape" = c(
    "recent_minus_mid_total_pkm", "mid_minus_early_total_pkm"
  ),
  "Pre-COVID injury level" = "log1p_mean_precovid_total_pkm",
  "Lockdown injury level" = "log1p_mean_lockdown_total_pkm",
  "Recovery injury level" = "log1p_mean_recovery_total_pkm",
  "COVID response change" = "covid_minus_precovid_total_pkm",
  "Recovery response change" = "recovery_minus_covid_total_pkm",
  "Overall pre-treatment injury level" = "log1p_mean_total_pkm",
  "Zero-quarter share" = "log1p_zero_quarter_share_pre"
)

balance_group_order <- names(balance_groups)

# =============================================================================
# LOAD DATA
# =============================================================================

matched_full <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))
matched_treated  <- readRDS(here("data", "processed", "OA_matched_treated_pooled.rds"))
matched_controls <- readRDS(here("data", "processed", "OA_matched_donors_pooled.rds"))
isolated_flags   <- readRDS(here("data", "processed", "OA_common_support_flags_pooled.rds"))
matching_pairs   <- readRDS(here("data", "processed", "OA_matching_pairs_pooled.rds"))
full_data        <- readRDS(here("data", "processed", "OA_matching_census.rds"))

cat("=== MATCHED DATASET (POOLED) ===\n")
cat("Treated:", nrow(matched_treated),
    "| Controls:", sum(matched_full$treat_indicator == 0), "\n")
cat("Unique control OAs:", n_distinct(matched_controls$OA), "\n\n")

# Prepare log-transformed columns
add_grouped_stage1_vars <- function(data) {
  if (!"ethnic_minority_pct" %in% names(data)) {
    other_eth_col <- intersect(
      c("Other_ethnicity_pct", "Other_ethnic_pct", "Other_ethnic_group_pct"),
      names(data)
    )
    eth_source <- c("Mixed_pct", "Asian_pct", "Black_pct", other_eth_col[1])
    if (all(!is.na(eth_source)) && all(eth_source %in% names(data))) {
      data[["ethnic_minority_pct"]] <- rowSums(data[eth_source], na.rm = FALSE)
    }
  }
  
  no_car_col <- intersect(
    c("cars_none_pct", "no_cars_pct", "cars_zero_pct"),
    names(data)
  )
  if (!"no_car_households_pct" %in% names(data) && length(no_car_col) > 0) {
    data[["no_car_households_pct"]] <- data[[no_car_col[1]]]
  }
  
  if (!"two_plus_car_households_pct" %in% names(data) &&
      "cars_twoPlus_pct" %in% names(data)) {
    data[["two_plus_car_households_pct"]] <- data[["cars_twoPlus_pct"]]
  }
  
  if (!"car_commute_pct" %in% names(data) &&
      all(c("Drive_Car_pct", "Passenger_Car_pct") %in% names(data))) {
    data[["car_commute_pct"]] <-
      data[["Drive_Car_pct"]] + data[["Passenger_Car_pct"]]
  }
  
  public_transport_cols <- c(
    "Underground_train_tram_pct", "Train_pct", "bus_Coach_pct", "Taxi_pct"
  )
  if (!"public_transport_to_work_pct" %in% names(data) &&
      all(public_transport_cols %in% names(data))) {
    data[["public_transport_to_work_pct"]] <-
      rowSums(data[public_transport_cols], na.rm = FALSE)
  }
  
  if (!"active_travel_to_work_pct" %in% names(data) &&
      all(c("Walk_pct", "Bicycle_pct") %in% names(data))) {
    data[["active_travel_to_work_pct"]] <-
      data[["Walk_pct"]] + data[["Bicycle_pct"]]
  }
  
  if (!"work_from_home_pct" %in% names(data) &&
      "workAthome_pct" %in% names(data)) {
    data[["work_from_home_pct"]] <- data[["workAthome_pct"]]
  }
  
  if (!"children_0_19_pct" %in% names(data) &&
      all(c("X4under_pct", "X5to19_pct") %in% names(data))) {
    data[["children_0_19_pct"]] <- data[["X4under_pct"]] + data[["X5to19_pct"]]
  }
  
  if (!"young_people_5_19_pct" %in% names(data) &&
      "X5to19_pct" %in% names(data)) {
    data[["young_people_5_19_pct"]] <- data[["X5to19_pct"]]
  }
  
  if (!"young_adults_20_34_pct" %in% names(data) &&
      all(c("X20to24_pct", "X25to29_pct", "X30to34_pct") %in% names(data))) {
    data[["young_adults_20_34_pct"]] <-
      data[["X20to24_pct"]] + data[["X25to29_pct"]] + data[["X30to34_pct"]]
  }
  
  if (!"working_age_20_64_pct" %in% names(data) &&
      "X20to64_pct" %in% names(data)) {
    data[["working_age_20_64_pct"]] <- data[["X20to64_pct"]]
  }
  
  if (!"older_65plus_pct" %in% names(data) &&
      "X65plus_pct" %in% names(data)) {
    data[["older_65plus_pct"]] <- data[["X65plus_pct"]]
  }
  
  data
}

add_log_vars <- function(data) {
  for (v in log_transform_s1) {
    if (v %in% names(data))
      data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  for (v in log_nozero_s1) {
    if (v %in% names(data))
      data[[paste0("log_", v)]] <- log(pmax(data[[v]], 1e-10))
  }
  for (v in stage2_levels) {
    if (v %in% names(data))
      data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
  }
  data
}

matched_full <- matched_full %>%
  add_grouped_stage1_vars() %>%
  add_log_vars()
matched_treated <- add_grouped_stage1_vars(matched_treated)
matched_controls <- add_grouped_stage1_vars(matched_controls)

required_stage2_diag_vars <- c(stage2_trends, paste0("log1p_", stage2_levels))
missing_stage2_diag_vars <- setdiff(required_stage2_diag_vars, names(matched_full))

if (length(missing_stage2_diag_vars) > 0) {
  stop(
    "Missing Stage 2 diagnostic variables in matched_full: ",
    paste(missing_stage2_diag_vars, collapse = ", "),
    "\nRegenerate the matching-variable dataset, merge census data, and rerun matching first."
  )
}

# Unmatched pool (for before-matching comparisons)
unmatched_pool <- full_data %>%
  filter(
    substr(LAD24CD, 1, 1) == "E",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA == 0,
    n_roads   > 0,
    !(treated_OA == 1 & zero_injury_OA == 1),
    !(control_group2_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1)) %>%
  add_grouped_stage1_vars() %>%
  add_log_vars()

missing_stage2_unmatched_vars <- setdiff(required_stage2_diag_vars, names(unmatched_pool))

if (length(missing_stage2_unmatched_vars) > 0) {
  stop(
    "Missing Stage 2 diagnostic variables in unmatched_pool: ",
    paste(missing_stage2_unmatched_vars, collapse = ", "),
    "\nRegenerate OA_matching_census.rds before running diagnostics."
  )
}

# =============================================================================
# HELPER FUNCTIONS
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

desc_stats <- function(data, var) {
  t <- data[[var]][data$treat_indicator == 1]
  c <- data[[var]][data$treat_indicator == 0]
  tibble(
    variable     = var,
    label        = coalesce(var_labels[var], var),
    treated_mean = round(mean(t, na.rm = TRUE), 4),
    treated_sd   = round(sd(t,   na.rm = TRUE), 4),
    control_mean = round(mean(c, na.rm = TRUE), 4),
    control_sd   = round(sd(c,   na.rm = TRUE), 4)
  )
}

balance_group_for_variable <- function(variable) {
  matched_group <- names(balance_groups)[
    map_lgl(balance_groups, ~ variable %in% .x)
  ]
  if (length(matched_group) == 0) return("Other")
  matched_group[1]
}

add_balance_group <- function(smd_df) {
  smd_df %>%
    mutate(
      balance_group = map_chr(variable, balance_group_for_variable),
      balance_group = factor(balance_group,
                             levels = c(balance_group_order, "Other"))
    )
}

max_or_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

worst_label_or_na <- function(labels, values) {
  if (all(is.na(values))) return(NA_character_)
  labels[which.max(replace_na(values, -Inf))][1]
}

summarise_balance_groups <- function(smd_df, group_vars = character()) {
  smd_df %>%
    add_balance_group() %>%
    mutate(
      abs_smd_before = abs(smd_before),
      abs_smd_after  = abs(smd_after)
    ) %>%
    group_by(across(all_of(group_vars)), balance_group) %>%
    summarise(
      n_variables = n(),
      mean_abs_smd_before = round(mean(abs_smd_before, na.rm = TRUE), 4),
      mean_abs_smd_after  = round(mean(abs_smd_after,  na.rm = TRUE), 4),
      max_abs_smd_after   = round(max_or_na(abs_smd_after), 4),
      n_imbalanced_0_10   = sum(abs_smd_after >= 0.10, na.rm = TRUE),
      n_imbalanced_0_05   = sum(abs_smd_after >= 0.05, na.rm = TRUE),
      worst_variable      = worst_label_or_na(label, abs_smd_after),
      all_balanced_0_10   = all(abs_smd_after < 0.10, na.rm = TRUE),
      all_balanced_0_05   = all(abs_smd_after < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(across(all_of(group_vars)), balance_group)
}

# â•”â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•—
# â•‘                     PART A â€” OVERALL DIAGNOSTICS                    â•‘
# â•šâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART A â€” OVERALL POST-MATCHING DIAGNOSTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
# =============================================================================
# SCHEME-LEVEL MATCHING SUMMARY
# =============================================================================

cat("--- 0. Scheme-level matching summary ---\n")

# Check matching_pairs column names
print(names(matching_pairs))

# Detect treated/control OA column names
pair_treated_col <- case_when(
  "treated_OA"  %in% names(matching_pairs) ~ "treated_OA",
  "treated_oa"  %in% names(matching_pairs) ~ "treated_oa",
  "treated"     %in% names(matching_pairs) ~ "treated",
  "OA_treated"  %in% names(matching_pairs) ~ "OA_treated",
  TRUE ~ NA_character_
)

pair_control_col <- case_when(
  "control_OA"  %in% names(matching_pairs) ~ "control_OA",
  "control_oa"  %in% names(matching_pairs) ~ "control_oa",
  "control"     %in% names(matching_pairs) ~ "control",
  "OA_control"  %in% names(matching_pairs) ~ "OA_control",
  TRUE ~ NA_character_
)

if (is.na(pair_treated_col) | is.na(pair_control_col)) {
  stop(
    "Could not find treated/control OA columns in matching_pairs. ",
    "Run names(matching_pairs) and update pair_treated_col / pair_control_col manually."
  )
}

scheme_summary_table <- matching_pairs %>%
  group_by(scheme) %>%
  summarise(
    treated         = n_distinct(.data[[pair_treated_col]]),
    unique_controls = n_distinct(.data[[pair_control_col]]),
    total_pairs     = n(),
    ratio           = paste0("1:", round(total_pairs / treated, 1)),
    .groups = "drop"
  )

scheme_summary_table <- bind_rows(
  scheme_summary_table,
  tibble(
    scheme = "Total",
    treated = sum(scheme_summary_table$treated, na.rm = TRUE),
    unique_controls = sum(scheme_summary_table$unique_controls, na.rm = TRUE),
    total_pairs = sum(scheme_summary_table$total_pairs, na.rm = TRUE),
    ratio = paste0("1:", round(
      sum(scheme_summary_table$total_pairs, na.rm = TRUE) /
        sum(scheme_summary_table$treated, na.rm = TRUE), 1
    ))
  )
)

print(scheme_summary_table)

write_csv(
  scheme_summary_table,
  file.path(outdir, "00_scheme_matching_summary.csv")
)

cat("  Saved: 00_scheme_matching_summary.csv\n\n")# 1. DESCRIPTIVE SUMMARY TABLE
# =============================================================================

cat("--- 1. Descriptive summary ---\n")

desc_raw_vars <- c(stage1_vars_raw, stage2_trends, stage2_levels)
desc_raw_vars <- intersect(desc_raw_vars, names(matched_full))

desc_table <- map_df(desc_raw_vars, ~ desc_stats(matched_full, .x))
write_csv(desc_table, file.path(outdir, "01_descriptive_table.csv"))
cat("  Saved: 01_descriptive_table.csv\n")



desc_table |> print(n=50)
# =============================================================================
# SMD TABLES
# =============================================================================

cat("--- 2. SMD tables ---\n")

smd_table <- function(vars, before_data, after_data, label) {
  vars_avail <- intersect(vars, intersect(names(before_data), names(after_data)))
  map_df(vars_avail, function(v) {
    tibble(
      stage         = label,
      variable      = v,
      label         = coalesce(var_labels[v], v),
      smd_before    = round(compute_smd(before_data, v), 4),
      smd_after     = round(compute_smd(after_data,  v), 4),
      balanced      = abs(compute_smd(after_data, v)) < 0.10
    )
  })
}

smd_s1 <- smd_table(stage1_vars_log, unmatched_pool, matched_full, "Stage 1")
smd_s2 <- smd_table(stage2_vars_log, unmatched_pool, matched_full, "Stage 2")

write_csv(smd_s1, file.path(outdir, "02_smd_stage1.csv"))
write_csv(smd_s2, file.path(outdir, "03_smd_stage2.csv"))

# I keep this full table for the appendix, because it shows every matching
# variable without forcing the main paper table to become unreadable.
smd_appendix <- bind_rows(smd_s1, smd_s2) %>%
  add_balance_group() %>%
  arrange(stage, balance_group, variable)

write_csv(smd_appendix, file.path(outdir, "03b_smd_full_appendix.csv"))

# I use this grouped table as the main balance summary to report.
balance_group_summary <- summarise_balance_groups(smd_appendix, group_vars = "stage")
write_csv(balance_group_summary,
          file.path(outdir, "03c_smd_compact_balance_by_group.csv"))

cat("  Stage 1: mean |SMD| before =", round(mean(abs(smd_s1$smd_before), na.rm = TRUE), 4),
    "| after =", round(mean(abs(smd_s1$smd_after), na.rm = TRUE), 4), "\n")
cat("  Stage 2: mean |SMD| before =", round(mean(abs(smd_s2$smd_before), na.rm = TRUE), 4),
    "| after =", round(mean(abs(smd_s2$smd_after), na.rm = TRUE), 4), "\n")
cat("  Stage 2 all balanced (<0.10):", all(smd_s2$balanced, na.rm = TRUE), "\n\n")
cat("  Saved: 03b_smd_full_appendix.csv\n")
cat("  Saved: 03c_smd_compact_balance_by_group.csv\n\n")

print(balance_group_summary, n = Inf)

# =============================================================================
#  LOVE PLOTS
# =============================================================================

cat("--- 3. Love plots ---\n")

make_love_plot <- function(smd_df, title, subtitle = "") {
  plot_data <- smd_df %>%
    filter(!is.na(smd_before), !is.na(smd_after)) %>%
    mutate(label = factor(label, levels = label[order(abs(smd_before))])) %>%
    pivot_longer(c(smd_before, smd_after),
                 names_to = "timing", values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_before", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  ggplot(plot_data,
         aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = 0.10, linetype = "dashed", colour = "#999999") +
    geom_vline(xintercept = 0, colour = "#DDDDDD") +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.3) +
    geom_point(size = POINT_SIZE) +
    scale_colour_manual(values = c("Before matching" = COL_BEFORE,
                                   "After matching"  = COL_AFTER)) +
    scale_shape_manual(values = c("Before matching" = 16,
                                  "After matching"  = 17)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, subtitle = subtitle,
         x = "|SMD|", y = NULL, colour = NULL, shape = NULL,
         caption = "Dashed = 0.10 threshold") +
    theme_diag() +
    theme(legend.position = "bottom")
}

p_love_s1 <- make_love_plot(
  smd_s1,
  "Stage 1 balance: structural & sociodemographic covariates",
  "Pooled"
)
save_fig(p_love_s1, "fig01_love_plot_stage1.png", width = 14, height = 16)

p_love_s2 <- make_love_plot(
  smd_s2,
  "Stage 2 balance - total injury trajectory + level",
  "Pooled matching - pre-COVID trend, trajectory-shape, and level variables"
)
save_fig(p_love_s2, "fig02_love_plot_stage2.png", width = 12, height = 6)

# =============================================================================
# 4. WEIGHT DISTRIBUTION
# =============================================================================

cat("--- 4. Weight distribution ---\n")

ctrl <- matched_full %>% filter(treat_indicator == 0)
w <- ctrl$weights

eff_n <- sum(w)^2 / sum(w^2)
efficiency <- eff_n / length(w)

wt_summary <- tibble(
  n_controls   = length(w),
  mean_wt      = round(mean(w), 3),
  median_wt    = round(median(w), 3),
  sd_wt        = round(sd(w), 3),
  p90_wt       = round(quantile(w, 0.90), 3),
  p95_wt       = round(quantile(w, 0.95), 3),
  max_wt       = round(max(w), 3),
  effective_n  = round(eff_n, 1),
  efficiency   = round(efficiency, 3),
  n_at_cap5    = sum(w >= 5),
  pct_at_cap5  = round(100 * mean(w >= 5), 2)
)

write_csv(wt_summary, file.path(outdir, "04_weight_distribution.csv"))
print(wt_summary)

# Control reuse frequency
reuse <- matching_pairs %>%
  count(control_OA, name = "times_used") %>%
  count(times_used, name = "n_controls")
write_csv(reuse, file.path(outdir, "04b_control_reuse.csv"))

# Top 20 controls by weight
top_ctrl <- ctrl %>%
  distinct(OA, weights) %>%
  arrange(desc(weights)) %>%
  head(20) %>%
  mutate(cum_weight = cumsum(weights),
         cum_pct    = round(100 * cum_weight / sum(ctrl$weights), 2))
write_csv(top_ctrl, file.path(outdir, "04c_top_weight_controls.csv"))

# Weight plots
p_wt_hist <- ggplot(ctrl, aes(x = weights)) +
  geom_histogram(bins = 50, fill = COL_CONTROL, alpha = 0.7) +
  labs(
    title = "Control weight distribution",
    subtitle = sprintf("N = %d | Effective N = %.0f | Efficiency = %.1f%%",
                       length(w), eff_n, 100 * efficiency),
    x = "Weight", y = "Count"
  ) +
  theme_diag()

p_wt_ecdf <- ggplot(ctrl, aes(x = weights)) +
  stat_ecdf(linewidth = 0.8, colour = COL_CONTROL) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = COL_BEFORE) +
  labs(title = "Weight ECDF", x = "Weight", y = "Cumulative proportion") +
  theme_diag()

save_fig(p_wt_hist + p_wt_ecdf, "fig06_weight_diagnostics.png",
         width = 14, height = 6)
cat("\n")

# =============================================================================
# STRATUM CHARACTERISTICS
# =============================================================================

cat("--- 5. Stratum characteristics ---\n")

stratum_vars <- intersect(
  c("mean_total_pkm", "road_length_km", "road_density_m_km2",
    "pop_density", "dist_BUA_centroid", "IMD",
    "pct_A_road", "pct_minor_road", "car_commute_pct",
    "public_transport_to_work_pct", "active_travel_to_work_pct"),
  names(matched_full)
)

stratum_stats <- matched_full %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  group_by(baseline_injury_stratum) %>%
  summarise(
    n = n(),
    across(all_of(stratum_vars),
           list(mean = ~ round(mean(., na.rm = TRUE), 3),
                sd   = ~ round(sd(.,   na.rm = TRUE), 3)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write_csv(stratum_stats, file.path(outdir, "05_stratum_characteristics.csv"))
cat("  Saved: 05_stratum_characteristics.csv\n\n")

# =============================================================================
# 6. COMMON SUPPORT / ISOLATED OAs
# =============================================================================

cat("--- 6. Common support ---\n")

n_isolated <- sum(isolated_flags$structurally_isolated)
cat("  Isolated OAs:", n_isolated, "/", nrow(isolated_flags), "\n")

if (n_isolated > 0) {
  isolated_oas <- isolated_flags %>%
    filter(structurally_isolated) %>%
    pull(treated_OA)
  
  iso_compare <- matched_full %>%
    filter(treat_indicator == 1) %>%
    mutate(isolated = OA %in% isolated_oas) %>%
    group_by(isolated) %>%
    summarise(
      n = n(),
      across(all_of(intersect(c("mean_total_pkm", "road_length_km",
                                "pop_density", "IMD"), names(.))),
             ~ round(mean(., na.rm = TRUE), 3)),
      .groups = "drop"
    )
  
  write_csv(iso_compare, file.path(outdir, "06a_isolated_OA_characteristics.csv"))
  print(iso_compare)
}
cat("\n")

# =============================================================================
# MAHALANOBIS DISTANCE
# =============================================================================


if ("mdist" %in% names(matched_full)) {
  dist_summary <- matched_full %>%
    filter(treat_indicator == 1, !is.na(mdist)) %>%
    summarise(
      n          = n(),
      median_d   = round(median(mdist), 2),
      mean_d     = round(mean(mdist), 2),
      pct_le5    = round(100 * mean(mdist <= 5), 1),
      pct_le10   = round(100 * mean(mdist <= 10), 1),
      pct_gt20   = round(100 * mean(mdist > 20), 1)
    )
  
  cat("  Distance summary (treated OAs):\n")
  print(dist_summary)
  
  p_mdist <- ggplot(matched_full %>% filter(!is.na(mdist)),
                    aes(x = mdist, colour = factor(treat_indicator))) +
    stat_ecdf(linewidth = 0.8) +
    geom_vline(xintercept = c(5, 10, 20),
               linetype = c("dashed", "dotted", "longdash"),
               colour = "#888888") +
    scale_colour_manual(values = c("0" = COL_CONTROL, "1" = COL_TREATED),
                        labels = c("Control", "Treated")) +
    labs(
      title = "Mahalanobis distance ECDF (Stage 2)",
      subtitle = "Stage 2 variables: trend, trajectory-shape, and level",
      x = "Mahalanobis distance", y = "Cumulative proportion",
      colour = NULL
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  save_fig(p_mdist, "fig07_mahalanobis_distance.png", width = 10, height = 7)
} else {
  cat("  No mdist column in matched data â€” skipping.\n")
}
cat("\n")

# =============================================================================
# 8. PARALLEL TRENDS DENSITY PLOTS
# =============================================================================

cat("--- 8. Parallel trends distributions ---\n")

trend_vars <- intersect(
  c(
    "trend_total_pkm"
    #,
    # "recent_minus_mid_total_pkm",
    #  "mid_minus_early_total_pkm"
  ),
  names(matched_full)
)

if (length(trend_vars) > 0) {
  trend_data <- matched_full %>%
    select(OA, treat_indicator, all_of(trend_vars)) %>%
    pivot_longer(all_of(trend_vars), names_to = "variable", values_to = "value") %>%
    mutate(
      group = if_else(treat_indicator == 1, "Treated", "Control"),
      label = coalesce(var_labels[variable], variable)
    )
  
  p_trends <- ggplot(trend_data, aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
    facet_wrap(~label, scales = "free") +
    scale_fill_manual(values = c("Treated" = COL_TREATED, "Control" = COL_CONTROL)) +
    scale_colour_manual(values = c("Treated" = COL_TREATED, "Control" = COL_CONTROL)) +
    labs(
      title    = "Pre-treatment trend distributions: treated vs matched controls",
      x = "Pre-treatment trend value", y = "Density",
      fill = NULL, colour = NULL
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  save_fig(p_trends, "fig14_parallel_trends.png", width = 10, height = 7)
  
  # Per-scheme version
  trend_scheme_data <- matched_full %>%
    select(OA, treat_indicator, scheme, all_of(trend_vars)) %>%
    pivot_longer(all_of(trend_vars), names_to = "variable", values_to = "value") %>%
    mutate(
      group = if_else(treat_indicator == 1, "Treated", "Control"),
      label = coalesce(var_labels[variable], variable)
    )
  
  p_trends_scheme <- ggplot(trend_scheme_data,
                            aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
    facet_wrap(~scheme, scales = "free_y") +
    scale_fill_manual(values = c("Treated" = COL_TREATED, "Control" = COL_CONTROL)) +
    scale_colour_manual(values = c("Treated" = COL_TREATED, "Control" = COL_CONTROL)) +
    labs(
      title    = "Pre-treatment trend distributions by scheme: treated vs matched controls",
      subtitle = "Pre-COVID slope and trajectory-shape variables",
      x = "Pre-treatment trend/shape value", y = "Density",
      fill = NULL, colour = NULL
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  save_fig(p_trends_scheme, "fig14b_parallel_trends_per_scheme.png",
           width = 16, height = 12)
  
  # Per-scheme pre-treatment time series
  library(zoo)
  
  oa_q <- readRDS(here("data", "processed", "OA_injuries_quarterly.rds")) %>%
    mutate(quarter_year = as.yearqtr(quarter_year))
  
  scheme_starts <- readRDS(here("data", "processed", "roads_caz_props.rds")) %>%
    distinct(scheme, caz_start_q) %>%
    filter(!is.na(scheme))
  
  oa_panel <- oa_q %>%
    inner_join(
      matched_full %>%
        select(OA, treat_indicator, weights, scheme, road_length_km),
      by = "OA"
    ) %>%
    left_join(scheme_starts, by = "scheme") %>%
    mutate(injuries_per_km = total_injuries / pmax(road_length_km, 0.001))
  
  oa_pre_quarterly <- oa_panel %>%
    filter(quarter_year < caz_start_q) %>%
    group_by(scheme, quarter_year, treat_indicator) %>%
    summarise(
      wtd_mean_inj = weighted.mean(injuries_per_km, w = weights, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      group  = if_else(treat_indicator == 1, "Treated", "Matched control"),
      q_date = as.Date(quarter_year)
    )
  
  vline_data <- scheme_starts %>%
    filter(scheme %in% unique(oa_pre_quarterly$scheme)) %>%
    mutate(start_date = as.Date(caz_start_q))
  
  p_ts_scheme <- ggplot(
    oa_pre_quarterly %>%
      mutate(group = factor(group, levels = c("Treated", "Matched control"))),
    aes(x = q_date, y = wtd_mean_inj, colour = group, linetype = group)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = POINT_SIZE - 1, alpha = 0.7) +
    geom_vline(
      data = vline_data,
      aes(xintercept = start_date),
      linetype = "dotted", colour = "#888888", linewidth = 0.5
    ) +
    scale_colour_manual(
      values = c("Treated" = COL_TREATED, "Matched control" = COL_CONTROL)
    ) +
    scale_linetype_manual(
      values = c("Treated" = "solid", "Matched control" = "longdash")
    ) +
    facet_wrap(~scheme, ncol = 3, scales = "free_y") +
    labs(
      title    = "Pre-treatment injury trends by scheme â€” pooled matching",
      subtitle = "Quarterly weighted mean injuries per km | dotted line = scheme start",
      x = NULL, y = "Weighted mean injuries per km",
      colour = NULL, linetype = NULL
    ) +
    theme_diag() +
    theme(
      legend.position = "bottom",
      axis.text.x     = element_text(size = AXIS_TEXT_SIZE - 2, angle = 45, hjust = 1),
      strip.text      = element_text(size = STRIP_TEXT_SIZE,     face = "bold"),
      panel.spacing   = unit(0.6, "lines")
    )
  
  save_fig(p_ts_scheme, "fig14c_parallel_trends_timeseries_per_scheme.png",
           width = 16, height = 14)
  
  # All schemes combined: treated vs control LOESS trajectories
  p_loess_all <- ggplot(
    oa_pre_quarterly %>%
      mutate(group = factor(group, levels = c("Treated", "Matched control"))),
    aes(x = q_date, y = wtd_mean_inj, colour = group, fill = group)
  ) +
    geom_point(size = POINT_SIZE - 1, alpha = 0.2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1.2) +
    scale_colour_manual(values = c("Treated" = COL_TREATED,
                                   "Matched control" = COL_CONTROL)) +
    scale_fill_manual(values = c("Treated" = COL_TREATED,
                                 "Matched control" = COL_CONTROL)) +
    labs(
      title    = "Pre-treatment injury trajectories: treated vs matched controls (all schemes)",
      subtitle = "LOESS smoother \u00b1 SE | parallel smoothers = parallel trends supported",
      x = NULL, y = "Weighted mean injuries per km",
      colour = NULL, fill = NULL
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  save_fig(p_loess_all, "fig14f_parallel_trends_loess_all_schemes.png",
           width = 12, height = 8)
  
  # Per-scheme LOESS: treated vs control trajectories
  p_loess_scheme <- ggplot(
    oa_pre_quarterly %>%
      mutate(group = factor(group, levels = c("Treated", "Matched control"))),
    aes(x = q_date, y = wtd_mean_inj, colour = group, fill = group)
  ) +
    geom_point(size = POINT_SIZE - 1, alpha = 0.3) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1.1) +
    facet_wrap(~scheme, ncol = 3, scales = "free_y") +
    scale_colour_manual(values = c("Treated" = COL_TREATED,
                                   "Matched control" = COL_CONTROL)) +
    scale_fill_manual(values = c("Treated" = COL_TREATED,
                                 "Matched control" = COL_CONTROL)) +
    labs(
      title    = "Pre-treatment injury trajectories: treated vs matched controls",
      subtitle = "LOESS smoother \u00b1 SE | parallel smoothers = parallel trends supported",
      x = NULL, y = "Weighted mean injuries per km",
      colour = NULL, fill = NULL
    ) +
    theme_diag() +
    theme(
      legend.position = "bottom",
      axis.text.x     = element_text(size = AXIS_TEXT_SIZE - 2, angle = 45, hjust = 1),
      strip.text      = element_text(size = STRIP_TEXT_SIZE,     face = "bold"),
      panel.spacing   = unit(0.6, "lines")
    )
  
  save_fig(p_loess_scheme, "fig14d_parallel_trends_loess_per_scheme.png",
           width = 16, height = 14)
}
cat("\n")

# â•”â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•—
# â•‘                  PART B â€” PER-SCHEME DIAGNOSTICS                    â•‘
# â•šâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART B â€” PER-SCHEME DIAGNOSTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Schemes
schemes <- matched_full %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("Schemes:", paste(schemes, collapse = ", "), "\n\n")

# Control â†’ scheme lookup from matching pairs (scheme column added by
# the per-scheme matching loop in script 16).
ctrl_scheme_lookup <- matching_pairs %>%
  select(OA = control_OA, scheme) %>%
  distinct()

# =============================================================================
# PER-SCHEME SMD COMPUTATION
# =============================================================================


smd_for_scheme <- function(scheme_name, vars) {
  treated_oas <- matched_full %>%
    filter(scheme == scheme_name, treat_indicator == 1) %>%
    pull(OA)
  
  control_oas <- ctrl_scheme_lookup %>%
    filter(scheme == scheme_name) %>%
    pull(OA)
  
  # Before: this scheme's treated vs full unmatched control pool
  before_df <- unmatched_pool %>%
    filter(treated_OA == 1 & OA %in% treated_oas | control_group2_OA == 1) %>%
    mutate(treat_indicator = as.integer(OA %in% treated_oas))
  
  # After: this scheme's treated + their matched controls
  after_df <- matched_full %>%
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
  smd_for_scheme(s, all_match_vars)
})

write_csv(smd_all_schemes, file.path(outdir, "09_smd_per_scheme.csv"))

# I also summarise the same balance groups within each scheme. This is the
# clearest way to show whether one scheme is driving a remaining imbalance.
smd_group_summary_by_scheme <- summarise_balance_groups(
  smd_all_schemes,
  group_vars = "scheme"
)

write_csv(
  smd_group_summary_by_scheme,
  file.path(outdir, "09b_smd_group_summary_per_scheme.csv")
)

# =============================================================================
# PER-SCHEME LOVE PLOTS (faceted)
# =============================================================================

cat("\n--- 10. Per-scheme love plots ---\n")

# Stage 2 variables + key Stage 1 variables
key_vars <- c("trend_total_pkm",
              "recent_minus_mid_total_pkm",
              "mid_minus_early_total_pkm",
              "covid_minus_precovid_total_pkm",
              "recovery_minus_covid_total_pkm",
              "log1p_mean_total_pkm",
              "log1p_mean_precovid_total_pkm",
              "log1p_mean_lockdown_total_pkm",
              "log1p_mean_recovery_total_pkm",
              "log1p_zero_quarter_share_pre",
              "log1p_road_density_m_km2", "log1p_road_length_km",
              "pct_A_road", "pct_minor_road",
              "log1p_dist_BUA_centroid", "log1p_pop_density", "IMD")

love_scheme_data <- smd_all_schemes %>%
  filter(variable %in% key_vars) %>%
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
  geom_point(size = POINT_SIZE) +
  scale_colour_manual(values = c("Before" = COL_BEFORE, "After" = COL_AFTER)) +
  scale_shape_manual(values = c("Before" = 16, "After" = 17)) +
  facet_wrap(~scheme, ncol = 4) +
  labs(
    title    = "Balance by scheme: key covariates",
    subtitle = "Before vs after matching",
    x = "|SMD|", y = NULL, colour = NULL, shape = NULL,
    caption = "Dashed = 0.10 threshold"
  ) +
  theme_diag() +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = AXIS_TEXT_SIZE - 1)
  )

save_fig(p_love_scheme, "fig10_love_plots_per_scheme.png",
         width = 16, height = 10)

# =============================================================================
#  SMD HEATMAP ACROSS SCHEMES (trends only)
# =============================================================================



smd_breaks  <- c(0, 0.05, 0.10, 0.15, 0.20, Inf)
smd_labels  <- c("<0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20", ">0.20")
smd_colours <- c("#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027")

s2_trend_heatmap <- smd_all_schemes %>%
  filter(variable %in% stage2_trends, !is.na(smd_after)) %>%
  mutate(
    scheme     = factor(scheme, levels = sort(unique(scheme))),
    smd_band   = cut(abs(smd_after), breaks = smd_breaks,
                     labels = smd_labels, right = FALSE, include.lowest = TRUE),
    cell_label = sprintf("%.3f", abs(smd_after))
  )

p_heatmap <- ggplot(s2_trend_heatmap,
                    aes(x = scheme, y = label, fill = smd_band)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label), size = CELL_LABEL_SIZE, colour = "#222222",
            fontface = "bold") +
  scale_fill_manual(
    values = setNames(smd_colours, smd_labels),
    name   = "|SMD| after matching",
    drop   = FALSE
  ) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(
    title    = "Injury trend balance \u2014 England",
    subtitle = "Each cell = |SMD| after matching",
    x = NULL, y = NULL,
    caption  = "Green < 0.05  \u00b7  yellow 0.05-0.10  \u00b7  orange 0.10-0.15  \u00b7  red > 0.20"
  ) +
  theme_diag() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.text.x      = element_text(size = AXIS_TEXT_SIZE, face = "bold"),
    axis.text.y      = element_text(size = AXIS_TEXT_SIZE),
    legend.position  = "bottom",
    legend.key.width = unit(1.4, "cm"),
    legend.text      = element_text(size = LEGEND_TEXT_SIZE)
  )

save_fig(p_heatmap, "fig11_smd_heatmap_per_scheme.png",
         width = max(8, length(schemes) * 1.6 + 4), height = 7)

# Full heatmap (all matching variables)
full_heatmap <- smd_all_schemes %>%
  filter(!is.na(smd_after)) %>%
  add_balance_group() %>%
  mutate(
    var_group = balance_group,
    scheme     = factor(scheme, levels = sort(unique(scheme))),
    smd_band   = cut(abs(smd_after), breaks = smd_breaks,
                     labels = smd_labels, right = FALSE, include.lowest = TRUE),
    cell_label = sprintf("%.2f", abs(smd_after))
  )

p_heatmap_full <- ggplot(full_heatmap,
                         aes(x = scheme, y = label, fill = smd_band)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label), size = CELL_LABEL_SIZE - 0.5,
            colour = "#222222", fontface = "bold") +
  scale_fill_manual(
    values = setNames(smd_colours, smd_labels),
    name   = "|SMD| after matching",
    drop   = FALSE
  ) +
  facet_wrap(~var_group, ncol = 1, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(
    title    = "Full balance heatmap by scheme (all matching variables)",
    subtitle = "Each cell = |SMD| after matching",
    x = NULL, y = NULL,
    caption  = "per-scheme |SMD| after matching"
  ) +
  theme_diag() +
  theme(
    axis.text.x      = element_text(size = AXIS_TEXT_SIZE, face = "bold"),
    axis.text.y      = element_text(size = AXIS_TEXT_SIZE - 2),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    strip.text       = element_text(size = STRIP_TEXT_SIZE, face = "bold"),
    legend.position  = "bottom",
    legend.key.width = unit(1.4, "cm"),
    legend.text      = element_text(size = LEGEND_TEXT_SIZE)
  )

save_fig(p_heatmap_full, "fig11b_smd_heatmap_full.png",
         width = 14, height = 22)

# =============================================================================
# 12. PER-SCHEME BALANCE SUMMARY TABLE
# =============================================================================

cat("--- 12. Per-scheme balance summary ---\n\n")

scheme_balance <- smd_all_schemes %>%
  filter(variable %in% stage2_vars_log) %>%
  group_by(scheme) %>%
  summarise(
    mean_stage2_smd_before = round(mean(abs(smd_before), na.rm = TRUE), 4),
    mean_stage2_smd_after  = round(mean(abs(smd_after),  na.rm = TRUE), 4),
    max_stage2_smd_after   = round(max_or_na(abs(smd_after)), 4),
    n_stage2_imbalanced    = sum(abs(smd_after) >= 0.10, na.rm = TRUE),
    all_stage2_balanced    = all(abs(smd_after) < 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    matched_full %>%
      filter(treat_indicator == 1) %>%
      count(scheme, name = "n_treated"),
    by = "scheme"
  ) %>%
  left_join(
    ctrl_scheme_lookup %>% count(scheme, name = "n_controls"),
    by = "scheme"
  )

write_csv(scheme_balance, file.path(outdir, "12_scheme_balance_summary.csv"))
write_csv(
  smd_group_summary_by_scheme,
  file.path(outdir, "12b_scheme_balance_by_reporting_group.csv")
)

cat("=== PER-SCHEME STAGE 2 BALANCE SUMMARY ===\n")
print(scheme_balance, n = Inf)

cat("\n=== PER-SCHEME BALANCE BY REPORTING GROUP ===\n")
print(smd_group_summary_by_scheme, n = Inf)



cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL DIAGNOSTICS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Output directory:", outdir, "\n")
cat("\nCSV outputs:\n")
list.files(outdir, pattern = "\\.csv$") %>%
  walk(~ cat("  ", .x, "\n"))
cat("\nFigure outputs:\n")
list.files(outdir, pattern = "\\.png$") %>%
  walk(~ cat("  ", .x, "\n"))
