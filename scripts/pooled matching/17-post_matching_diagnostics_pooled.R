# =============================================================================
# POST-MATCHING DIAGNOSTICS & DESCRIPTIVES — POOLED MATCHING
# =============================================================================
#
## matching pipeline (total injuries only, England only).
#
# PART A — Overall diagnostics (script 17 equivalent):
#   1. Descriptive summary tables (treated vs control)
#   2. SMD tables (Stage 1 + Stage 2)
#   3. Love plots (Stage 1 + Stage 2)
#   4. Weight distribution diagnostics + plots
#   5. Stratum characteristics
#   6. Common support / isolated OAs
#   7. Mahalanobis distance diagnostics
#   8. Parallel trends density plots
#
# PART B — Per-scheme diagnostics (script 18b equivalent):
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

# Central text-size constants — change these to rescale everything uniformly
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

# Stage 1 (same as original script 16)
stage1_vars_raw <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_BUA_centroid", "pop_density", "area_km2",
  "business_retail_per_km2", "IMD",
  "cars_one_pct", "cars_twoPlus_pct",
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct",
  "age_under15_pct", "age_15to24_pct", "age_25to44_pct",
  "age_45to64_pct", "age_65to84_pct"
)

log_transform_s1 <- c("road_length_km", "pop_density", "dist_BUA_centroid",
                      "road_density_m_km2", "business_retail_per_km2")
log_nozero_s1    <- c("area_km2")

stage1_vars_log <- c(
  paste0("log1p_", log_transform_s1),
  paste0("log_", log_nozero_s1),
  setdiff(stage1_vars_raw, c(log_transform_s1, log_nozero_s1))
)

# Stage 2 — pooled (only 2 variables)
stage2_trends <- "trend_total_pkm"
stage2_levels <- "mean_total_pkm"
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
  "workAthome_pct"             = "% Work from home",
  "Other_pct"                  = "% Other commute",
  "White_pct"                  = "% White",
  "Mixed_pct"                  = "% Mixed ethnicity",
  "Asian_pct"                  = "% Asian",
  "Black_pct"                  = "% Black",
  "age_under15_pct"            = "% Under 15",
  "age_15to24_pct"             = "% 15-24",
  "age_25to44_pct"             = "% 25-44",
  "age_45to64_pct"             = "% 45-64",
  "age_65to84_pct"             = "% 65-84",
  "trend_total_pkm"            = "Trend: total injuries/km",
  "log1p_mean_total_pkm"       = "Level: mean total injuries/km (log)"
)

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

matched_full <- add_log_vars(matched_full)

# Unmatched pool (for before-matching comparisons)
unmatched_pool <- full_data %>%
  filter(
    substr(LAD24CD, 1, 1) == "E",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA == 0,
    n_roads   > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1)) %>%
  add_log_vars()

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

# ╔═══════════════════════════════════════════════════════════════════════╗
# ║                     PART A — OVERALL DIAGNOSTICS                    ║
# ╚═══════════════════════════════════════════════════════════════════════╝

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART A — OVERALL POST-MATCHING DIAGNOSTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# 0. SCHEME-LEVEL MATCHING SUMMARY
# =============================================================================

cat("--- 0. Scheme-level matching summary ---\n")

scheme_summary_table <- map_df(seq_along(all_results), function(i) {
  md <- all_results[[i]]$matched_data
  m  <- all_results[[i]]$matchit_s2
  n_t <- sum(md$treat_indicator == 1)
  tibble(
    scheme          = unique(md$scheme[!is.na(md$scheme)]),
    treated         = n_t,
    unique_controls = n_distinct(md$OA[md$treat_indicator == 0]),
    total_pairs     = n_t * m$info$ratio,
    ratio           = paste0("1:", m$info$ratio)
  )
})

scheme_summary_table <- bind_rows(
  scheme_summary_table,
  tibble(scheme = "Total",
         treated = sum(scheme_summary_table$treated),
         unique_controls = sum(scheme_summary_table$unique_controls),
         total_pairs = sum(scheme_summary_table$total_pairs),
         ratio = "")
)

print(scheme_summary_table)
write_csv(scheme_summary_table, file.path(outdir, "00_scheme_matching_summary.csv"))
cat("  Saved: 00_scheme_matching_summary.csv\n\n")

# =============================================================================
# 1. DESCRIPTIVE SUMMARY TABLE
# =============================================================================

cat("--- 1. Descriptive summary ---\n")

desc_raw_vars <- c(stage1_vars_raw, stage2_trends, stage2_levels)
desc_raw_vars <- intersect(desc_raw_vars, names(matched_full))

desc_table <- map_df(desc_raw_vars, ~ desc_stats(matched_full, .x))
write_csv(desc_table, file.path(outdir, "01_descriptive_table.csv"))
cat("  Saved: 01_descriptive_table.csv\n")

# =============================================================================
# 2. SMD TABLES
# =============================================================================

cat("--- 2. SMD tables ---\n")

smd_table <- function(vars, before_data, after_data, label) {
  vars_avail <- intersect(vars, intersect(names(before_data), names(after_data)))
  map_df(vars_avail, function(v) {
    tibble(
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

cat("  Stage 1: mean |SMD| before =", round(mean(abs(smd_s1$smd_before), na.rm = TRUE), 4),
    "| after =", round(mean(abs(smd_s1$smd_after), na.rm = TRUE), 4), "\n")
cat("  Stage 2: mean |SMD| before =", round(mean(abs(smd_s2$smd_before), na.rm = TRUE), 4),
    "| after =", round(mean(abs(smd_s2$smd_after), na.rm = TRUE), 4), "\n")
cat("  Stage 2 all balanced (<0.10):", all(smd_s2$balanced, na.rm = TRUE), "\n\n")

# =============================================================================
# 3. LOVE PLOTS
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
  "Stage 1 balance — structural & sociodemographic covariates",
  "Pooled matching (England, other-city controls only)"
)
save_fig(p_love_s1, "fig01_love_plot_stage1.png", width = 14, height = 16)

p_love_s2 <- make_love_plot(
  smd_s2,
  "Stage 2 balance — total injury trend + level",
  "Pooled matching — 2 variables (trend_total_pkm + log1p_mean_total_pkm)"
)
save_fig(p_love_s2, "fig02_love_plot_stage2.png", width = 10, height = 5)

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
# 5. STRATUM CHARACTERISTICS
# =============================================================================

cat("--- 5. Stratum characteristics ---\n")

stratum_vars <- intersect(
  c("mean_total_pkm", "road_length_km", "road_density_m_km2",
    "pop_density", "dist_BUA_centroid", "IMD",
    "pct_A_road", "pct_minor_road", "Drive_Car_pct", "Walk_pct"),
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
# 7. MAHALANOBIS DISTANCE
# =============================================================================

cat("--- 7. Mahalanobis distance ---\n")

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
      subtitle = "2 Stage 2 variables",
      x = "Mahalanobis distance", y = "Cumulative proportion",
      colour = NULL
    ) +
    theme_diag() +
    theme(legend.position = "bottom")
  
  save_fig(p_mdist, "fig07_mahalanobis_distance.png", width = 10, height = 7)
} else {
  cat("  No mdist column in matched data — skipping.\n")
}
cat("\n")

# =============================================================================
# 8. PARALLEL TRENDS DENSITY PLOTS
# =============================================================================

cat("--- 8. Parallel trends distributions ---\n")

trend_vars <- intersect(c("trend_total_pkm"), names(matched_full))

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
      subtitle = "total injuries per road-km",
      x = "Pre-treatment slope", y = "Density",
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
      subtitle = "Total injuries per road-km",
      x = "Pre-treatment slope", y = "Density",
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
      title    = "Pre-treatment injury trends by scheme — pooled matching",
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

# ╔═══════════════════════════════════════════════════════════════════════╗
# ║                  PART B — PER-SCHEME DIAGNOSTICS                    ║
# ╚═══════════════════════════════════════════════════════════════════════╝

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PART B — PER-SCHEME DIAGNOSTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Schemes
schemes <- matched_full %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("Schemes:", paste(schemes, collapse = ", "), "\n\n")

# Control → scheme lookup from matching pairs (scheme column added by
# the per-scheme matching loop in script 16).
ctrl_scheme_lookup <- matching_pairs %>%
  select(OA = control_OA, scheme) %>%
  distinct()

# =============================================================================
# 9. PER-SCHEME SMD COMPUTATION
# =============================================================================

cat("--- 9. Per-scheme SMD ---\n")

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

# =============================================================================
# 10. PER-SCHEME LOVE PLOTS (faceted)
# =============================================================================

cat("\n--- 10. Per-scheme love plots ---\n")

# Stage 2 variables + key Stage 1 variables
key_vars <- c("trend_total_pkm", "log1p_mean_total_pkm",
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
# 11. SMD HEATMAP ACROSS SCHEMES (trends only)
# =============================================================================

cat("--- 11. SMD heatmap ---\n")

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
  mutate(
    var_group = case_when(
      variable %in% stage2_vars_log ~ "Stage 2: Injury",
      variable %in% paste0("log1p_", log_transform_s1) |
        variable %in% paste0("log_", log_nozero_s1) |
        variable %in% c("pct_A_road", "pct_B_road", "pct_minor_road") ~ "Road network",
      variable %in% c("log1p_dist_BUA_centroid", "log1p_pop_density",
                      "log1p_business_retail_per_km2") ~ "Urban form",
      TRUE ~ "Sociodemographic"
    ),
    var_group = factor(var_group,
                       levels = c("Stage 2: Injury", "Road network",
                                  "Urban form", "Sociodemographic")),
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

cat("=== PER-SCHEME TREND BALANCE SUMMARY ===\n")
print(scheme_balance, n = Inf)

# =============================================================================
# DONE
# =============================================================================

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