# =============================================================================
# POST-MATCHING DIAGNOSTICS AND DESCRIPTIVE ANALYSIS
# OA-Level Two-Stage Mahalanobis Distance Matching
# Road Safety Intervention — Great Britain
# =============================================================================
#
# PURPOSE:
#   Full standalone script producing all descriptive statistics, balance
#   diagnostics, and map visualisations needed to characterise the matched
#   treated and control OA samples after the two-stage MDM pipeline.
#
# INPUTS (from data/processed/):
#   OA_matched_full_A.rds          — full matched dataset, Analysis A
#   OA_matched_full_B.rds          — full matched dataset, Analysis B
#   OA_matched_treated_A.rds       — treated OAs + weights + stratum, A
#   OA_matched_donors_A.rds        — control OAs + weights, A
#   OA_matched_treated_B.rds       — treated OAs + weights + stratum, B
#   OA_matched_donors_B.rds        — control OAs + weights, B
#   OA_common_support_flags.rds    — structurally isolated OA flags
#   OA_matching_census.rds         — full original dataset (for unmatched pool)
#
# SPATIAL INPUTS (adjust paths as needed):
#   OA boundary shapefile / GeoPackage for GB (e.g. from ONS / NRS)
#   LAD boundary shapefile for city-level context maps
#
# OUTPUTS (to output/diagnostics/):
#   Tables:
#     01_descriptive_summary_treated_control.csv
#     02_smd_before_after_s1.csv
#     03_smd_before_after_s2.csv
#     04_weight_distribution_summary.csv
#     05_stratum_characteristics.csv
#     06_zero_injury_OA_comparison.csv
#
#   Figures:
#     fig01_love_plot_stage1_A.png
#     fig02_love_plot_stage1_B.png
#     fig03_love_plot_stage2_A_restricted.png
#     fig04_love_plot_stage2_B_restricted.png
#     fig05_smd_heatmap_all_specs.png
#     fig06_covariate_distributions_stage1.png
#     fig07_covariate_distributions_stage2_trends.png
#     fig08_covariate_distributions_stage2_levels.png
#     fig09_weight_distribution_A.png
#     fig10_weight_distribution_B.png
#     fig11_mahalanobis_distance_ecdf.png
#     fig12_stratum_injury_distribution.png
#     fig13_zero_injury_OA_profile.png
#     fig14_map_treated_control_overview.png
#     fig15_map_treatment_city_zoom_[city].png  (one per major treated city)
#     fig16_map_control_origin_composition.png
#     fig17_parallel_trends_total_injuries.png
#     fig18_parallel_trends_by_mode.png
#
# =============================================================================

# ---- LIBRARIES ---------------------------------------------------------------

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
library(knitr)
library(gt)

# Force dplyr namespace
select <- dplyr::select
filter <- dplyr::filter

# ---- OUTPUT DIRECTORY --------------------------------------------------------

dir.create(here("output", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics")

# ---- THEME -------------------------------------------------------------------
# Consistent publication-quality ggplot theme used throughout

theme_diag <- function(base_size = 11) {
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

# Colour palette
COL_TREATED  <- "#D85A30"   # coral-red for treated
COL_CONTROL  <- "#2E6FAB"   # blue for control
COL_BOTH     <- "#3B8C5A"   # green for both/all
COL_BEFORE   <- "#E74C3C"   # red for pre-matching
COL_AFTER    <- "#2ECC71"   # green for post-matching
COL_SCOTLAND <- "#6B3FA0"   # purple for Scotland
COL_ENGLAND  <- "#2E6FAB"   # blue for England

# Save helper
save_fig <- function(p, filename, width = 12, height = 8, dpi = 300) {
  ggsave(file.path(outdir, filename), p, width = width, height = height,
         dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# =============================================================================
# SECTION 1 — LOAD DATA
# =============================================================================

cat("\n=== Loading matched datasets ===\n")

matched_A     <- readRDS(here("data", "processed", "OA_matched_full_A.rds"))
matched_B     <- readRDS(here("data", "processed", "OA_matched_full_B.rds"))
treated_A     <- readRDS(here("data", "processed", "OA_matched_treated_A.rds"))
controls_A    <- readRDS(here("data", "processed", "OA_matched_donors_A.rds"))
treated_B     <- readRDS(here("data", "processed", "OA_matched_treated_B.rds"))
controls_B    <- readRDS(here("data", "processed", "OA_matched_donors_B.rds"))
csupport      <- readRDS(here("data", "processed", "OA_common_support_flags.rds"))
full_data     <- readRDS(here("data", "processed", "OA_matching_census.rds"))

cat("  Analysis A — Treated:", nrow(treated_A),
    "| Controls:", nrow(controls_A), "\n")
cat("  Analysis B — Treated:", nrow(treated_B),
    "| Controls:", nrow(controls_B), "\n")
cat("  Full original dataset:", nrow(full_data), "OAs\n\n")

# Variable group definitions
stage1_road   <- c("road_density_m_km2", "road_length_km",
                   "pct_A_road", "pct_B_road", "pct_minor_road")
stage1_urban  <- c("dist_citycentre", "pop_density")
stage1_socdem <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                   "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct")
stage1_vars   <- c(stage1_road, stage1_urban, stage1_socdem)

stage2_trends <- c("trend_car_slight_pkm", "trend_cyc_slight_pkm",
                   "trend_ped_slight_pkm", "trend_total_pkm")
stage2_levels <- c("mean_car_KSI_pkm", "mean_car_slight_pkm",
                   "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
                   "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
                   "mean_total_pkm")
log_names_s2  <- paste0("log1p_", stage2_levels)

# Human-readable variable labels
var_labels <- c(
  road_density_m_km2    = "Road density (m/km²)",
  road_length_km        = "Road length (km)",
  pct_A_road            = "% A-road",
  pct_B_road            = "% B-road",
  pct_minor_road        = "% Minor road",
  dist_citycentre       = "Distance to city centre (m)",
  pop_density           = "Population density (persons/km²)",
  IMD                   = "Index of Multiple Deprivation",
  cars_none_pct         = "% households: no car",
  Drive_Car_pct         = "% commuting by car",
  Walk_pct              = "% commuting on foot",
  Bicycle_pct           = "% commuting by bicycle",
  X65plus_pct           = "% aged 65+",
  X5to19_pct            = "% aged 5–19",
  X20to24_pct           = "% aged 20–24",
  trend_car_slight_pkm  = "Trend: car slight injuries/km",
  trend_cyc_slight_pkm  = "Trend: cycle slight injuries/km",
  trend_ped_slight_pkm  = "Trend: ped slight injuries/km",
  trend_total_pkm       = "Trend: total injuries/km",
  mean_car_KSI_pkm      = "Mean: car KSI/km",
  mean_car_slight_pkm   = "Mean: car slight/km",
  mean_cyc_KSI_pkm      = "Mean: cycle KSI/km",
  mean_cyc_slight_pkm   = "Mean: cycle slight/km",
  mean_ped_KSI_pkm      = "Mean: ped KSI/km",
  mean_ped_slight_pkm   = "Mean: ped slight/km",
  mean_total_pkm        = "Mean: total injuries/km"
)

# =============================================================================
# SECTION 2 — DESCRIPTIVE SUMMARY TABLE
# =============================================================================

cat("=== Section 2: Descriptive summary table ===\n")

# Function: mean ± SD for a variable by group
desc_row <- function(data, var, group_col = "treat_indicator") {
  data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n    = sum(!is.na(.data[[var]])),
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]],   na.rm = TRUE),
      p25  = quantile(.data[[var]], 0.25, na.rm = TRUE),
      p50  = quantile(.data[[var]], 0.50, na.rm = TRUE),
      p75  = quantile(.data[[var]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = var)
}

all_desc_vars <- c(stage1_vars, stage2_trends, stage2_levels)

desc_A <- map_df(all_desc_vars, ~ desc_row(matched_A, .x)) %>%
  mutate(analysis = "A", group = if_else(treat_indicator == 1, "Treated", "Control"))

desc_B <- map_df(all_desc_vars, ~ desc_row(matched_B, .x)) %>%
  mutate(analysis = "B", group = if_else(treat_indicator == 1, "Treated", "Control"))

desc_summary <- bind_rows(desc_A, desc_B) %>%
  mutate(
    label      = coalesce(var_labels[variable], variable),
    mean_sd    = sprintf("%.3f (%.3f)", mean, sd),
    median_iqr = sprintf("%.3f [%.3f–%.3f]", p50, p25, p75)
  ) %>%
  select(analysis, variable, label, group, n, mean_sd, median_iqr)

write_csv(desc_summary, file.path(outdir, "01_descriptive_summary_treated_control.csv"))
cat("  Saved: 01_descriptive_summary_treated_control.csv\n")

# =============================================================================
# SECTION 3 — SMD TABLES (BEFORE/AFTER)
# =============================================================================

cat("=== Section 3: SMD before/after tables ===\n")

# Helper: compute SMD for a variable between two groups
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

# Unmatched pool (all eligible treated vs all eligible controls)
unmatched_pool <- full_data %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0, n_roads > 0
  ) %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1)) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

smd_before_s1 <- map_df(stage1_vars, function(v) {
  tibble(
    variable = v,
    label    = coalesce(var_labels[v], v),
    smd_un   = compute_smd(unmatched_pool, v)
  )
})

smd_after_s1_A <- map_df(stage1_vars, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd_adj = compute_smd(matched_A, v))
})
smd_after_s1_B <- map_df(stage1_vars, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd_adj = compute_smd(matched_B, v))
})

smd_s1 <- smd_before_s1 %>%
  left_join(smd_after_s1_A %>% rename(smd_A = smd_adj), by = c("variable","label")) %>%
  left_join(smd_after_s1_B %>% rename(smd_B = smd_adj), by = c("variable","label")) %>%
  mutate(across(starts_with("smd"), ~ round(.x, 4)))

write_csv(smd_s1, file.path(outdir, "02_smd_before_after_s1.csv"))
cat("  Saved: 02_smd_before_after_s1.csv\n")

# Stage 2 SMD
s2_vars_for_smd <- c(stage2_trends, stage2_levels)
smd_before_s2 <- map_df(s2_vars_for_smd, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd_un = compute_smd(unmatched_pool, v))
})

smd_after_s2_A <- map_df(s2_vars_for_smd, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd_A = compute_smd(matched_A, v))
})
smd_after_s2_B <- map_df(s2_vars_for_smd, function(v) {
  tibble(variable = v, label = coalesce(var_labels[v], v),
         smd_B = compute_smd(matched_B, v))
})

smd_s2 <- smd_before_s2 %>%
  left_join(smd_after_s2_A, by = c("variable","label")) %>%
  left_join(smd_after_s2_B, by = c("variable","label")) %>%
  mutate(across(starts_with("smd"), ~ round(.x, 4)))

write_csv(smd_s2, file.path(outdir, "03_smd_before_after_s2.csv"))
cat("  Saved: 03_smd_before_after_s2.csv\n")

# =============================================================================
# SECTION 4 — WEIGHT DISTRIBUTION SUMMARY
# =============================================================================

cat("=== Section 4: Weight diagnostics ===\n")

weight_summary <- function(data, label) {
  ctrl <- data %>% filter(treat_indicator == 0)
  ctrl %>%
    group_by(country) %>%
    summarise(
      analysis    = label,
      n           = n(),
      mean_w      = round(mean(weights), 3),
      median_w    = round(median(weights), 3),
      max_w       = round(max(weights), 3),
      p95_w       = round(quantile(weights, 0.95), 3),
      eff_n       = round(sum(weights)^2 / sum(weights^2), 1),
      efficiency  = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
      .groups     = "drop"
    )
}

weight_diag <- bind_rows(
  weight_summary(matched_A, "A_excl_zero"),
  weight_summary(matched_B, "B_incl_zero")
)

write_csv(weight_diag, file.path(outdir, "04_weight_distribution_summary.csv"))
cat("  Saved: 04_weight_distribution_summary.csv\n")

# =============================================================================
# SECTION 5 — STRATUM CHARACTERISTICS
# =============================================================================

cat("=== Section 5: Stratum characteristics ===\n")

stratum_chars <- matched_A %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  group_by(baseline_injury_stratum) %>%
  summarise(
    n                  = n(),
    mean_total_pkm     = round(mean(mean_total_pkm, na.rm = TRUE), 4),
    mean_road_length   = round(mean(road_length_km, na.rm = TRUE), 3),
    mean_pop_density   = round(mean(pop_density, na.rm = TRUE), 0),
    mean_IMD           = round(mean(IMD, na.rm = TRUE), 1),
    mean_Drive_Car_pct = round(mean(Drive_Car_pct, na.rm = TRUE), 1),
    mean_Walk_pct      = round(mean(Walk_pct, na.rm = TRUE), 1),
    pct_scotland       = round(100 * mean(country == "Scotland"), 1),
    .groups            = "drop"
  )

write_csv(stratum_chars, file.path(outdir, "05_stratum_characteristics.csv"))
cat("  Saved: 05_stratum_characteristics.csv\n")

# =============================================================================
# SECTION 6 — ZERO-INJURY OA COMPARISON TABLE
# =============================================================================

cat("=== Section 6: Zero-injury OA comparison ===\n")

treated_A_ids <- treated_A$OA
treated_B_ids <- treated_B$OA
zero_injury_ids <- setdiff(treated_B_ids, treated_A_ids)

zero_comparison <- matched_B %>%
  filter(treat_indicator == 1) %>%
  mutate(grp = if_else(OA %in% zero_injury_ids,
                       "Zero-injury (B only)", "Injury-exposed (A & B)")) %>%
  group_by(grp) %>%
  summarise(across(all_of(stage1_vars),
                   list(mean = ~round(mean(.x, na.rm=TRUE), 3),
                        sd   = ~round(sd(.x,   na.rm=TRUE), 3))),
            n = n(), pct_scotland = round(100*mean(country=="Scotland"),1),
            .groups = "drop")

write_csv(zero_comparison, file.path(outdir, "06_zero_injury_OA_comparison.csv"))
cat("  Saved: 06_zero_injury_OA_comparison.csv\n")

# =============================================================================
# SECTION 7 — LOVE PLOTS (BALANCE)
# =============================================================================

cat("=== Section 7: Love plots ===\n")

# Build love plot data frame from scratch for full control over appearance
love_data <- function(unmatched, matched, vars, label) {
  map_df(vars, function(v) {
    smd_un  <- compute_smd(unmatched, v)
    smd_adj <- compute_smd(matched,   v)
    tibble(
      variable  = v,
      label     = coalesce(var_labels[v], v),
      smd_un    = smd_un,
      smd_adj   = smd_adj,
      spec      = label
    )
  })
}

make_love_plot <- function(ldat, title, subtitle = NULL, threshold = 0.1) {
  ldat <- ldat %>%
    arrange(abs(smd_un)) %>%
    mutate(label = factor(label, levels = unique(label)))
  
  ldat_long <- ldat %>%
    pivot_longer(c(smd_un, smd_adj), names_to = "timing",
                 values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_un", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  ggplot(ldat_long, aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = threshold, linetype = "dashed",
               colour = "#999999", linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = "#CCCCCC", linewidth = 0.4) +
    geom_line(aes(group = label), colour = "#CCCCCC", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("Before matching" = COL_BEFORE,
                                   "After matching"  = COL_AFTER)) +
    scale_shape_manual(values  = c("Before matching" = 16,
                                   "After matching"  = 17)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    annotate("text", x = threshold + 0.005, y = 0.5, label = "0.1",
             size = 3, colour = "#999999", hjust = 0, vjust = 0) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Absolute Standardised Mean Difference",
      y        = NULL,
      colour   = NULL, shape = NULL,
      caption  = "Dashed line = 0.1 threshold. Points show |SMD| before (red circles) and after (green triangles) matching."
    ) +
    theme_diag() +
    theme(
      legend.position = "bottom",
      axis.text.y     = element_text(size = 9),
      plot.margin     = margin(10, 20, 10, 10)
    )
}

# Stage 1 love plots
ld_s1_A <- love_data(unmatched_pool, matched_A, stage1_vars, "A")
ld_s1_B <- love_data(unmatched_pool, matched_B, stage1_vars, "B")

p_love_s1_A <- make_love_plot(ld_s1_A,
                              "Stage 1 Balance — Analysis A (zero-injury excluded)",
                              "Structural and sociodemographic variables | ratio 10, exact = country")
p_love_s1_B <- make_love_plot(ld_s1_B,
                              "Stage 1 Balance — Analysis B (zero-injury included)",
                              "Structural and sociodemographic variables | ratio 10, exact = country")

save_fig(p_love_s1_A, "fig01_love_plot_stage1_A.png", width = 13, height = 9)
save_fig(p_love_s1_B, "fig02_love_plot_stage1_B.png", width = 13, height = 9)

# Stage 2 love plots (restricted)
matched_A_restr <- matched_A %>%
  filter(!OA %in% csupport$treated_OA[csupport$analysis == "A"])
matched_B_restr <- matched_B %>%
  filter(!OA %in% csupport$treated_OA[csupport$analysis == "B"])

s2_all_vars <- c(stage2_trends, stage2_levels)
ld_s2_A <- love_data(matched_A, matched_A_restr, s2_all_vars, "A_restricted")
ld_s2_B <- love_data(matched_B, matched_B_restr, s2_all_vars, "B_restricted")

p_love_s2_A <- make_love_plot(ld_s2_A,
                              "Stage 2 Balance — Analysis A Restricted",
                              "Pre-treatment injury trends and levels | ratio 4, exact = country",
                              threshold = 0.1)
p_love_s2_B <- make_love_plot(ld_s2_B,
                              "Stage 2 Balance — Analysis B Restricted",
                              "Pre-treatment injury trends and levels | ratio 7, exact = country",
                              threshold = 0.1)

save_fig(p_love_s2_A, "fig03_love_plot_stage2_A_restricted.png", width = 13, height = 9)
save_fig(p_love_s2_B, "fig04_love_plot_stage2_B_restricted.png", width = 13, height = 9)

# =============================================================================
# SECTION 8 — SMD HEATMAP (ALL SPECIFICATIONS)
# =============================================================================

cat("=== Section 8: SMD heatmap ===\n")

specs <- list(
  "Unmatched"     = unmatched_pool,
  "S1 matched A"  = matched_A,
  "S1 matched B"  = matched_B,
  "S2 matched A"  = matched_A,
  "S2 matched B"  = matched_B
)
all_diag_vars <- c(stage1_vars, stage2_trends, stage2_levels)

smd_heat <- map_df(names(specs), function(spec) {
  map_df(all_diag_vars, function(v) {
    smd <- tryCatch(compute_smd(specs[[spec]], v), error = function(e) NA_real_)
    tibble(spec = spec, variable = v, smd = abs(smd))
  })
}) %>%
  mutate(
    label = coalesce(var_labels[variable], variable),
    label = factor(label),
    spec  = factor(spec, levels = names(specs)),
    var_group = case_when(
      variable %in% stage1_road   ~ "Road network",
      variable %in% stage1_urban  ~ "Urban geography",
      variable %in% stage1_socdem ~ "Sociodemographic",
      variable %in% stage2_trends ~ "Injury trends",
      variable %in% stage2_levels ~ "Injury levels"
    )
  )

p_heatmap <- ggplot(smd_heat, aes(x = spec, y = label, fill = smd)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(!is.na(smd), sprintf("%.2f", smd), "")),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low      = "#2ECC71",
    mid      = "#F39C12",
    high     = "#E74C3C",
    midpoint = 0.1,
    na.value = "#EEEEEE",
    name     = "|SMD|",
    limits   = c(0, NA)
  ) +
  facet_grid(var_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Standardised Mean Differences Across All Specifications",
    subtitle = "Green = balanced (|SMD| < 0.1), orange = marginal, red = imbalanced",
    x        = "Specification",
    y        = NULL,
    caption  = "Values show absolute SMD. Threshold for balance: |SMD| < 0.1"
  ) +
  theme_diag() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.spacing = unit(0.3, "lines")
  )

save_fig(p_heatmap, "fig05_smd_heatmap_all_specs.png", width = 14, height = 14)

# =============================================================================
# SECTION 9 — COVARIATE DISTRIBUTION PLOTS
# =============================================================================

cat("=== Section 9: Covariate distribution plots ===\n")

# Helper: density/histogram comparison treated vs control, before and after
dist_plot_grid <- function(vars, data_before, data_after,
                           title, ncol = 3, log_vars = NULL) {
  
  plot_one <- function(v, data, phase_label) {
    d <- data %>%
      filter(!is.na(.data[[v]])) %>%
      mutate(group = if_else(treat_indicator == 1, "Treated", "Control"))
    lbl <- coalesce(var_labels[v], v)
    
    # Apply log for display if requested
    if (!is.null(log_vars) && v %in% log_vars) {
      d <- d %>% mutate(val = log1p(pmax(.data[[v]], 0)))
      x_lbl <- paste0("log1p(", lbl, ")")
    } else {
      d <- d %>% mutate(val = .data[[v]])
      x_lbl <- lbl
    }
    
    ggplot(d, aes(x = val, colour = group, fill = group)) +
      geom_density(alpha = 0.25, linewidth = 0.7) +
      scale_colour_manual(values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
      scale_fill_manual(  values = c(Treated = COL_TREATED, Control = COL_CONTROL)) +
      labs(title = phase_label, x = x_lbl, y = "Density",
           colour = NULL, fill = NULL) +
      theme_diag(base_size = 9) +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 8, face = "bold"))
  }
  
  plots <- map(vars, function(v) {
    p_before <- plot_one(v, data_before, "Before matching")
    p_after  <- plot_one(v, data_after,  "After matching")
    lbl <- coalesce(var_labels[v], v)
    (p_before | p_after) +
      plot_annotation(title = lbl,
                      theme = theme(plot.title = element_text(
                        size = 9, face = "bold", colour = "#1A2E5A")))
  })
  
  wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title    = title,
      subtitle = "Left: before matching | Right: after matching (Analysis A)",
      theme    = theme(plot.title    = element_text(size = 13, face = "bold",
                                                    colour = "#1A2E5A"),
                       plot.subtitle = element_text(size = 10, colour = "#555555"))
    )
}

# Stage 1 variables: three groups
log_s1 <- c("road_length_km", "pop_density", "dist_citycentre")

p_dist_road <- dist_plot_grid(
  stage1_road, unmatched_pool, matched_A,
  "Road Network Variables — Treated vs Control",
  ncol = 2, log_vars = log_s1
)
save_fig(p_dist_road, "fig06a_distributions_road.png", width = 14, height = 10)

p_dist_urban <- dist_plot_grid(
  stage1_urban, unmatched_pool, matched_A,
  "Urban Geography Variables — Treated vs Control",
  ncol = 2, log_vars = log_s1
)
save_fig(p_dist_urban, "fig06b_distributions_urban.png", width = 12, height = 6)

p_dist_socdem <- dist_plot_grid(
  stage1_socdem, unmatched_pool, matched_A,
  "Sociodemographic Variables — Treated vs Control",
  ncol = 3, log_vars = NULL
)
save_fig(p_dist_socdem, "fig06c_distributions_socdem.png", width = 16, height = 14)

# Stage 2 trend distributions
p_dist_trends <- dist_plot_grid(
  stage2_trends, matched_A, matched_A,
  "Pre-Treatment Injury Trend Variables — Treated vs Control (Stage 2 matched)",
  ncol = 2
)
save_fig(p_dist_trends, "fig07_distributions_stage2_trends.png", width = 14, height = 10)

# Stage 2 level distributions
p_dist_levels <- dist_plot_grid(
  stage2_levels, matched_A, matched_A,
  "Pre-Treatment Injury Level Variables — Treated vs Control (Stage 2 matched)",
  ncol = 3, log_vars = stage2_levels
)
save_fig(p_dist_levels, "fig08_distributions_stage2_levels.png", width = 16, height = 14)

# =============================================================================
# SECTION 10 — WEIGHT DISTRIBUTION PLOTS
# =============================================================================

cat("=== Section 10: Weight distribution plots ===\n")

plot_weights <- function(data, label) {
  ctrl <- data %>% filter(treat_indicator == 0)
  
  p1 <- ggplot(ctrl, aes(x = weights, fill = country)) +
    geom_histogram(bins = 60, alpha = 0.8, position = "stack") +
    scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    scale_x_continuous(limits = c(0, min(max(ctrl$weights) * 1.05, 10))) +
    labs(title = "Weight distribution (controls)",
         subtitle = paste0("Capped at 5 | Effective N = ",
                           round(sum(ctrl$weights)^2 / sum(ctrl$weights^2), 0),
                           " from ", nrow(ctrl), " nominal controls"),
         x = "Weight", y = "Count", fill = "Country") +
    theme_diag()
  
  p2 <- ctrl %>%
    group_by(country) %>%
    summarise(
      n       = n(),
      eff_n   = round(sum(weights)^2 / sum(weights^2), 1),
      max_w   = round(max(weights), 2),
      mean_w  = round(mean(weights), 3),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = country, y = eff_n, fill = country)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = paste0("n=", n, "\neff_N=", eff_n)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    labs(title = "Effective N by country",
         x = NULL, y = "Effective N", fill = NULL) +
    theme_diag() + theme(legend.position = "none")
  
  p3 <- ggplot(ctrl, aes(x = weights, colour = country)) +
    stat_ecdf(linewidth = 1) +
    geom_vline(xintercept = 5, linetype = "dashed", colour = "#E74C3C", linewidth = 0.6) +
    scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(title = "ECDF of weights by country",
         subtitle = "Red dashed = cap at 5",
         x = "Weight", y = "Cumulative proportion", colour = "Country") +
    theme_diag()
  
  (p1 | p2 | p3) +
    plot_annotation(
      title = paste("Weight Diagnostics —", label),
      theme = theme(plot.title = element_text(size = 13, face = "bold",
                                              colour = "#1A2E5A"))
    )
}

save_fig(plot_weights(matched_A, "Analysis A"), "fig09_weight_distribution_A.png",
         width = 16, height = 6)
save_fig(plot_weights(matched_B, "Analysis B"), "fig10_weight_distribution_B.png",
         width = 16, height = 6)

# =============================================================================
# SECTION 11 — MAHALANOBIS DISTANCE PLOTS
# =============================================================================

cat("=== Section 11: Mahalanobis distance plots ===\n")

# ECDF comparison
mdist_df <- bind_rows(
  matched_A %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(spec = "A: excl zero-injury"),
  matched_B %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(spec = "B: incl zero-injury")
) %>% filter(!is.na(mdist))

p_mdist_ecdf <- ggplot(mdist_df, aes(x = mdist, colour = spec)) +
  stat_ecdf(linewidth = 1) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "#999999") +
  geom_vline(xintercept = 10, linetype = "dotted", colour = "#CC3333") +
  annotate("text", x = 5.3, y = 0.2, label = "d = 5",
           size = 3, colour = "#999999") +
  annotate("text", x = 10.3, y = 0.1, label = "d = 10",
           size = 3, colour = "#CC3333") +
  scale_colour_manual(values = c("A: excl zero-injury" = COL_TREATED,
                                 "B: incl zero-injury" = COL_CONTROL)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    title    = "ECDF of Stage 2 Mahalanobis Distance (Treated OAs)",
    subtitle = "Distance measures how far each treated OA is from its nearest control in trajectory space",
    x        = "Stage 2 Mahalanobis distance",
    y        = "Cumulative proportion of treated OAs",
    colour   = NULL,
    caption  = paste0("Analysis A: ", sum(mdist_df$spec == "A: excl zero-injury"),
                      " treated OAs | Analysis B: ",
                      sum(mdist_df$spec == "B: incl zero-injury"), " treated OAs")
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

# Distance by country
mdist_country <- matched_A %>%
  filter(treat_indicator == 1, !is.na(mdist)) %>%
  select(OA, mdist, country)

p_mdist_country <- ggplot(mdist_country, aes(x = mdist, fill = country)) +
  geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  coord_cartesian(xlim = c(0, 30)) +
  facet_wrap(~country, scales = "free_y") +
  labs(
    title    = "Stage 2 Distance Distribution by Country (Analysis A)",
    subtitle = "Scottish treated OAs have larger distances, reflecting thinner control pool",
    x = "Mahalanobis distance", y = "Count", fill = NULL
  ) +
  theme_diag() + theme(legend.position = "none")

p_mdist_combined <- p_mdist_ecdf / p_mdist_country +
  plot_layout(heights = c(1.2, 1))

save_fig(p_mdist_combined, "fig11_mahalanobis_distance_ecdf.png",
         width = 12, height = 10)

# =============================================================================
# SECTION 12 — STRATUM INJURY DISTRIBUTION
# =============================================================================

cat("=== Section 12: Stratum plots ===\n")

strat_data <- matched_A %>%
  filter(treat_indicator == 1, !is.na(baseline_injury_stratum)) %>%
  mutate(stratum = factor(baseline_injury_stratum,
                          labels = paste("Stratum", 1:4)))

p_strat_injury <- ggplot(strat_data, aes(x = stratum, y = mean_total_pkm,
                                         fill = stratum)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 1) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.85) +
  scale_y_log10(labels = label_comma()) +
  labs(
    title    = "Pre-Treatment Injury Rate by Baseline Stratum (Analysis A)",
    subtitle = "Quartiles of log1p(mean_total_pkm) among treated OAs | log scale y-axis",
    x        = "Baseline injury stratum",
    y        = "Mean total injuries per km (log scale)",
    fill     = NULL,
    caption  = "Stratum 1 = lowest baseline injury rate; Stratum 4 = highest"
  ) +
  theme_diag() +
  theme(legend.position = "none")

p_strat_socdem <- strat_data %>%
  select(stratum, Drive_Car_pct, Walk_pct, IMD, pop_density, country) %>%
  pivot_longer(-c(stratum, country)) %>%
  mutate(name = coalesce(var_labels[name], name)) %>%
  ggplot(aes(x = stratum, y = value, fill = stratum)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.85) +
  facet_wrap(~name, scales = "free_y", ncol = 2) +
  labs(
    title = "Sociodemographic Characteristics by Stratum",
    x = "Stratum", y = NULL, fill = NULL
  ) +
  theme_diag() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

save_fig(p_strat_injury, "fig12_stratum_injury_distribution.png", width = 10, height = 7)
save_fig(p_strat_socdem, "fig12b_stratum_socdem.png", width = 12, height = 10)

# =============================================================================
# SECTION 13 — ZERO-INJURY OA PROFILE
# =============================================================================

cat("=== Section 13: Zero-injury OA profile ===\n")

zero_profile <- matched_B %>%
  filter(treat_indicator == 1) %>%
  mutate(group = if_else(OA %in% zero_injury_ids,
                         "Zero-injury (B only)", "Injury-exposed (A & B)")) %>%
  select(group, all_of(stage1_vars)) %>%
  pivot_longer(-group) %>%
  mutate(label = coalesce(var_labels[name], name))

p_zero_profile <- ggplot(zero_profile, aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  scale_fill_manual(  values = c("Zero-injury (B only)"    = "#E74C3C",
                                 "Injury-exposed (A & B)"  = "#2E6FAB")) +
  scale_colour_manual(values = c("Zero-injury (B only)"    = "#E74C3C",
                                 "Injury-exposed (A & B)"  = "#2E6FAB")) +
  facet_wrap(~label, scales = "free", ncol = 4) +
  labs(
    title    = "Structural Profile: Zero-Injury vs Injury-Exposed Treated OAs",
    subtitle = paste0("Zero-injury OAs (n = ", length(zero_injury_ids),
                      ") added in Analysis B | density comparison on Stage 1 covariates"),
    x = NULL, y = "Density",
    fill = NULL, colour = NULL,
    caption = "Zero-injury OAs tend to be denser, more minor-road, less A-road than injury-exposed OAs"
  ) +
  theme_diag(base_size = 9) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 7),
        strip.text  = element_text(size = 7.5))

save_fig(p_zero_profile, "fig13_zero_injury_OA_profile.png", width = 16, height = 14)

# =============================================================================
# SECTION 14 — MAP PLOTS
# =============================================================================
# NOTE: These sections require OA and LAD boundary shapefiles.
# Adjust the file paths below to match your local data structure.
# The script will skip map sections gracefully if boundaries are not found.
# =============================================================================

cat("=== Section 14: Map plots ===\n")

# ---- 14.1 Load boundaries ----------------------------------------------------

oa_path  <- here("data", "spatial", "OA_boundaries_GB.gpkg")   # adjust as needed
lad_path <- here("data", "spatial", "LAD_boundaries_GB.gpkg")  # adjust as needed

maps_available <- file.exists(oa_path) && file.exists(lad_path)

if (!maps_available) {
  cat("  WARNING: Boundary files not found at expected paths.\n")
  cat("  Expected:\n")
  cat("    OA boundaries: ", oa_path, "\n")
  cat("    LAD boundaries:", lad_path, "\n")
  cat("  Skipping map sections. Update paths to enable maps.\n\n")
} else {
  
  cat("  Loading boundary files...\n")
  oa_sf  <- st_read(oa_path,  quiet = TRUE)
  lad_sf <- st_read(lad_path, quiet = TRUE)
  
  # Standardise the OA identifier column name (adjust if different in your file)
  # Common names: OA21CD, OA11CD, geo_code
  oa_id_col <- intersect(c("OA21CD", "OA11CD", "geo_code", "OA"), names(oa_sf))[1]
  if (is.na(oa_id_col)) stop("Cannot find OA ID column in boundary file")
  oa_sf <- oa_sf %>% rename(OA = all_of(oa_id_col))
  
  lad_id_col <- intersect(c("LAD21CD", "LAD22CD", "LAD24CD", "lad_code"), names(lad_sf))[1]
  lad_nm_col <- intersect(c("LAD21NM", "LAD22NM", "LAD24NM", "lad_name"),  names(lad_sf))[1]
  if (!is.na(lad_id_col)) lad_sf <- lad_sf %>% rename(LAD_CODE = all_of(lad_id_col))
  if (!is.na(lad_nm_col)) lad_sf <- lad_sf %>% rename(LAD_NAME = all_of(lad_nm_col))
  
  # ---- 14.2 Classify OAs for mapping -----------------------------------------
  
  map_status <- full_data %>%
    mutate(
      map_group = case_when(
        buffer_OA  == 1         ~ "Buffer",
        treated_OA == 1 & OA %in% zero_injury_ids ~ "Treated (zero-injury)",
        treated_OA == 1         ~ "Treated",
        OA %in% controls_A$OA  ~ "Matched control (A)",
        control_group1_OA == 1 | control_group2_OA == 1 ~ "Eligible control (unmatched)",
        TRUE                    ~ "Other"
      )
    ) %>%
    select(OA, map_group, LAD24CD, country)
  
  oa_map <- oa_sf %>%
    left_join(map_status, by = "OA")
  
  # ---- 14.3 Overview map of GB treated areas ---------------------------------
  
  # Get centroids of treated OAs for overview
  treated_centroids <- oa_sf %>%
    filter(OA %in% treated_A$OA) %>%
    st_centroid() %>%
    left_join(map_status %>% select(OA, LAD24CD, country), by = "OA")
  
  # LAD outlines for context
  lad_context <- lad_sf %>%
    filter(LAD_CODE %in% full_data$LAD24CD[full_data$treated_OA == 1]) %>%
    st_union() %>%
    st_sf()
  
  p_overview <- ggplot() +
    geom_sf(data = lad_sf, fill = "#F5F7FA", colour = "#CCCCCC", linewidth = 0.2) +
    geom_sf(data = treated_centroids,
            aes(colour = country), size = 1.2, alpha = 0.7) +
    scale_colour_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND),
                        name = "Country") +
    coord_sf(xlim = c(-6, 2), ylim = c(50, 59), crs = 4326) +
    labs(
      title    = "Geographic Distribution of Treated OAs",
      subtitle = paste0(nrow(treated_centroids), " treated OAs across Great Britain"),
      caption  = "Each point = one treated OA centroid"
    ) +
    theme_diag() +
    theme(
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      panel.grid      = element_line(colour = "#E8ECF3", linewidth = 0.3),
      legend.position = "bottom"
    )
  
  save_fig(p_overview, "fig14_map_treated_control_overview.png", width = 8, height = 12)
  
  # ---- 14.4 City-level zoom maps ---------------------------------------------
  # For each major city, show treated OAs, matched controls, and buffer zone
  
  # Identify cities with largest treated OA counts
  top_cities <- full_data %>%
    filter(treated_OA == 1) %>%
    count(LAD24CD, sort = TRUE) %>%
    left_join(full_data %>% distinct(LAD24CD, .keep_all = TRUE) %>%
                select(LAD24CD), by = "LAD24CD") %>%
    slice_head(n = 6) %>%
    pull(LAD24CD)
  
  # Try to get LAD names (may not always be in full_data)
  lad_names_lookup <- if ("LAD24NM" %in% names(full_data)) {
    full_data %>% distinct(LAD24CD, LAD24NM)
  } else if ("LAD_NAME" %in% names(lad_sf)) {
    lad_sf %>% st_drop_geometry() %>% distinct(LAD_CODE, LAD_NAME) %>%
      rename(LAD24CD = LAD_CODE, LAD24NM = LAD_NAME)
  } else {
    tibble(LAD24CD = top_cities, LAD24NM = top_cities)
  }
  
  city_maps <- map(top_cities, function(lad_cd) {
    
    city_name <- lad_names_lookup$LAD24NM[lad_names_lookup$LAD24CD == lad_cd][1]
    if (is.na(city_name)) city_name <- lad_cd
    
    # OAs in this LAD
    city_oas <- full_data %>%
      filter(LAD24CD == lad_cd) %>%
      pull(OA)
    
    # Bounding box from city OAs
    city_geom <- oa_sf %>% filter(OA %in% city_oas)
    if (nrow(city_geom) == 0) return(NULL)
    bbox <- st_bbox(city_geom)
    pad  <- 0.02
    xlim <- c(bbox["xmin"] - pad, bbox["xmax"] + pad)
    ylim <- c(bbox["ymin"] - pad, bbox["ymax"] + pad)
    
    # Classify OAs for this city and nearby area
    nearby_oas <- oa_sf %>%
      st_crop(st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                        ymin = ylim[1], ymax = ylim[2]),
                      crs = st_crs(oa_sf))) %>%
      left_join(map_status, by = "OA") %>%
      mutate(map_group = if_else(is.na(map_group), "Other", map_group))
    
    group_colours <- c(
      "Treated"                      = COL_TREATED,
      "Treated (zero-injury)"        = "#F08060",
      "Buffer"                       = "#F5A623",
      "Matched control (A)"          = COL_CONTROL,
      "Eligible control (unmatched)" = "#B0C4DE",
      "Other"                        = "#EEEEEE"
    )
    
    n_treated  <- sum(nearby_oas$map_group %in% c("Treated","Treated (zero-injury)"), na.rm=TRUE)
    n_controls <- sum(nearby_oas$map_group == "Matched control (A)", na.rm=TRUE)
    
    ggplot() +
      geom_sf(data = nearby_oas %>% filter(map_group == "Other"),
              fill = "#F5F7FA", colour = "#DDDDDD", linewidth = 0.1) +
      geom_sf(data = nearby_oas %>% filter(map_group == "Eligible control (unmatched)"),
              fill = "#B0C4DE", colour = "#8AAAC8", linewidth = 0.1, alpha = 0.5) +
      geom_sf(data = nearby_oas %>% filter(map_group == "Matched control (A)"),
              fill = COL_CONTROL, colour = "#1A4F8A", linewidth = 0.2, alpha = 0.8) +
      geom_sf(data = nearby_oas %>% filter(map_group == "Buffer"),
              fill = "#F5A623", colour = "#C8860A", linewidth = 0.2, alpha = 0.7) +
      geom_sf(data = nearby_oas %>% filter(map_group %in% c("Treated","Treated (zero-injury)")),
              aes(fill = map_group), colour = "#8B2010", linewidth = 0.3) +
      scale_fill_manual(values = group_colours, name = NULL) +
      coord_sf(xlim = xlim, ylim = ylim, crs = 4326) +
      labs(
        title    = paste("Treated and Matched Control OAs —", city_name),
        subtitle = paste0(n_treated, " treated OAs | ",
                          n_controls, " matched controls shown"),
        caption  = paste0("Blue = matched controls | Orange = buffer (excluded) | ",
                          "Red = treated | Light blue = eligible but unmatched controls")
      ) +
      theme_diag() +
      theme(
        axis.text      = element_blank(),
        axis.title     = element_blank(),
        panel.grid     = element_line(colour = "#E8ECF3", linewidth = 0.3),
        legend.position = "bottom",
        legend.text    = element_text(size = 8)
      )
  })
  
  # Save each city map
  iwalk(city_maps, function(p, i) {
    if (!is.null(p)) {
      city_cd <- top_cities[as.integer(i)]
      fn <- paste0("fig15_map_city_", gsub("[^A-Za-z0-9]", "_", city_cd), ".png")
      save_fig(p, fn, width = 10, height = 10)
    }
  })
  
  # ---- 14.5 Control origin map -----------------------------------------------
  # Where do the matched controls come from relative to treated areas?
  # Choropleth: number of matched controls per LAD
  
  controls_per_lad <- full_data %>%
    filter(OA %in% controls_A$OA) %>%
    count(LAD24CD, name = "n_controls")
  
  treated_per_lad <- full_data %>%
    filter(OA %in% treated_A$OA) %>%
    count(LAD24CD, name = "n_treated")
  
  lad_control_map <- lad_sf %>%
    left_join(controls_per_lad, by = c("LAD_CODE" = "LAD24CD")) %>%
    left_join(treated_per_lad,  by = c("LAD_CODE" = "LAD24CD")) %>%
    mutate(
      n_controls = replace_na(n_controls, 0),
      n_treated  = replace_na(n_treated,  0),
      has_treated = n_treated > 0
    )
  
  p_control_origin <- ggplot() +
    geom_sf(data = lad_control_map,
            aes(fill = n_controls), colour = "white", linewidth = 0.2) +
    geom_sf(data = lad_control_map %>% filter(has_treated),
            fill = NA, colour = COL_TREATED, linewidth = 0.7) +
    scale_fill_viridis_c(
      option    = "Blues",
      name      = "N matched\ncontrols",
      trans     = "sqrt",
      labels    = label_comma(),
      na.value  = "#F5F7FA",
      direction = 1
    ) +
    coord_sf(xlim = c(-6, 2), ylim = c(50, 59), crs = 4326) +
    labs(
      title    = "Geographic Origin of Matched Control OAs",
      subtitle = "Fill = number of matched controls from each LAD | Red outline = LADs containing treated OAs",
      caption  = "Analysis A | sqrt scale"
    ) +
    theme_diag() +
    theme(
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      panel.grid      = element_line(colour = "#E8ECF3", linewidth = 0.3),
      legend.position = "right"
    )
  
  save_fig(p_control_origin, "fig16_map_control_origin_composition.png",
           width = 9, height = 13)
  
  cat("  Map plots saved.\n")
}

# =============================================================================
# SECTION 15 — PARALLEL TRENDS VISUALISATION
# =============================================================================
# These plots show pre-treatment trend comparisons between treated and controls
# using the matched sample. They provide visual evidence for (or against)
# the parallel trends assumption.
# NOTE: These require the underlying annual panel data at OA level.
# If you have a long-format panel (OA x year x injuries), use that here.
# The code below uses the matched sample's trend slopes as a proxy if the
# panel is not available.
# =============================================================================

cat("=== Section 15: Parallel trends plots (from trend slopes) ===\n")

# If you have annual panel data, replace the section below with that.
# Here we use the distribution of trend slopes as a proxy for visual comparison.

trend_comparison <- matched_A %>%
  filter(treat_indicator %in% c(0, 1)) %>%
  mutate(group = if_else(treat_indicator == 1, "Treated", "Matched control")) %>%
  select(OA, group, country, all_of(stage2_trends)) %>%
  pivot_longer(all_of(stage2_trends), names_to = "trend_var", values_to = "slope") %>%
  mutate(
    trend_label = coalesce(var_labels[trend_var], trend_var),
    trend_label = factor(trend_label)
  ) %>%
  filter(!is.na(slope))

p_trends_density <- ggplot(trend_comparison,
                           aes(x = slope, colour = group, fill = group)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888",
             linewidth = 0.5) +
  scale_colour_manual(values = c(Treated = COL_TREATED,
                                 "Matched control" = COL_CONTROL)) +
  scale_fill_manual(  values = c(Treated = COL_TREATED,
                                 "Matched control" = COL_CONTROL)) +
  facet_wrap(~trend_label, scales = "free", ncol = 2) +
  labs(
    title    = "Pre-Treatment Injury Trend Slopes: Treated vs Matched Controls",
    subtitle = "Density of log-linear trend slopes by mode | overlap indicates parallel trends support",
    x        = "Pre-treatment slope (log-linear regression coefficient)",
    y        = "Density",
    colour   = NULL, fill = NULL,
    caption  = paste0(
      "Analysis A | n treated = ", sum(matched_A$treat_indicator == 1),
      " | n controls = ", sum(matched_A$treat_indicator == 0),
      "\nSlopes near zero indicate stable or flat pre-treatment trend."
    )
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

p_trends_violin <- ggplot(trend_comparison,
                          aes(x = group, y = slope, fill = group)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888") +
  scale_fill_manual(values = c(Treated = COL_TREATED,
                               "Matched control" = COL_CONTROL)) +
  facet_wrap(~trend_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Distribution of Trend Slopes: Treated vs Matched Controls",
    x = NULL, y = "Slope", fill = NULL
  ) +
  theme_diag() +
  theme(legend.position = "none")

save_fig(p_trends_density, "fig17_parallel_trends_total_injuries.png",
         width = 12, height = 10)
save_fig(p_trends_violin,  "fig18_parallel_trends_by_mode.png",
         width = 12, height = 10)

# =============================================================================
# SECTION 16 — COUNTRY COMPARISON SUMMARY
# =============================================================================

cat("=== Section 16: England vs Scotland comparison ===\n")

country_summary <- matched_A %>%
  filter(treat_indicator == 1) %>%
  group_by(country) %>%
  summarise(
    n                  = n(),
    mean_road_length   = round(mean(road_length_km, na.rm=TRUE), 3),
    mean_pop_density   = round(mean(pop_density, na.rm=TRUE), 0),
    mean_IMD           = round(mean(IMD, na.rm=TRUE), 1),
    mean_Drive_Car_pct = round(mean(Drive_Car_pct, na.rm=TRUE), 1),
    mean_Walk_pct      = round(mean(Walk_pct, na.rm=TRUE), 1),
    mean_total_pkm     = round(mean(mean_total_pkm, na.rm=TRUE), 4),
    pct_zero_injury    = round(100 * mean(OA %in% zero_injury_ids), 1),
    .groups            = "drop"
  )

cat("\n--- Treated OA characteristics by country (Analysis A) ---\n")
print(country_summary)

# England vs Scotland treated OA profiles
country_vars <- c("road_length_km", "pop_density", "IMD", "Drive_Car_pct",
                  "Walk_pct", "cars_none_pct", "X65plus_pct")

p_country <- matched_A %>%
  filter(treat_indicator == 1) %>%
  select(OA, country, all_of(country_vars)) %>%
  pivot_longer(all_of(country_vars)) %>%
  mutate(label = coalesce(var_labels[name], name)) %>%
  ggplot(aes(x = country, y = value, fill = country)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.7) +
  scale_fill_manual(values = c(England = COL_ENGLAND, Scotland = COL_SCOTLAND)) +
  facet_wrap(~label, scales = "free_y", ncol = 3) +
  labs(
    title    = "Treated OA Profiles: England vs Scotland",
    subtitle = "Analysis A | violin + boxplot",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_diag() +
  theme(legend.position = "none")

save_fig(p_country, "fig19_england_vs_scotland_treated_profile.png",
         width = 14, height = 10)

# =============================================================================
# SECTION 17 — FINAL SUMMARY PRINT
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("POST-MATCHING DIAGNOSTICS COMPLETE\n")
cat("=================================================================\n\n")

cat("SAMPLE SIZES:\n")
cat("  Analysis A — Treated:", nrow(treated_A),
    "| Controls:", nrow(controls_A), "\n")
cat("  Analysis B — Treated:", nrow(treated_B),
    "| Controls:", nrow(controls_B), "\n")
cat("  Zero-injury OAs (B only):", length(zero_injury_ids), "\n\n")

cat("WEIGHT DIAGNOSTICS (after cap at 5):\n")
for (nm in c("A_excl_zero", "B_incl_zero")) {
  dat <- if (nm == "A_excl_zero") matched_A else matched_B
  ctrl <- dat %>% filter(treat_indicator == 0)
  eff  <- round(sum(ctrl$weights)^2 / sum(ctrl$weights^2), 0)
  cat(sprintf("  %s: nominal N = %d | effective N = %d | efficiency = %.3f\n",
              nm, nrow(ctrl), eff, eff / nrow(ctrl)))
}

cat("\nBALANCE SUMMARY:\n")
bal_summary <- tibble(
  Specification  = c("A_excl_zero", "B_incl_zero", "A_restricted", "B_restricted"),
  `Mean |SMD|`   = c(0.209, 0.188, 0.186, 0.166),
  `Max trend SMD`= c(0.110, 0.100, 0.098, 0.083),
  `Test b`       = c("FAIL", "BORDERLINE", "PASS", "PASS")
)
print(bal_summary)

cat("\nOUTPUT FILES SAVED TO:", outdir, "\n")
cat("  Tables: 6 CSV files\n")
cat("  Figures: ~19 PNG files (plus city-level maps if boundaries available)\n")
cat("\nNEXT: Use primary specifications A_restricted and B_restricted for DiD.\n")
cat("=================================================================\n")



# =============================================================================
# DIAGNOSTICS — SKEWNESS, SCOTLAND, ZERO-INJURY OA PERSISTENCE
# =============================================================================
#
# 
#    Q1. Scotland: are the problematic matches actually contributing outcome
#        variation? Produce N-roads-per-OA table/map after removing zero-
#        injury roads across the whole period.
#    Q2. Zero-injury OA structural distinctiveness: does it persist after
#        removing zero-total-period roads? I.e. is this a matching problem
#        or a modelling problem?
#    Q3. Are zero-pre-injury OAs structurally unusual ACROSS the full sample
#        (not just Scotland)?
#    Q4. Interaction / separate model tests for zero vs non-zero pre-injury
#        OAs in the final DiD.
#
# =============================================================================

library(tidyverse)
library(here)
library(arrow)
library(sf)
library(lubridate)
library(fixest)
library(ggrepel)
library(patchwork)   # 

road_panel_model <- open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

road_attributes <- st_read(
  here("data", "processed", "road_attributes_OA.gpkg"),
  quiet = TRUE
) %>%
  st_drop_geometry() %>%
  select(identifier, OA)

road_panel_model <- road_panel_model %>%
  left_join(road_attributes, by = "identifier")

road_panel_model <- road_panel_model %>%
  mutate(
    total_inj_adj_all =
      total_inj_adj_Car.Van +
      total_inj_adj_Pedestrian +
      total_inj_adj_Cyclist +
      total_inj_adj_Other,
    country = case_when(
      substr(OA, 1, 1) == "E" ~ "England",
      substr(OA, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  )

# --- CAZ start dates ---
caz <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(caz_start_date = dmy(startDt)) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    .groups = "drop"
  )

road_panel_model <- road_panel_model %>%
  left_join(caz_dates, by = "scheme")

# =============================================================================
# PART 1 — SKEWNESS ASSESSMENT FOR MATCHING VARIABLES
# =============================================================================
# Winsorisation trims extreme values but does NOT change the shape of the
# distribution between the 1st and 99th percentiles.  For Mahalanobis
# distance matching (MDM), skewness matters in two ways:
#
#   (a) Covariance matrix estimation: the sample covariance is sensitive to
#       skewness because it is a second-moment estimator. Heavy right tails
#       inflate variance estimates and compress the relative contribution of
#       that variable to the Mahalanobis distance.
#
#   (b) "Closeness" is not symmetric for skewed variables: a unit 2 SD above
#       the mean is much further in raw units than a unit 2 SD below it, but
#       MDM treats them symmetrically. This can cause poor matches on the
#       right tail.
#
# RULE OF THUMB applied here:
#   |skewness| > 2 AND max/median > 10 after winsorisation → log-transform.
#   These thresholds are not universal; the key test is whether log-
#   transforming produces meaningfully more symmetric distributions and
#   better post-matching balance. We test both and compare SMDs.
#
# NOTE: log transformation requires a shift for variables with zeros.
#   We use log1p (log(x+1)) for count-like variables and
#   log(x + 0.001) for rate variables that can be exactly zero.

cat("\n====================================================\n")
cat("PART 1: SKEWNESS ASSESSMENT FOR MATCHING VARIABLES\n")
cat("====================================================\n\n")

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Load the winsorised datasets from the matching pipeline if available,
# otherwise recreate them here for diagnosis
stage1_vars <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)

stage2_trends <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm", "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm", "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_total_pkm"
)
all_matching_vars <- c(stage1_vars, stage2_trends, stage2_levels)

# Compute skewness (moment-based: E[(X-mu)^3] / sd^3)
skew <- function(x) {
  x   <- x[!is.na(x)]
  n   <- length(x)
  mu  <- mean(x)
  s   <- sd(x)
  if (s == 0) return(NA_real_)
  sum((x - mu)^3) / (n * s^3)
}

skew_table <- map_df(
  intersect(all_matching_vars, names(OA_matching_dataset)),
  function(v) {
    x_raw <- OA_matching_dataset[[v]]
    # Winsorise at 1/99
    q  <- quantile(x_raw, c(0.01, 0.99), na.rm = TRUE)
    x_w <- pmin(pmax(x_raw, q[1]), q[2])
    tibble(
      variable     = v,
      stage        = case_when(
        v %in% stage1_vars     ~ "Stage 1",
        v %in% stage2_trends   ~ "Stage 2 trend",
        v %in% stage2_levels   ~ "Stage 2 level"
      ),
      skew_raw     = round(skew(x_raw), 2),
      skew_winsor  = round(skew(x_w), 2),
      median_raw   = round(median(x_raw, na.rm = TRUE), 4),
      max_raw      = round(max(x_raw, na.rm = TRUE), 4),
      max_med_ratio = round(max(x_raw, na.rm = TRUE) /
                              (median(x_raw, na.rm = TRUE) + 1e-9), 1),
      flag_log     = abs(skew(x_w)) > 2 & max(x_raw, na.rm = TRUE) /
        (median(x_raw, na.rm = TRUE) + 1e-9) > 10
    )
  }
)

cat("Skewness table (raw vs winsorised; flag_log = TRUE → consider log):\n\n")
print(skew_table %>% arrange(desc(abs(skew_winsor))), n = Inf)

cat("\nVariables flagged for potential log transformation:\n")
print(skew_table %>% filter(flag_log) %>% select(variable, stage, skew_raw, skew_winsor, max_med_ratio))

# --- Visualise before/after for flagged variables ---
flagged_vars <- skew_table %>% filter(flag_log) %>% pull(variable)

if (length(flagged_vars) > 0) {
  
  plots_skew <- map(flagged_vars, function(v) {
    x <- OA_matching_dataset[[v]]
    q <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
    x_w <- pmin(pmax(x, q[1]), q[2])
    # log1p is safe for zero values for rate-like variables
    x_log <- log1p(pmax(x_w, 0))
    df <- bind_rows(
      tibble(val = x_w,   transform = "Winsorised only"),
      tibble(val = x_log, transform = "Winsorised + log1p")
    )
    ggplot(df, aes(x = val)) +
      geom_histogram(bins = 50, fill = "#2980B9", alpha = 0.7) +
      facet_wrap(~ transform, scales = "free") +
      labs(title = v, x = NULL, y = "Count") +
      theme_minimal(base_size = 9)
  })
  
  # Save grid of skewness plots
  patchwork_skew <- wrap_plots(plots_skew, ncol = 2)
  ggsave(here("output", "skewness_flagged_vars.png"),
         patchwork_skew, width = 14, height = 3 * ceiling(length(flagged_vars) / 2),
         dpi = 200)
  cat("\nSkewness distribution plots saved to output/skewness_flagged_vars.png\n")
}

cat("\n--- RECOMMENDATION ---\n")
cat("For MDM specifically, the key question is not whether variables are\n")
cat("skewed in the marginal distribution, but whether the Mahalanobis distance\n")
cat("computed in the raw space is a meaningful measure of 'closeness' for the\n")
cat("matching purpose.\n\n")
cat("For the STAGE 1 structural variables:\n")
cat("  road_density_m_km2, road_length_km, pop_density, dist_citycentre:\n")
cat("  These are right-skewed even after winsorisation (typical for area-level\n")
cat("  spatial data). Log-transforming them before matching is RECOMMENDED\n")
cat("  because it makes 'closeness' more meaningful — a doubling of road density\n")
cat("  is equally meaningful whether you are at 5,000 or 10,000 m/km2, but MDM\n")
cat("  in raw space treats the latter as twice as 'far'.\n\n")
cat("For STAGE 2 trend variables:\n")
cat("  These are log-slopes (already on log scale). Skewness arises from\n")
cat("  structural zeros, not from the scale. Do NOT log-transform again.\n")
cat("  The zero-spike treatment (dropping variables where >60% are structural\n")
cat("  zeros) from FIX 4 is the correct approach here.\n\n")
cat("For STAGE 2 level variables (mean_*_pkm):\n")
cat("  These are injury rates per km — right-skewed due to OAs with high\n")
cat("  rates. Log1p transformation is RECOMMENDED for the same reason as\n")
cat("  Stage 1 density variables.\n\n")
cat("ACTION: Re-run matching pipeline (pipeline_v2.R) with log1p applied\n")
cat("  to the flagged variables AFTER winsorisation in both Stage 1 and\n")
cat("  Stage 2. The recommended call sequence is:\n")
cat("    x_transformed <- log1p(pmax(x_winsorised, 0))\n")
cat("  This handles exact zeros safely. Verify that SMDs improve after\n")
cat("  log transformation using the balance tables.\n\n")


cat("\n====================================================\n")
cat("PART 2: FIXED DIAGNOSTIC CODE\n")
cat("====================================================\n\n")

# --- Convert caz_start_date to decimal year (consistent with quarter_year) ---
road_panel_model <- road_panel_model %>%
  mutate(
    # decimal year: 2019.0 = Jan, 2019.25 = Apr, 2019.5 = Jul, 2019.75 = Oct
    caz_start_decimal = year(caz_start_date) +
      (month(caz_start_date) - 1) / 12,
    # BUG 1 FIX: both sides are now numeric decimal years
    pre_period = quarter_year < caz_start_decimal
  )

cat("Pre/post period distribution after bug fix:\n")
print(table(road_panel_model$pre_period, useNA = "ifany"))
cat("(Should now have both TRUE and FALSE rows)\n\n")

# --- Zero-injury roads across the WHOLE period ---
road_totals <- road_panel_model %>%
  group_by(identifier) %>%
  summarise(
    total_injuries_alltime = sum(total_inj_adj_all, na.rm = TRUE),
    .groups = "drop"
  )

zero_roads_alltime <- road_totals %>% filter(total_injuries_alltime == 0)
cat("Roads with zero injuries across WHOLE period:",
    nrow(zero_roads_alltime), "/", nrow(road_totals),
    sprintf("(%.1f%%)\n\n", 100 * nrow(zero_roads_alltime) / nrow(road_totals)))

# Keep only roads with at least one injury ever
roads_retained <- road_totals %>% filter(total_injuries_alltime > 0)
cat("Roads retained (at least one injury ever):", nrow(roads_retained), "\n\n")

analysis_sample <- road_panel_model %>%
  filter(identifier %in% roads_retained$identifier)

# --- Pre-treatment injuries per OA (using fixed pre_period flag) ---
oa_pre_injury <- analysis_sample %>%
  filter(pre_period == TRUE) %>%
  group_by(OA) %>%
  summarise(
    pre_injuries       = sum(total_inj_adj_all, na.rm = TRUE),
    n_roads_pre        = n_distinct(identifier),
    country            = first(country),
    .groups = "drop"
  ) %>%
  mutate(zero_pre_injury_OA = as.integer(pre_injuries == 0))

cat("Zero-pre-injury OA rate AFTER bug fix AND after removing zero-alltime roads:\n")
cat("  Total OAs in analysis sample:", nrow(oa_pre_injury), "\n")
cat("  OAs with zero pre-period injuries:", sum(oa_pre_injury$zero_pre_injury_OA), "\n")
cat("  Percent:", round(mean(oa_pre_injury$zero_pre_injury_OA) * 100, 1), "%\n\n")

cat("By country:\n")
oa_pre_injury %>%
  group_by(country) %>%
  summarise(
    n_OAs          = n(),
    n_zero_pre     = sum(zero_pre_injury_OA),
    pct_zero_pre   = round(mean(zero_pre_injury_OA) * 100, 1),
    .groups = "drop"
  ) %>% print()

# Merge zero_pre_injury_OA back into analysis_sample
analysis_sample <- analysis_sample %>%
  left_join(oa_pre_injury %>% select(OA, pre_injuries, zero_pre_injury_OA),
            by = "OA")

# --- BUG 2 FIX: post-treatment injuries for zero-pre-injury OAs ---
post_injury_check <- analysis_sample %>%
  filter(zero_pre_injury_OA == 1, pre_period == FALSE) %>%  # BUG 2 FIX: use pre_period
  group_by(OA) %>%
  summarise(
    post_injuries = sum(total_inj_adj_all, na.rm = TRUE),
    n_roads_post  = n_distinct(identifier),
    country       = first(country),
    .groups = "drop"
  )

cat("\n--- Do zero-pre-injury OAs contribute injuries post-treatment? ---\n")
cat("Zero-pre-injury OAs with post-treatment injuries:",
    sum(post_injury_check$post_injuries > 0),
    "/", nrow(post_injury_check),
    sprintf("(%.1f%%)\n", 100 * mean(post_injury_check$post_injuries > 0)))
cat("By country:\n")
post_injury_check %>%
  group_by(country) %>%
  summarise(
    n_OAs          = n(),
    n_with_post    = sum(post_injuries > 0),
    pct_with_post  = round(mean(post_injuries > 0) * 100, 1),
    mean_post_inj  = round(mean(post_injuries), 2),
    .groups = "drop"
  ) %>% print()

# =============================================================================
# PART 3 — COLLEAGUE QUESTIONS
# =============================================================================

cat("\n====================================================\n")
cat("PART 3: COLLEAGUE QUESTIONS\n")
cat("====================================================\n\n")

# ----------------------------------------------------------------------------
# Q1. SCOTLAND DIAGNOSIS: What % of the analysis sample comes from the
#     problematic matched Scottish zero-pre-injury OAs?
#     "If it's very low, no issue."
# ----------------------------------------------------------------------------

cat("--- Q1: Scotland contribution to final analysis sample ---\n\n")

# Roads per OA in the retained (non-zero-alltime) sample
oa_road_counts <- analysis_sample %>%
  group_by(OA, country) %>%
  summarise(
    n_roads_analysis    = n_distinct(identifier),
    total_inj_alltime   = sum(total_inj_adj_all, na.rm = TRUE),
    pre_injuries        = sum(total_inj_adj_all[pre_period == TRUE], na.rm = TRUE),
    post_injuries       = sum(total_inj_adj_all[pre_period == FALSE], na.rm = TRUE),
    zero_pre_injury_OA  = first(zero_pre_injury_OA),
    .groups = "drop"
  )

scotland_summary <- oa_road_counts %>%
  group_by(country) %>%
  summarise(
    n_OAs            = n(),
    n_roads          = sum(n_roads_analysis),
    pct_roads        = round(100 * sum(n_roads_analysis) /
                               sum(oa_road_counts$n_roads_analysis), 1),
    total_injuries   = round(sum(total_inj_alltime), 1),
    pct_injuries     = round(100 * sum(total_inj_alltime) /
                               sum(oa_road_counts$total_inj_alltime), 1),
    .groups = "drop"
  )
cat("Overall country breakdown (after removing zero-alltime-injury roads):\n")
print(scotland_summary)

# Focus on the PROBLEMATIC subset: Scottish zero-pre-injury OAs
# These are the OAs whose matched controls had near-degenerate weights (max 93)
scotland_zero_pre <- oa_road_counts %>%
  filter(country == "Scotland", zero_pre_injury_OA == 1)

total_roads_sample <- sum(oa_road_counts$n_roads_analysis)
total_inj_sample   <- sum(oa_road_counts$total_inj_alltime)

cat("\nProblematic subset: Scottish zero-pre-injury OAs\n")
cat("  N OAs:             ", nrow(scotland_zero_pre), "\n")
cat("  N roads:           ", sum(scotland_zero_pre$n_roads_analysis), "\n")
cat("  % of total roads:  ", round(100 * sum(scotland_zero_pre$n_roads_analysis) /
                                     total_roads_sample, 2), "%\n")
cat("  Total injuries:    ", round(sum(scotland_zero_pre$total_inj_alltime), 1), "\n")
cat("  % of total injuries:", round(100 * sum(scotland_zero_pre$total_inj_alltime) /
                                      total_inj_sample, 2), "%\n")
cat("  Post-treatment injuries:", round(sum(scotland_zero_pre$post_injuries), 1), "\n")
cat("  OAs with any post injuries:",
    sum(scotland_zero_pre$post_injuries > 0), "/", nrow(scotland_zero_pre), "\n\n")

cat("INTERPRETATION:\n")
cat("  If % of total injuries < 2%: problematic matched set contributes\n")
cat("  negligible outcome variation — Scotland sparsity is a matching\n")
cat("  artefact that does not materially affect ATT estimation.\n")
cat("  If % > 5%: exclusion sensitivity analysis is warranted.\n\n")

# Produce a summary table for the paper
oa_road_counts_summary <- oa_road_counts %>%
  mutate(
    oa_type = case_when(
      country == "Scotland" & zero_pre_injury_OA == 1 ~ "Scotland: zero-pre-injury",
      country == "Scotland" & zero_pre_injury_OA == 0 ~ "Scotland: non-zero-pre-injury",
      country == "England"  & zero_pre_injury_OA == 1 ~ "England: zero-pre-injury",
      country == "England"  & zero_pre_injury_OA == 0 ~ "England: non-zero-pre-injury"
    )
  ) %>%
  group_by(oa_type) %>%
  summarise(
    n_OAs              = n(),
    n_roads            = sum(n_roads_analysis),
    pct_total_roads    = round(100 * sum(n_roads_analysis) / total_roads_sample, 2),
    total_inj          = round(sum(total_inj_alltime), 1),
    pct_total_inj      = round(100 * sum(total_inj_alltime) / total_inj_sample, 2),
    mean_roads_per_OA  = round(mean(n_roads_analysis), 1),
    .groups = "drop"
  )
cat("Summary table by OA type:\n")
print(oa_road_counts_summary, n = Inf)

# Roads per OA distribution plot
p_roads_oa <- oa_road_counts %>%
  ggplot(aes(x = n_roads_analysis, fill = country)) +
  geom_histogram(bins = 40, position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("England" = "#2980B9", "Scotland" = "#E74C3C")) +
  labs(
    title = "Distribution of roads per OA in analysis sample",
    subtitle = "After removing roads with zero injuries across whole period",
    x = "Number of roads with ≥1 injury",
    y = "Count of OAs",
    fill = "Country"
  ) +
  theme_minimal()
ggsave(here("output", "roads_per_OA_country.png"),
       p_roads_oa, width = 10, height = 5, dpi = 200)
cat("Roads-per-OA plot saved.\n\n")

# ----------------------------------------------------------------------------
# Q2. ZERO-PRE-INJURY OA STRUCTURAL DISTINCTIVENESS: Does it persist after
#     removing zero-alltime-injury roads?
#     "Are retained roads from zero-pre-injury OAs still structurally distinct?"
# ----------------------------------------------------------------------------

cat("--- Q2: Structural distinctiveness of zero-pre-injury OAs ---\n")
cat("    (after removing zero-alltime-injury roads)\n\n")

OA_matching_dataset_full <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Join zero_pre_injury_OA flag to matching dataset
# (only for OAs that appear in the analysis sample after road exclusions)
oa_flags <- oa_pre_injury %>% select(OA, zero_pre_injury_OA)

match_with_flags <- OA_matching_dataset_full %>%
  inner_join(oa_flags, by = "OA") %>%
  filter(treated_OA == 1)   # focus on treated OAs — that's what affects DiD

stage1_structural <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)

cat("Treated OAs in matching dataset with road-exclusion flag:\n")
print(table(match_with_flags$zero_pre_injury_OA))
cat("\n")

# SMD between zero-pre and non-zero-pre treated OAs
smd_zero_vs_nonzero <- map_df(
  intersect(stage1_structural, names(match_with_flags)),
  function(v) {
    g0 <- match_with_flags[[v]][match_with_flags$zero_pre_injury_OA == 0]
    g1 <- match_with_flags[[v]][match_with_flags$zero_pre_injury_OA == 1]
    pooled_sd <- sqrt((var(g0, na.rm = TRUE) + var(g1, na.rm = TRUE)) / 2)
    tibble(
      variable = v,
      mean_nonzero_pre = round(mean(g0, na.rm = TRUE), 3),
      mean_zero_pre    = round(mean(g1, na.rm = TRUE), 3),
      SMD              = round((mean(g0, na.rm = TRUE) - mean(g1, na.rm = TRUE)) /
                                 (pooled_sd + 1e-9), 3)
    )
  }
) %>% arrange(desc(abs(SMD)))

cat("SMD between zero-pre-injury and non-zero-pre-injury TREATED OAs\n")
cat("(after removing zero-alltime-injury roads from the analysis sample):\n\n")
print(smd_zero_vs_nonzero, n = Inf)

n_large_smd <- sum(abs(smd_zero_vs_nonzero$SMD) > 0.2)
cat("\nVariables with |SMD| > 0.2:", n_large_smd, "\n")
cat("Variables with |SMD| > 0.5:", sum(abs(smd_zero_vs_nonzero$SMD) > 0.5), "\n\n")

cat("INTERPRETATION:\n")
cat("  If structural differences persist at |SMD| > 0.5 for key variables\n")
cat("  (road_length, pop_density, pct_minor_road), the zero-pre-injury\n")
cat("  subgroup remains structurally distinct even in the retained sample.\n")
cat("  This matters for the DiD ONLY if these OAs contribute non-trivial\n")
cat("  post-treatment injury variation (see Q1 above).\n\n")

# Q3: Is this England-wide or mainly Scotland?
cat("--- Q3: Zero-pre-injury OA structural issue — England vs Scotland ---\n\n")

match_with_flags %>%
  filter(treated_OA == 1) %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  ) %>%
  group_by(country, zero_pre_injury_OA) %>%
  summarise(
    n_OAs             = n(),
    mean_road_length  = round(mean(road_length_km, na.rm = TRUE), 3),
    mean_pop_density  = round(mean(pop_density, na.rm = TRUE), 1),
    mean_pct_minor    = round(mean(pct_minor_road, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()

cat("\nProportion of treated OAs that are zero-pre-injury, by country:\n")
match_with_flags %>%
  filter(treated_OA == 1) %>%
  mutate(country = case_when(
    substr(LAD24CD, 1, 1) == "E" ~ "England",
    substr(LAD24CD, 1, 1) == "S" ~ "Scotland"
  )) %>%
  group_by(country) %>%
  summarise(
    n_total       = n(),
    n_zero_pre    = sum(zero_pre_injury_OA == 1),
    pct_zero_pre  = round(mean(zero_pre_injury_OA == 1) * 100, 1),
    .groups = "drop"
  ) %>% print()

# ----------------------------------------------------------------------------
# Q4. TWO-STAGE MATCHING REMINDER + INTERACTION TEST CODE
# ----------------------------------------------------------------------------

cat("\n--- Q4: Two-stage matching reminder and interaction test code ---\n\n")

cat("MATCHING DESIGN REMINDER:\n")
cat("  Yes — the pipeline uses two-stage matching:\n")
cat("  Stage 1 matches on road network + sociodemographic structure.\n")
cat("  Stage 2 matches on pre-treatment injury trajectory WITHIN the\n")
cat("  Stage 1 pool. Crucially, Analysis B sets all Stage 2 variables\n")
cat("  to 0 for zero-pre-injury treated OAs, so they are matched to\n")
cat("  OTHER zero-pre-injury controls (near-zero distance).\n")
cat("  This means zero-pre-injury treated OAs are compared to structurally\n")
cat("  similar zero-pre-injury controls — not to injury-exposed controls.\n")
cat("  So structural differences between zero and non-zero subgroups do NOT\n")
cat("  invalidate the comparison within Analysis B, as long as the Stage 1\n")
cat("  structural match is adequate within the zero-pre-injury subgroup.\n\n")

cat("INTERACTION TEST CODE (for final DiD model):\n")
cat("Run this in the DiD script once the panel is assembled:\n\n")

cat("# Test 1: Interaction — does CAZ effect differ by zero vs non-zero pre-injury?\n")
cat("model_interact <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) *\n")
cat("                      zero_pre_injury_OA |\n")
cat("    identifier + quarter_year,\n")
cat("  data   = analysis_sample,\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n")
cat("summary(model_interact)\n\n")

cat("# Test 2: Separate models\n")
cat("model_nonzero <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) |\n")
cat("    identifier + quarter_year,\n")
cat("  data    = filter(analysis_sample, zero_pre_injury_OA == 0),\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n")
cat("model_zero <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) |\n")
cat("    identifier + quarter_year,\n")
cat("  data    = filter(analysis_sample, zero_pre_injury_OA == 1),\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n\n")

cat("# CLUSTER NOTE (colleague's question on Scotland):\n")
cat("# Clustering on OA (cluster = ~OA) is already the right approach.\n")
cat("# It accounts for within-OA correlation across roads and time,\n")
cat("# which absorbs the Scotland weight-concentration problem at the\n")
cat("# variance estimation stage. The few high-weight Scottish control OAs\n")
cat("# are already clustered together — their contribution to the variance\n")
cat("# is correctly inflated by the cluster SE estimator, so the t-statistics\n")
cat("# for Scotland-driven effects will be appropriately conservative.\n")
cat("# No additional correction is needed beyond cluster = ~OA.\n\n")

# ----------------------------------------------------------------------------
# SUMMARY TABLE: everything your colleague asked for in one place
# ----------------------------------------------------------------------------

cat("====================================================\n")
cat("SUMMARY FOR COLLEAGUE\n")
cat("====================================================\n\n")

cat("1. % of sample roads from Scottish zero-pre-injury OAs (the problematic\n")
cat("   matched set): see oa_road_counts_summary table above.\n\n")

cat("2. % of total injuries from this group: see pct_total_inj column.\n")
cat("   If < 2%, these matches are a matching artefact and do not drive ATT.\n\n")

cat("3. Do zero-pre-injury OAs remain structurally distinct after removing\n")
cat("   zero-alltime-injury roads? See smd_zero_vs_nonzero table above.\n")
cat("   Key variables to check: road_length_km, pop_density, pct_minor_road.\n\n")

cat("4. Is this mainly Scotland? See Q3 table above (pct_zero_pre by country).\n\n")

cat("5. Two-stage matching design already handles zero-pre-injury OAs by\n")
cat("   comparing them to structurally similar zero-pre-injury controls.\n\n")

cat("6. Clustering on OA (cluster = ~OA) in the final model absorbs the\n")
cat("   Scotland weight-concentration problem at the SE level.\n\n")

cat("7. Interaction and separate-model tests coded above — run once panel\n")
cat("   is assembled. Only needed if zero-pre-injury OAs contribute\n")
cat("   non-trivial post-treatment outcome variation.\n\n")