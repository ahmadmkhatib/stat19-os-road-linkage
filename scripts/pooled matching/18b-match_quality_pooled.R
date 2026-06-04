# =============================================================================
# MATCHING QUALITY BY SCHEME — POOLED MATCHING
# =============================================================================
#
# Per-scheme balance diagnostics for the pooled matching (total injuries only).
# Simplified: Stage 2 has only 2 variables (trend_total_pkm, log1p_mean_total_pkm).
#
# INPUTS:
#   OA_matched_full_pooled.rds
#   OA_matching_census.rds
#
# OUTPUTS:
#   output/diagnostics/pooled/scheme_comparison/
#
# =============================================================================

library(tidyverse)
library(here)
library(patchwork)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "pooled", "scheme_comparison")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

matched_full <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))
full_data    <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Stage 2 variables used in pooled matching
stage2_trends <- "trend_total_pkm"
stage2_levels <- "mean_total_pkm"

# Stage 1 variables (same as original)
stage1_vars <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density", "area_km2",
  "business_retail_per_km2", "IMD",
  "cars_one_pct", "cars_twoPlus_pct",
  "Drive_Car_pct", "Passenger_Car_pct", "Walk_pct", "Bicycle_pct",
  "bus_Coach_pct", "Train_pct", "Underground_train_tram_pct",
  "Taxi_pct", "workAthome_pct", "Other_pct",
  "White_pct", "Mixed_pct", "Asian_pct", "Black_pct"
)

all_vars <- c(stage1_vars, stage2_trends, stage2_levels)

# Prepare unmatched pool
add_log_vars <- function(data) {
  for (v in c(stage2_levels)) {
    if (v %in% names(data)) {
      data[[paste0("log1p_", v)]] <- log1p(pmax(data[[v]], 0))
    }
  }
  data
}

unmatched_pool <- full_data %>%
  filter(
    (treated_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  ) %>%
  mutate(
    treat_indicator = as.integer(treated_OA == 1),
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(country == "England") %>%
  add_log_vars()

matched_full <- matched_full %>% add_log_vars()

# Scheme assignment for matched controls
schemes <- sort(unique(matched_full$scheme[matched_full$treat_indicator == 1]))
cat("Schemes:", paste(schemes, collapse = ", "), "\n\n")

# =============================================================================
# COMPUTE SMD PER SCHEME
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

smd_all <- map_df(schemes, function(s) {
  cat("  Computing SMD:", s, "\n")

  treated_oas <- matched_full %>%
    filter(scheme == s, treat_indicator == 1) %>%
    pull(OA)

  # Before: this scheme's treated vs full unmatched control pool
  before_df <- unmatched_pool %>%
    filter(treated_OA == 1 & OA %in% treated_oas | control_group2_OA == 1) %>%
    mutate(treat_indicator = as.integer(OA %in% treated_oas))

  after_df <- matched_full %>% filter(scheme == s)

  check_vars <- c(stage2_trends, paste0("log1p_", stage2_levels), stage1_vars)
  check_vars <- intersect(check_vars, intersect(names(before_df), names(after_df)))

  map_df(check_vars, function(v) {
    tibble(
      scheme     = s,
      variable   = v,
      smd_before = round(compute_smd(before_df, v), 4),
      smd_after  = round(compute_smd(after_df,  v), 4)
    )
  })
})

# =============================================================================
# SUMMARY TABLE
# =============================================================================

smd_summary <- smd_all %>%
  group_by(scheme) %>%
  summarise(
    mean_abs_smd_before = round(mean(abs(smd_before), na.rm = TRUE), 4),
    mean_abs_smd_after  = round(mean(abs(smd_after),  na.rm = TRUE), 4),
    max_abs_smd_after   = round(max(abs(smd_after),   na.rm = TRUE), 4),
    trend_smd_after     = round(
      abs(smd_after[variable == "trend_total_pkm"]), 4),
    level_smd_after     = round(
      abs(smd_after[variable == "log1p_mean_total_pkm"]), 4),
    n_treated = sum(matched_full$scheme == first(scheme) &
                      matched_full$treat_indicator == 1),
    n_control = sum(matched_full$scheme == first(scheme) &
                      matched_full$treat_indicator == 0),
    .groups = "drop"
  ) %>%
  mutate(
    trend_pass = trend_smd_after < 0.10,
    improvement = mean_abs_smd_before - mean_abs_smd_after
  )

cat("\n=== BALANCE SUMMARY BY SCHEME (POOLED MATCHING) ===\n")
print(smd_summary, n = Inf)

write_csv(smd_summary, file.path(outdir, "balance_summary_pooled.csv"))

# =============================================================================
# LOVE PLOT — ALL SCHEMES FACETED
# =============================================================================

# Focus on Stage 2 + key Stage 1 variables
key_vars <- c("trend_total_pkm", "log1p_mean_total_pkm",
              "road_density_m_km2", "road_length_km",
              "pct_A_road", "pct_minor_road",
              "dist_citycentre", "pop_density", "IMD")

love_data <- smd_all %>%
  filter(variable %in% key_vars) %>%
  pivot_longer(c(smd_before, smd_after),
               names_to = "timing", values_to = "smd") %>%
  mutate(
    timing = if_else(timing == "smd_before", "Before", "After"),
    timing = factor(timing, levels = c("Before", "After"))
  )

p_love <- ggplot(love_data,
                 aes(x = abs(smd), y = variable, colour = timing, shape = timing)) +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "#999999") +
  geom_vline(xintercept = 0, colour = "#DDDDDD") +
  geom_line(aes(group = variable), colour = "#DDDDDD", linewidth = 0.3) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Before" = "#E74C3C", "After" = "#2ECC71")) +
  scale_shape_manual(values = c("Before" = 16, "After" = 17)) +
  facet_wrap(~scheme, ncol = 4) +
  labs(
    title = "Balance diagnostics: pooled matching (total injuries)",
    subtitle = "Before vs after matching — key covariates",
    x = "|SMD|", y = NULL, colour = NULL, shape = NULL,
    caption = "Dashed line = 0.10 threshold"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(outdir, "fig_love_plots_pooled.png"),
       p_love, width = 16, height = 10, dpi = 300, bg = "white")

# =============================================================================
# SMD HEATMAP — STAGE 2 VARIABLES BY SCHEME
# =============================================================================

s2_smd <- smd_all %>%
  filter(variable %in% c("trend_total_pkm", "log1p_mean_total_pkm"))

p_heatmap <- ggplot(s2_smd, aes(x = scheme, y = variable, fill = abs(smd_after))) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", smd_after)), size = 3.5) +
  scale_fill_gradient2(low = "#2ECC71", mid = "#F1C40F", high = "#E74C3C",
                       midpoint = 0.10, name = "|SMD|") +
  labs(
    title = "Stage 2 balance after matching — pooled total injuries",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(outdir, "fig_smd_heatmap_pooled.png"),
       p_heatmap, width = 12, height = 4, dpi = 300, bg = "white")

cat("\n=== ALL DIAGNOSTICS SAVED ===\n")
cat("Output dir:", outdir, "\n")
