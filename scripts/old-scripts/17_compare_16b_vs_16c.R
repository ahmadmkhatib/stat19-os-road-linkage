# =============================================================================
# COMPARISON: 16b vs 16c MATCHING APPROACHES
#
# 16b: England only — other-city controls, Scotland excluded entirely
# 16c: England (other-city controls) + Scotland (same-city controls),
#      matched separately and combined
#
# PURPOSE:
#   Decide which approach is more appropriate for the primary analysis.
#   Three conditions are evaluated:
#     A. 16b  — England only
#     B. 16c  — England subset (should mirror 16b if design is consistent)
#     C. 16c  — Scotland subset (same-city controls, potential attenuation bias)
#
# KEY QUESTION:
#   Does Scotland's post-matching balance on injury TREND variables meet the
#   |SMD| < 0.10 threshold? If not, including Scotland adds noisy comparisons.
#
# INPUTS (from data/processed/):
#   16b outputs: OA_matched_full_England.rds, OA_balance_tests_England.rds
#   16c outputs: OA_matched_full_mixed_England.rds, OA_matched_full_mixed_Scotland.rds,
#                OA_balance_tests_mixed.rds
#   Shared:      OA_matching_census.rds
#
# OUTPUTS (to output/diagnostics/):
#   compare_16b_16c_decision_table.csv
#   compare_16b_16c_love_plot.png
#   compare_16b_16c_smd_bars.png
# =============================================================================

library(tidyverse)
library(here)
library(patchwork)
library(scales)

select <- dplyr::select
filter <- dplyr::filter

dir.create(here("output", "diagnostics"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics")

# =============================================================================
# VARIABLE DEFINITIONS — must match 16b / 16c
# =============================================================================

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

stage1_key <- c(
  "road_length_km", "road_density_m_km2", "pct_minor_road",
  "dist_citycentre", "pop_density", "IMD",
  "Drive_Car_pct", "Walk_pct", "Bicycle_pct"
)

var_labels <- c(
  trend_car_KSI_pkm       = "Car KSI trend",
  trend_car_slight_pkm    = "Car slight trend",
  trend_cyc_KSI_pkm       = "Cycle KSI trend",
  trend_cyc_slight_pkm    = "Cycle slight trend",
  trend_ped_KSI_pkm       = "Ped KSI trend",
  trend_ped_slight_pkm    = "Ped slight trend",
  trend_other_KSI_pkm     = "Other KSI trend",
  trend_other_slight_pkm  = "Other slight trend",
  trend_total_pkm         = "Total trend",
  mean_car_KSI_pkm        = "Car KSI level",
  mean_car_slight_pkm     = "Car slight level",
  mean_cyc_KSI_pkm        = "Cycle KSI level",
  mean_cyc_slight_pkm     = "Cycle slight level",
  mean_ped_KSI_pkm        = "Ped KSI level",
  mean_ped_slight_pkm     = "Ped slight level",
  mean_other_KSI_pkm      = "Other KSI level",
  mean_other_slight_pkm   = "Other slight level",
  mean_total_pkm          = "Total level",
  road_length_km          = "Road length (km)",
  road_density_m_km2      = "Road density",
  pct_minor_road          = "% minor road",
  dist_citycentre         = "Dist. city centre",
  pop_density             = "Pop. density",
  IMD                     = "IMD",
  Drive_Car_pct           = "Drive car %",
  Walk_pct                = "Walk %",
  Bicycle_pct             = "Bicycle %"
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading matched datasets...\n")

full_16b_eng <- readRDS(here("data", "processed", "OA_matched_full_England.rds"))
bal_16b      <- readRDS(here("data", "processed", "OA_balance_tests_England.rds"))

full_16c_eng <- readRDS(here("data", "processed", "OA_matched_full_mixed_England.rds"))
full_16c_sco <- readRDS(here("data", "processed", "OA_matched_full_mixed_Scotland.rds"))
bal_16c      <- readRDS(here("data", "processed", "OA_balance_tests_mixed.rds"))

full_data    <- readRDS(here("data", "processed", "OA_matching_census.rds")) %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    ),
    treat_indicator = as.integer(treated_OA == 1)
  )

cat("Loaded.\n\n")

# =============================================================================
# WEIGHTED SMD HELPER
# =============================================================================

# Compute absolute weighted SMD for one variable in a matched dataset.
# Pooled SD is taken from the unmatched pool to keep it comparable across
# approaches (using matched SD would conflate balance with sample reduction).

compute_smd <- function(matched_data, var, unmatched_pool) {
  if (!var %in% names(matched_data) || !var %in% names(unmatched_pool)) return(NA_real_)

  treated  <- matched_data %>% filter(treat_indicator == 1)
  controls <- matched_data %>% filter(treat_indicator == 0)

  mean_t <- mean(treated[[var]], na.rm = TRUE)
  mean_c <- weighted.mean(controls[[var]], w = controls$weights, na.rm = TRUE)

  sd_t <- var(unmatched_pool[[var]][unmatched_pool$treat_indicator == 1], na.rm = TRUE)
  sd_c <- var(unmatched_pool[[var]][unmatched_pool$treat_indicator == 0], na.rm = TRUE)
  sd_pool <- sqrt((sd_t + sd_c) / 2)

  if (is.na(sd_pool) || sd_pool == 0) return(NA_real_)
  abs(mean_t - mean_c) / sd_pool
}

# Unweighted SMD for unmatched pool (before matching reference)
compute_smd_unmatched <- function(pool, var) {
  if (!var %in% names(pool)) return(NA_real_)
  t_vals <- pool[[var]][pool$treat_indicator == 1]
  c_vals <- pool[[var]][pool$treat_indicator == 0]
  sd_pool <- sqrt((var(t_vals, na.rm = TRUE) + var(c_vals, na.rm = TRUE)) / 2)
  if (is.na(sd_pool) || sd_pool == 0) return(NA_real_)
  abs(mean(t_vals, na.rm = TRUE) - mean(c_vals, na.rm = TRUE)) / sd_pool
}

# =============================================================================
# UNMATCHED POOLS (per country, for reference SMDs and pooled SDs)
# =============================================================================

treated_lads_england <- full_data %>%
  filter(treated_OA == 1, country == "England") %>%
  distinct(LAD24CD) %>% pull(LAD24CD)

unmatched_eng <- full_data %>%
  filter(
    country == "England",
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0, n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1),
    !(control_group1_OA == 1 & LAD24CD %in% treated_lads_england)
  )

unmatched_sco <- full_data %>%
  filter(
    country == "Scotland",
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0, n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  )

# =============================================================================
# COMPUTE SMD TABLES FOR ALL CONDITIONS
# =============================================================================

all_vars <- c(stage2_trends, stage2_levels, stage1_key)

conditions <- list(
  "Unmatched — England"    = list(data = unmatched_eng, pool = unmatched_eng,
                                   unmatched = TRUE),
  "Unmatched — Scotland"   = list(data = unmatched_sco, pool = unmatched_sco,
                                   unmatched = TRUE),
  "16b — England (matched)" = list(data = full_16b_eng,  pool = unmatched_eng,
                                   unmatched = FALSE),
  "16c — England (matched)" = list(data = full_16c_eng,  pool = unmatched_eng,
                                   unmatched = FALSE),
  "16c — Scotland (matched)"= list(data = full_16c_sco,  pool = unmatched_sco,
                                   unmatched = FALSE)
)

smd_long <- map_df(names(conditions), function(cond_name) {
  cond <- conditions[[cond_name]]
  map_df(all_vars, function(v) {
    smd_val <- if (cond$unmatched) {
      compute_smd_unmatched(cond$data, v)
    } else {
      compute_smd(cond$data, v, cond$pool)
    }
    tibble(
      condition = cond_name,
      variable  = v,
      label     = coalesce(var_labels[v], v),
      smd       = smd_val,
      var_type  = case_when(
        v %in% stage2_trends ~ "Stage 2: Injury trends",
        v %in% stage2_levels ~ "Stage 2: Injury levels",
        TRUE                 ~ "Stage 1: Structural"
      )
    )
  })
})

# =============================================================================
# SECTION 1 — DECISION SUMMARY TABLE
# =============================================================================

cat("=== DECISION SUMMARY TABLE ===\n\n")

effective_n <- function(d) {
  w <- d$weights[d$treat_indicator == 0]
  round(sum(w)^2 / sum(w^2), 0)
}

summary_tbl <- tribble(
  ~Approach,                        ~Data,               ~n_treated, ~n_controls, ~eff_N,
  "16b — England only",             full_16b_eng,        NA, NA, NA,
  "16c — England subset",           full_16c_eng,        NA, NA, NA,
  "16c — Scotland subset",          full_16c_sco,        NA, NA, NA,
  "16c — Combined (Eng + Scot)",
    bind_rows(full_16c_eng %>% mutate(country = "England"),
              full_16c_sco %>% mutate(country = "Scotland")),
    NA, NA, NA
) %>% select(-Data)

datasets <- list(
  "16b — England only"          = full_16b_eng,
  "16c — England subset"        = full_16c_eng,
  "16c — Scotland subset"       = full_16c_sco,
  "16c — Combined (Eng + Scot)" = bind_rows(full_16c_eng, full_16c_sco)
)

smd_stats <- smd_long %>%
  filter(condition %in% c("16b — England (matched)",
                           "16c — England (matched)",
                           "16c — Scotland (matched)")) %>%
  group_by(condition) %>%
  summarise(
    mean_smd_trends  = round(mean(smd[var_type == "Stage 2: Injury trends"],
                                  na.rm = TRUE), 4),
    max_smd_trends   = round(max(smd[var_type  == "Stage 2: Injury trends"],
                                  na.rm = TRUE), 4),
    mean_smd_all     = round(mean(smd, na.rm = TRUE), 4),
    n_trends_over_10 = sum(smd[var_type == "Stage 2: Injury trends"] > 0.10,
                            na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trend_balance_pass = if_else(max_smd_trends < 0.10, "PASS ✓", "FAIL ✗")
  )

size_tbl <- map_df(names(datasets), function(nm) {
  d <- datasets[[nm]]
  tibble(
    Approach     = nm,
    n_treated    = sum(d$treat_indicator == 1),
    n_controls   = sum(d$treat_indicator == 0),
    effective_N  = effective_n(d)
  )
})

decision_tbl <- size_tbl %>%
  left_join(
    smd_stats %>%
      rename(Approach = condition) %>%
      mutate(Approach = recode(Approach,
        "16b — England (matched)"  = "16b — England only",
        "16c — England (matched)"  = "16c — England subset",
        "16c — Scotland (matched)" = "16c — Scotland subset"
      )),
    by = "Approach"
  )

print(decision_tbl)
cat("\n")
cat("Threshold: max trend |SMD| < 0.10 for acceptable balance\n")
cat("NOTE: Scotland uses same-city controls — estimates may be attenuated.\n\n")

write_csv(decision_tbl,
          file.path(outdir, "compare_16b_16c_decision_table.csv"))
cat("Saved: compare_16b_16c_decision_table.csv\n\n")

# =============================================================================
# SECTION 2 — LOVE PLOT: MATCHED CONDITIONS
# =============================================================================

plot_vars <- c(stage2_trends, stage2_levels, stage1_key)

love_data <- smd_long %>%
  filter(variable %in% plot_vars) %>%
  mutate(
    condition = factor(condition, levels = c(
      "Unmatched — England",
      "Unmatched — Scotland",
      "16b — England (matched)",
      "16c — England (matched)",
      "16c — Scotland (matched)"
    )),
    label = factor(label, levels = rev(unique(label[order(var_type, label)])))
  )

pal <- c(
  "Unmatched — England"     = "#BBBBBB",
  "Unmatched — Scotland"    = "#DDAAAA",
  "16b — England (matched)" = "#2E6FAB",
  "16c — England (matched)" = "#27A060",
  "16c — Scotland (matched)"= "#9B3FAF"
)

shapes <- c(
  "Unmatched — England"     = 1,
  "Unmatched — Scotland"    = 1,
  "16b — England (matched)" = 16,
  "16c — England (matched)" = 17,
  "16c — Scotland (matched)"= 15
)

p_love <- ggplot(love_data, aes(x = smd, y = label,
                                colour = condition, shape = condition)) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dotted",
             colour = "#CCCCCC", linewidth = 0.4) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ var_type, ncol = 1, scales = "free_y") +
  labs(
    title    = "Covariate Balance: 16b vs 16c",
    subtitle = paste0(
      "Dashed = |SMD| 0.10 threshold | Dotted = 0.05\n",
      "16b: England only (other-city controls)\n",
      "16c England: same design as 16b for England\n",
      "16c Scotland: same-city controls — potential downward bias"
    ),
    x = "Absolute Standardised Mean Difference",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1A2E5A", size = 14),
    plot.subtitle   = element_text(colour = "#555555", size = 10),
    strip.text      = element_text(face = "bold", colour = "#1A2E5A"),
    strip.background= element_rect(fill = "#EEF2F8", colour = NA),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    panel.grid.major.y = element_line(colour = "#F0F0F0"),
    panel.grid.major.x = element_line(colour = "#E5E5E5")
  )

ggsave(file.path(outdir, "compare_16b_16c_love_plot.png"),
       p_love, width = 14, height = 18, dpi = 300, bg = "white")
cat("Saved: compare_16b_16c_love_plot.png\n\n")

# =============================================================================
# SECTION 3 — BAR CHART: MAX TREND |SMD| BY APPROACH
# =============================================================================

bar_data <- smd_long %>%
  filter(
    var_type  == "Stage 2: Injury trends",
    condition %in% c("16b — England (matched)",
                     "16c — England (matched)",
                     "16c — Scotland (matched)")
  ) %>%
  mutate(
    condition = recode(condition,
      "16b — England (matched)"  = "16b\nEngland only",
      "16c — England (matched)"  = "16c\nEngland",
      "16c — Scotland (matched)" = "16c\nScotland"
    ),
    label = factor(coalesce(var_labels[variable], variable))
  )

pal_bar <- c(
  "16b\nEngland only" = "#2E6FAB",
  "16c\nEngland"      = "#27A060",
  "16c\nScotland"     = "#9B3FAF"
)

p_bars <- ggplot(bar_data, aes(x = condition, y = smd, fill = condition)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#E74C3C", linewidth = 0.6) +
  scale_fill_manual(values = pal_bar, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ label, nrow = 3) +
  labs(
    title    = "Post-Matching |SMD| for Injury Trend Variables",
    subtitle = "Red dashed = 0.10 threshold | Trend variables drive the parallel-trends assumption",
    x        = NULL,
    y        = "Absolute SMD"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", colour = "#1A2E5A"),
    plot.subtitle   = element_text(colour = "#555555", size = 9),
    strip.text      = element_text(face = "bold", size = 9),
    strip.background= element_rect(fill = "#EEF2F8", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(outdir, "compare_16b_16c_smd_bars.png"),
       p_bars, width = 16, height = 10, dpi = 300, bg = "white")
cat("Saved: compare_16b_16c_smd_bars.png\n\n")

# =============================================================================
# SECTION 4 — PRINTED RECOMMENDATION GUIDE
# =============================================================================

cat(paste(rep("=", 65), collapse = ""), "\n")
cat("RECOMMENDATION GUIDE\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n")

scot_max  <- decision_tbl$max_smd_trends[decision_tbl$Approach == "16c — Scotland subset"]
scot_pass <- !is.na(scot_max) && scot_max < 0.10

eng_16b   <- decision_tbl$max_smd_trends[decision_tbl$Approach == "16b — England only"]
eng_16c   <- decision_tbl$max_smd_trends[decision_tbl$Approach == "16c — England subset"]
eng_consistent <- !is.na(eng_16b) && !is.na(eng_16c) &&
                  abs(eng_16b - eng_16c) < 0.01

cat(sprintf("England balance consistent between 16b and 16c: %s (|diff| = %.4f)\n",
            if (eng_consistent) "YES ✓" else "NO — investigate",
            abs(eng_16b - eng_16c)))
cat(sprintf("Scotland max trend |SMD|: %.4f  — Balance threshold (< 0.10): %s\n\n",
            scot_max, if (scot_pass) "PASS ✓" else "FAIL ✗"))

if (!eng_consistent) {
  cat("⚠  England results differ between 16b and 16c — check variable definitions.\n\n")
}

if (scot_pass) {
  cat("Scotland achieves acceptable trend balance.\n")
  cat("16c (England + Scotland combined) can be considered for primary analysis.\n")
  cat("Caveat: Scottish controls are same-city (potential attenuation bias).\n")
  cat("Recommendation: use 16c with Scotland reported separately and limitations flagged.\n")
} else {
  cat("Scotland does NOT achieve acceptable trend balance.\n")
  cat("Same-city controls share injury trends with treated OAs — matching cannot\n")
  cat("separate treatment effect from shared time trends.\n")
  cat("Recommendation: use 16b (England only) as the primary analysis.\n")
  cat("Scotland can be reported as a supplementary sensitivity analysis with\n")
  cat("strong caveats about potential attenuation bias.\n")
}

cat("\n", paste(rep("=", 65), collapse = ""), "\n")
