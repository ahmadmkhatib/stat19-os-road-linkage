
# =============================================================================
# MATCHING QUALITY BY SCHEME — 11-CITY COMPARISON
#
# INPUTS (from data/processed/):
#   OA_matched_full_mixed.rds         — full matched dataset (all schemes)
#   OA_matching_census.rds            — full unmatched pool
#
# Scheme identifier: column `scheme` in matched data
#
# OUTPUTS (to output/diagnostics/scheme_comparison/):
#   Per-scheme love plots  : fig_scheme_<name>_love.png     (one per scheme)
#   SMD summary heatmap    : fig_scheme_smd_heatmap.png
#   Weight/efficiency panel: fig_scheme_weights.png
#   Mahalanobis distance   : fig_scheme_mdist.png
#
# =============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(patchwork)
library(scales)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "scheme_comparison")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# THEME
# =============================================================================

theme_diag <- function(base_size = 13) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      colour = "#1A2E5A", margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base_size - 2, colour = "#555555",
                                      margin = margin(b = 10)),
      plot.caption     = element_text(size = base_size - 3, colour = "#888888",
                                      hjust = 0, margin = margin(t = 6)),
      axis.title       = element_text(size = base_size - 2, colour = "#333333"),
      axis.text        = element_text(size = base_size - 3, colour = "#444444"),
      strip.text       = element_text(size = base_size - 2, face = "bold",
                                      colour = "#1A2E5A"),
      strip.background = element_rect(fill = "#EEF2F8", colour = NA),
      panel.grid.major = element_line(colour = "#E5E9F0", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "#CCCCCC", fill = NA,
                                      linewidth = 0.35),
      legend.title     = element_text(size = base_size - 2, face = "bold"),
      legend.text      = element_text(size = base_size - 3),
      legend.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(10, 14, 10, 10)
    )
}

save_fig <- function(p, filename, width = 14, height = 10, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# Scheme colour palette — 11 distinct colours
SCHEME_COLS <- c(
  "#2E6FAB", "#D85A30", "#2ECC71", "#9B59B6", "#E67E22",
  "#1ABC9C", "#E74C3C", "#3498DB", "#F39C12", "#16A085", "#8E44AD"
)

# =============================================================================
# VARIABLE DEFINITIONS  (copy from your diagnostics script)
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

stage2_trends <- c(
  "trend_car_KSI_pkm",   "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",   "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",   "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm",
  "trend_total_pkm"
)
stage2_levels_log <- c(
  "log1p_mean_car_KSI_pkm",   "log1p_mean_car_slight_pkm",
  "log1p_mean_cyc_KSI_pkm",   "log1p_mean_cyc_slight_pkm",
  "log1p_mean_ped_KSI_pkm",   "log1p_mean_ped_slight_pkm",
  "log1p_mean_other_KSI_pkm", "log1p_mean_other_slight_pkm",
  "log1p_mean_total_pkm"
)
all_match_vars <- c(stage1_vars, stage2_trends, stage2_levels_log)

var_group_lookup <- c(
  setNames(rep("1. Road network",    length(stage1_road)),   stage1_road),
  setNames(rep("2. Urban geography", length(stage1_urban)),  stage1_urban),
  setNames(rep("3. Sociodemographic",length(stage1_socdem)), stage1_socdem),
  setNames(rep("4. Injury trends",   length(stage2_trends)), stage2_trends),
  setNames(rep("5. Injury levels",   length(stage2_levels_log)), stage2_levels_log)
)

var_labels <- c(
  log1p_road_length_km          = "Road length (km)",
  log1p_road_density_m_km2      = "Road density (m/km²)",
  log_area_km2                  = "Area (km²)",
  pct_A_road                    = "% A-road",
  pct_B_road                    = "% B-road",
  pct_minor_road                = "% Minor road",
  log1p_dist_citycentre         = "Distance to city centre",
  log1p_pop_density             = "Population density",
  log1p_business_retail_per_km2 = "Retail businesses/km²",
  IMD                           = "IMD",
  cars_one_pct                  = "% 1 car",
  cars_twoPlus_pct              = "% 2+ cars",
  Drive_Car_pct                 = "% drive car",
  Passenger_Car_pct             = "% car passenger",
  Walk_pct                      = "% walk",
  Bicycle_pct                   = "% bicycle",
  bus_Coach_pct                 = "% bus/coach",
  Train_pct                     = "% train",
  Underground_train_tram_pct    = "% underground/tram",
  Taxi_pct                      = "% taxi",
  workAthome_pct                = "% work from home",
  Other_pct                     = "% other mode",
  White_pct                     = "% White",
  Mixed_pct                     = "% Mixed",
  Asian_pct                     = "% Asian",
  Black_pct                     = "% Black",
  age_under15_pct               = "% aged <15",
  age_15to24_pct                = "% aged 15–24",
  age_25to44_pct                = "% aged 25–44",
  age_45to64_pct                = "% aged 45–64",
  age_65to84_pct                = "% aged 65–84",
  trend_car_KSI_pkm             = "Trend: car KSI",
  trend_car_slight_pkm          = "Trend: car slight",
  trend_cyc_KSI_pkm             = "Trend: cycling KSI",
  trend_cyc_slight_pkm          = "Trend: cycling slight",
  trend_ped_KSI_pkm             = "Trend: ped KSI",
  trend_ped_slight_pkm          = "Trend: ped slight",
  trend_other_KSI_pkm           = "Trend: other KSI",
  trend_other_slight_pkm        = "Trend: other slight",
  trend_total_pkm               = "Trend: total",
  log1p_mean_car_KSI_pkm        = "Mean: car KSI (log)",
  log1p_mean_car_slight_pkm     = "Mean: car slight (log)",
  log1p_mean_cyc_KSI_pkm        = "Mean: cycling KSI (log)",
  log1p_mean_cyc_slight_pkm     = "Mean: cycling slight (log)",
  log1p_mean_ped_KSI_pkm        = "Mean: ped KSI (log)",
  log1p_mean_ped_slight_pkm     = "Mean: ped slight (log)",
  log1p_mean_other_KSI_pkm      = "Mean: other KSI (log)",
  log1p_mean_other_slight_pkm   = "Mean: other slight (log)",
  log1p_mean_total_pkm          = "Mean: total (log)"
)

matched_full <- readRDS(here("data", "processed", "OA_matched_full_mixed.rds"))
full_data    <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Controls have no scheme — propagate from their matched treated OA via pairs file.
# OA_matching_pairs_mixed.rds links each control OA back to its treated OA.
if (any(matched_full$scheme %in% c("Control", NA))) {
  
  pairs <- readRDS(here("data", "processed", "OA_matching_pairs_mixed.rds"))
  cat("Pairs file columns:", paste(names(pairs), collapse = ", "), "\n")
  
  # Detect treated/control OA columns robustly
  treated_col <- names(pairs)[str_detect(tolower(names(pairs)), "treat")][1]
  control_col <- names(pairs)[str_detect(tolower(names(pairs)), "control|donor")][1]
  cat("Using pairs columns — treated:", treated_col, "| control:", control_col, "\n\n")
  
  # Build lookup: control OA -> scheme via treated OA -> scheme
  treated_scheme <- matched_full %>%
    filter(treat_indicator == 1) %>%
    select(OA, scheme)
  
  ctrl_scheme_lookup <- pairs %>%
    select(treated_OA = all_of(treated_col),
           control_OA = all_of(control_col)) %>%
    left_join(treated_scheme, by = c("treated_OA" = "OA")) %>%
    select(OA = control_OA, scheme_from_treated = scheme) %>%
    distinct(OA, .keep_all = TRUE)
  
  matched_full <- matched_full %>%
    left_join(ctrl_scheme_lookup, by = "OA") %>%
    mutate(scheme = case_when(
      treat_indicator == 1        ~ scheme,
      !is.na(scheme_from_treated) ~ scheme_from_treated,
      TRUE                        ~ scheme
    )) %>%
    select(-scheme_from_treated)
  
  n_still_bad <- sum(matched_full$scheme %in% c("Control", NA), na.rm = TRUE)
  if (n_still_bad > 0)
    warning(n_still_bad, " rows still labelled Control/NA — ",
            "check the column names printed above match your pairs file.")
  else
    cat("Scheme propagated to all control OAs via pairs file.\n\n")
}

# Derive log vars + age bands on unmatched pool (same as diagnostics script)
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
    log1p_mean_car_KSI_pkm        = log1p(pmax(mean_car_KSI_pkm,        0)),
    log1p_mean_car_slight_pkm     = log1p(pmax(mean_car_slight_pkm,     0)),
    log1p_mean_cyc_KSI_pkm        = log1p(pmax(mean_cyc_KSI_pkm,        0)),
    log1p_mean_cyc_slight_pkm     = log1p(pmax(mean_cyc_slight_pkm,     0)),
    log1p_mean_ped_KSI_pkm        = log1p(pmax(mean_ped_KSI_pkm,        0)),
    log1p_mean_ped_slight_pkm     = log1p(pmax(mean_ped_slight_pkm,     0)),
    log1p_mean_other_KSI_pkm      = log1p(pmax(mean_other_KSI_pkm,      0)),
    log1p_mean_other_slight_pkm   = log1p(pmax(mean_other_slight_pkm,   0)),
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
  ) %>%
  add_log_vars()

schemes <- sort(unique(matched_full$scheme))
cat("Schemes found:", paste(schemes, collapse = ", "), "\n")
cat("N schemes:", length(schemes), "\n\n")

# Assign colours
scheme_col_map <- setNames(SCHEME_COLS[seq_along(schemes)], schemes)

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

# Build SMD table for one scheme:
#   before = unmatched pool filtered to that scheme's treated OAs + full control pool
#   after  = matched data for that scheme
smd_for_scheme <- function(scheme_name, vars = all_match_vars) {
  treated_oas <- matched_full %>%
    filter(scheme == scheme_name, treat_indicator == 1) %>%
    pull(OA)
  
  # "Before": treated OAs for this scheme vs full unmatched control pool
  before_df <- unmatched_pool %>%
    filter(treated_OA == 1 & OA %in% treated_oas |
             control_group1_OA == 1 | control_group2_OA == 1) %>%
    mutate(treat_indicator = as.integer(OA %in% treated_oas))
  
  after_df <- matched_full %>%
    filter(scheme == scheme_name)
  
  map_df(vars, function(v) {
    tibble(
      scheme     = scheme_name,
      variable   = v,
      label      = coalesce(var_labels[v], v),
      var_group  = coalesce(var_group_lookup[v], "Other"),
      smd_before = round(compute_smd(before_df, v), 4),
      smd_after  = round(compute_smd(after_df,  v), 4)
    )
  })
}

# =============================================================================
# COMPUTE SMD FOR ALL SCHEMES
# =============================================================================

cat("Computing SMDs for all schemes...\n")
smd_all <- map_df(schemes, function(s) {
  cat("  →", s, "\n")
  smd_for_scheme(s, vars = stage2_trends)
})
cat("\n")

# =============================================================================
# PLOT 1 — PER-SCHEME LOVE PLOTS (one PNG per scheme)
# =============================================================================

cat("=== Per-scheme love plots ===\n")

make_love_plot_scheme <- function(smd_df, scheme_name) {
  plot_data <- smd_df %>%
    filter(!is.na(smd_before), !is.na(smd_after)) %>%
    mutate(label = factor(label, levels = label[order(abs(smd_before))])) %>%
    pivot_longer(c(smd_before, smd_after),
                 names_to = "timing", values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_before", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  n_treated  <- sum(matched_full$scheme == scheme_name &
                      matched_full$treat_indicator == 1)
  n_controls <- sum(matched_full$scheme == scheme_name &
                      matched_full$treat_indicator == 0)
  
  ggplot(plot_data,
         aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = 0.10, linetype = "dashed",
               colour = "#999999", linewidth = 0.5) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.4) +
    geom_point(size = 3.5) +
    scale_colour_manual(
      values = c("Before matching" = "#E74C3C",
                 "After matching"  = "#2ECC71")) +
    scale_shape_manual(
      values = c("Before matching" = 16, "After matching" = 17)) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("Injury trend balance — ", scheme_name),
      subtitle = paste0("Treated OAs: ", n_treated,
                        " | Matched controls: ", n_controls),
      x        = "|SMD|", y = NULL,
      colour   = NULL, shape = NULL,
      caption  = "Dashed = 0.10 threshold  |  ● Before matching  ▲ After matching"
    ) +
    theme_diag() +
    theme(legend.position = "bottom",
          axis.text.y     = element_text(size = 11),
          axis.text.x     = element_text(size = 11))
}

# Diagnose NA coverage before plotting
na_check <- smd_all %>%
  group_by(scheme) %>%
  summarise(
    n_vars        = n(),
    n_na_before   = sum(is.na(smd_before)),
    n_na_after    = sum(is.na(smd_after)),
    n_plottable   = sum(!is.na(smd_before) & !is.na(smd_after)),
    .groups       = "drop"
  )
cat("SMD NA diagnostic:\n")
print(na_check, n = Inf)
cat("\n")

# Schemes where treated OAs were not found in unmatched_pool will have
# smd_before = NA for all variables — flag them clearly
problem_schemes <- na_check %>% filter(n_plottable == 0) %>% pull(scheme)
if (length(problem_schemes) > 0) {
  cat("WARNING: the following schemes have zero plottable variables\n")
  cat("  (treated OAs not matched to unmatched_pool — check OA identifiers):\n")
  cat(" ", paste(problem_schemes, collapse = ", "), "\n\n")
}

# Derive country per scheme from matched_full (treated rows only, no NAs)
scheme_country <- matched_full %>%
  filter(treat_indicator == 1, !is.na(scheme), !is.na(country)) %>%
  distinct(scheme, country) %>%
  mutate(country = case_when(
    str_detect(tolower(country), "scot")     ~ "Scotland",
    str_detect(tolower(country), "eng|wal")  ~ "England",
    TRUE                                     ~ country
  )) %>%
  filter(!is.na(country))

cat("Scheme-country mapping:\n")
print(as.data.frame(scheme_country))
cat("\n")

missing_country <- setdiff(schemes, scheme_country$scheme)
if (length(missing_country) > 0)
  cat("WARNING: no country found for schemes:", paste(missing_country, collapse = ", "), "\n\n")

# Shared y-axis order: variables sorted by mean |SMD| before, across all schemes
var_order <- smd_all %>%
  filter(!is.na(smd_before)) %>%
  mutate(label = coalesce(var_labels[variable], variable)) %>%
  group_by(label) %>%
  summarise(mean_before = mean(abs(smd_before), na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_before) %>%
  pull(label)

make_love_grid <- function(country_name) {
  country_schemes <- scheme_country %>%
    filter(country == country_name) %>%
    pull(scheme)
  
  # Sample size label for each scheme strip
  size_labels <- matched_full %>%
    filter(scheme %in% country_schemes) %>%
    group_by(scheme) %>%
    summarise(n_t = sum(treat_indicator == 1),
              n_c = sum(treat_indicator == 0), .groups = "drop") %>%
    mutate(scheme_lab = paste0(scheme, "\nT=", n_t, "  C=", n_c))
  
  plot_data <- smd_all %>%
    filter(scheme %in% country_schemes,
           !is.na(smd_before), !is.na(smd_after)) %>%
    mutate(label = coalesce(var_labels[variable], variable),
           label = factor(label, levels = var_order)) %>%
    left_join(size_labels, by = "scheme") %>%
    pivot_longer(c(smd_before, smd_after),
                 names_to = "timing", values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_before", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  n_schemes <- length(country_schemes)
  
  ggplot(plot_data,
         aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = 0.10, linetype = "dashed",
               colour = "#999999", linewidth = 0.5) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.4) +
    geom_point(size = 3.5) +
    scale_colour_manual(
      values = c("Before matching" = "#E74C3C",
                 "After matching"  = "#2ECC71")) +
    scale_shape_manual(
      values = c("Before matching" = 16, "After matching" = 17)) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06))) +
    facet_wrap(~scheme_lab, nrow = 1) +
    labs(
      title    = paste0("Injury trend balance \u2014 ", country_name, " schemes"),
      subtitle = "T = treated OAs  |  C = matched controls  |  shared y-axis across all panels",
      x = "|SMD|", y = NULL,
      colour = NULL, shape = NULL,
      caption = "Dashed = 0.10 threshold  \u25cf Before matching  \u25b2 After matching"
    ) +
    theme_diag() +
    theme(
      legend.position = "bottom",
      axis.text.y     = element_text(size = 10),
      axis.text.x     = element_text(size = 9),
      strip.text      = element_text(size = 10, face = "bold"),
      panel.spacing.x = unit(0.8, "lines")
    )
}

p_eng <- make_love_grid("England")
p_sco <- make_love_grid("Scotland")

n_eng <- sum(scheme_country$country == "England")
n_sco <- sum(scheme_country$country == "Scotland")

save_fig(p_eng, "fig_love_england_schemes.png",
         width = max(10, n_eng * 2.8), height = 9)
save_fig(p_sco, "fig_love_scotland_schemes.png",
         width = max(10, n_sco * 2.8), height = 9)
cat("\n")

# --- Individual per-scheme love plots ----------------------------------------
cat("=== Individual per-scheme love plots ===\n")

for (s in schemes) {
  smd_scheme <- smd_all %>% filter(scheme == s)
  
  if (nrow(smd_scheme %>% filter(!is.na(smd_before), !is.na(smd_after))) == 0) {
    cat("  SKIP (no plottable data):", s, "\n")
    next
  }
  
  ctry <- scheme_country %>% filter(scheme == s) %>% pull(country)
  ctry <- if (length(ctry) == 0) "" else ctry
  
  plot_data <- smd_scheme %>%
    filter(!is.na(smd_before), !is.na(smd_after)) %>%
    mutate(label = coalesce(var_labels[variable], variable),
           label = factor(label, levels = var_order)) %>%
    pivot_longer(c(smd_before, smd_after),
                 names_to = "timing", values_to = "smd") %>%
    mutate(
      timing = if_else(timing == "smd_before", "Before matching", "After matching"),
      timing = factor(timing, levels = c("Before matching", "After matching"))
    )
  
  n_t <- sum(matched_full$scheme == s & matched_full$treat_indicator == 1)
  n_c <- sum(matched_full$scheme == s & matched_full$treat_indicator == 0)
  
  p <- ggplot(plot_data,
              aes(x = abs(smd), y = label, colour = timing, shape = timing)) +
    geom_vline(xintercept = 0.10, linetype = "dashed",
               colour = "#999999", linewidth = 0.5) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.3) +
    geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.4) +
    geom_point(size = 3.5) +
    scale_colour_manual(
      values = c("Before matching" = "#E74C3C",
                 "After matching"  = "#2ECC71")) +
    scale_shape_manual(
      values = c("Before matching" = 16, "After matching" = 17)) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(
      title    = paste0("Injury trend balance \u2014 ", s,
                        if (nchar(ctry) > 0) paste0(" (", ctry, ")") else ""),
      subtitle = paste0("Treated OAs: ", n_t, "  |  Matched controls: ", n_c),
      x = "|SMD|", y = NULL,
      colour = NULL, shape = NULL,
      caption = "Dashed = 0.10 threshold  \u25cf Before matching  \u25b2 After matching"
    ) +
    theme_diag() +
    theme(legend.position = "bottom",
          axis.text.y     = element_text(size = 11),
          axis.text.x     = element_text(size = 10))
  
  fname <- paste0("fig_love_", gsub("[^a-zA-Z0-9]", "_", s), ".png")
  save_fig(p, fname, width = 8, height = 7)
}
cat("\n")

# =============================================================================
# =============================================================================
# PLOT 2 — SMD HEATMAP: separate plots for England and Scotland
# Rows = trend variables, Columns = schemes, Fill = |SMD| after matching
# =============================================================================

cat("=== SMD summary heatmaps (England + Scotland) ===\n")

smd_breaks  <- c(0, 0.05, 0.10, 0.15, 0.20, Inf)
smd_labels  <- c("<0.05", "0.05-0.10", "0.10-0.15", "0.15-0.20", ">0.20")
smd_colours <- c("#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027")

make_heatmap <- function(country_name) {
  country_schemes <- scheme_country %>%
    filter(country == country_name) %>%
    pull(scheme)
  
  hdat <- smd_all %>%
    filter(variable %in% stage2_trends,
           scheme %in% country_schemes) %>%
    mutate(
      label  = coalesce(var_labels[variable], variable),
      scheme = factor(scheme, levels = sort(country_schemes)),
      label  = factor(label, levels = rev(var_order))
    ) %>%
    group_by(scheme, label) %>%
    summarise(
      smd_after  = mean(abs(smd_after),  na.rm = TRUE),
      smd_before = mean(abs(smd_before), na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      smd_band   = cut(smd_after, breaks = smd_breaks,
                       labels = smd_labels, right = FALSE, include.lowest = TRUE),
      cell_label = sprintf("%.3f", smd_after)
    )
  
  ggplot(hdat, aes(x = scheme, y = label, fill = smd_band)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3, colour = "#222222") +
    scale_fill_manual(
      values = setNames(smd_colours, smd_labels),
      name   = "|SMD| after matching",
      drop   = FALSE
    ) +
    scale_x_discrete(guide = guide_axis(angle = 35)) +
    labs(
      title    = paste0("Injury trend balance \u2014 ", country_name),
      subtitle = "Each cell = |SMD| after matching",
      x = NULL, y = NULL,
      caption  = "Green < 0.05  \u00b7  yellow 0.05-0.10  \u00b7  orange 0.10-0.15  \u00b7  red > 0.20"
    ) +
    theme_diag(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.text.x      = element_text(size = 11, face = "bold"),
      axis.text.y      = element_text(size = 10),
      legend.position  = "bottom",
      legend.key.width = unit(1.2, "cm")
    )
}

p_heat_eng <- make_heatmap("England")
p_heat_sco <- make_heatmap("Scotland")

n_eng <- sum(scheme_country$country == "England")
n_sco <- sum(scheme_country$country == "Scotland")

save_fig(p_heat_eng, "fig_heatmap_england.png",
         width = max(8, n_eng * 1.6 + 4), height = 7)
save_fig(p_heat_sco, "fig_heatmap_scotland.png",
         width = max(8, n_sco * 1.6 + 4), height = 7)
cat("\n")

# --- Individual per-scheme heatmaps ------------------------------------------
cat("=== Individual per-scheme heatmaps ===\n")

for (s in schemes) {
  hdat_s <- smd_all %>%
    filter(variable %in% stage2_trends, scheme == s,
           !is.na(smd_after)) %>%
    mutate(
      label  = coalesce(var_labels[variable], variable),
      label  = factor(label, levels = rev(var_order))
    ) %>%
    mutate(
      smd_band   = cut(abs(smd_after), breaks = smd_breaks,
                       labels = smd_labels, right = FALSE, include.lowest = TRUE),
      cell_label = sprintf("%.3f", abs(smd_after))
    )
  
  if (nrow(hdat_s) == 0) {
    cat("  SKIP (no data):", s, "\n")
    next
  }
  
  ctry <- scheme_country %>% filter(scheme == s) %>% pull(country)
  ctry <- if (length(ctry) == 0) "" else ctry
  
  p_h <- ggplot(hdat_s, aes(x = "After matching", y = label, fill = smd_band)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.5, colour = "#222222") +
    scale_fill_manual(
      values = setNames(smd_colours, smd_labels),
      name   = "|SMD| after matching",
      drop   = FALSE
    ) +
    labs(
      title    = paste0("Injury trend balance \u2014 ", s,
                        if (nchar(ctry) > 0) paste0(" (", ctry, ")") else ""),
      subtitle = "Each cell = |SMD| after matching for one trend variable",
      x = NULL, y = NULL,
      caption  = "Green < 0.05  \u00b7  yellow 0.05-0.10  \u00b7  orange 0.10-0.15  \u00b7  red > 0.20"
    ) +
    theme_diag(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(size = 11),
      legend.position  = "bottom",
      legend.key.width = unit(1.2, "cm")
    )
  
  fname <- paste0("fig_heatmap_", gsub("[^a-zA-Z0-9]", "_", s), ".png")
  save_fig(p_h, fname, width = 5, height = 7)
}
cat("\n")

# =============================================================================
# PLOT 3 — WEIGHT / EFFICIENCY PANEL (all schemes, 3 sub-panels)
#
# Panel A: effective N as % of nominal N (efficiency)
# Panel B: max weight per scheme
# Panel C: % controls at weight cap (≥5)
# =============================================================================

cat("=== Weight diagnostics panel ===\n")

weight_stats <- matched_full %>%
  filter(treat_indicator == 0) %>%
  group_by(scheme) %>%
  summarise(
    n_controls    = n(),
    n_treated     = sum(matched_full$scheme == first(scheme) &
                          matched_full$treat_indicator == 1),
    eff_N         = sum(weights)^2 / sum(weights^2),
    efficiency    = eff_N / n_controls,
    max_weight    = max(weights),
    median_weight = median(weights),
    pct_at_cap    = mean(weights >= 5) * 100,
    .groups       = "drop"
  ) %>%
  mutate(scheme = factor(scheme, levels = schemes))

# Panel A: efficiency
pA <- ggplot(weight_stats, aes(x = scheme, y = efficiency,
                               fill = scheme, colour = scheme)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_hline(yintercept = c(0.5, 0.75), linetype = "dashed",
             colour = c("#E74C3C", "#F39C12"), linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", efficiency)),
            vjust = -0.4, size = 3, fontface = "bold",
            colour = "#333333") +
  scale_fill_manual(values   = scheme_col_map, guide = "none") +
  scale_colour_manual(values = scheme_col_map, guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format()) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(title = "A — Weight efficiency (Eff N / Nominal N)",
       x = NULL, y = "Efficiency",
       caption = "Orange dashed = 0.75 · Red dashed = 0.50") +
  theme_diag(base_size = 11)

# Panel B: max weight
pB <- ggplot(weight_stats, aes(x = scheme, y = max_weight,
                               fill = scheme, colour = scheme)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_hline(yintercept = 5, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", max_weight)),
            vjust = -0.4, size = 3, fontface = "bold",
            colour = "#333333") +
  scale_fill_manual(values   = scheme_col_map, guide = "none") +
  scale_colour_manual(values = scheme_col_map, guide = "none") +
  scale_y_continuous(limits = c(0, max(weight_stats$max_weight) * 1.15)) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(title = "B — Max control weight",
       x = NULL, y = "Max weight",
       caption = "Dashed = cap at 5") +
  theme_diag(base_size = 11)

# Panel C: % at cap
pC <- ggplot(weight_stats, aes(x = scheme, y = pct_at_cap,
                               fill = scheme, colour = scheme)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct_at_cap)),
            vjust = -0.4, size = 3, fontface = "bold",
            colour = "#333333") +
  scale_fill_manual(values   = scheme_col_map, guide = "none") +
  scale_colour_manual(values = scheme_col_map, guide = "none") +
  scale_y_continuous(limits = c(0, max(weight_stats$pct_at_cap) * 1.2 + 1)) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(title = "C — % controls at weight cap (≥5)",
       x = NULL, y = "% at cap") +
  theme_diag(base_size = 11)

# Sample size inset table — treated + controls per scheme
size_tbl <- matched_full %>%
  count(scheme, treat_indicator) %>%
  mutate(group = if_else(treat_indicator == 1, "Treated", "Controls")) %>%
  select(scheme, group, n) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  mutate(ratio = sprintf("1:%.1f", Controls / Treated))

pD <- ggplot(size_tbl %>% pivot_longer(c(Treated, Controls),
                                       names_to = "group", values_to = "n") %>%
               mutate(group   = factor(group, levels = c("Treated","Controls")),
                      scheme  = factor(scheme, levels = schemes)),
             aes(x = scheme, y = n, fill = group)) +
  geom_col(position = "dodge", width = 0.65, alpha = 0.85) +
  geom_text(
    data = size_tbl %>% mutate(scheme = factor(scheme, levels = schemes)),
    aes(x = scheme, y = Controls + max(size_tbl$Controls) * 0.04,
        label = paste0("1:", round(Controls/Treated,1))),
    inherit.aes = FALSE, size = 2.8, colour = "#333333", vjust = 0
  ) +
  scale_fill_manual(values = c(Treated = "#D85A30", Controls = "#2E6FAB")) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(title = "D — Sample sizes (ratio label above each pair)",
       x = NULL, y = "N OAs", fill = NULL) +
  theme_diag(base_size = 11) +
  theme(legend.position = "bottom")

p_weights <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Weight & sample-size diagnostics by scheme",
    subtitle = "All figures based on post-capping weights (cap = 5)",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", colour = "#1A2E5A"),
      plot.subtitle = element_text(size = 11, colour = "#555555"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

save_fig(p_weights, "fig_scheme_weights.png", width = 18, height = 14)
cat("\n")

# =============================================================================
# PLOT 4 — MAHALANOBIS DISTANCE (if mdist column present)
# Four sub-panels:
#   A: median distance per scheme (bar)
#   B: P90 distance per scheme (bar)
#   C: % treated OAs with distance < 5
#   D: overlapping ECDF curves, coloured by scheme
# =============================================================================

cat("=== Mahalanobis distance panel ===\n")

if ("mdist" %in% names(matched_full)) {
  
  mdist_data <- matched_full %>%
    filter(treat_indicator == 1, !is.na(mdist)) %>%
    mutate(scheme = factor(scheme, levels = schemes))
  
  mdist_stats <- mdist_data %>%
    group_by(scheme) %>%
    summarise(
      n         = n(),
      median_d  = median(mdist),
      p90_d     = quantile(mdist, 0.90),
      pct_lt5   = mean(mdist < 5) * 100,
      pct_lt10  = mean(mdist < 10) * 100,
      .groups   = "drop"
    )
  
  pM_A <- ggplot(mdist_stats, aes(x = fct_reorder(scheme, median_d),
                                  y = median_d, fill = scheme)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_text(aes(label = sprintf("%.1f", median_d)),
              hjust = -0.15, size = 3, fontface = "bold", colour = "#333333") +
    scale_fill_manual(values = scheme_col_map, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    coord_flip() +
    labs(title = "A — Median Mahalanobis distance",
         x = NULL, y = "Median distance",
         caption = "Schemes ordered by median distance") +
    theme_diag(base_size = 11)
  
  pM_B <- ggplot(mdist_stats, aes(x = fct_reorder(scheme, p90_d),
                                  y = p90_d, fill = scheme)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_hline(yintercept = c(10, 20),
               linetype   = c("dashed", "dotted"),
               colour     = c("#E74C3C", "#8B0000"),
               linewidth  = 0.5) +
    geom_text(aes(label = sprintf("%.1f", p90_d)),
              hjust = -0.15, size = 3, fontface = "bold", colour = "#333333") +
    scale_fill_manual(values = scheme_col_map, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    coord_flip() +
    labs(title = "B — P90 Mahalanobis distance",
         x = NULL, y = "P90 distance",
         caption = "Red dashed = 10 · Dark red dotted = 20") +
    theme_diag(base_size = 11)
  
  pM_C <- ggplot(mdist_stats %>%
                   pivot_longer(c(pct_lt5, pct_lt10),
                                names_to = "threshold", values_to = "pct") %>%
                   mutate(threshold = if_else(threshold == "pct_lt5",
                                              "< 5", "< 10")),
                 aes(x = fct_reorder(scheme, pct), y = pct,
                     fill = scheme, alpha = threshold)) +
    geom_col(position = "dodge", width = 0.65) +
    scale_fill_manual(values = scheme_col_map, guide = "none") +
    scale_alpha_manual(values = c("< 5" = 0.95, "< 10" = 0.55),
                       name = "Distance threshold") +
    scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
    coord_flip() +
    labs(title = "C — % treated OAs with distance below threshold",
         x = NULL, y = "% treated OAs") +
    theme_diag(base_size = 11) +
    theme(legend.position = "bottom")
  
  # ECDF by scheme
  pM_D <- ggplot(mdist_data,
                 aes(x = mdist, colour = scheme)) +
    stat_ecdf(linewidth = 0.9, geom = "step") +
    geom_vline(xintercept = c(5, 10),
               linetype   = c("dashed", "dotted"),
               colour     = c("#888888", "#CC3333"),
               linewidth  = 0.5) +
    scale_colour_manual(values = scheme_col_map, name = "Scheme") +
    coord_cartesian(xlim = c(0, 35)) +
    labs(title    = "D — ECDF of Mahalanobis distance by scheme",
         x        = "Stage 2 Mahalanobis distance",
         y        = "Cumulative proportion",
         caption  = "Grey dashed = 5 · Red dotted = 10") +
    theme_diag(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text      = element_text(size = 8),
          legend.key.width = unit(0.6, "cm"))
  
  p_mdist <- (pM_A | pM_B) / (pM_C | pM_D) +
    plot_annotation(
      title    = "Mahalanobis distance diagnostics by scheme",
      subtitle = "Lower distance = better match quality; based on treated OAs only",
      theme    = theme(
        plot.title      = element_text(size = 14, face = "bold", colour = "#1A2E5A"),
        plot.subtitle   = element_text(size = 11, colour = "#555555"),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    )
  
  save_fig(p_mdist, "fig_scheme_mdist.png", width = 18, height = 14)
  
} else {
  cat("  WARNING: mdist column not found — skipping distance panel.\n")
  cat("  Re-run matching script to include mdist in matched data.\n")
}

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("\n================================================================\n")
cat("SCHEME COMPARISON SUMMARY\n")
cat("================================================================\n\n")

summary_tbl <- smd_all %>%
  group_by(scheme) %>%
  summarise(
    mean_smd_before = round(mean(abs(smd_before), na.rm = TRUE), 3),
    mean_smd_after  = round(mean(abs(smd_after),  na.rm = TRUE), 3),
    pct_balanced    = round(mean(abs(smd_after) < 0.10, na.rm = TRUE) * 100, 1),
    max_smd_after   = round(max(abs(smd_after),  na.rm = TRUE), 3),
    .groups         = "drop"
  ) %>%
  left_join(
    matched_full %>%
      filter(treat_indicator == 0) %>%
      group_by(scheme) %>%
      summarise(
        n_controls = n(),
        eff_N      = round(sum(weights)^2 / sum(weights^2), 0),
        efficiency = round((sum(weights)^2 / sum(weights^2)) / n(), 3),
        .groups    = "drop"
      ),
    by = "scheme"
  ) %>%
  left_join(
    matched_full %>%
      filter(treat_indicator == 1) %>%
      count(scheme, name = "n_treated"),
    by = "scheme"
  ) %>%
  arrange(mean_smd_after)

print(summary_tbl, n = Inf)

cat("\nOUTPUTS SAVED TO:", outdir, "\n")
cat("  Per-scheme love plots : fig_scheme_<name>_love.png (one per scheme)\n")
cat("  SMD heatmap           : fig_scheme_smd_heatmap.png\n")
cat("  Weight diagnostics    : fig_scheme_weights.png\n")
cat("  Mahalanobis distance  : fig_scheme_mdist.png\n")
cat("================================================================\n")


#### GALSGOW 

# How many unique controls does Glasgow actually have?
matched_full %>%
  filter(scheme == "Glasgow", treat_indicator == 0) %>%
  summarise(
    n_controls      = n(),
    n_unique        = n_distinct(OA),
    mean_weight     = mean(weights),
    max_weight      = max(weights),
    pct_at_cap      = mean(weights >= 5) * 100
  )

# Which control OAs are doing all the heavy lifting?
matched_full %>%
  filter(scheme == "Glasgow", treat_indicator == 0) %>%
  arrange(desc(weights)) %>%
  select(OA, weights, trend_total_pkm, trend_car_KSI_pkm) %>%
  head(10)

# What is the Mahalanobis distance distribution for Glasgow treated OAs?
matched_full %>%
  filter(scheme == "Glasgow", treat_indicator == 1) %>%
  summarise(
    median_mdist = median(mdist, na.rm = TRUE),
    p90_mdist    = quantile(mdist, 0.9, na.rm = TRUE),
    max_mdist    = max(mdist, na.rm = TRUE)
  )

###  10% of Glasgow's treated OAs have extremely poor matches

# Which treated OAs have terrible matches?
glasgow_treated <- matched_full %>%
  filter(scheme == "Glasgow", treat_indicator == 1) %>%
  arrange(desc(mdist)) %>%
  select(OA, mdist, trend_total_pkm, trend_car_KSI_pkm, 
         trend_cyc_KSI_pkm, trend_cyc_slight_pkm, trend_ped_slight_pkm)

head(glasgow_treated, 15)

# How many are above common distance thresholds?
glasgow_treated %>%
  summarise(
    pct_gt5  = mean(mdist > 5)  * 100,
    pct_gt10 = mean(mdist > 10) * 100,
    pct_gt20 = mean(mdist > 20) * 100,
    n_gt5    = sum(mdist > 5),
    n_gt10   = sum(mdist > 10),
    n_gt20   = sum(mdist > 20)
  )

# Compare trend distributions: well-matched vs poorly-matched treated OAs
matched_full %>%
  filter(scheme == "Glasgow", treat_indicator == 1) %>%
  mutate(match_quality = if_else(mdist > 10, "Poor (mdist > 10)", "Good (mdist <= 10)")) %>%
  group_by(match_quality) %>%
  summarise(
    n                    = n(),
    mean_trend_total     = mean(trend_total_pkm,     na.rm = TRUE),
    mean_trend_cyc_KSI   = mean(trend_cyc_KSI_pkm,   na.rm = TRUE),
    mean_trend_ped_slight = mean(trend_ped_slight_pkm, na.rm = TRUE),
    mean_trend_car_KSI   = mean(trend_car_KSI_pkm,   na.rm = TRUE),
    .groups = "drop"
  )

############################

# Schemes flagged from heatmap: Dundee, Glasgow, Sheffield, Bristol
# (Aberdeen and Newcastle have scattered orange cells — included for completeness)

# Scheme colours 
SCHEME_COLS <- c(
  Dundee    = "#D85A30",
  Glasgow   = "#9B59B6",
  Sheffield = "#E67E22",
  Bristol   = "#2E6FAB",
  Aberdeen  = "#1ABC9C",
  Newcastle = "#E74C3C"
)

# Propagate scheme to controls 
if (any(matched_full$scheme %in% c("Control", NA), na.rm = TRUE)) {
  pairs <- readRDS(here("data", "processed", "OA_matching_pairs_mixed.rds"))
  treated_col <- names(pairs)[str_detect(tolower(names(pairs)), "treat")][1]
  control_col <- names(pairs)[str_detect(tolower(names(pairs)), "control|donor")][1]
  treated_scheme <- matched_full %>%
    filter(treat_indicator == 1) %>%
    select(OA, scheme)
  ctrl_lookup <- pairs %>%
    select(treated_OA = all_of(treated_col), control_OA = all_of(control_col)) %>%
    left_join(treated_scheme, by = c("treated_OA" = "OA")) %>%
    select(OA = control_OA, scheme_from_treated = scheme) %>%
    distinct(OA, .keep_all = TRUE)
  matched_full <- matched_full %>%
    left_join(ctrl_lookup, by = "OA") %>%
    mutate(scheme = case_when(
      treat_indicator == 1        ~ scheme,
      !is.na(scheme_from_treated) ~ scheme_from_treated,
      TRUE                        ~ scheme
    )) %>%
    select(-scheme_from_treated)
}


# Flagged schemes — edit this vector if you want to add/remove
flagged <- c("Dundee", "Glasgow", "Sheffield", "Bristol", "Aberdeen", "Newcastle")

compute_smd <- function(data, var) {
  t <- data[[var]][data$treat_indicator == 1]
  c <- data[[var]][data$treat_indicator == 0]
  t <- t[!is.na(t)]; c <- c[!is.na(c)]
  if (length(t) < 2 || length(c) < 2) return(NA_real_)
  pooled_sd <- sqrt((var(t) + var(c)) / 2)
  if (pooled_sd == 0) return(0)
  (mean(t) - mean(c)) / pooled_sd
}

# =============================================================================
# TABLE 1 — MDIST SUMMARY PER FLAGGED SCHEME
# =============================================================================

cat("=== Mahalanobis distance summary ===\n")

mdist_summary <- matched_full %>%
  filter(scheme %in% flagged, treat_indicator == 1, !is.na(mdist)) %>%
  group_by(scheme) %>%
  summarise(
    n_treated    = n(),
    median_mdist = round(median(mdist),              2),
    mean_mdist   = round(mean(mdist),                2),
    p75_mdist    = round(quantile(mdist, 0.75),      2),
    p90_mdist    = round(quantile(mdist, 0.90),      2),
    p95_mdist    = round(quantile(mdist, 0.95),      2),
    max_mdist    = round(max(mdist),                 2),
    n_gt5        = sum(mdist > 5),
    n_gt10       = sum(mdist > 10),
    n_gt15       = sum(mdist > 15),
    n_gt20       = sum(mdist > 20),
    pct_gt5      = round(mean(mdist > 5)  * 100, 1),
    pct_gt10     = round(mean(mdist > 10) * 100, 1),
    pct_gt15     = round(mean(mdist > 15) * 100, 1),
    pct_gt20     = round(mean(mdist > 20) * 100, 1),
    .groups      = "drop"
  )

print(as.data.frame(mdist_summary))
write_csv(mdist_summary, file.path(outdir, "19_mdist_summary_table.csv"))
cat("  Saved: 19_mdist_summary_table.csv\n\n")

# =============================================================================
# TABLE 2 — CALIPER IMPACT: N DROPPED + MEAN |SMD| AT EACH THRESHOLD
# =============================================================================

cat("=== Caliper impact table ===\n")

calipers <- c(5, 10, 15, 20, 25)

caliper_impact <- map_df(flagged, function(s) {
  base_df <- matched_full %>% filter(scheme == s)
  map_df(calipers, function(cap) {
    calipered <- base_df %>%
      filter(treat_indicator == 0 |
               (treat_indicator == 1 & (is.na(mdist) | mdist <= cap)))
    n_dropped  <- sum(base_df$treat_indicator == 1) -
      sum(calipered$treat_indicator == 1)
    pct_dropped <- round(100 * n_dropped / sum(base_df$treat_indicator == 1), 1)
    mean_smd   <- mean(map_dbl(stage2_trends, ~abs(compute_smd(calipered, .x))),
                       na.rm = TRUE)
    max_smd    <- max(map_dbl(stage2_trends,  ~abs(compute_smd(calipered, .x))),
                      na.rm = TRUE)
    tibble(scheme      = s,
           caliper     = cap,
           n_treated   = sum(calipered$treat_indicator == 1),
           n_dropped   = n_dropped,
           pct_dropped = pct_dropped,
           mean_smd_after = round(mean_smd, 3),
           max_smd_after  = round(max_smd,  3))
  })
})

# Add uncalipered baseline (caliper = Inf)
baseline <- map_df(flagged, function(s) {
  df <- matched_full %>% filter(scheme == s)
  tibble(
    scheme      = s,
    caliper     = Inf,
    n_treated   = sum(df$treat_indicator == 1),
    n_dropped   = 0L,
    pct_dropped = 0,
    mean_smd_after = round(mean(map_dbl(stage2_trends,
                                        ~abs(compute_smd(df, .x))), na.rm = TRUE), 3),
    max_smd_after  = round(max(map_dbl(stage2_trends,
                                       ~abs(compute_smd(df, .x))), na.rm = TRUE), 3)
  )
})

caliper_impact_full <- bind_rows(baseline, caliper_impact) %>%
  arrange(scheme, caliper)

print(as.data.frame(caliper_impact_full))
write_csv(caliper_impact_full,
          file.path(outdir, "19_caliper_impact_table.csv"))
cat("  Saved: 19_caliper_impact_table.csv\n\n")

# =============================================================================
# PLOT 1 — ECDF OF MDIST, ALL FLAGGED SCHEMES
# =============================================================================

cat("=== Plot 1: ECDF ===\n")

mdist_data <- matched_full %>%
  filter(scheme %in% flagged, treat_indicator == 1, !is.na(mdist)) %>%
  mutate(scheme = factor(scheme, levels = flagged))

# Add reference lines at key calipers
ref_lines <- tibble(
  xint    = c(5, 10, 15, 20),
  ltype   = c("dotted", "dashed", "dashed", "dotdash"),
  colour  = c("#AAAAAA", "#E67E22", "#E74C3C", "#8B0000"),
  label   = c("5", "10", "15", "20")
)

p_ecdf <- ggplot(mdist_data, aes(x = mdist, colour = scheme)) +
  stat_ecdf(linewidth = 1.1, geom = "step") +
  geom_vline(data = ref_lines,
             aes(xintercept = xint),
             linetype  = ref_lines$ltype,
             colour    = ref_lines$colour,
             linewidth = 0.6) +
  annotate("text", x = ref_lines$xint + 0.4, y = 0.08,
           label = ref_lines$label, size = 3.2,
           colour = ref_lines$colour, hjust = 0) +
  scale_colour_manual(values = SCHEME_COLS, name = "Scheme") +
  coord_cartesian(xlim = c(0, 40)) +
  labs(
    title    = "ECDF of Stage 2 Mahalanobis distance — flagged schemes",
    subtitle = "Treated OAs only  |  steeper = more OAs with poor matches",
    x        = "Mahalanobis distance",
    y        = "Cumulative proportion of treated OAs",
    caption  = "Vertical lines: candidate caliper thresholds at 5, 10, 15, 20"
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p_ecdf, "fig_mdist_ecdf_problematic.png", width = 12, height = 7)

# =============================================================================
# PLOT 2 — CALIPER IMPACT: MEAN |SMD| vs % TREATED OAs DROPPED
#
# Each scheme gets a curve: x = % dropped, y = mean |SMD| after
# Ideal = bottom-left (low drop, low SMD)
# =============================================================================

cat("=== Plot 2: Caliper trade-off ===\n")

caliper_plot_data <- caliper_impact_full %>%
  filter(is.finite(caliper)) %>%
  mutate(
    scheme  = factor(scheme, levels = flagged),
    caliper = factor(caliper)
  )

# Also add the no-caliper point
no_caliper <- caliper_impact_full %>%
  filter(!is.finite(caliper)) %>%
  mutate(scheme = factor(scheme, levels = flagged))

p_tradeoff <- ggplot(caliper_impact_full %>%
                       mutate(caliper_lab = if_else(is.finite(caliper),
                                                    as.character(caliper), "None"),
                              scheme = factor(scheme, levels = flagged)),
                     aes(x = pct_dropped, y = mean_smd_after,
                         colour = scheme, group = scheme)) +
  geom_path(linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(shape = caliper_lab), size = 3.5) +
  geom_hline(yintercept = 0.10, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +
  scale_colour_manual(values = SCHEME_COLS, name = "Scheme") +
  scale_shape_manual(
    values = c("None" = 4, "5" = 15, "10" = 16, "15" = 17, "20" = 18, "25" = 19),
    name   = "Caliper"
  ) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Caliper threshold trade-off: balance vs sample retention",
    subtitle = "Each point = one caliper value  |  path moves left to right as caliper tightens",
    x        = "% treated OAs dropped",
    y        = "Mean |SMD| after caliper (across 9 trend variables)",
    caption  = "Dashed line = 0.10 balance threshold  |  X = no caliper applied"
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

save_fig(p_tradeoff, "fig_mdist_caliper_tradeoff.png", width = 12, height = 8)

# =============================================================================
# PLOT 3 — PER-VARIABLE SMD BEFORE/AFTER CALIPER (facet by scheme)
#
# Shows which specific trend variables are fixed by calipers
# Uses caliper = 15 as the illustrative threshold — change if needed
# =============================================================================

cat("=== Plot 3: Per-variable SMD at caliper = 15 ===\n")

ILLUSTRATIVE_CALIPER <- 15

var_labels_trend <- c(
  trend_car_KSI_pkm      = "Car KSI",
  trend_car_slight_pkm   = "Car slight",
  trend_cyc_KSI_pkm      = "Cycling KSI",
  trend_cyc_slight_pkm   = "Cycling slight",
  trend_ped_KSI_pkm      = "Ped KSI",
  trend_ped_slight_pkm   = "Ped slight",
  trend_other_KSI_pkm    = "Other KSI",
  trend_other_slight_pkm = "Other slight",
  trend_total_pkm        = "Total"
)

smd_comparison <- map_df(flagged, function(s) {
  full_df      <- matched_full %>% filter(scheme == s)
  calipered_df <- full_df %>%
    filter(treat_indicator == 0 |
             (treat_indicator == 1 & (is.na(mdist) | mdist <= ILLUSTRATIVE_CALIPER)))
  
  map_df(stage2_trends, function(v) {
    tibble(
      scheme     = s,
      variable   = v,
      label      = coalesce(var_labels_trend[v], v),
      smd_full   = abs(compute_smd(full_df,      v)),
      smd_caliper = abs(compute_smd(calipered_df, v)),
      n_full      = sum(full_df$treat_indicator      == 1),
      n_calipered = sum(calipered_df$treat_indicator == 1)
    )
  })
}) %>%
  mutate(scheme = factor(scheme, levels = flagged),
         label  = factor(label,  levels = rev(var_labels_trend)))

smd_long <- smd_comparison %>%
  pivot_longer(c(smd_full, smd_caliper),
               names_to = "version", values_to = "smd") %>%
  mutate(
    version = if_else(version == "smd_full",
                      "No caliper", paste0("Caliper \u2264 ", ILLUSTRATIVE_CALIPER)),
    version = factor(version,
                     levels = c("No caliper",
                                paste0("Caliper \u2264 ", ILLUSTRATIVE_CALIPER)))
  )

# Strip label with sample sizes
strip_labels <- smd_comparison %>%
  distinct(scheme, n_full, n_calipered) %>%
  mutate(strip = paste0(scheme, "\n(", n_calipered, "/", n_full, " retained)"))
strip_map <- setNames(strip_labels$strip, strip_labels$scheme)

smd_long <- smd_long %>%
  mutate(scheme_lab = strip_map[as.character(scheme)],
         scheme_lab = factor(scheme_lab, levels = strip_map))

p_var_smd <- ggplot(smd_long,
                    aes(x = smd, y = label,
                        colour = version, shape = version)) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "#EEEEEE", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "#DDDDDD", linewidth = 0.5) +
  geom_point(size = 3.5) +
  scale_colour_manual(
    values = setNames(
      c("#E74C3C", "#2ECC71"),
      c("No caliper", paste0("Caliper \u2264 ", ILLUSTRATIVE_CALIPER))
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = setNames(
      c(16, 17),
      c("No caliper", paste0("Caliper \u2264 ", ILLUSTRATIVE_CALIPER))
    ),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.06))) +
  facet_wrap(~scheme_lab, nrow = 1) +
  labs(
    title    = paste0("Per-variable |SMD|: no caliper vs caliper \u2264 ",
                      ILLUSTRATIVE_CALIPER, " — flagged schemes"),
    subtitle = "Strip label shows OAs retained after caliper  |  shared y-axis",
    x = "|SMD|", y = NULL,
    caption  = "Dashed = 0.10 threshold  \u25cf No caliper  \u25b2 Calipered"
  ) +
  theme_diag() +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 10),
    axis.text.x     = element_text(size = 9),
    strip.text      = element_text(size = 9, face = "bold"),
    panel.spacing.x = unit(0.8, "lines")
  )

save_fig(p_var_smd, "fig_mdist_caliper_smd_comparison.png",
         width = max(14, length(flagged) * 2.8), height = 8)

# =============================================================================
# PLOT 4 — MDIST VS WORST TREND VARIABLE (scatter, facet by scheme)
# Identifies structurally unusual OAs driving imbalance
# =============================================================================

cat("=== Plot 4: Distance vs trend scatter ===\n")

# Use trend_total_pkm as x — most consistently problematic across flagged schemes
scatter_data <- matched_full %>%
  filter(scheme %in% flagged, treat_indicator == 1, !is.na(mdist)) %>%
  mutate(
    scheme      = factor(scheme, levels = flagged),
    above_15    = mdist > ILLUSTRATIVE_CALIPER,
    label_point = if_else(mdist > 20, OA, NA_character_)
  )

p_scatter <- ggplot(scatter_data,
                    aes(x = trend_total_pkm, y = mdist,
                        colour = above_15)) +
  geom_hline(yintercept = ILLUSTRATIVE_CALIPER, linetype = "dashed",
             colour = "#E74C3C", linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dotted",
             colour = "#E67E22", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 2) +
  ggrepel::geom_text_repel(aes(label = label_point),
                           size = 2.5, max.overlaps = 10,
                           colour = "#333333") +
  scale_colour_manual(
    values = c("FALSE" = "#2E6FAB", "TRUE" = "#E74C3C"),
    labels = c(paste0("mdist \u2264 ", ILLUSTRATIVE_CALIPER),
               paste0("mdist > ",  ILLUSTRATIVE_CALIPER)),
    name   = NULL
  ) +
  facet_wrap(~scheme, nrow = 2, scales = "free_x") +
  labs(
    title    = "Stage 2 Mahalanobis distance vs total injury trend — flagged schemes",
    subtitle = "OAs with mdist > 20 labelled by OA code",
    x        = "Pre-treatment total injury trend (slope)",
    y        = "Mahalanobis distance",
    caption  = paste0("Red dashed = caliper at ", ILLUSTRATIVE_CALIPER,
                      "  |  Orange dotted = caliper at 10")
  ) +
  theme_diag() +
  theme(legend.position = "bottom",
        strip.text      = element_text(size = 11))

save_fig(p_scatter, "fig_mdist_scatter_flagged.png", width = 14, height = 10)

# =============================================================================
# CONSOLE SUMMARY
# =============================================================================

cat("\n================================================================\n")
cat("DISTANCE DIAGNOSTIC SUMMARY — FLAGGED SCHEMES\n")
cat("================================================================\n\n")

cat("At caliper = 15:\n")
caliper_impact_full %>%
  filter(caliper %in% c(Inf, 15)) %>%
  mutate(caliper = if_else(is.infinite(caliper), "None", as.character(caliper))) %>%
  select(scheme, caliper, n_treated, pct_dropped, mean_smd_after, max_smd_after) %>%
  arrange(scheme, caliper) %>%
  print(n = Inf)

cat("\nOutputs saved to:", outdir, "\n")
cat("  19_mdist_summary_table.csv\n")
cat("  19_caliper_impact_table.csv\n")
cat("  fig_mdist_ecdf_problematic.png\n")
cat("  fig_mdist_caliper_tradeoff.png\n")
cat("  fig_mdist_caliper_smd_comparison.png\n")
cat("  fig_mdist_scatter_flagged.png\n")
cat("================================================================\n")

# =============================================================================
# PLOT 5 — PARALLEL TRENDS BY SCHEME
#
# Panel A: Semi-annual pre-treatment time series (treated vs matched control),
#          faceted by scheme, ordered by mean |SMD| across 9 trend variables.
# Panel B: Ranking bar chart of mean |SMD| with quality colour bands.
#
# OUTPUTS:
#   fig_parallel_trends_timeseries.png
#   fig_parallel_trends_ranking.png
# =============================================================================

cat("\n================================================================\n")
cat("PARALLEL TRENDS DIAGNOSTIC BY SCHEME\n")
cat("================================================================\n\n")

library(zoo)

# --- 1. Build pre-treatment OA-quarter panel for matched sample ---------------

oa_q <- readRDS(here("data", "processed", "OA_injuries_quarterly.rds")) %>%
  mutate(quarter_year = as.yearqtr(quarter_year))

scheme_starts <- readRDS(here("data", "processed", "roads_caz_props.rds")) %>%
  distinct(scheme, caz_start_q) %>%
  filter(!is.na(scheme))

oa_panel <- oa_q %>%
  inner_join(
    matched_full %>%
      select(OA, treat_indicator, weights, scheme, country),
    by = "OA"
  ) %>%
  left_join(scheme_starts, by = "scheme")

oa_pre <- oa_panel %>%
  filter(quarter_year < caz_start_q)

cat("Pre-treatment OA-quarter rows:", nrow(oa_pre), "\n\n")

# --- 2. Ranking metric: mean |SMD| across 9 trend variables ------------------

trend_rank <- smd_all %>%
  filter(variable %in% stage2_trends) %>%
  group_by(scheme) %>%
  summarise(
    mean_abs_smd = mean(abs(smd_after), na.rm = TRUE),
    max_abs_smd  = max(abs(smd_after), na.rm = TRUE),
    n_above_010  = sum(abs(smd_after) >= 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    matched_full %>%
      filter(treat_indicator == 1) %>%
      distinct(scheme, country),
    by = "scheme"
  ) %>%
  arrange(mean_abs_smd) %>%
  mutate(rank = row_number())

cat("Parallel trends ranking (by mean |SMD| of trend variables):\n")
print(as.data.frame(trend_rank), row.names = FALSE)
cat("\n")

# --- 3. Semi-annual aggregation for time-series plot --------------------------

oa_pre_semi <- oa_pre %>%
  mutate(
    half_year = paste0(
      lubridate::year(as.Date(quarter_year)),
      if_else(as.numeric(format(as.Date(quarter_year), "%m")) <= 6, " H1", " H2")
    )
  ) %>%
  group_by(scheme, half_year, treat_indicator) %>%
  summarise(
    wtd_mean_inj = weighted.mean(total_injuries, w = weights, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = if_else(treat_indicator == 1, "Treated", "Matched control"),
    half_date = as.Date(paste0(
      str_extract(half_year, "\\d{4}"), "-",
      if_else(str_detect(half_year, "H1"), "03", "09"), "-01"
    ))
  )

facet_info <- trend_rank %>%
  select(scheme, rank, country, mean_abs_smd) %>%
  mutate(
    facet_label = paste0(rank, ". ", scheme, " (", country, ")"),
    facet_label = fct_reorder(facet_label, rank)
  )

plot_semi <- oa_pre_semi %>%
  left_join(facet_info, by = "scheme") %>%
  left_join(scheme_starts, by = "scheme") %>%
  mutate(group = factor(group, levels = c("Treated", "Matched control")))

vline_semi <- facet_info %>%
  left_join(scheme_starts, by = "scheme") %>%
  mutate(start_date = as.Date(caz_start_q))

# --- 4. Figure A: pre-treatment time series -----------------------------------

p_pt_timeseries <- ggplot(
  plot_semi,
  aes(x = half_date, y = wtd_mean_inj, colour = group, linetype = group)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(
    data = vline_semi,
    aes(xintercept = start_date),
    linetype = "dotted", colour = "#888888", linewidth = 0.5
  ) +
  scale_colour_manual(
    values = c("Treated" = "#D85A30", "Matched control" = "#2E6FAB")
  ) +
  scale_linetype_manual(
    values = c("Treated" = "solid", "Matched control" = "longdash")
  ) +
  facet_wrap(~facet_label, ncol = 3, scales = "free_y") +
  labs(
    title    = "Pre-treatment injury trends by scheme",
    subtitle = paste0(
      "Semi-annual weighted mean total injuries per OA | ",
      "ranked by mean |SMD| across 9 trend variables (1 = best)"
    ),
    x = NULL, y = "Weighted mean injuries per OA",
    colour = NULL, linetype = NULL,
    caption = "Dotted vertical line = scheme start date | Ranking based on post-matching mean |SMD| of 9 pre-treatment injury trend slopes"
  ) +
  theme_diag() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(size = 8, angle = 45, hjust = 1),
    strip.text      = element_text(size = 9, face = "bold"),
    panel.spacing   = unit(0.6, "lines")
  )

save_fig(p_pt_timeseries, "fig_parallel_trends_timeseries.png",
         width = 16, height = 14)

# --- 5. Figure B: ranking bar chart -------------------------------------------

quality_breaks <- c(0, 0.05, 0.10, 0.15, Inf)
quality_labels <- c("Excellent (<0.05)", "Acceptable (0.05\u20130.10)",
                     "Marginal (0.10\u20130.15)", "Poor (>0.15)")
quality_colours <- c(
  "Excellent (<0.05)"            = "#1a9850",
  "Acceptable (0.05\u20130.10)"  = "#91cf60",
  "Marginal (0.10\u20130.15)"    = "#fc8d59",
  "Poor (>0.15)"                 = "#d73027"
)

rank_bar_data <- trend_rank %>%
  mutate(
    quality = cut(mean_abs_smd, breaks = quality_breaks,
                  labels = quality_labels, right = FALSE, include.lowest = TRUE),
    scheme_lab = paste0(scheme, " (", country, ")"),
    scheme_lab = fct_reorder(scheme_lab, -mean_abs_smd)
  )

p_pt_ranking <- ggplot(
  rank_bar_data,
  aes(y = scheme_lab, x = mean_abs_smd, fill = quality)
) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_vline(xintercept = 0.10, linetype = "dashed",
             colour = "#555555", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean_abs_smd)),
            hjust = -0.15, size = 3.2, colour = "#333333", fontface = "bold") +
  scale_fill_manual(values = quality_colours, name = "Balance quality", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title    = "Parallel trends assumption: scheme ranking",
    subtitle = "Mean |SMD| across 9 pre-treatment injury trend variables (lower = better)",
    x = "Mean |SMD| after matching", y = NULL,
    caption  = "Dashed line = |SMD| = 0.10 conventional threshold"
  ) +
  theme_diag() +
  theme(
    legend.position    = "bottom",
    axis.text.y        = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_blank()
  )

save_fig(p_pt_ranking, "fig_parallel_trends_ranking.png",
         width = 10, height = 7)

cat("\nParallel trends outputs saved:\n")
cat("  fig_parallel_trends_timeseries.png\n")
cat("  fig_parallel_trends_ranking.png\n")
cat("================================================================\n")

