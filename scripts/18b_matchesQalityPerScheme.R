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

# =============================================================================
# LOAD DATA
# =============================================================================

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
# COUNTRY-LEVEL LOVE PLOTS (SEPARATE GRIDS)
# =============================================================================

cat("=== Country-level love plots (separate grids) ===\n")

# join country info safely
smd_all_country <- smd_all %>%
  left_join(
    matched_full %>% distinct(scheme, country),
    by = "scheme"
  )

make_country_love_plot <- function(df, country_name) {
  
  plot_data <- df %>%
    filter(country == country_name) %>%
    filter(!is.na(smd_before), !is.na(smd_after)) %>%
    mutate(
      label = factor(label, levels = unique(label[order(abs(smd_before))]))
    ) %>%
    pivot_longer(
      cols = c(smd_before, smd_after),
      names_to = "timing",
      values_to = "smd"
    ) %>%
    mutate(
      timing = if_else(timing == "smd_before",
                       "Before matching",
                       "After matching"),
      timing = factor(timing,
                      levels = c("Before matching", "After matching"))
    )
  
  ggplot(plot_data,
         aes(x = abs(smd), y = label,
             colour = timing, shape = timing)) +
    
    geom_vline(xintercept = 0.10, linetype = "dashed",
               colour = "#999999", linewidth = 0.5) +
    
    geom_line(aes(group = label),
              colour = "#DDDDDD", linewidth = 0.4) +
    
    geom_point(size = 2.8) +
    
    facet_wrap(~scheme, scales = "free_y", ncol = 3) +
    
    scale_colour_manual(values = c(
      "Before matching" = "#E74C3C",
      "After matching"  = "#2ECC71"
    )) +
    
    scale_shape_manual(values = c(16, 17)) +
    
    labs(
      title = paste0("Injury trend balance — ", country_name),
      subtitle = "Each panel = one scheme",
      x = "|SMD|",
      y = NULL,
      colour = NULL,
      shape = NULL,
      caption = "Dashed line = 0.10 threshold"
    ) +
    
    theme_diag() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold")
    )
}

# =========================
# RUN FOR EACH COUNTRY
# =========================

countries <- c("England", "Scotland")

for (ctry in countries) {
  
  cat("  →", ctry, "\n")
  
  p <- make_country_love_plot(smd_all_country, ctry)
  
  fname <- paste0(
    "fig_love_", tolower(ctry), "_schemes_grid.png"
  )
  
  save_fig(p, fname, width = 14, height = 9)
}
# =============================================================================
# PLOT 2 — SMD SUMMARY HEATMAP (schemes × variable groups)
#
# Rows    = variable groups (5)
# Columns = schemes (11)
# Fill    = mean |SMD| after matching
# Outline = mean |SMD| before (shown as text, top-right corner of cell)
# =============================================================================

cat("=== SMD summary heatmap ===\n")

heatmap_data <- smd_all %>%
  filter(variable %in% stage2_trends) %>%
  mutate(
    label = coalesce(var_labels[variable], variable),
    scheme = factor(scheme, levels = schemes)
  ) %>%
  group_by(scheme, label) %>%
  summarise(
    mean_smd_after  = mean(abs(smd_after),  na.rm = TRUE),
    mean_smd_before = mean(abs(smd_before), na.rm = TRUE),
    .groups = "drop"
  )

# Colour breaks for fill
smd_breaks  <- c(0, 0.05, 0.10, 0.15, 0.20, Inf)
smd_labels  <- c("<0.05", "0.05–0.10", "0.10–0.15", "0.15–0.20", ">0.20")
smd_colours <- c("#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027")

heatmap_data <- heatmap_data %>%
  mutate(
    smd_band   = cut(mean_smd_after, breaks = smd_breaks,
                     labels = smd_labels, right = FALSE, include.lowest = TRUE),
    cell_label = sprintf("%.3f", mean_smd_after)
  )

p_heatmap <- ggplot(heatmap_data,
                    aes(x = scheme, y = fct_rev(label), fill = smd_band)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label), size = 2.8, colour = "#222222") +
  scale_fill_manual(
    values = setNames(smd_colours, smd_labels),
    name   = "|SMD| after matching",
    drop   = FALSE
  ) +
  scale_x_discrete(guide = guide_axis(angle = 35)) +
  labs(
    title    = "Injury trend balance across schemes — |SMD| after matching",
    subtitle = "Each cell = one trend variable for one scheme",
    x = NULL, y = NULL,
    caption  = "Green < 0.05 · yellow 0.05–0.10 · orange 0.10–0.15 · red > 0.20"
  ) +
  theme_diag(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.text.x     = element_text(size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm")
  )

save_fig(p_heatmap, "fig_scheme_smd_heatmap.png", width = 16, height = 7)
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

