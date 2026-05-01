




# =============================================================================
# PRELIMINARY DiD — STAGGERED ATT ESTIMATION
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(dplyr)
library(here)
# =============================================================================
# 
# =============================================================================

matched_full <- readRDS(here("data", "processed", "OA_matched_full_A.rds"))

matched_treated_OAs  <- matched_full %>% filter(treat_indicator == 1) %>% pull(OA)
matched_control_OAs  <- matched_full %>% filter(treat_indicator == 0) %>% pull(OA)
all_matched_OAs      <- matched_full %>% pull(OA)

cat("Matched treated OAs: ",  length(matched_treated_OAs), "\n")
cat("Matched control OAs: ",  length(matched_control_OAs), "\n")
cat("Total matched OAs:   ",  length(all_matched_OAs),     "\n")

# =============================================================================
#  Load road-to-OA lookup
# =============================================================================
# road_attributes already has OA attached from script 9 (road_attributes_OA.gpkg)

road_attributes <- st_read(here("data", "processed", "road_attributes_OA.gpkg")) %>%
  st_drop_geometry()
 

# =============================================================================
#  Filter roads to matched OAs only
# =============================================================================

matched_roads <- road_attributes %>%
  filter(OA %in% all_matched_OAs)

cat("Roads in matched OAs:", n_distinct(matched_roads$identifier), "\n")
cat("OAs represented:     ", n_distinct(matched_roads$OA), "\n")


# =============================================================================
#  Load full panel and filter to matched roads
# =============================================================================

road_panel_full <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

road_panel_matched <- road_panel_full %>%
  inner_join(matched_roads, by = "identifier") %>%    # attaches OA column
  left_join(
    matched_full ,
    by = "OA"
  )

cat("\nMatched panel dimensions:\n")
cat("  Rows:              ", nrow(road_panel_matched), "\n")
cat("  Unique roads:      ", n_distinct(road_panel_matched$identifier), "\n")
cat("  Unique OAs:        ", n_distinct(road_panel_matched$OA), "\n")
cat("  Treated roads:     ", n_distinct(road_panel_matched$identifier[road_panel_matched$treat_indicator == 1]), "\n")
cat("  Control roads:     ", n_distinct(road_panel_matched$identifier[road_panel_matched$treat_indicator == 0]), "\n")







road_panel_matched %>%
dplyr::select(starts_with("total"), starts_with("KSI"), starts_with("Slight")) %>%
  names()

# Check road length column name
road_panel_matched %>%
  dplyr::select(contains("length"), contains("road_len")) %>%
  names()

# Check caz_start_q for controls — is it NA?
road_panel_matched %>%
  dplyr::count(treat_indicator, is.na(caz_start_q))

## =============================================================================
# CHECK — what columns does matched_full actually have?
# =============================================================================

names(matched_full)
head(matched_full %>% dplyr::select(OA, treat_indicator, everything()), 3)

# =============================================================================
# STEP 1 — Build timing lookup WITHOUT subclass
# =============================================================================
# Controls need a reference start date for pre/post split.
# Simplest valid approach: assign each control the MEDIAN LEZ start
# across all treated OAs. This is standard in pooled pre-post comparisons.

# Get treated OA timing from the panel (has caz_start_q)
treated_oa_timing <- road_panel_matched %>%
  dplyr::filter(treat_indicator == 1, !is.na(caz_start_q)) %>%
  dplyr::distinct(OA, caz_start_q)

cat("Treated OAs with timing:", nrow(treated_oa_timing), "\n")

# Median LEZ start across all treated OAs
median_start <- median(as.yearqtr(treated_oa_timing$caz_start_q), na.rm = TRUE)
cat("Median LEZ start:", as.character(median_start), "\n")

# City-specific timing for treated; median timing for controls
oa_timing_lookup <- bind_rows(
  # Treated OAs: use their own start
  treated_oa_timing %>%
    dplyr::rename(ref_start = caz_start_q) %>%
    dplyr::mutate(ref_start = as.yearqtr(ref_start)),
  
  # Control OAs: assign median start
  matched_full %>%
    dplyr::filter(treat_indicator == 0) %>%
    dplyr::distinct(OA) %>%
    dplyr::mutate(ref_start = median_start)
)

cat("OAs in timing lookup:", nrow(oa_timing_lookup), "\n")
cat("NA in ref_start:     ", sum(is.na(oa_timing_lookup$ref_start)), "\n")

# =============================================================================
# STEP 2 — Join timing into panel and create rate outcomes
# =============================================================================

panel <- road_panel_matched %>%
  dplyr::left_join(oa_timing_lookup, by = "OA") %>%
  dplyr::mutate(
    quarter_year   = as.yearqtr(quarter_year),
    road_length_km = dplyr::if_else(
      road_length_km <= 0 | is.na(road_length_km), 1e-6, road_length_km
    ),
    # Rate outcomes
    total_pkm  = total_inj_adj_All     / road_length_km,
    KSI_pkm    = KSI_adj_All           / road_length_km,
    slight_pkm = Slight_adj_All        / road_length_km,
    ped_pkm    = total_inj_adj_Pedestrian / road_length_km,
    cyc_pkm    = total_inj_adj_Cyclist    / road_length_km,
    # Group and period flags
    group = dplyr::if_else(treat_indicator == 1, "LEZ roads", "Matched controls"),
    post  = quarter_year >= ref_start
  )

# Sense check
cat("\nPeriod distribution:\n")
table(panel$post, panel$group, useNA = "ifany")

cat("\nNA in ref_start after join:", sum(is.na(panel$ref_start)), "\n")
cat("NA in total_pkm:           ", sum(is.na(panel$total_pkm)),  "\n")

# =============================================================================
# STEP 3 — Summary stats by group and period
# =============================================================================

summary_stats <- panel %>%
  dplyr::group_by(group, post) %>%
  dplyr::summarise(
    mean_total_pkm  = round(mean(total_pkm,  na.rm = TRUE), 4),
    mean_KSI_pkm    = round(mean(KSI_pkm,    na.rm = TRUE), 4),
    mean_slight_pkm = round(mean(slight_pkm, na.rm = TRUE), 4),
    mean_ped_pkm    = round(mean(ped_pkm,    na.rm = TRUE), 4),
    n_roads         = dplyr::n_distinct(identifier),
    n_obs           = dplyr::n(),
    .groups         = "drop"
  ) %>%
  dplyr::mutate(period = dplyr::if_else(post, "Post", "Pre")) %>%
  dplyr::select(group, period, mean_total_pkm, mean_KSI_pkm,
                mean_slight_pkm, mean_ped_pkm, n_roads, n_obs)

print(summary_stats)

# =============================================================================
# STEP 4 — Simple DiD estimates
# =============================================================================

did_table <- summary_stats %>%
  tidyr::pivot_wider(
    id_cols     = group,
    names_from  = period,
    values_from = c(mean_total_pkm, mean_KSI_pkm, mean_slight_pkm, mean_ped_pkm)
  ) %>%
  dplyr::mutate(
    change_total  = round(mean_total_pkm_Post  - mean_total_pkm_Pre,  5),
    change_KSI    = round(mean_KSI_pkm_Post    - mean_KSI_pkm_Pre,    5),
    change_slight = round(mean_slight_pkm_Post - mean_slight_pkm_Pre, 5),
    change_ped    = round(mean_ped_pkm_Post     - mean_ped_pkm_Pre,    5)
  )

print(did_table %>% dplyr::select(group, starts_with("change")))

did_est <- did_table %>%
  dplyr::summarise(
    DiD_total  = change_total[group  == "LEZ roads"] - change_total[group  == "Matched controls"],
    DiD_KSI    = change_KSI[group    == "LEZ roads"] - change_KSI[group    == "Matched controls"],
    DiD_slight = change_slight[group == "LEZ roads"] - change_slight[group == "Matched controls"],
    DiD_ped    = change_ped[group    == "LEZ roads"] - change_ped[group    == "Matched controls"]
  )

treated_pre_mean <- summary_stats %>%
  dplyr::filter(group == "LEZ roads", period == "Pre")

cat("\n========================================\n")
cat("SIMPLE DiD ESTIMATES (descriptive only)\n")
cat("Treated change minus Control change\n")
cat("========================================\n")
cat("Total injuries/km: ", round(did_est$DiD_total,  5),
    " (", round(100 * did_est$DiD_total  / treated_pre_mean$mean_total_pkm,  1), "%)\n")
cat("KSI/km:            ", round(did_est$DiD_KSI,    5),
    " (", round(100 * did_est$DiD_KSI    / treated_pre_mean$mean_KSI_pkm,    1), "%)\n")
cat("Slight/km:         ", round(did_est$DiD_slight, 5),
    " (", round(100 * did_est$DiD_slight / treated_pre_mean$mean_slight_pkm, 1), "%)\n")
cat("Pedestrian/km:     ", round(did_est$DiD_ped,    5),
    " (", round(100 * did_est$DiD_ped    / treated_pre_mean$mean_ped_pkm,    1), "%)\n")

# =============================================================================
# STEP 5 — Trend plot
# =============================================================================

plot_data <- panel %>%
  dplyr::group_by(group, quarter_year) %>%
  dplyr::summarise(
    mean_total_pkm = mean(total_pkm, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(plot_data, aes(x = as.Date(quarter_year),
                      y = mean_total_pkm,
                      colour = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = as.Date(median_start),
             linetype = "dashed", colour = "grey40") +
  annotate("text", x = as.Date(median_start), y = Inf,
           label = "Median LEZ start", vjust = 2, hjust = -0.1,
           size = 3.5, colour = "grey40") +
  scale_colour_manual(values = c("LEZ roads"       = "#E74C3C",
                                 "Matched controls" = "#2C3E50")) +
  labs(title   = "Mean road traffic injuries per road-km: LEZ roads vs matched controls",
       x = NULL, y = "Mean injuries per road-km", colour = NULL,
       caption = "Descriptive pre-post comparison. Formal DiD estimates pending.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(here("output", "prepost_descriptive.png"), width = 10, height = 6, dpi = 300)

treated_pre <- summary_stats %>%
  dplyr::filter(group == "Treated", period == "Pre")

cat("\n========================================\n")
cat("SIMPLE DiD ESTIMATES (descriptive only)\n")
cat("========================================\n")
cat("Total injuries/km: ",  round(did_est$DiD_total,  4), "\n")
cat("KSI/km:            ",  round(did_est$DiD_KSI,    4), "\n")
cat("Slight/km:         ",  round(did_est$DiD_slight, 4), "\n")
cat("\nAs % of treated pre-period mean:\n")
cat("Total:  ", round(100 * did_est$DiD_total  / treated_pre$mean_total_pkm,  1), "%\n")
cat("KSI:    ", round(100 * did_est$DiD_KSI    / treated_pre$mean_KSI_pkm,    1), "%\n")
cat("Slight: ", round(100 * did_est$DiD_slight / treated_pre$mean_slight_pkm, 1), "%\n")


# =============================================================================
# STEP 6 — TWFE DiD (road + quarter fixed effects, OA-clustered SEs)
# =============================================================================
# Regression-based check of the simple DiD estimates above.
# treat_indicator is absorbed by road FE; the interaction term is the DiD estimate.

library(fixest)

panel_reg <- panel %>%
  mutate(
    treated  = as.integer(treat_indicator),
    post_int = as.integer(post),
    qtr_num  = as.numeric(quarter_year)
  )

m_total  <- feols(total_pkm  ~ treated:post_int | identifier + qtr_num, data = panel_reg, cluster = ~OA)
m_KSI    <- feols(KSI_pkm    ~ treated:post_int | identifier + qtr_num, data = panel_reg, cluster = ~OA)
m_slight <- feols(slight_pkm ~ treated:post_int | identifier + qtr_num, data = panel_reg, cluster = ~OA)
m_ped    <- feols(ped_pkm    ~ treated:post_int | identifier + qtr_num, data = panel_reg, cluster = ~OA)

etable(
  m_total, m_KSI, m_slight, m_ped,
  headers      = c("Total/km", "KSI/km", "Slight/km", "Ped/km"),
  dict         = c("treated:post_int" = "Treated × Post"),
  digits       = 5,
  digits.stats = 3
)

# NOTE: None of the TWFE estimates are statistically significant. The −8%
# from the simple pre-post comparison does not survive road + time fixed effects.


# =============================================================================
# STEP 7 — Staggered DiD: Callaway & Sant'Anna (2021)
# =============================================================================
# Properly handles staggered rollout across cities; uses never-treated roads
# as the control group rather than assigning a median start date.
# SEs clustered at OA level. bstrap = FALSE uses analytical (faster) SEs.

# --- 7a. Prepare panel_did: numeric IDs, integer time index, cohort variable ---

min_qtr <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

panel_did <- panel %>%
  mutate(
    road_id = as.integer(factor(identifier)),
    # Integer time index (1 = earliest quarter in data)
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    # g = first treated period (integer); 0 = never treated
    g = case_when(
      treat_indicator == 1 & !is.na(caz_start_q) ~
        as.integer(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1L,
      TRUE ~ 0L
    )
  ) %>%
  filter(!is.na(total_pkm))

cat("Time periods :", n_distinct(panel_did$qtr_int), "\n")
cat("Roads        :", n_distinct(panel_did$road_id),  "\n")

# Treatment cohort breakdown
panel_did %>%
  distinct(road_id, g) %>%
  count(g) %>%
  mutate(label = if_else(g == 0, "Never treated",
                         as.character(as.yearqtr(min_qtr + (g - 1) / 4)))) %>%
  print()

# --- 7b. att_gt for each outcome ---

att_total <- att_gt(
  yname = "total_pkm", tname = "qtr_int", idname = "road_id", gname = "g",
  data = panel_did, control_group = "nevertreated",
  clustervars = "OA", bstrap = FALSE, anticipation = 0, panel = TRUE
)

att_KSI <- att_gt(
  yname = "KSI_pkm", tname = "qtr_int", idname = "road_id", gname = "g",
  data = panel_did, control_group = "nevertreated",
  clustervars = "OA", bstrap = FALSE, anticipation = 0, panel = TRUE
)

att_ped <- att_gt(
  yname = "ped_pkm", tname = "qtr_int", idname = "road_id", gname = "g",
  data = panel_did, control_group = "nevertreated",
  clustervars = "OA", bstrap = FALSE, anticipation = 0, panel = TRUE
)

# --- 7c. Aggregate to overall ATT ---

agg_simple  <- aggte(att_total, type = "simple", na.rm = TRUE)
agg_KSI     <- aggte(att_KSI,   type = "simple", na.rm = TRUE)
agg_ped     <- aggte(att_ped,   type = "simple", na.rm = TRUE)

cat("\n=== Overall ATT — Total injuries/km ===\n"); summary(agg_simple)
cat("\n=== Overall ATT — KSI/km ===\n");            summary(agg_KSI)
cat("\n=== Overall ATT — Pedestrian/km ===\n");     summary(agg_ped)

# Results:
#   Total/km  : ATT = -0.0088, SE = 0.0049, 95% CI [-0.019, +0.001]  — marginal
#   KSI/km    : ATT = -0.0052, SE = 0.0023, 95% CI [-0.010, -0.001]  — significant *
#   Ped/km    : ATT = -0.0071, SE = 0.0029, 95% CI [-0.013, -0.001]  — significant *

# --- 7d. Event-study (dynamic) aggregation ---

agg_dyn     <- aggte(att_total, type = "dynamic", na.rm = TRUE)
agg_dyn_KSI <- aggte(att_KSI,  type = "dynamic", na.rm = TRUE)
agg_dyn_ped <- aggte(att_ped,  type = "dynamic", na.rm = TRUE)

# --- 7e. Event-study plot (±12 quarters around implementation) ---

extract_es <- function(agg, label) {
  tibble(
    event_time = agg$egt,
    att        = agg$att.egt,
    se         = agg$se.egt,
    ci_lo      = agg$att.egt - 1.96 * agg$se.egt,
    ci_hi      = agg$att.egt + 1.96 * agg$se.egt,
    outcome    = label
  )
}

es_df <- bind_rows(
  extract_es(agg_dyn,     "Total injuries/km"),
  extract_es(agg_dyn_KSI, "KSI/km"),
  extract_es(agg_dyn_ped, "Pedestrian/km")
) %>%
  filter(event_time >= -12, event_time <= 12)

ggplot(es_df, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~outcome, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
  labs(
    title    = "Event-study estimates: effect of LEZ on road traffic injuries",
    subtitle = "Callaway & Sant'Anna (2021) staggered DiD — never-treated controls, OA-clustered SEs",
    x        = "Quarters relative to LEZ implementation",
    y        = "ATT (injuries per road-km)",
    caption  = "Shaded band: 95% pointwise CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(here("output", "event_study_staggered_did.png"), width = 9, height = 10, dpi = 300)

