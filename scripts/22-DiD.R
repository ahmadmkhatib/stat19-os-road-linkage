# =============================================================================
# PRELIMINARY DiD — STAGGERED ATT ESTIMATION
# =============================================================================
#
# Reads the matched road × quarter panel from script 20 and runs:
#   - Simple pre/post descriptive DiD
#   - TWFE DiD (road + quarter FE, OA-clustered SEs)
#   - Callaway & Sant'Anna (2021) staggered DiD with never-treated controls
#
# Input:
#   road_panel_matched.parquet  (from script 20)
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(fixest)

# =============================================================================
# Load matched panel from script 20
# =============================================================================

road_panel_matched <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched.parquet")
)

cat("Rows:", nrow(road_panel_matched),
    "| Panel units:", n_distinct(road_panel_matched$panel_id),
    "| OAs:", n_distinct(road_panel_matched$OA),
    "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n")

road_panel_matched <- road_panel_matched %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

# =============================================================================
# Create rate outcomes (injuries per road-km)
# =============================================================================

panel <- road_panel_matched %>%
  mutate(
    road_length_km = length / 1000,
    road_length_km = if_else(
      road_length_km <= 0 | is.na(road_length_km), 1e-6, road_length_km
    ),
    total_pkm  = total_inj_adj_All     / road_length_km,
    KSI_pkm    = KSI_adj_All           / road_length_km,
    slight_pkm = Slight_adj_All        / road_length_km,
    ped_pkm    = total_inj_adj_Pedestrian / road_length_km,
    cyc_pkm    = total_inj_adj_Cyclist    / road_length_km,
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  )

# =============================================================================
# STEP 1 — Reference timing for pre/post split
# =============================================================================
# Each road uses its scheme's caz_start_q.
# Controls inherit scheme from their matched treated OA (done in script 20),
# so they get the scheme-specific start date rather than a pooled median.

scheme_timing <- panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

cat("Scheme start dates:\n")
print(scheme_timing %>% arrange(caz_start_q))

panel <- panel %>%
  left_join(
    scheme_timing %>% rename(ref_start = caz_start_q),
    by = "scheme"
  ) %>%
  mutate(
    post_flag = as.integer(quarter_year >= ref_start)
  )

cat("\nPeriod distribution:\n")
table(panel$post_flag, panel$group, useNA = "ifany")

cat("\nNA in ref_start:", sum(is.na(panel$ref_start)), "\n")

# =============================================================================
# STEP 2 — Summary stats by group and period
# =============================================================================

summary_stats <- panel %>%
  group_by(group, post_flag) %>%
  summarise(
    mean_total_pkm  = round(mean(total_pkm,  na.rm = TRUE), 4),
    mean_KSI_pkm    = round(mean(KSI_pkm,    na.rm = TRUE), 4),
    mean_slight_pkm = round(mean(slight_pkm, na.rm = TRUE), 4),
    mean_ped_pkm    = round(mean(ped_pkm,    na.rm = TRUE), 4),
    n_roads         = n_distinct(identifier),
    n_obs           = n(),
    .groups         = "drop"
  ) %>%
  mutate(period = if_else(post_flag == 1L, "Post", "Pre")) %>%
  select(group, period, mean_total_pkm, mean_KSI_pkm,
         mean_slight_pkm, mean_ped_pkm, n_roads, n_obs)

print(summary_stats)

# =============================================================================
# STEP 3 — Simple DiD estimates (descriptive)
# =============================================================================

did_table <- summary_stats %>%
  pivot_wider(
    id_cols     = group,
    names_from  = period,
    values_from = c(mean_total_pkm, mean_KSI_pkm, mean_slight_pkm, mean_ped_pkm)
  ) %>%
  mutate(
    change_total  = round(mean_total_pkm_Post  - mean_total_pkm_Pre,  5),
    change_KSI    = round(mean_KSI_pkm_Post    - mean_KSI_pkm_Pre,    5),
    change_slight = round(mean_slight_pkm_Post - mean_slight_pkm_Pre, 5),
    change_ped    = round(mean_ped_pkm_Post     - mean_ped_pkm_Pre,    5)
  )

print(did_table %>% select(group, starts_with("change")))

did_est <- did_table %>%
  summarise(
    DiD_total  = change_total[group  == "CAZ roads"] - change_total[group  == "Matched controls"],
    DiD_KSI    = change_KSI[group    == "CAZ roads"] - change_KSI[group    == "Matched controls"],
    DiD_slight = change_slight[group == "CAZ roads"] - change_slight[group == "Matched controls"],
    DiD_ped    = change_ped[group    == "CAZ roads"] - change_ped[group    == "Matched controls"]
  )

treated_pre_mean <- summary_stats %>%
  filter(group == "CAZ roads", period == "Pre")

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
# STEP 4 — Trend plot
# =============================================================================

plot_data <- panel %>%
  group_by(group, quarter_year) %>%
  summarise(
    mean_total_pkm = mean(total_pkm, na.rm = TRUE),
    .groups = "drop"
  )

# Median start for the reference line on the plot
median_start <- median(scheme_timing$caz_start_q, na.rm = TRUE)

ggplot(plot_data, aes(x = as.Date(quarter_year),
                      y = mean_total_pkm,
                      colour = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = as.Date(median_start),
             linetype = "dashed", colour = "grey40") +
  annotate("text", x = as.Date(median_start), y = Inf,
           label = "Median CAZ start", vjust = 2, hjust = -0.1,
           size = 3.5, colour = "grey40") +
  scale_colour_manual(values = c("CAZ roads"         = "#E74C3C",
                                 "Matched controls"  = "#2C3E50")) +
  labs(title   = "Mean road traffic injuries per road-km: CAZ roads vs matched controls",
       x = NULL, y = "Mean injuries per road-km", colour = NULL,
       caption = "Descriptive pre-post comparison. Formal DiD estimates pending.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(here("output", "prepost_descriptive.png"), width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 5 — TWFE DiD (road + quarter FE, OA-clustered SEs)
# =============================================================================

panel_reg <- panel %>%
  mutate(
    treated  = as.integer(treat_group),
    post_int = post_flag,
    qtr_num  = as.numeric(quarter_year)
  )

m_total  <- feols(total_pkm  ~ treated:post_int | panel_id + qtr_num, data = panel_reg, cluster = ~OA)
m_KSI    <- feols(KSI_pkm    ~ treated:post_int | panel_id + qtr_num, data = panel_reg, cluster = ~OA)
m_slight <- feols(slight_pkm ~ treated:post_int | panel_id + qtr_num, data = panel_reg, cluster = ~OA)
m_ped    <- feols(ped_pkm    ~ treated:post_int | panel_id + qtr_num, data = panel_reg, cluster = ~OA)

etable(
  m_total, m_KSI, m_slight, m_ped,
  headers      = c("Total/km", "KSI/km", "Slight/km", "Ped/km"),
  dict         = c("treated:post_int" = "Treated x Post"),
  digits       = 5,
  digits.stats = 3
)

# =============================================================================
# STEP 6 — Staggered DiD: Callaway & Sant'Anna (2021)
# =============================================================================
# Uses never-treated roads as the control group.
# SEs clustered at OA level. bstrap = FALSE uses analytical SEs.

# --- 6a. Prepare panel_did: numeric IDs, integer time index, cohort variable ---

min_qtr <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

panel_did <- panel %>%
  mutate(
    road_id = as.integer(factor(panel_id)),
    # Integer time index (1 = earliest quarter in data)
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    # g = first treated period (numeric); 0 = never treated
    # Must be numeric (not integer) — att_gt internally sets g = Inf for controls
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    )
  ) %>%
  filter(!is.na(total_pkm))

cat("\nStaggered DiD setup:\n")
cat("Time periods  :", n_distinct(panel_did$qtr_int), "\n")
cat("Panel units   :", n_distinct(panel_did$road_id),  "\n")

# Treatment cohort breakdown
panel_did %>%
  distinct(road_id, g) %>%
  count(g) %>%
  mutate(label = if_else(g == 0, "Never treated",
                         as.character(as.yearqtr(min_qtr + (g - 1) / 4)))) %>%
  print()

# --- 6b. att_gt for each outcome ---

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

# --- 6c. Aggregate to overall ATT ---

agg_simple  <- aggte(att_total, type = "simple", na.rm = TRUE)
agg_KSI     <- aggte(att_KSI,   type = "simple", na.rm = TRUE)
agg_ped     <- aggte(att_ped,   type = "simple", na.rm = TRUE)

cat("\n=== Overall ATT — Total injuries/km ===\n"); summary(agg_simple)
cat("\n=== Overall ATT — KSI/km ===\n");            summary(agg_KSI)
cat("\n=== Overall ATT — Pedestrian/km ===\n");     summary(agg_ped)

# --- 6d. Event-study (dynamic) aggregation ---

agg_dyn     <- aggte(att_total, type = "dynamic", na.rm = TRUE)
agg_dyn_KSI <- aggte(att_KSI,  type = "dynamic", na.rm = TRUE)
agg_dyn_ped <- aggte(att_ped,  type = "dynamic", na.rm = TRUE)

# --- 6e. Event-study plot (+-12 quarters around implementation) ---

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
    title    = "Event-study estimates: effect of CAZ on road traffic injuries",
    subtitle = "Callaway & Sant'Anna (2021) staggered DiD — never-treated controls, OA-clustered SEs",
    x        = "Quarters relative to CAZ implementation",
    y        = "ATT (injuries per road-km)",
    caption  = "Shaded band: 95% pointwise CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(here("output", "event_study_staggered_did.png"), width = 9, height = 10, dpi = 300)

# =============================================================================
# STEP 7 — Staggered DiD by severity × mode (8 estimates)
# =============================================================================
# KSI and Slight, each split by Pedestrian / Car-Van / Cyclist / Other

panel_did8 <- panel %>%
  mutate(
    KSI_ped_pkm    = KSI_adj_Pedestrian    / road_length_km,
    KSI_car_pkm    = KSI_adj_Car.Van       / road_length_km,
    KSI_cyc_pkm    = KSI_adj_Cyclist       / road_length_km,
    KSI_oth_pkm    = KSI_adj_Other         / road_length_km,
    Slight_ped_pkm = Slight_adj_Pedestrian / road_length_km,
    Slight_car_pkm = Slight_adj_Car.Van    / road_length_km,
    Slight_cyc_pkm = Slight_adj_Cyclist    / road_length_km,
    Slight_oth_pkm = Slight_adj_Other      / road_length_km,
    road_id = as.integer(factor(panel_id)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    )
  ) %>%
  filter(
    !is.na(KSI_ped_pkm),  !is.na(KSI_car_pkm),
    !is.na(KSI_cyc_pkm),  !is.na(KSI_oth_pkm),
    !is.na(Slight_ped_pkm), !is.na(Slight_car_pkm),
    !is.na(Slight_cyc_pkm), !is.na(Slight_oth_pkm)
  )

outcomes_8 <- c(
  "KSI_ped_pkm", "KSI_car_pkm", "KSI_cyc_pkm", "KSI_oth_pkm",
  "Slight_ped_pkm", "Slight_car_pkm", "Slight_cyc_pkm", "Slight_oth_pkm"
)
labels_8 <- c(
  "KSI — Pedestrian", "KSI — Car/Van", "KSI — Cyclist", "KSI — Other",
  "Slight — Pedestrian", "Slight — Car/Van", "Slight — Cyclist", "Slight — Other"
)

results_8 <- map2_df(outcomes_8, labels_8, function(yvar, lab) {
  cat("Running:", lab, "\n")
  att <- att_gt(
    yname = yvar, tname = "qtr_int", idname = "road_id", gname = "g",
    data = panel_did8, control_group = "nevertreated",
    clustervars = "OA", bstrap = FALSE, anticipation = 0, panel = TRUE
  )
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  tibble(
    outcome = lab,
    att     = agg$overall.att,
    se      = agg$overall.se,
    ci_lo   = agg$overall.att - 1.96 * agg$overall.se,
    ci_hi   = agg$overall.att + 1.96 * agg$overall.se,
    pval    = 2 * pnorm(-abs(agg$overall.att / agg$overall.se))
  )
})

results_8 <- results_8 %>%
  mutate(
    sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                    pval < 0.1 ~ "*", TRUE ~ ""),
    severity = if_else(grepl("^KSI", outcome), "KSI", "Slight"),
    outcome = factor(outcome, levels = rev(labels_8)),
    sig_label = paste0(sprintf("%.4f", att), sig)
  )

print(results_8)

# --- Coefficient plot ---

ggplot(results_8, aes(x = att, y = outcome, colour = severity)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = sig_label), hjust = -0.15, size = 3,
            show.legend = FALSE) +
  scale_colour_manual(values = c("KSI" = "#c0392b", "Slight" = "#2980b9")) +
  labs(
    title    = "Staggered DiD ATT estimates by severity and mode",
    subtitle = "Callaway & Sant'Anna (2021) — never-treated controls, OA-clustered SEs",
    x        = "ATT (injuries per road-km)",
    y        = NULL,
    colour   = "Severity",
    caption  = "Error bars: 95% CI  |  * p<0.1  ** p<0.05  *** p<0.01"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave(here("output", "did_severity_mode_coefficients.png"),
       width = 10, height = 7, dpi = 300)

# =============================================================================
# STEP 8 — Scheme-specific DiD by severity × mode
# =============================================================================
# Same 8 outcomes as Step 7, estimated separately for each scheme.

schemes_all <- sort(unique(panel_did8$scheme))

results_by_scheme <- map_df(schemes_all, function(s) {
  cat("=== Scheme:", s, "===\n")

  d <- panel_did8 %>%
    filter(scheme == s) %>%
    mutate(
      road_id = as.integer(factor(identifier)),
      g = as.numeric(g)
    )

  map2_df(outcomes_8, labels_8, function(yvar, lab) {
    att <- tryCatch(
      att_gt(
        yname = yvar, tname = "qtr_int", idname = "road_id", gname = "g",
        data = d, control_group = "nevertreated",
        bstrap = FALSE, anticipation = 0, panel = TRUE
      ),
      error = function(e) { cat("  FAILED:", lab, "\n"); NULL }
    )
    if (is.null(att)) return(tibble(
      scheme = s, outcome = lab, att = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_
    ))
    agg <- aggte(att, type = "simple", na.rm = TRUE)
    tibble(
      scheme  = s,
      outcome = lab,
      att     = agg$overall.att,
      se      = agg$overall.se,
      ci_lo   = agg$overall.att - 1.96 * agg$overall.se,
      ci_hi   = agg$overall.att + 1.96 * agg$overall.se,
      pval    = 2 * pnorm(-abs(agg$overall.att / agg$overall.se))
    )
  })
})

results_by_scheme <- results_by_scheme %>%
  mutate(
    sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                    pval < 0.1 ~ "*", TRUE ~ ""),
    severity = if_else(grepl("^KSI", outcome), "KSI", "Slight"),
    outcome  = factor(outcome, levels = rev(labels_8))
  )

print(results_by_scheme %>% filter(sig != ""), n = Inf)

# --- Scheme-specific coefficient plot ---

ggplot(results_by_scheme, aes(x = att, y = outcome, colour = severity)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.25, linewidth = 0.5, orientation = "y") +
  geom_point(aes(shape = sig != ""), size = 2) +
  scale_colour_manual(values = c("KSI" = "#c0392b", "Slight" = "#2980b9")) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("p >= 0.1", "p < 0.1"), name = "Significance") +
  facet_wrap(~scheme, ncol = 4, scales = "free_x") +
  labs(
    title    = "Staggered DiD ATT estimates by scheme, severity, and mode",
    subtitle = "Callaway & Sant'Anna (2021) — never-treated controls | filled = p < 0.1",
    x        = "ATT (injuries per road-km)",
    y        = NULL,
    colour   = "Severity",
    caption  = "Error bars: 95% CI  |  x-axis scales differ across panels"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 8)
  )

ggsave(here("output", "did_severity_mode_by_scheme.png"),
       width = 14, height = 10, dpi = 300)
