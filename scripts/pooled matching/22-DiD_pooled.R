# =============================================================================
# DiD — POOLED MATCHING (TOTAL INJURIES)
# =============================================================================
#
# Reads the pooled matched panel and estimates:
#   1. Simple pre/post descriptive DiD
#   2. TWFE DiD (road + quarter FE, OA-clustered SEs)
#   3. Callaway & Sant'Anna (2021) staggered DiD — never-treated controls
#   4. Scheme-specific staggered DiD estimates
#
# Single outcome: total injuries per road-km (any injury whatsoever).
#
# Input:  road_panel_matched_pooled.parquet
# Output: output/pooled/
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(fixest)

dir.create(here("output", "pooled"), showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "pooled")

# =============================================================================
# LOAD DATA
# =============================================================================

road_panel_matched <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

cat("Rows:", nrow(road_panel_matched),
    "| Roads:", n_distinct(road_panel_matched$identifier),
    "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n")

# =============================================================================
# CREATE OUTCOME
# =============================================================================

scheme_timing <- road_panel_matched %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

panel <- road_panel_matched %>%
  mutate(
    road_length_km = length / 1000,
    road_length_km = if_else(
      road_length_km <= 0 | is.na(road_length_km), 1e-6, road_length_km
    ),
    total_pkm = total_inj_adj_All / road_length_km,
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  ) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  mutate(post_flag = as.integer(quarter_year >= ref_start))

# =============================================================================
# STEP 1 — DESCRIPTIVE DiD
# =============================================================================

cat("\n================================================================\n")
cat("STEP 1 — DESCRIPTIVE DiD\n")
cat("================================================================\n\n")

summary_stats <- panel %>%
  group_by(group, post_flag) %>%
  summarise(
    mean_total_pkm = round(mean(total_pkm, na.rm = TRUE), 4),
    n_roads = n_distinct(identifier),
    n_obs   = n(),
    .groups = "drop"
  ) %>%
  mutate(period = if_else(post_flag == 1L, "Post", "Pre")) %>%
  select(group, period, mean_total_pkm, n_roads, n_obs)

print(summary_stats)

did_table <- summary_stats %>%
  pivot_wider(id_cols = group, names_from = period,
              values_from = mean_total_pkm) %>%
  mutate(change = round(Post - Pre, 5))

did_est <- did_table$change[did_table$group == "CAZ roads"] -
           did_table$change[did_table$group == "Matched controls"]

treated_pre <- summary_stats %>%
  filter(group == "CAZ roads", period == "Pre") %>%
  pull(mean_total_pkm)

cat("\nSimple DiD estimate (total injuries/km):", round(did_est, 5),
    "(", round(100 * did_est / treated_pre, 1), "%)\n")

# =============================================================================
# STEP 2 — TREND PLOT
# =============================================================================

plot_data <- panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_total_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop")

median_start <- median(scheme_timing$caz_start_q, na.rm = TRUE)

ggplot(plot_data, aes(x = as.Date(quarter_year), y = mean_total_pkm,
                      colour = group)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = as.Date(median_start),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("CAZ roads" = "#E74C3C",
                                 "Matched controls" = "#2C3E50")) +
  labs(
    title = "Total injuries per road-km: CAZ roads vs matched controls (pooled matching)",
    x = NULL, y = "Mean injuries/road-km", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "prepost_descriptive_pooled.png"),
       width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 3 — TWFE DiD
# =============================================================================

cat("\n================================================================\n")
cat("STEP 3 — TWFE DiD\n")
cat("================================================================\n\n")

panel_reg <- panel %>%
  mutate(
    treated  = as.integer(treat_group),
    post_int = post_flag,
    qtr_num  = as.numeric(quarter_year)
  )

m_total <- feols(total_pkm ~ treated:post_int | identifier + qtr_num,
                 data = panel_reg, cluster = ~OA)

etable(m_total,
       headers = "Total injuries/km",
       dict = c("treated:post_int" = "Treated x Post"),
       digits = 5, digits.stats = 3)

# =============================================================================
# STEP 4 — STAGGERED DiD: Callaway & Sant'Anna (2021)
# =============================================================================

cat("\n================================================================\n")
cat("STEP 4 — STAGGERED DiD (Callaway & Sant'Anna)\n")
cat("================================================================\n\n")

min_qtr <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

panel_did <- panel %>%
  mutate(
    road_id = as.integer(factor(identifier)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    # g must be numeric (att_gt sets g = Inf for controls internally)
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    )
  ) %>%
  filter(!is.na(total_pkm))

cat("Time periods:", n_distinct(panel_did$qtr_int),
    "| Roads:", n_distinct(panel_did$road_id), "\n")

panel_did %>%
  distinct(road_id, g) %>%
  count(g) %>%
  mutate(label = if_else(g == 0, "Never treated",
                         as.character(as.yearqtr(min_qtr + (g - 1) / 4)))) %>%
  print()

# --- Overall ATT ---

att_total <- att_gt(
  yname = "total_pkm", tname = "qtr_int", idname = "road_id", gname = "g",
  data = panel_did, control_group = "nevertreated",
  clustervars = "OA", bstrap = FALSE, anticipation = 0, panel = TRUE
)

agg_simple <- aggte(att_total, type = "simple", na.rm = TRUE)
cat("\n=== Overall ATT — Total injuries/km ===\n")
summary(agg_simple)

# --- Event-study ---

agg_dyn <- aggte(att_total, type = "dynamic", na.rm = TRUE)

es_df <- tibble(
  event_time = agg_dyn$egt,
  att        = agg_dyn$att.egt,
  se         = agg_dyn$se.egt,
  ci_lo      = agg_dyn$att.egt - 1.96 * agg_dyn$se.egt,
  ci_hi      = agg_dyn$att.egt + 1.96 * agg_dyn$se.egt
) %>%
  filter(event_time >= -12, event_time <= 12)

ggplot(es_df, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
  labs(
    title = "Event-study: effect of CAZ on total injuries (pooled matching)",
    subtitle = "Callaway & Sant'Anna (2021) — never-treated controls, OA-clustered SEs",
    x = "Quarters relative to CAZ implementation",
    y = "ATT (injuries/road-km)",
    caption = "Shaded: 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "event_study_pooled.png"),
       width = 10, height = 7, dpi = 300)

# =============================================================================
# STEP 5 — SCHEME-SPECIFIC STAGGERED DiD
# =============================================================================

cat("\n================================================================\n")
cat("STEP 5 — SCHEME-SPECIFIC ESTIMATES\n")
cat("================================================================\n\n")

schemes_all <- sort(unique(panel_did$scheme))

results_by_scheme <- map_df(schemes_all, function(s) {
  cat("  Scheme:", s, "\n")
  d <- panel_did %>%
    filter(scheme == s) %>%
    mutate(road_id = as.integer(factor(identifier)),
           g = as.numeric(g))

  att <- tryCatch(
    att_gt(
      yname = "total_pkm", tname = "qtr_int", idname = "road_id", gname = "g",
      data = d, control_group = "nevertreated",
      bstrap = FALSE, anticipation = 0, panel = TRUE
    ),
    error = function(e) { cat("    FAILED\n"); NULL }
  )
  if (is.null(att)) return(tibble(
    scheme = s, att = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_
  ))
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  tibble(
    scheme = s,
    att    = agg$overall.att,
    se     = agg$overall.se,
    ci_lo  = agg$overall.att - 1.96 * agg$overall.se,
    ci_hi  = agg$overall.att + 1.96 * agg$overall.se,
    pval   = 2 * pnorm(-abs(agg$overall.att / agg$overall.se))
  )
})

results_by_scheme <- results_by_scheme %>%
  mutate(
    sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                    pval < 0.1 ~ "*", TRUE ~ ""),
    sig_label = paste0(sprintf("%.4f", att), sig)
  )

print(results_by_scheme)

# Coefficient plot
ggplot(results_by_scheme,
       aes(x = att, y = fct_reorder(scheme, att))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.25, linewidth = 0.7, orientation = "y") +
  geom_point(size = 3, colour = "#2E6FAB") +
  geom_text(aes(label = sig_label), hjust = -0.15, size = 3) +
  labs(
    title = "Staggered DiD ATT by scheme — total injuries/km (pooled matching)",
    subtitle = "Callaway & Sant'Anna (2021) — never-treated controls, OA-clustered SEs",
    x = "ATT (injuries/road-km)", y = NULL,
    caption = "Error bars: 95% CI  |  * p<0.1  ** p<0.05  *** p<0.01"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "did_by_scheme_pooled.png"),
       width = 10, height = 7, dpi = 300)

cat("\n=== ALL OUTPUTS SAVED ===\n")
cat("Output dir:", outdir, "\n")
