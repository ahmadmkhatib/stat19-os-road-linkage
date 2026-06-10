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
# STEP 0 — INJURY SUMMARY STATISTICS
# =============================================================================

cat("\n================================================================\n")
cat("STEP 0 — INJURY SUMMARY STATISTICS\n")
cat("================================================================\n\n")

cat("--- Quarterly injury counts per road link ---\n")
panel %>%
  group_by(group) %>%
  summarise(
    n_obs    = format(n(), big.mark = ","),
    mean     = round(mean(total_inj_adj_All, na.rm = TRUE), 4),
    sd       = round(sd(total_inj_adj_All, na.rm = TRUE), 4),
    median   = median(total_inj_adj_All, na.rm = TRUE),
    q95      = round(quantile(total_inj_adj_All, 0.95, na.rm = TRUE), 4),
    max      = round(max(total_inj_adj_All, na.rm = TRUE), 2),
    pct_zero = round(100 * mean(total_inj_adj_All == 0, na.rm = TRUE), 1),
    .groups  = "drop"
  ) %>% print()

cat("\n--- Injuries per road-km ---\n")
panel %>%
  group_by(group) %>%
  summarise(
    mean     = round(mean(total_pkm, na.rm = TRUE), 4),
    sd       = round(sd(total_pkm, na.rm = TRUE), 4),
    median   = round(median(total_pkm, na.rm = TRUE), 4),
    q95      = round(quantile(total_pkm, 0.95, na.rm = TRUE), 4),
    max      = round(max(total_pkm, na.rm = TRUE), 2),
    pct_zero = round(100 * mean(total_pkm == 0, na.rm = TRUE), 1),
    .groups  = "drop"
  ) %>% print()

cat("\n--- Pre vs Post (treated roads) ---\n")
panel %>%
  filter(treat_group == 1) %>%
  group_by(period = if_else(post_flag == 1, "Post", "Pre")) %>%
  summarise(
    n_obs      = format(n(), big.mark = ","),
    mean_count = round(mean(total_inj_adj_All, na.rm = TRUE), 4),
    mean_pkm   = round(mean(total_pkm, na.rm = TRUE), 4),
    pct_zero   = round(100 * mean(total_inj_adj_All == 0, na.rm = TRUE), 1),
    .groups    = "drop"
  ) %>% print()

cat("\n--- Road-level total injuries (summed across all quarters) ---\n")
panel %>%
  group_by(panel_id, group) %>%
  summarise(total = sum(total_inj_adj_All, na.rm = TRUE), .groups = "drop") %>%
  group_by(group) %>%
  summarise(
    n_roads      = format(n(), big.mark = ","),
    mean         = round(mean(total), 3),
    median       = median(total),
    pct_always_0 = round(100 * mean(total == 0), 1),
    pct_1_plus   = round(100 * mean(total >= 1), 1),
    .groups      = "drop"
  ) %>% print()

cat("\n")

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
  dplyr::select(group, period, mean_total_pkm, n_roads, n_obs)

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

# Unit ID must be unique per id × time. Because control roads can appear
# in multiple schemes (from replacement matching), use panel_id × scheme
# as the unit identifier to avoid duplicate id × time rows.
panel_did <- panel %>%
  mutate(
    uid     = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    )
  ) %>%
  filter(!is.na(total_pkm))

cat("Time periods:", n_distinct(panel_did$qtr_int),
    "| Panel units:", n_distinct(panel_did$uid_int), "\n")

panel_did %>%
  distinct(uid_int, g) %>%
  count(g) %>%
  mutate(label = if_else(g == 0, "Never treated",
                         as.character(as.yearqtr(min_qtr + (g - 1) / 4)))) %>%
  print()

# --- Overall ATT ---

att_total <- att_gt(
  yname = "total_pkm", tname = "qtr_int", idname = "uid_int", gname = "g",
  data = panel_did, control_group = "nevertreated",
  clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE
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
    mutate(uid_int = as.integer(factor(uid)),
           g = as.numeric(g))

  att <- tryCatch(
    att_gt(
      yname = "total_pkm", tname = "qtr_int", idname = "uid_int", gname = "g",
      data = d, control_group = "nevertreated",
      clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE
    ),
    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL }
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

# --- Event study per scheme ---

cat("\n--- Scheme-specific event studies ---\n")

es_by_scheme <- map_df(schemes_all, function(s) {
  cat("  Event study:", s, "\n")
  d <- panel_did %>%
    filter(scheme == s) %>%
    mutate(uid_int = as.integer(factor(uid)),
           g = as.numeric(g))

  att <- tryCatch(
    att_gt(yname = "total_pkm", tname = "qtr_int", idname = "uid_int", gname = "g",
           data = d, control_group = "nevertreated",
           clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE),
    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL }
  )
  if (is.null(att)) return(NULL)

  agg_dyn <- tryCatch(
    aggte(att, type = "dynamic", na.rm = TRUE),
    error = function(e) { cat("    aggte FAILED:", e$message, "\n"); NULL }
  )
  if (is.null(agg_dyn)) return(NULL)

  tibble(
    scheme     = s,
    event_time = agg_dyn$egt,
    att        = agg_dyn$att.egt,
    se         = agg_dyn$se.egt,
    ci_lo      = agg_dyn$att.egt - 1.96 * agg_dyn$se.egt,
    ci_hi      = agg_dyn$att.egt + 1.96 * agg_dyn$se.egt
  )
})

es_by_scheme_plot <- es_by_scheme %>%
  filter(event_time >= -12, event_time <= 12)

ggplot(es_by_scheme_plot, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
  facet_wrap(~scheme, scales = "free_y", ncol = 2) +
  labs(
    title = "Event study by CAZ scheme — total injuries/km",
    subtitle = "Callaway & Sant'Anna (2021), OA-clustered SEs",
    x = "Quarters relative to CAZ implementation",
    y = "ATT (injuries / road-km)",
    caption = "Shaded: 95% CI. Dashed vertical = treatment onset."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(outdir, "event_study_by_scheme_pooled.png"),
       width = 12, height = 10, dpi = 300)

# =============================================================================
# STEP 6 — SENSITIVITY: RESTRICT TO INJURY-PRONE ROADS
# =============================================================================
# Roads with zero injuries in every pre-treatment quarter are structurally
# unable to show a treatment effect. Restricting to roads with ≥1 pre-treatment
# injury gives the ATT on the policy-relevant subpopulation.

cat("\n================================================================\n")
cat("STEP 6 — SENSITIVITY: INJURY-PRONE ROADS ONLY\n")
cat("================================================================\n\n")

# Identify units with ≥1 injury before their scheme's start
scheme_starts <- panel_did %>%
  filter(g > 0) %>%
  distinct(scheme, g) %>%
  rename(scheme_start = g)

pre_injury_units <- panel_did %>%
  left_join(scheme_starts, by = "scheme") %>%
  mutate(ref_start = if_else(g > 0, g, scheme_start)) %>%
  filter(qtr_int < ref_start) %>%
  group_by(uid_int) %>%
  summarise(any_pre_injury = any(total_pkm > 0, na.rm = TRUE), .groups = "drop")

keep_uids <- pre_injury_units$uid_int[pre_injury_units$any_pre_injury]

panel_filtered <- panel_did %>%
  filter(uid_int %in% keep_uids) %>%
  mutate(uid_filt = as.integer(factor(uid_int)))

cat("Retained", n_distinct(panel_filtered$uid_filt), "of",
    n_distinct(panel_did$uid_int), "panel units",
    sprintf("(%.1f%%)\n", 100 * n_distinct(panel_filtered$uid_filt) /
              n_distinct(panel_did$uid_int)))

panel_filtered %>%
  distinct(uid_filt, treat_group) %>%
  count(treat_group) %>%
  mutate(label = if_else(treat_group == 1, "Treated", "Control")) %>%
  print()

cat("Zero rate:", round(100 * mean(panel_filtered$total_pkm == 0), 1), "%\n")
cat("Mean outcome:", round(mean(panel_filtered$total_pkm, na.rm = TRUE), 4), "\n\n")

# --- Overall ATT ---

att_filtered <- att_gt(
  yname = "total_pkm", tname = "qtr_int", idname = "uid_filt", gname = "g",
  data = panel_filtered, control_group = "nevertreated",
  clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE
)

agg_filt <- aggte(att_filtered, type = "simple", na.rm = TRUE)
cat("=== Overall ATT (injury-prone roads) ===\n")
summary(agg_filt)

cat("\nComparison:\n")
cat("  Full sample: ATT =", round(agg_simple$overall.att, 4),
    "  SE =", round(agg_simple$overall.se, 4), "\n")
cat("  Filtered:    ATT =", round(agg_filt$overall.att, 4),
    "  SE =", round(agg_filt$overall.se, 4), "\n\n")

# --- Event study ---

agg_dyn_filt <- aggte(att_filtered, type = "dynamic", na.rm = TRUE)

es_filt <- tibble(
  event_time = agg_dyn_filt$egt,
  att        = agg_dyn_filt$att.egt,
  se         = agg_dyn_filt$se.egt
) %>%
  mutate(ci_lo = att - 1.96 * se, ci_hi = att + 1.96 * se) %>%
  filter(event_time >= -12, event_time <= 15)

ggplot(es_filt, aes(x = event_time, y = att)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq(-12, 15, by = 2)) +
  labs(
    title = "Event study: injury-prone roads only (pre-treatment injury > 0)",
    subtitle = "Callaway & Sant'Anna (2021), OA-clustered SEs",
    x = "Quarters relative to CAZ implementation",
    y = "ATT (injuries / road-km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "event_study_filtered_pooled.png"),
       width = 10, height = 7, dpi = 300)

# --- Scheme-specific ---

results_filt <- map_df(schemes_all, function(s) {
  cat("  Scheme:", s, "\n")
  d <- panel_filtered %>%
    filter(scheme == s) %>%
    mutate(uid_filt = as.integer(factor(uid_filt)), g = as.numeric(g))

  att <- tryCatch(
    att_gt(yname = "total_pkm", tname = "qtr_int", idname = "uid_filt", gname = "g",
           data = d, control_group = "nevertreated",
           clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE),
    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL }
  )
  if (is.null(att)) return(tibble(
    scheme = s, att = NA_real_, se = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_
  ))
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  tibble(scheme = s, att = agg$overall.att, se = agg$overall.se,
         ci_lo = agg$overall.att - 1.96 * agg$overall.se,
         ci_hi = agg$overall.att + 1.96 * agg$overall.se,
         pval = 2 * pnorm(-abs(agg$overall.att / agg$overall.se)))
}) %>%
  mutate(
    sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                    pval < 0.1 ~ "*", TRUE ~ ""),
    sig_label = paste0(sprintf("%.4f", att), sig)
  )

cat("\n=== Scheme-specific ATT (injury-prone roads) ===\n")
print(results_filt)

ggplot(results_filt, aes(x = att, y = fct_reorder(scheme, att))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.25, linewidth = 0.7, orientation = "y") +
  geom_point(size = 3, colour = "#2E6FAB") +
  geom_text(aes(label = sig_label), hjust = -0.15, size = 3) +
  labs(
    title = "DiD ATT by scheme — injury-prone roads only (pooled matching)",
    subtitle = "Callaway & Sant'Anna (2021), OA-clustered SEs",
    x = "ATT (injuries / road-km)", y = NULL,
    caption = "Error bars: 95% CI  |  * p<0.1  ** p<0.05  *** p<0.01"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "did_by_scheme_filtered_pooled.png"),
       width = 10, height = 7, dpi = 300)

cat("\n=== ALL OUTPUTS SAVED ===\n")
cat("Output dir:", outdir, "\n")
