# =============================================================================
# PANEL DIAGNOSTICS — POOLED MATCHING (TOTAL INJURIES)
# =============================================================================
#
# Same diagnostic suite as main script 21, adapted for the pooled matching
# pipeline. Single outcome focus: total injuries per road-km.
#
# Input:  road_panel_matched_pooled.parquet
# Output: output/diagnostics/panel_pooled/
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(fixest)
library(patchwork)

dir.create(here("output", "diagnostics", "panel_pooled"),
           showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics", "panel_pooled")

save_fig <- function(p, filename, width = 14, height = 10, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

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

cat("Panel loaded:", nrow(road_panel_matched), "rows |",
    n_distinct(road_panel_matched$identifier), "roads |",
    n_distinct(road_panel_matched$quarter_year), "quarters\n\n")

panel <- road_panel_matched %>%
  mutate(
    road_length_km = length / 1000,
    road_length_km = if_else(
      road_length_km <= 0 | is.na(road_length_km), 1e-6, road_length_km
    ),
    total_pkm = total_inj_adj_All / road_length_km,
    KSI_pkm   = KSI_adj_All      / road_length_km,
    slight_pkm = Slight_adj_All   / road_length_km,
    group = if_else(treat_group == 1, "Treated", "Control")
  )

scheme_timing <- panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

panel <- panel %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme")



# =============================================================================
# OA COVARIATES
# =============================================================================
# Time-invariant OA covariates are retained in the model panel so they are
# available for robustness analyses, especially C&S checks in the separate script.

matched_covars <- readRDS(
  here("data", "processed", "OA_matched_full_pooled.rds")
) %>%
  mutate(
    log1p_pop_density = log1p(pmax(pop_density, 0)),
    log1p_road_density_m_km2 = log1p(pmax(road_density_m_km2, 0))
  ) %>%
  select(OA, log1p_pop_density, IMD, log1p_road_density_m_km2) %>%
  distinct(OA, .keep_all = TRUE)



# =============================================================================
# 1. PANEL STRUCTURE
# =============================================================================


n_roads <- n_distinct(panel$identifier)
n_qtrs  <- n_distinct(panel$quarter_year)
cat("Roads:", n_roads, "| Quarters:", n_qtrs, "| Rows:", nrow(panel), "\n")
cat("Balanced:", nrow(panel) == n_roads * n_qtrs, "\n\n")

na_summary <- panel %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_na") %>%
  filter(n_na > 0) %>%
  mutate(pct_na = round(100 * n_na / nrow(panel), 2))
if (nrow(na_summary) > 0) { cat("Missing values:\n"); print(na_summary) }

cat("\nTreatment breakdown:\n")
panel %>%
  distinct(identifier, treat_group, scheme) %>%
  count(treat_group, scheme) %>%
  pivot_wider(names_from = treat_group, values_from = n,
              names_prefix = "group_") %>%
  rename(control = group_0, treated = group_1) %>%
  mutate(ratio = round(control / treated, 1)) %>%
  print()
cat("\n")

# =============================================================================
# PRE-TREATMENT DISTRIBUTIONS
# =============================================================================

cat("================================================================\n")
cat("2. PRE-TREATMENT DISTRIBUTIONS\n")
cat("================================================================\n\n")

panel_pre <- panel %>% filter(quarter_year < ref_start)

panel_pre %>%
  group_by(group) %>%
  summarise(
    n_obs      = n(),
    n_roads    = n_distinct(identifier),
    mean_total = round(mean(total_pkm, na.rm = TRUE), 4),
    sd_total   = round(sd(total_pkm,   na.rm = TRUE), 4),
    pct_zero   = round(mean(total_inj_adj_All == 0) * 100, 1),
    .groups    = "drop"
  ) %>%
  print()
cat("\n")

# =============================================================================
# 3. PARALLEL TRENDS — VISUAL
# =============================================================================

cat("================================================================\n")
cat("3. PARALLEL TRENDS — VISUAL\n")
cat("================================================================\n\n")

# Pooled trend
trend_pooled <- panel %>%
  group_by(group, quarter_year) %>%
  summarise(total_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop")

median_start <- median(scheme_timing$caz_start_q, na.rm = TRUE)

p_pooled <- ggplot(trend_pooled,
                   aes(x = as.Date(quarter_year), y = total_pkm, colour = group)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = as.Date(median_start),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Parallel trends: total injuries per road-km (pooled)",
    subtitle = "Dashed = median CAZ start",
    x = NULL, y = "Mean injuries/road-km", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_fig(p_pooled, "fig_parallel_trends_pooled.png", width = 12, height = 6)

# Per-scheme trends
trend_scheme <- panel %>%
  group_by(scheme, group, quarter_year) %>%
  summarise(total_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop")

p_scheme <- ggplot(trend_scheme,
                   aes(x = as.Date(quarter_year), y = total_pkm, colour = group)) +
  geom_line(linewidth = 0.7) +
  geom_vline(data = scheme_timing,
             aes(xintercept = as.Date(caz_start_q)),
             linetype = "dashed", colour = "grey40") +
  facet_wrap(~scheme, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Parallel trends by scheme: total injuries per road-km",
    x = NULL, y = "Mean injuries/road-km", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_scheme, "fig_parallel_trends_by_scheme.png", width = 16, height = 10)

# --- 3b. LOESS trajectory plots (road level) ---

# All schemes — LOESS
trend_pooled_loess <- panel %>%
  group_by(group, quarter_year) %>%
  summarise(total_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    group  = factor(group, levels = c("Treated", "Control")),
    q_date = as.Date(quarter_year)
  )

p_loess_all <- ggplot(
  trend_pooled_loess,
  aes(x = q_date, y = total_pkm, colour = group, fill = group)
) +
  geom_point(size = 1, alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  scale_colour_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  scale_fill_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  labs(
    title    = "Injury trajectories: treated vs matched controls (road level, all schemes)",
    subtitle = "LOESS smoother \u00b1 SE | parallel smoothers = parallel trends supported",
    x = NULL, y = "Mean injuries per road-km",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

save_fig(p_loess_all, "fig_parallel_trends_loess_all.png", width = 12, height = 8)

# Per scheme — LOESS
trend_scheme_loess <- panel %>%
  group_by(scheme, group, quarter_year) %>%
  summarise(total_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    group  = factor(group, levels = c("Treated", "Control")),
    q_date = as.Date(quarter_year)
  )

p_loess_scheme <- ggplot(
  trend_scheme_loess,
  aes(x = q_date, y = total_pkm, colour = group, fill = group)
) +
  geom_point(size = 1, alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1.1) +
  geom_vline(data = scheme_timing,
             aes(xintercept = as.Date(caz_start_q)),
             linetype = "dotted", colour = "#888888", linewidth = 0.5) +
  facet_wrap(~scheme, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  scale_fill_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  labs(
    title    = "Injury trajectories by scheme (road level)",
    subtitle = "LOESS smoother \u00b1 SE | dotted = scheme start | parallel smoothers = parallel trends",
    x = NULL, y = "Mean injuries per road-km",
    colour = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(size = 8, angle = 45, hjust = 1),
    strip.text      = element_text(size = 10, face = "bold"),
    panel.spacing   = unit(0.6, "lines")
  )

save_fig(p_loess_scheme, "fig_parallel_trends_loess_by_scheme.png",
         width = 16, height = 14)

# --- 3c. Pre-treatment trend slope densities (road level) ---

cat("Computing per-road pre-treatment trend slopes...\n")

road_trends <- panel_pre %>%
  mutate(t = as.numeric(quarter_year)) %>%
  group_by(identifier, group, scheme) %>%
  summarise(
    trend_slope = coef(lm(total_pkm ~ t))[2],
    .groups = "drop"
  )

# Trim axis to 1st–99th percentile for readability
trend_xlim <- quantile(road_trends$trend_slope, c(0.01, 0.99), na.rm = TRUE)

# All schemes
p_trend_dens_all <- ggplot(
  road_trends,
  aes(x = trend_slope, fill = group, colour = group)
) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
  coord_cartesian(xlim = trend_xlim) +
  scale_fill_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  scale_colour_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  labs(
    title    = "Pre-treatment trend distributions: treated vs matched controls (road level)",
    subtitle = "Per-road-link slope of total injuries/km over pre-treatment quarters",
    x = "Pre-treatment slope", y = "Density",
    fill = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

save_fig(p_trend_dens_all, "fig_trend_density_all.png", width = 12, height = 7)

# Per scheme
p_trend_dens_scheme <- ggplot(
  road_trends,
  aes(x = trend_slope, fill = group, colour = group)
) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#888888") +
  coord_cartesian(xlim = trend_xlim) +
  facet_wrap(~scheme, scales = "free_y") +
  scale_fill_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  scale_colour_manual(values = c("Treated" = "#D85A30", "Control" = "#2E6FAB")) +
  labs(
    title    = "Pre-treatment trend distributions by scheme (road level)",
    subtitle = "Per-road-link slope of total injuries/km over pre-treatment quarters",
    x = "Pre-treatment slope", y = "Density",
    fill = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_fig(p_trend_dens_scheme, "fig_trend_density_by_scheme.png",
         width = 16, height = 12)

# =============================================================================
# 4. EVENT-STUDY — sunab (integer quarters)
# =============================================================================

cat("================================================================\n")
cat("4. EVENT-STUDY (sunab)\n")
cat("================================================================\n\n")

min_qtr_num <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

panel_es <- panel %>%
  mutate(
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr_num) * 4)) + 1L,
    cohort_q = if_else(
      treat_group == 1 & !is.na(caz_start_q),
      as.numeric(round((as.numeric(caz_start_q) - min_qtr_num) * 4)) + 1,
      Inf
    )
  )

es_total <- feols(
  total_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
  data = panel_es, cluster = ~OA
)

extract_sunab <- function(mod, label) {
  broom::tidy(mod, conf.int = TRUE) %>%
    filter(grepl("^qtr_int::", term), !grepl("cohort", term)) %>%
    mutate(
      rel_time = as.numeric(gsub("^qtr_int::", "", term)),
      outcome  = label
    )
}

es_coefs <- extract_sunab(es_total, "Total injuries/km") %>%
  filter(rel_time >= -12, rel_time <= 12)

p_es <- ggplot(es_coefs, aes(x = rel_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
  labs(
    title = "Event-study: total injuries per road-km (pooled matching)",
    subtitle = "sunab() — road + quarter FE, OA-clustered SEs",
    x = "Quarters relative to CAZ implementation",
    y = "Coefficient (injuries/road-km)",
    caption = "Reference: t = -1 quarter  |  Shaded: 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_fig(p_es, "fig_event_study_pooled.png", width = 10, height = 7)

# Pre-trend F-test
all_coefs <- names(coef(es_total))
agg_coefs <- all_coefs[grepl("^qtr_int::", all_coefs) & !grepl("cohort", all_coefs)]
rel_times <- as.numeric(gsub("^qtr_int::", "", agg_coefs))
pre_coefs <- agg_coefs[rel_times < -1]

if (length(pre_coefs) > 0) {
  wt <- tryCatch(wald(es_total, keep = pre_coefs), error = function(e) NULL)
  if (!is.null(wt)) {
    cat("Joint F-test (all pre-treatment coefficients = 0):\n")
    cat("  F =", round(wt$stat, 3), "| p =", round(wt$p, 4), "\n")
    cat("  Interpretation: p > 0.05 → fail to reject parallel trends\n\n")
  }
}

# =============================================================================
# 5. SCHEME-SPECIFIC EVENT STUDIES
# =============================================================================

cat("================================================================\n")
cat("5. SCHEME-SPECIFIC EVENT STUDIES\n")
cat("================================================================\n\n")

schemes_all <- sort(unique(panel_es$scheme))

scheme_es <- map_df(schemes_all, function(s) {
  cat("  Event study:", s, "\n")
  d <- panel_es %>% filter(scheme == s)
  mod <- tryCatch(
    feols(total_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
          data = d, cluster = ~OA),
    error = function(e) { cat("    FAILED\n"); NULL }
  )
  if (is.null(mod)) return(tibble())

  cf <- extract_sunab(mod, "Total") %>% mutate(scheme = s)

  all_c <- names(coef(mod))
  agg_c <- all_c[grepl("^qtr_int::", all_c) & !grepl("cohort", all_c)]
  rt <- as.numeric(gsub("^qtr_int::", "", agg_c))
  pre_c <- agg_c[rt < -1]
  ft <- tryCatch(wald(mod, keep = pre_c), error = function(e) NULL)
  if (!is.null(ft)) cf <- cf %>% mutate(pretrend_p = round(ft$p, 4))
  cf
})

scheme_es_plot <- scheme_es %>% filter(rel_time >= -12, rel_time <= 12)

scheme_labels <- scheme_es %>%
  distinct(scheme, pretrend_p) %>%
  mutate(label = paste0(scheme, "\n(pre-trend p=", sprintf("%.3f", pretrend_p), ")"))

scheme_es_plot <- scheme_es_plot %>%
  left_join(scheme_labels %>% select(scheme, label), by = "scheme")

p_es_scheme <- ggplot(scheme_es_plot, aes(x = rel_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_wrap(~label, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
  labs(
    title = "Event-study by scheme: total injuries/road-km (pooled matching)",
    subtitle = "sunab() — pre-trend p-value in label",
    x = "Quarters relative to CAZ", y = "Coefficient",
    caption = "Reference: t=-1  |  p > 0.05 = fail to reject parallel trends"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

save_fig(p_es_scheme, "fig_event_study_by_scheme.png", width = 16, height = 10)

# =============================================================================
# 6. POWER ANALYSIS PER SCHEME
# =============================================================================

cat("================================================================\n")
cat("6. POWER ANALYSIS PER SCHEME\n")
cat("================================================================\n\n")

power_table <- map_df(schemes_all, function(s) {
  d <- panel %>% filter(scheme == s)
  d_pre <- d %>% filter(quarter_year < ref_start)

  n_treated  <- n_distinct(d$identifier[d$treat_group == 1])
  n_control  <- n_distinct(d$identifier[d$treat_group == 0])
  n_pre_qtrs <- n_distinct(d_pre$quarter_year[d_pre$treat_group == 1])
  n_post_qtrs <- n_distinct(d$quarter_year[d$treat_group == 1 &
                                             d$quarter_year >= d$ref_start])

  treated_pre_mean <- mean(d_pre$total_pkm[d_pre$treat_group == 1], na.rm = TRUE)

  mod_pre <- tryCatch(
    feols(total_pkm ~ 1 | identifier + quarter_year, data = d_pre, cluster = ~OA),
    error = function(e) NULL
  )
  resid_sd <- if (!is.null(mod_pre)) sd(resid(mod_pre)) else NA_real_

  d_ctrl_pre <- d_pre %>% filter(treat_group == 0)
  oa_means   <- d_ctrl_pre %>%
    group_by(OA) %>%
    summarise(oa_mean = mean(total_pkm, na.rm = TRUE), .groups = "drop")
  between_var <- var(oa_means$oa_mean, na.rm = TRUE)
  total_var   <- var(d_ctrl_pre$total_pkm, na.rm = TRUE)
  icc <- if (!is.na(total_var) && total_var > 0) between_var / total_var else NA_real_

  n_oa_treated <- n_distinct(d$OA[d$treat_group == 1])
  roads_per_oa <- n_treated / max(n_oa_treated, 1)
  deff <- 1 + (roads_per_oa - 1) * icc

  n_eff_treated <- n_treated / deff
  n_eff_control <- n_control / deff
  n_eff <- (n_eff_treated * n_eff_control) / (n_eff_treated + n_eff_control)

  se_did <- if (!is.na(resid_sd))
    resid_sd * sqrt(1 / (n_eff * n_post_qtrs) + 1 / (n_eff * n_pre_qtrs))
  else NA_real_

  mde <- se_did * (qnorm(0.975) + qnorm(0.80))

  tibble(
    scheme = s, n_treated = n_treated, n_control = n_control,
    n_pre_qtrs = n_pre_qtrs, n_post_qtrs = n_post_qtrs,
    treated_pre_mean = round(treated_pre_mean, 4),
    resid_sd = round(resid_sd, 4), ICC = round(icc, 4),
    deff = round(deff, 2), MDE = round(mde, 4),
    MDE_pct = round(100 * mde / treated_pre_mean, 1)
  )
})

cat("Power analysis — MDE (80% power, alpha = 0.05):\n")
print(power_table %>% select(scheme, n_treated, n_control, resid_sd, ICC, MDE, MDE_pct))

p_mde <- ggplot(power_table, aes(y = fct_reorder(scheme, MDE_pct), x = MDE_pct)) +
  geom_col(fill = "#2E6FAB", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = paste0(MDE_pct, "%")), hjust = -0.15, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "MDE by scheme (pooled matching, 80% power)",
       x = "MDE (% of pre-treatment mean)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

save_fig(p_mde, "fig_MDE_by_scheme.png", width = 10, height = 7)

# =============================================================================
# 7. SPARSITY & SERIAL CORRELATION
# =============================================================================

cat("================================================================\n")
cat("7. SPARSITY & SERIAL CORRELATION\n")
cat("================================================================\n\n")

# Zero rates
zero_rates <- panel %>%
  group_by(scheme, group) %>%
  summarise(
    pct_zero = round(mean(total_inj_adj_All == 0) * 100, 1),
    .groups = "drop"
  )
cat("Zero-inflation (% zero total injuries):\n")
print(zero_rates, n = Inf)

# Serial correlation
mod_pre_all <- feols(
  total_pkm ~ treat_group | identifier + quarter_year,
  data = panel_pre
)

resid_df <- panel_pre %>%
  select(identifier, quarter_year) %>%
  mutate(resid = resid(mod_pre_all)) %>%
  arrange(identifier, quarter_year) %>%
  group_by(identifier) %>%
  mutate(resid_lag = dplyr::lag(resid)) %>%
  ungroup() %>%
  filter(!is.na(resid_lag))

ar1_mod <- lm(resid ~ resid_lag, data = resid_df)
ar1_coef <- coef(ar1_mod)["resid_lag"]
ar1_pval <- summary(ar1_mod)$coefficients["resid_lag", "Pr(>|t|)"]

cat("\nAR(1) on FE residuals:", round(ar1_coef, 4),
    "| p =", format.pval(ar1_pval, digits = 4), "\n")
if (ar1_pval < 0.05) {
  cat(">> Serial correlation detected — cluster-robust SEs essential.\n")
} else {
  cat(">> No significant serial correlation.\n")
}

# =============================================================================
# SAVE COMBINED SUMMARY
# =============================================================================

pretrend_summary <- scheme_es %>%
  distinct(scheme, pretrend_p)

diagnostics_summary <- power_table %>%
  left_join(pretrend_summary, by = "scheme") %>%
  left_join(
    zero_rates %>% filter(group == "Treated") %>%
      select(scheme, pct_zero),
    by = "scheme"
  )

write_csv(diagnostics_summary,
          file.path(outdir, "panel_diagnostics_summary_pooled.csv"))

cat("\n=== ALL DIAGNOSTICS COMPLETE ===\n")
cat("Figures saved to:", outdir, "\n")
