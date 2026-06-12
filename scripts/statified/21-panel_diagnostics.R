# =============================================================================
# COMPREHENSIVE PANEL DIAGNOSTICS — PRE-DiD ASSESSMENT
# =============================================================================
#
# Purpose:
#   Reads the matched road × quarter panel from script 20 and runs all
#   diagnostics needed before estimating DiD:
#     1. Panel structure & descriptives
#     2. Pre-treatment outcome distributions
#     3. Parallel trends — visual (pooled + per-scheme)
#     4. Formal pre-trend tests (event-study with sunab)
#     5. Scheme-specific event-study plots
#     6. Power analysis per scheme
#     7. Outcome sparsity & distributional checks
#     8. Serial correlation check
#
# Input:
#   road_panel_matched.parquet  (from script 20)
#
# Output:
#   Figures saved to output/diagnostics/panel/
#   Summary CSV: output/diagnostics/panel/panel_diagnostics_summary.csv
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(fixest)
library(patchwork)

dir.create(here("output", "diagnostics", "panel"),
           showWarnings = FALSE, recursive = TRUE)
outdir <- here("output", "diagnostics", "panel")

save_fig <- function(p, filename, width = 14, height = 10, dpi = 300) {
  ggsave(file.path(outdir, filename), p,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

# =============================================================================
# LOAD DATA
# =============================================================================

road_panel_matched <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

cat("Panel loaded:", nrow(road_panel_matched), "rows |",
    n_distinct(road_panel_matched$identifier), "roads |",
    n_distinct(road_panel_matched$quarter_year), "quarters\n\n")

# Create per-km rates
panel <- road_panel_matched %>%
  mutate(
    road_length_km = length / 1000,
    road_length_km = if_else(
      road_length_km <= 0 | is.na(road_length_km), 1e-6, road_length_km
    ),
    total_pkm      = total_inj_adj_All       / road_length_km,
    KSI_pkm        = KSI_adj_All             / road_length_km,
    slight_pkm     = Slight_adj_All          / road_length_km,
    ped_pkm        = total_inj_adj_Pedestrian / road_length_km,
    cyc_pkm        = total_inj_adj_Cyclist    / road_length_km,
    car_pkm        = total_inj_adj_Car.Van    / road_length_km,
    group = if_else(treat_group == 1, "Treated", "Control")
  )

# Scheme-specific start dates
scheme_timing <- panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

panel <- panel %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme")

# =============================================================================
# 1. PANEL STRUCTURE & DESCRIPTIVES
# =============================================================================

cat("================================================================\n")
cat("1. PANEL STRUCTURE & DESCRIPTIVES\n")
cat("================================================================\n\n")

n_roads <- n_distinct(panel$identifier)
n_qtrs  <- n_distinct(panel$quarter_year)

cat("Roads:    ", n_roads, "\n")
cat("Quarters: ", n_qtrs, "\n")
cat("Rows:     ", nrow(panel), "\n")
cat("Expected: ", n_roads * n_qtrs, "\n")
cat("Balanced: ", nrow(panel) == n_roads * n_qtrs, "\n\n")

# Missingness
na_summary <- panel %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_na") %>%
  filter(n_na > 0) %>%
  mutate(pct_na = round(100 * n_na / nrow(panel), 2))

cat("Columns with missing values:\n")
if (nrow(na_summary) > 0) print(na_summary) else cat("  None\n")
cat("\n")

# Treatment breakdown
cat("Treatment breakdown (road level):\n")
panel %>%
  distinct(identifier, treat_group, scheme) %>%
  count(treat_group, scheme) %>%
  pivot_wider(names_from = treat_group, values_from = n,
              names_prefix = "group_") %>%
  rename(control = group_0, treated = group_1) %>%
  mutate(ratio = round(control / treated, 1)) %>%
  print()
cat("\n")

# Scheme × timing
cat("Scheme start dates:\n")
print(scheme_timing %>% arrange(caz_start_q))
cat("\n")

# Road class distribution by group
cat("Road class by treatment group:\n")
panel %>%
  distinct(identifier, treat_group, road_class) %>%
  count(treat_group, road_class) %>%
  group_by(treat_group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print(n = Inf)
cat("\n")

# =============================================================================
# 2. PRE-TREATMENT OUTCOME DISTRIBUTIONS
# =============================================================================

cat("================================================================\n")
cat("2. PRE-TREATMENT OUTCOME DISTRIBUTIONS\n")
cat("================================================================\n\n")

panel_pre <- panel %>% filter(quarter_year < ref_start)

outcome_summary <- panel_pre %>%
  group_by(group) %>%
  summarise(
    n_obs      = n(),
    n_roads    = n_distinct(identifier),
    mean_total = round(mean(total_pkm,  na.rm = TRUE), 4),
    sd_total   = round(sd(total_pkm,    na.rm = TRUE), 4),
    mean_KSI   = round(mean(KSI_pkm,    na.rm = TRUE), 4),
    mean_slight = round(mean(slight_pkm, na.rm = TRUE), 4),
    mean_ped   = round(mean(ped_pkm,    na.rm = TRUE), 4),
    mean_cyc   = round(mean(cyc_pkm,    na.rm = TRUE), 4),
    pct_zero_total = round(mean(total_inj_adj_All == 0) * 100, 1),
    .groups = "drop"
  )

print(outcome_summary)
cat("\n")

# Pre-treatment distribution histograms — road-level means
road_pre_means <- panel_pre %>%
  group_by(identifier, group) %>%
  summarise(
    mean_total = mean(total_pkm, na.rm = TRUE),
    mean_KSI   = mean(KSI_pkm,  na.rm = TRUE),
    .groups = "drop"
  )

p_dist <- road_pre_means %>%
  pivot_longer(c(mean_total, mean_KSI), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = if_else(outcome == "mean_total",
                           "Total injuries/km", "KSI/km")) %>%
  ggplot(aes(x = value, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  facet_wrap(~outcome, scales = "free") +
  scale_fill_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Pre-treatment outcome distributions (road-level means)",
    x = "Mean injuries per road-km", y = "Count", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_fig(p_dist, "fig_pre_treatment_distributions.png", width = 12, height = 6)

# =============================================================================
# 3. PARALLEL TRENDS — VISUAL
# =============================================================================

cat("================================================================\n")
cat("3. PARALLEL TRENDS — VISUAL\n")
cat("================================================================\n\n")

# --- 3a. Pooled trend plot (multiple outcomes) ---

trend_pooled <- panel %>%
  group_by(group, quarter_year) %>%
  summarise(
    total_pkm = mean(total_pkm, na.rm = TRUE),
    KSI_pkm   = mean(KSI_pkm,   na.rm = TRUE),
    slight_pkm = mean(slight_pkm, na.rm = TRUE),
    ped_pkm   = mean(ped_pkm,   na.rm = TRUE),
    cyc_pkm   = mean(cyc_pkm,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(total_pkm, KSI_pkm, slight_pkm, ped_pkm, cyc_pkm),
               names_to = "outcome", values_to = "value") %>%
  mutate(outcome = case_when(
    outcome == "total_pkm"  ~ "Total",
    outcome == "KSI_pkm"    ~ "KSI",
    outcome == "slight_pkm" ~ "Slight",
    outcome == "ped_pkm"    ~ "Pedestrian",
    outcome == "cyc_pkm"    ~ "Cyclist"
  ))

median_start <- median(scheme_timing$caz_start_q, na.rm = TRUE)

p_trends_pooled <- ggplot(
  trend_pooled,
  aes(x = as.Date(quarter_year), y = value, colour = group)
) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date(median_start),
             linetype = "dashed", colour = "grey40") +
  facet_wrap(~outcome, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title    = "Parallel trends: treated vs matched controls (pooled)",
    subtitle = "Dashed line = median CAZ start",
    x = NULL, y = "Mean injuries per road-km", colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_fig(p_trends_pooled, "fig_parallel_trends_pooled.png", width = 14, height = 10)

# --- 3b. Per-scheme trend plots (total injuries) ---

trend_scheme <- panel %>%
  group_by(scheme, group, quarter_year) %>%
  summarise(
    total_pkm = mean(total_pkm, na.rm = TRUE),
    KSI_pkm   = mean(KSI_pkm,  na.rm = TRUE),
    .groups = "drop"
  )

p_trends_scheme <- ggplot(
  trend_scheme,
  aes(x = as.Date(quarter_year), y = total_pkm, colour = group)
) +
  geom_line(linewidth = 0.7) +
  geom_vline(
    data = scheme_timing,
    aes(xintercept = as.Date(caz_start_q)),
    linetype = "dashed", colour = "grey40"
  ) +
  facet_wrap(~scheme, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Parallel trends by scheme: total injuries per road-km",
    subtitle = "Dashed line = scheme-specific CAZ start",
    x = NULL, y = "Mean injuries per road-km", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_trends_scheme, "fig_parallel_trends_by_scheme.png", width = 16, height = 10)

# Same for KSI
p_trends_scheme_ksi <- ggplot(
  trend_scheme,
  aes(x = as.Date(quarter_year), y = KSI_pkm, colour = group)
) +
  geom_line(linewidth = 0.7) +
  geom_vline(
    data = scheme_timing,
    aes(xintercept = as.Date(caz_start_q)),
    linetype = "dashed", colour = "grey40"
  ) +
  facet_wrap(~scheme, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Parallel trends by scheme: KSI per road-km",
    subtitle = "Dashed line = scheme-specific CAZ start",
    x = NULL, y = "Mean KSI per road-km", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(p_trends_scheme_ksi, "fig_parallel_trends_by_scheme_KSI.png",
         width = 16, height = 10)

# =============================================================================
# 4. FORMAL PRE-TREND TESTS — EVENT-STUDY WITH sunab()
# =============================================================================

cat("================================================================\n")
cat("4. FORMAL PRE-TREND TESTS — EVENT-STUDY (sunab)\n")
cat("================================================================\n\n")

# Prepare integer quarter index and cohort for sunab
# Using integer quarters (not numeric yearqtr) so the event-study x-axis
# is genuinely in quarters. sunab(cohort, period) computes relative time
# as period - cohort, so both must be on the same scale.

min_qtr_num <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

panel_es <- panel %>%
  mutate(
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr_num) * 4)) + 1L,
    # Cohort = first treated quarter (integer); Inf = never treated
    cohort_q = if_else(
      treat_group == 1 & !is.na(caz_start_q),
      as.numeric(round((as.numeric(caz_start_q) - min_qtr_num) * 4)) + 1,
      Inf
    )
  )

# --- Pooled event-study models ---

es_total <- feols(
  total_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
  data = panel_es, cluster = ~OA
)
es_KSI <- feols(
  KSI_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
  data = panel_es, cluster = ~OA
)
es_slight <- feols(
  slight_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
  data = panel_es, cluster = ~OA
)

# Extract aggregated ATT coefficients (qtr_int::X without cohort_q::)
extract_sunab <- function(mod, label) {
  cf <- broom::tidy(mod, conf.int = TRUE) %>%
    filter(grepl("^qtr_int::", term), !grepl("cohort", term)) %>%
    mutate(
      rel_time = as.numeric(gsub("^qtr_int::", "", term)),
      outcome  = label
    )
  cf
}

es_coefs <- bind_rows(
  extract_sunab(es_total,  "Total"),
  extract_sunab(es_KSI,    "KSI"),
  extract_sunab(es_slight, "Slight")
) %>%
  filter(rel_time >= -12, rel_time <= 12)

# Event-study plot
p_es_pooled <- ggplot(es_coefs, aes(x = rel_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~outcome, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
  labs(
    title    = "Event-study: pre- and post-treatment coefficients (pooled)",
    subtitle = "Sun & Abraham (2021) via sunab() — road + quarter FE, OA-clustered SEs",
    x = "Quarters relative to CAZ implementation",
    y = "Coefficient (injuries per road-km)",
    caption = "Reference period: t = -1  |  Shaded band: 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

save_fig(p_es_pooled, "fig_event_study_pooled.png", width = 10, height = 12)

# --- Joint F-test: all pre-treatment coefficients = 0 ---

pre_test <- function(mod, label) {
  # Aggregated coefficients: qtr_int::X (no cohort)
  all_coefs <- names(coef(mod))
  agg_coefs <- all_coefs[grepl("^qtr_int::", all_coefs) &
                           !grepl("cohort", all_coefs)]
  rel_times <- as.numeric(gsub("^qtr_int::", "", agg_coefs))
  pre_coefs <- agg_coefs[rel_times < -1]
  if (length(pre_coefs) == 0) return(tibble(outcome = label, F_stat = NA,
                                             p_value = NA, n_pre_coefs = 0L))
  wt <- tryCatch(
    wald(mod, keep = pre_coefs),
    error = function(e) NULL
  )
  if (is.null(wt)) return(tibble(outcome = label, F_stat = NA,
                                   p_value = NA, n_pre_coefs = length(pre_coefs)))
  tibble(
    outcome     = label,
    F_stat      = round(wt$stat, 3),
    p_value     = round(wt$p, 4),
    n_pre_coefs = length(pre_coefs)
  )
}

pretrend_tests <- bind_rows(
  pre_test(es_total,  "Total"),
  pre_test(es_KSI,    "KSI"),
  pre_test(es_slight, "Slight")
)

cat("Joint F-test: all pre-treatment coefficients = 0\n")
print(pretrend_tests)
cat("Interpretation: p > 0.05 → fail to reject parallel trends (good)\n\n")

# =============================================================================
# 5. SCHEME-SPECIFIC EVENT-STUDY PLOTS
# =============================================================================

cat("================================================================\n")
cat("5. SCHEME-SPECIFIC EVENT-STUDY PLOTS\n")
cat("================================================================\n\n")

schemes_all <- sort(unique(panel_es$scheme))

scheme_es_results <- map_df(schemes_all, function(s) {
  cat("  Event study:", s, "\n")
  d <- panel_es %>% filter(scheme == s)

  mod <- tryCatch(
    feols(
      total_pkm ~ sunab(cohort_q, qtr_int, ref.p = -1) | identifier + qtr_int,
      data = d, cluster = ~OA
    ),
    error = function(e) { cat("    FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(mod)) return(tibble())

  # Coefficients
  cf <- extract_sunab(mod, "Total") %>%
    mutate(scheme = s)

  # Pre-trend F-test on aggregated coefficients
  all_coefs <- names(coef(mod))
  agg_coefs <- all_coefs[grepl("^qtr_int::", all_coefs) &
                           !grepl("cohort", all_coefs)]
  rel_times <- as.numeric(gsub("^qtr_int::", "", agg_coefs))
  pre_only  <- agg_coefs[rel_times < -1]
  ftest <- tryCatch(wald(mod, keep = pre_only), error = function(e) NULL)

  if (!is.null(ftest)) {
    cf <- cf %>%
      mutate(
        pretrend_F = round(ftest$stat, 3),
        pretrend_p = round(ftest$p, 4)
      )
  }
  cf
})

# Plot: scheme-specific event studies (total injuries)
scheme_es_plot <- scheme_es_results %>%
  filter(rel_time >= -12, rel_time <= 12)

# Add pre-trend test result to facet label
scheme_labels <- scheme_es_results %>%
  distinct(scheme, pretrend_p) %>%
  mutate(
    label = paste0(scheme, "\n(pre-trend p=",
                   sprintf("%.3f", pretrend_p), ")")
  )

scheme_es_plot <- scheme_es_plot %>%
  left_join(scheme_labels %>% select(scheme, label), by = "scheme")

p_es_scheme <- ggplot(scheme_es_plot, aes(x = rel_time, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~label, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
  labs(
    title    = "Event-study by scheme: total injuries per road-km",
    subtitle = "sunab() — road + quarter FE, OA-clustered SEs  |  pre-trend p-value in label",
    x = "Quarters relative to CAZ implementation",
    y = "Coefficient",
    caption = "Reference: t = -1  |  p > 0.05 = fail to reject parallel trends"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

save_fig(p_es_scheme, "fig_event_study_by_scheme.png", width = 16, height = 10)

# Pre-trend test summary table
pretrend_scheme_summary <- scheme_es_results %>%
  distinct(scheme, pretrend_F, pretrend_p) %>%
  mutate(
    parallel_trends = if_else(pretrend_p > 0.05, "PASS", "FAIL"),
    parallel_trends_10 = if_else(pretrend_p > 0.10, "PASS", "FAIL (marginal)")
  ) %>%
  arrange(pretrend_p)

cat("\nScheme-level pre-trend tests (total injuries):\n")
print(pretrend_scheme_summary)
cat("\n")

# =============================================================================
# 6. POWER ANALYSIS PER SCHEME
# =============================================================================

cat("================================================================\n")
cat("6. POWER ANALYSIS PER SCHEME\n")
cat("================================================================\n\n")

# Strategy: for each scheme, compute the residual SD from a road + quarter FE
# regression on control roads in the pre-treatment period. Then compute MDE
# using the DiD SE formula with clustering adjustment.

power_table <- map_df(schemes_all, function(s) {
  d <- panel %>% filter(scheme == s)
  d_pre <- d %>% filter(quarter_year < ref_start)

  n_treated  <- n_distinct(d$identifier[d$treat_group == 1])
  n_control  <- n_distinct(d$identifier[d$treat_group == 0])
  n_pre_qtrs <- n_distinct(d_pre$quarter_year[d_pre$treat_group == 1])
  n_post_qtrs <- n_distinct(d$quarter_year[d$treat_group == 1 &
                                             d$quarter_year >= d$ref_start])
  n_total_qtrs <- n_distinct(d$quarter_year)

  # Pre-treatment mean for treated group
  treated_pre_mean <- mean(
    d_pre$total_pkm[d_pre$treat_group == 1], na.rm = TRUE
  )

  # Residual SD from FE regression on pre-treatment data
  mod_pre <- tryCatch(
    feols(total_pkm ~ 1 | identifier + quarter_year,
          data = d_pre, cluster = ~OA),
    error = function(e) NULL
  )

  resid_sd <- if (!is.null(mod_pre)) sd(resid(mod_pre)) else NA_real_

  # ICC at OA level (pre-treatment, controls)
  d_ctrl_pre <- d_pre %>% filter(treat_group == 0)
  oa_means   <- d_ctrl_pre %>%
    group_by(OA) %>%
    summarise(oa_mean = mean(total_pkm, na.rm = TRUE), .groups = "drop")
  between_var <- var(oa_means$oa_mean, na.rm = TRUE)
  total_var   <- var(d_ctrl_pre$total_pkm, na.rm = TRUE)
  icc <- if (!is.na(total_var) && total_var > 0) between_var / total_var else NA_real_

  n_oa_treated <- n_distinct(d$OA[d$treat_group == 1])
  n_oa_control <- n_distinct(d$OA[d$treat_group == 0])

  # Design effect for clustering
  roads_per_oa <- n_treated / max(n_oa_treated, 1)
  deff <- 1 + (roads_per_oa - 1) * icc

  # Effective N
  n_eff_treated <- n_treated / deff
  n_eff_control <- n_control / deff
  n_eff <- (n_eff_treated * n_eff_control) / (n_eff_treated + n_eff_control)

  # SE of DiD estimator (simplified)
  se_did <- if (!is.na(resid_sd))
    resid_sd * sqrt(1 / (n_eff * n_post_qtrs) + 1 / (n_eff * n_pre_qtrs))
  else NA_real_

  # MDE at 80% power, alpha = 0.05 (two-sided)
  mde <- se_did * (qnorm(0.975) + qnorm(0.80))

  tibble(
    scheme          = s,
    n_treated       = n_treated,
    n_control       = n_control,
    n_oa_treated    = n_oa_treated,
    n_oa_control    = n_oa_control,
    n_pre_qtrs      = n_pre_qtrs,
    n_post_qtrs     = n_post_qtrs,
    treated_pre_mean = round(treated_pre_mean, 4),
    resid_sd        = round(resid_sd, 4),
    ICC             = round(icc, 4),
    deff            = round(deff, 2),
    se_did          = round(se_did, 5),
    MDE             = round(mde, 4),
    MDE_pct         = round(100 * mde / treated_pre_mean, 1)
  )
})

cat("Power analysis — minimum detectable effect (80% power, alpha = 0.05):\n")
print(power_table %>% select(scheme, n_treated, n_control, n_oa_treated,
                              resid_sd, ICC, MDE, MDE_pct))
cat("\n")

# MDE bar chart
p_mde <- ggplot(
  power_table,
  aes(y = fct_reorder(scheme, MDE_pct), x = MDE_pct)
) +
  geom_col(fill = "#2E6FAB", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = paste0(MDE_pct, "%")),
            hjust = -0.15, size = 3.5, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Minimum detectable effect by scheme (80% power, alpha = 0.05)",
    subtitle = "MDE as % of pre-treatment treated mean total injuries/km",
    x = "MDE (% of pre-treatment mean)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

save_fig(p_mde, "fig_MDE_by_scheme.png", width = 10, height = 7)

# Sample size summary plot
size_data <- power_table %>%
  select(scheme, Treated = n_treated, Control = n_control) %>%
  pivot_longer(-scheme, names_to = "group", values_to = "n") %>%
  mutate(scheme = fct_reorder(scheme, n, .fun = sum))

p_size <- ggplot(size_data, aes(x = n, y = scheme, fill = group)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Treated" = "#c0392b", "Control" = "#2980b9")) +
  labs(
    title = "Sample sizes per scheme (road level)",
    x = "Number of roads", y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

save_fig(p_size, "fig_sample_sizes_by_scheme.png", width = 10, height = 7)

# =============================================================================
# 7. OUTCOME SPARSITY & DISTRIBUTIONAL CHECKS
# =============================================================================

cat("================================================================\n")
cat("7. OUTCOME SPARSITY & DISTRIBUTIONAL CHECKS\n")
cat("================================================================\n\n")

# Zero-inflation by scheme × group
zero_rates <- panel %>%
  group_by(scheme, group) %>%
  summarise(
    n_obs          = n(),
    pct_zero_total = round(mean(total_inj_adj_All == 0) * 100, 1),
    pct_zero_KSI   = round(mean(KSI_adj_All == 0) * 100, 1),
    pct_zero_slight = round(mean(Slight_adj_All == 0) * 100, 1),
    .groups = "drop"
  )

cat("Zero-inflation rates (% of road-quarter obs with zero injuries):\n")
print(zero_rates, n = Inf)
cat("\n")

# Variance-to-mean ratio (overdispersion check)
vmt <- panel_pre %>%
  group_by(scheme, group) %>%
  summarise(
    mean_total = round(mean(total_pkm, na.rm = TRUE), 4),
    var_total  = round(var(total_pkm,  na.rm = TRUE), 4),
    vmr_total  = round(var(total_pkm,  na.rm = TRUE) /
                          mean(total_pkm, na.rm = TRUE), 2),
    mean_KSI   = round(mean(KSI_pkm, na.rm = TRUE), 4),
    var_KSI    = round(var(KSI_pkm,  na.rm = TRUE), 4),
    vmr_KSI    = round(var(KSI_pkm,  na.rm = TRUE) /
                          mean(KSI_pkm, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Variance-to-mean ratio (VMR) — pre-treatment:\n")
cat("  VMR >> 1 indicates overdispersion\n")
print(vmt %>% select(scheme, group, mean_total, vmr_total, mean_KSI, vmr_KSI),
      n = Inf)
cat("\n")

# Zero rate heatmap
p_zero <- zero_rates %>%
  pivot_longer(starts_with("pct_zero"), names_to = "outcome",
               values_to = "pct_zero") %>%
  mutate(
    outcome = case_when(
      outcome == "pct_zero_total"  ~ "Total",
      outcome == "pct_zero_KSI"   ~ "KSI",
      outcome == "pct_zero_slight" ~ "Slight"
    ),
    outcome = factor(outcome, levels = c("Total", "KSI", "Slight"))
  ) %>%
  ggplot(aes(x = group, y = scheme, fill = pct_zero)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(pct_zero, "%")), size = 3.5) +
  facet_wrap(~outcome) +
  scale_fill_gradient(low = "#eef7ee", high = "#d73027",
                      name = "% zero") +
  labs(
    title = "Zero-inflation rates by scheme, group, and outcome",
    subtitle = "% of road-quarter observations with zero injuries",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

save_fig(p_zero, "fig_zero_inflation_heatmap.png", width = 14, height = 7)

# =============================================================================
# 8. SERIAL CORRELATION CHECK
# =============================================================================

cat("================================================================\n")
cat("8. SERIAL CORRELATION CHECK\n")
cat("================================================================\n\n")

# Estimate TWFE on pre-treatment data, then check residual autocorrelation
# using a regression of residuals on lagged residuals (panel AR(1) test)

mod_pre_all <- feols(
  total_pkm ~ treat_group | identifier + quarter_year,
  data = panel_pre
)

panel_pre_resid <- panel_pre %>%
  select(identifier, quarter_year) %>%
  mutate(resid = resid(mod_pre_all)) %>%
  arrange(identifier, quarter_year) %>%
  group_by(identifier) %>%
  mutate(resid_lag = dplyr::lag(resid)) %>%
  ungroup() %>%
  filter(!is.na(resid_lag))

ar1_mod <- lm(resid ~ resid_lag, data = panel_pre_resid)
ar1_coef <- coef(ar1_mod)["resid_lag"]
ar1_pval <- summary(ar1_mod)$coefficients["resid_lag", "Pr(>|t|)"]

cat("AR(1) coefficient on FE residuals:", round(ar1_coef, 4), "\n")
cat("p-value:", format.pval(ar1_pval, digits = 4), "\n")
if (ar1_pval < 0.05) {
  cat(">> Serial correlation DETECTED — use cluster-robust SEs at minimum.\n")
  cat("   Consider Newey-West SEs or block bootstrap.\n")
} else {
  cat(">> No significant serial correlation at 5% level.\n")
}
cat("\n")

# Per-scheme AR(1)
ar1_by_scheme <- map_df(schemes_all, function(s) {
  d <- panel_pre %>% filter(scheme == s)
  mod <- tryCatch(
    feols(total_pkm ~ treat_group | identifier + quarter_year, data = d),
    error = function(e) NULL
  )
  if (is.null(mod)) return(tibble(scheme = s, ar1 = NA, p = NA))
  r <- d %>%
    select(identifier, quarter_year) %>%
    mutate(resid = resid(mod)) %>%
    arrange(identifier, quarter_year) %>%
    group_by(identifier) %>%
    mutate(resid_lag = dplyr::lag(resid)) %>%
    ungroup() %>%
    filter(!is.na(resid_lag))
  m <- lm(resid ~ resid_lag, data = r)
  tibble(
    scheme = s,
    ar1    = round(coef(m)["resid_lag"], 4),
    p      = round(summary(m)$coefficients["resid_lag", "Pr(>|t|)"], 4)
  )
})

cat("AR(1) by scheme:\n")
print(ar1_by_scheme)
cat("\n")

# =============================================================================
# SAVE COMBINED DIAGNOSTICS TABLE
# =============================================================================

cat("================================================================\n")
cat("SAVING DIAGNOSTICS SUMMARY\n")
cat("================================================================\n\n")

diagnostics_summary <- power_table %>%
  left_join(
    pretrend_scheme_summary %>% select(scheme, pretrend_F, pretrend_p),
    by = "scheme"
  ) %>%
  left_join(ar1_by_scheme %>% rename(ar1_coef = ar1, ar1_p = p), by = "scheme") %>%
  left_join(
    zero_rates %>%
      filter(group == "Treated") %>%
      select(scheme, pct_zero_total, pct_zero_KSI),
    by = "scheme"
  )

write_csv(diagnostics_summary,
          file.path(outdir, "panel_diagnostics_summary.csv"))

cat("Summary table saved to:\n  ", file.path(outdir, "panel_diagnostics_summary.csv"), "\n\n")

cat("Diagnostics table:\n")
diagnostics_summary %>%
  select(scheme, n_treated, n_control, MDE_pct,
         pretrend_p, ar1_coef, pct_zero_total) %>%
  print()

cat("\n================================================================\n")
cat("ALL PANEL DIAGNOSTICS COMPLETE\n")
cat("================================================================\n")
cat("Figures saved to:", outdir, "\n")
