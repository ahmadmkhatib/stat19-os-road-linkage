# =============================================================================
# DiD — POOLED MATCHING — STRATIFIED BY INJURY MODE
# =============================================================================
#
# Runs the full analysis pipeline for three outcomes:
#   1. All injuries (total_inj_adj_All)
#   2. Vehicle injuries — car/van + other motor (total_inj_adj_Vehicle)
#   3. Active travel injuries — cyclist + pedestrian (total_inj_adj_ActiveTravel)
#
# Primary specification: Callaway & Sant'Anna (2021) doubly-robust estimator
# with covariate adjustment for residual post-matching imbalance.
#
# For each outcome:
#   0  Injury summary statistics
#   1  Descriptive DiD
#   2  Trend plot
#   3  TWFE DiD (road + quarter FE)
#   4  Main analysis: covariate-adjusted C&S (overall, per-scheme, event studies)
#   5  Sensitivity: unadjusted vs adjusted comparison
#   6  Sensitivity: injury-prone roads only (adjusted)
#   7  Sensitivity: PPML TWFE (pre-filtered sample)
#   8  Sensitivity: minimum road length >= 10m (adjusted C&S + PPML)
#   9  Sensitivity: exclude 2024 / STATS19 specification change (adj C&S + PPML)
#  10  Robustness summary table
#
# Input:  road_panel_matched_pooled.parquet
# Output: output/pooled/{All,Vehicle,ActiveTravel}/
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(fixest)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

run_cs <- function(data, id = "uid_int", xformla = NULL, label = "") {
  if (nchar(label) > 0) cat("  ", label, "\n")
  att <- att_gt(
    yname = "outcome_pkm", tname = "qtr_int", idname = id, gname = "g",
    data = data, control_group = "nevertreated",
    xformla = xformla,
    clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE
  )
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  list(
    att_obj = att, agg_obj = agg,
    result = tibble(
      att = agg$overall.att, se = agg$overall.se,
      ci_lo = agg$overall.att - 1.96 * agg$overall.se,
      ci_hi = agg$overall.att + 1.96 * agg$overall.se,
      pval = 2 * pnorm(-abs(agg$overall.att / agg$overall.se))
    )
  )
}

run_cs_by_scheme <- function(data, schemes, id = "uid_int", xformla = NULL) {
  map_df(schemes, function(s) {
    cat("  Scheme:", s, "\n")
    d <- data %>% filter(scheme == s) %>%
      mutate(!!id := as.integer(factor(.data[[id]])), g = as.numeric(g))
    fit <- tryCatch(run_cs(d, id = id, xformla = xformla),
                    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL })
    if (is.null(fit)) return(tibble(scheme = s, att = NA_real_, se = NA_real_,
                                    ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_))
    bind_cols(tibble(scheme = s), fit$result)
  }) %>%
    mutate(sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                           pval < 0.1 ~ "*", TRUE ~ ""),
           sig_label = paste0(sprintf("%.4f", att), sig))
}

es_from_att <- function(att_obj, trim = c(-12, 12)) {
  agg_dyn <- aggte(att_obj, type = "dynamic", na.rm = TRUE)
  tibble(event_time = agg_dyn$egt, att = agg_dyn$att.egt, se = agg_dyn$se.egt) %>%
    mutate(ci_lo = att - 1.96 * se, ci_hi = att + 1.96 * se) %>%
    filter(event_time >= trim[1], event_time <= trim[2])
}

plot_es <- function(es_df, title, subtitle = "") {
  ggplot(es_df, aes(x = event_time, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
    labs(title = title, subtitle = subtitle,
         x = "Quarters relative to CAZ implementation",
         y = "ATT (injuries / road-km)") +
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
}

plot_coefs <- function(res_df, title, subtitle = "") {
  ggplot(res_df, aes(x = att, y = fct_reorder(scheme, att))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                  width = 0.25, linewidth = 0.7, orientation = "y") +
    geom_point(size = 3, colour = "#2E6FAB") +
    geom_text(aes(label = sig_label), hjust = -0.15, size = 3) +
    labs(title = title, subtitle = subtitle,
         x = "ATT (injuries / road-km)", y = NULL,
         caption = "Error bars: 95% CI  |  * p<0.1  ** p<0.05  *** p<0.01") +
    theme_minimal(base_size = 12) + theme(panel.grid.major.y = element_blank())
}

run_ppml_one <- function(data) {
  m <- fepois(outcome_count ~ treated:post_int | ppml_id + qtr_num,
              data = data, cluster = ~OA, offset = ~log_length)
  b  <- coef(m)["treated:post_int"]
  se <- sqrt(vcov(m)["treated:post_int", "treated:post_int"])
  list(model = m, result = tibble(
    pct_change = (exp(b)-1)*100, ci_lo_pct = (exp(b-1.96*se)-1)*100,
    ci_hi_pct = (exp(b+1.96*se)-1)*100, pval = 2*pnorm(-abs(b/se)),
    n_obs = nobs(m)
  ))
}

run_ppml <- function(data, schemes) {
  overall_fit <- run_ppml_one(data)
  overall <- overall_fit$result %>% mutate(scheme = "Overall", .before = 1)
  per_scheme <- map_df(schemes, function(s) {
    cat("  PPML:", s, "\n")
    d <- data %>% filter(scheme == s) %>%
      mutate(ppml_id = as.integer(factor(ppml_id)), qtr_num = as.numeric(quarter_year))
    fit <- tryCatch(run_ppml_one(d),
                    error = function(e) { cat("    FAILED\n"); NULL })
    if (is.null(fit)) return(tibble(scheme = s, pct_change = NA_real_,
      ci_lo_pct = NA_real_, ci_hi_pct = NA_real_, pval = NA_real_, n_obs = NA_integer_))
    fit$result %>% mutate(scheme = s, .before = 1)
  })
  result <- bind_rows(overall, per_scheme) %>%
    mutate(sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                           pval < 0.1 ~ "*", TRUE ~ ""),
           sig_label = paste0(sprintf("%+.1f%%", pct_change), sig))
  list(model = overall_fit$model, result = result)
}

fmt_cell <- function(pct, pval) {
  stars <- case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                     pval < 0.1 ~ "*", TRUE ~ "")
  sprintf("%+.1f%s", pct, stars)
}

# =============================================================================
# LOAD DATA
# =============================================================================

road_panel_matched <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(quarter_year = as.yearqtr(quarter_year),
         caz_start_q  = as.yearqtr(caz_start_q))

cat("Rows:", nrow(road_panel_matched),
    "| Roads:", n_distinct(road_panel_matched$identifier),
    "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n")

# Scheme timing
scheme_timing <- road_panel_matched %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

# OA-level covariates for doubly-robust adjustment
matched_full_cov <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))
matched_covars <- matched_full_cov %>%
  mutate(log1p_pop_density = log1p(pmax(pop_density, 0)),
         log1p_road_density_m_km2 = log1p(pmax(road_density_m_km2, 0))) %>%
  select(OA, log1p_pop_density, IMD, log1p_road_density_m_km2) %>%
  distinct(OA, .keep_all = TRUE)
rm(matched_full_cov)

xformla_adj <- ~ log1p_pop_density + IMD + log1p_road_density_m_km2

# Outcome columns to loop over
outcomes <- c("All", "Vehicle", "ActiveTravel")
outcome_labels <- c(All = "Total injuries", Vehicle = "Vehicle injuries (car/van + other)",
                    ActiveTravel = "Active travel injuries (cyclist + pedestrian)")

# =============================================================================
# MAIN ANALYSIS FUNCTION — runs Steps 0–10 for one outcome
# =============================================================================

run_full_analysis <- function(outcome_name, base_panel, covars, xformla,
                              scheme_timing_df, outdir_base) {

  outcome_col   <- paste0("total_inj_adj_", outcome_name)
  outcome_label <- outcome_labels[outcome_name]
  outdir <- file.path(outdir_base, outcome_name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  cat("\n")
  cat(paste(rep("#", 72), collapse = ""), "\n")
  cat("##  OUTCOME:", outcome_label, "\n")
  cat(paste(rep("#", 72), collapse = ""), "\n\n")

  # --- Build panel ---

  panel <- base_panel %>%
    mutate(
      road_length_km = length / 1000,
      road_length_km = if_else(road_length_km <= 0 | is.na(road_length_km),
                               1e-6, road_length_km),
      outcome_count = .data[[outcome_col]],
      outcome_pkm   = outcome_count / road_length_km,
      group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
    ) %>%
    left_join(scheme_timing_df %>% rename(ref_start = caz_start_q), by = "scheme") %>%
    mutate(
      post_flag = as.integer(quarter_year >= ref_start),
      # COVID period flag: 2020 Q2 – 2021 Q2
      covid_period = as.integer(quarter_year >= as.yearqtr("2020 Q2") &
                                quarter_year <= as.yearqtr("2021 Q2"))
    )

  min_qtr <- min(as.numeric(panel$quarter_year), na.rm = TRUE)

  panel_did <- panel %>%
    mutate(
      uid = paste0(panel_id, "_", scheme),
      uid_int = as.integer(factor(uid)),
      qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
      g = case_when(
        treat_group == 1 & !is.na(caz_start_q) ~
          as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
        TRUE ~ 0
      )
    ) %>%
    filter(!is.na(outcome_pkm)) %>%
    left_join(covars, by = "OA")

  schemes_all <- sort(unique(panel_did$scheme))

  # COVID qtr_int range for exclusion sensitivity
  covid_start_int <- as.integer(round((as.numeric(as.yearqtr("2020 Q2")) - min_qtr) * 4)) + 1L
  covid_end_int   <- as.integer(round((as.numeric(as.yearqtr("2021 Q2")) - min_qtr) * 4)) + 1L

  # =========================================================================
  # — SUMMARY STATISTICS
  # =========================================================================

  cat("\n== STEP 0 — SUMMARY STATISTICS ==\n\n")

  panel %>%
    group_by(group) %>%
    summarise(
      n_obs    = format(n(), big.mark = ","),
      mean     = round(mean(outcome_count, na.rm = TRUE), 4),
      pct_zero = round(100 * mean(outcome_count == 0, na.rm = TRUE), 1),
      .groups  = "drop"
    ) %>% print()

  # =========================================================================
  # — DESCRIPTIVE DiD
  # =========================================================================

  cat("\n== STEP 1 — DESCRIPTIVE DiD ==\n\n")

  ss <- panel %>%
    group_by(group, post_flag) %>%
    summarise(mean_pkm = round(mean(outcome_pkm, na.rm = TRUE), 4),
              .groups = "drop") %>%
    mutate(period = if_else(post_flag == 1, "Post", "Pre"))
  print(ss)

  dt <- ss %>% pivot_wider(id_cols = group, names_from = period, values_from = mean_pkm) %>%
    mutate(change = Post - Pre)
  did_est <- dt$change[dt$group == "CAZ roads"] - dt$change[dt$group == "Matched controls"]
  cat("\nDescriptive DiD:", round(did_est, 5), "\n")

  # =========================================================================
  # — TREND PLOT
  # =========================================================================

  pd <- panel %>%
    group_by(group, quarter_year) %>%
    summarise(mean_pkm = mean(outcome_pkm, na.rm = TRUE), .groups = "drop")

  ggplot(pd, aes(x = as.Date(quarter_year), y = mean_pkm, colour = group)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = as.Date(median(scheme_timing_df$caz_start_q)),
               linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("CAZ roads" = "#E74C3C", "Matched controls" = "#2C3E50")) +
    labs(title = paste0("Trend: ", outcome_label),
         x = NULL, y = "Mean injuries/road-km", colour = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(outdir, "trend_plot.png"), width = 10, height = 6, dpi = 300)

  # =========================================================================
  # — COVID DIFFERENTIAL CHECK
  # =========================================================================
  # COVID lockdowns (2020 Q2 – 2021 Q2) may have hit CAZ city-centre roads
  # harder than dispersed control areas. Check differential impact.

  cat("\n== STEP 2b — COVID DIFFERENTIAL CHECK ==\n\n")

  covid_check <- panel %>%
    filter(quarter_year >= as.yearqtr("2019 Q1"),
           quarter_year <= as.yearqtr("2022 Q4")) %>%
    group_by(group, quarter_year) %>%
    summarise(total_inj = sum(outcome_count, na.rm = TRUE), .groups = "drop")

  # Index to 2019 Q1
  covid_base <- covid_check %>%
    filter(quarter_year == as.yearqtr("2019 Q1")) %>%
    select(group, base = total_inj)

  covid_idx <- covid_check %>%
    left_join(covid_base, by = "group") %>%
    mutate(index = round(total_inj / base * 100, 1))

  covid_wide <- covid_idx %>%
    select(quarter_year, group, index) %>%
    pivot_wider(names_from = group, values_from = index) %>%
    mutate(gap = `CAZ roads` - `Matched controls`)

  cat("COVID differential impact (2019 Q1 = 100):\n")
  print(covid_wide, n = 16)

  # Plot
  ggplot(covid_idx, aes(x = as.Date(quarter_year), y = index, colour = group)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    geom_hline(yintercept = 100, linetype = "dotted", colour = "grey60") +
    annotate("rect", xmin = as.Date("2020-04-01"), xmax = as.Date("2021-07-01"),
             ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "red") +
    annotate("text", x = as.Date("2020-10-01"), y = max(covid_idx$index) * 0.95,
             label = "COVID", size = 3.5, colour = "grey40") +
    scale_colour_manual(values = c("CAZ roads" = "#E74C3C", "Matched controls" = "#2C3E50")) +
    labs(title = paste0("COVID differential: ", outcome_label),
         x = NULL, y = "Index (2019 Q1 = 100)", colour = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(outdir, "covid_differential_check.png"),
         width = 10, height = 6, dpi = 300)

  # =========================================================================
  # STEP 3 — TWFE (with COVID interaction)
  # =========================================================================
  # The treated × covid_period interaction absorbs the differential COVID
  # impact as a nuisance parameter, so the treated × post coefficient
  # captures the CAZ effect net of COVID disruption.

  cat("\n== STEP 3 — TWFE (with COVID interaction) ==\n\n")

  panel_reg <- panel %>%
    mutate(treated = as.integer(treat_group), post_int = post_flag,
           qtr_num = as.numeric(quarter_year))

  m_twfe <- feols(outcome_pkm ~ treated:post_int + treated:covid_period |
                    identifier + qtr_num,
                  data = panel_reg, cluster = ~OA)

  etable(m_twfe, headers = outcome_label,
         dict = c("treated:post_int" = "Treated x Post",
                   "treated:covid_period" = "Treated x COVID"),
         digits = 5, digits.stats = 3)

  # Also run without COVID interaction for comparison
  m_twfe_nocov <- feols(outcome_pkm ~ treated:post_int | identifier + qtr_num,
                        data = panel_reg, cluster = ~OA)

  cat("\n--- TWFE comparison: with vs without COVID interaction ---\n")
  etable(m_twfe_nocov, m_twfe,
         headers = c("Without COVID", "With COVID interaction"),
         dict = c("treated:post_int" = "Treated x Post",
                   "treated:covid_period" = "Treated x COVID"),
         digits = 5, digits.stats = 3)

  # =========================================================================
  # STEP 4 — MAIN: COVARIATE-ADJUSTED C&S
  # =========================================================================

  cat("\n== STEP 4 — MAIN: ADJUSTED C&S ==\n\n")

  main_fit <- run_cs(panel_did, xformla = xformla, label = "Overall (adjusted)")
  agg_main <- main_fit$agg_obj
  summary(agg_main)

  es_main <- es_from_att(main_fit$att_obj)
  plot_es(es_main, paste0("Event study: ", outcome_label, " (adjusted)"))
  ggsave(file.path(outdir, "event_study_adjusted.png"), width = 10, height = 7, dpi = 300)

  results_main <- run_cs_by_scheme(panel_did, schemes_all, xformla = xformla)
  print(results_main)

  plot_coefs(results_main, paste0("ATT by scheme: ", outcome_label, " (adjusted)"))
  ggsave(file.path(outdir, "coef_plot_adjusted.png"), width = 10, height = 7, dpi = 300)

  # Per-scheme event studies
  es_by_scheme <- map_df(schemes_all, function(s) {
    d <- panel_did %>% filter(scheme == s) %>%
      mutate(uid_int = as.integer(factor(uid)), g = as.numeric(g))
    att <- tryCatch(
      att_gt(yname = "outcome_pkm", tname = "qtr_int", idname = "uid_int", gname = "g",
             data = d, control_group = "nevertreated", xformla = xformla,
             clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = TRUE),
      error = function(e) NULL)
    if (is.null(att)) return(NULL)
    tryCatch(es_from_att(att) %>% mutate(scheme = s), error = function(e) NULL)
  })

  if (nrow(es_by_scheme) > 0) {
    ggplot(es_by_scheme, aes(x = event_time, y = att)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
      scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
      facet_wrap(~scheme, scales = "free_y", ncol = 2) +
      labs(title = paste0("Event study by scheme: ", outcome_label),
           x = "Quarters relative to CAZ", y = "ATT (injuries/road-km)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
    ggsave(file.path(outdir, "event_study_by_scheme.png"), width = 12, height = 10, dpi = 300)
  }

  # =========================================================================
  # STEP 5 — SENSITIVITY: UNADJUSTED VS ADJUSTED
  # =========================================================================

  cat("\n== STEP 5 — UNADJUSTED VS ADJUSTED ==\n\n")

  unadj_fit <- run_cs(panel_did, xformla = NULL, label = "Overall (unadjusted)")
  agg_unadj <- unadj_fit$agg_obj
  results_unadj <- run_cs_by_scheme(panel_did, schemes_all, xformla = NULL)

  cat(sprintf("  Unadjusted: %.4f (SE %.4f)\n  Adjusted:   %.4f (SE %.4f)\n",
              agg_unadj$overall.att, agg_unadj$overall.se,
              agg_main$overall.att, agg_main$overall.se))

  # =========================================================================
  # STEP 5b — SENSITIVITY: EXCLUDE COVID (2020 Q2 – 2021 Q2)
  # =========================================================================
  # COVID lockdowns hit CAZ city-centre roads differentially (deeper drop,
  # faster rebound). Bath (2021 Q1) and Birmingham (2021 Q2) started
  # treatment during COVID, contaminating their post-treatment estimates.
  # We drop 2020 Q2 – 2021 Q2 and re-estimate. panel = FALSE because the
  # temporal gap creates an unbalanced panel.

  cat("\n== STEP 5b — EXCLUDE COVID (2020 Q2 – 2021 Q2) ==\n\n")

  panel_nocovid <- panel_did %>%
    filter(!(qtr_int >= covid_start_int & qtr_int <= covid_end_int)) %>%
    mutate(uid_nocov = as.integer(factor(uid)))

  cat("Excluded qtr_int", covid_start_int, "to", covid_end_int,
      "| Retained", n_distinct(panel_nocovid$qtr_int), "of",
      n_distinct(panel_did$qtr_int), "quarters\n")

  fit_nocovid <- tryCatch({
    att_nc <- att_gt(
      yname = "outcome_pkm", tname = "qtr_int", idname = "uid_nocov", gname = "g",
      data = panel_nocovid, control_group = "nevertreated",
      xformla = xformla,
      clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = FALSE
    )
    agg_nc <- aggte(att_nc, type = "simple", na.rm = TRUE)
    list(att_obj = att_nc, agg_obj = agg_nc,
         result = tibble(att = agg_nc$overall.att, se = agg_nc$overall.se,
                         pval = 2 * pnorm(-abs(agg_nc$overall.att / agg_nc$overall.se))))
  }, error = function(e) { cat("  Overall FAILED:", e$message, "\n"); NULL })

  if (!is.null(fit_nocovid)) {
    agg_nocovid <- fit_nocovid$agg_obj
    cat("\n=== Overall ATT (excl COVID, adjusted) ===\n")
    summary(agg_nocovid)
  } else {
    agg_nocovid <- NULL
  }

  results_nocovid_safe <- map_df(schemes_all, function(s) {
    cat("  Scheme:", s, "\n")
    d <- panel_nocovid %>% filter(scheme == s) %>%
      mutate(uid_nocov = as.integer(factor(uid_nocov)), g = as.numeric(g))
    att <- tryCatch(
      att_gt(yname = "outcome_pkm", tname = "qtr_int", idname = "uid_nocov", gname = "g",
             data = d, control_group = "nevertreated", xformla = xformla,
             clustervars = "OA", bstrap = TRUE, anticipation = 0, panel = FALSE),
      error = function(e) { cat("    FAILED:", e$message, "\n"); NULL })
    if (is.null(att)) return(tibble(scheme = s, att = NA_real_, se = NA_real_,
                                    ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_))
    agg <- aggte(att, type = "simple", na.rm = TRUE)
    tibble(scheme = s, att = agg$overall.att, se = agg$overall.se,
           ci_lo = agg$overall.att - 1.96 * agg$overall.se,
           ci_hi = agg$overall.att + 1.96 * agg$overall.se,
           pval = 2 * pnorm(-abs(agg$overall.att / agg$overall.se)))
  }) %>%
    mutate(sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                           pval < 0.1 ~ "*", TRUE ~ ""),
           sig_label = paste0(sprintf("%.4f", att), sig))

  cat("\n=== Per-scheme (excl COVID, adjusted) ===\n")
  print(results_nocovid_safe)

  # =========================================================================
  # STEP 6 — SENSITIVITY: INJURY-PRONE ROADS (ADJUSTED)
  # =========================================================================

  cat("\n== STEP 6 — INJURY-PRONE ROADS (ADJUSTED) ==\n\n")

  scheme_starts <- panel_did %>% filter(g > 0) %>% distinct(scheme, g) %>% rename(scheme_start = g)

  pre_inj_ids <- panel_did %>%
    left_join(scheme_starts, by = "scheme") %>%
    mutate(ref_g = if_else(g > 0, g, scheme_start)) %>%
    filter(qtr_int < ref_g) %>%
    group_by(uid_int) %>%
    summarise(any_pre = any(outcome_pkm > 0, na.rm = TRUE), .groups = "drop") %>%
    filter(any_pre) %>% pull(uid_int)

  panel_filt <- panel_did %>% filter(uid_int %in% pre_inj_ids) %>%
    mutate(uid_filt = as.integer(factor(uid_int)))

  cat("Retained", n_distinct(panel_filt$uid_filt), "of",
      n_distinct(panel_did$uid_int), "units\n")

  filt_fit <- run_cs(panel_filt, id = "uid_filt", xformla = xformla,
                     label = "Filtered (adjusted)")
  agg_filt <- filt_fit$agg_obj
  summary(agg_filt)

  es_filt <- es_from_att(filt_fit$att_obj, trim = c(-12, 15))
  plot_es(es_filt, paste0("Event study: injury-prone roads — ", outcome_label))
  ggsave(file.path(outdir, "event_study_filtered.png"), width = 10, height = 7, dpi = 300)

  results_filt <- run_cs_by_scheme(panel_filt, schemes_all,
                                   id = "uid_filt", xformla = xformla)

  # =========================================================================
  # STEP 7 — SENSITIVITY: PPML TWFE (PRE-FILTERED)
  # =========================================================================

  cat("\n== STEP 7 — PPML (PRE-FILTERED) ==\n\n")

  panel_ppml <- panel_filt %>%
    mutate(treated = as.integer(treat_group),
           post_int = as.integer(quarter_year >= ref_start),
           qtr_num = as.numeric(quarter_year),
           log_length = log(road_length_km),
           ppml_id = as.integer(factor(uid_filt)))

  ppml_fit <- run_ppml(panel_ppml, schemes_all)
  print(ppml_fit$result)

  # =========================================================================
  # STEP 8 — SENSITIVITY: >= 10m ROADS (ADJUSTED C&S + PPML)
  # =========================================================================

  cat("\n== STEP 8 — ROADS >= 10m ==\n\n")

  panel_10m <- panel_did %>% filter(road_length_km >= 0.01) %>%
    mutate(uid_10m = as.integer(factor(uid)))

  fit_10m <- run_cs(panel_10m, id = "uid_10m", xformla = xformla,
                    label = "Overall >= 10m (adjusted)")
  agg_10m <- fit_10m$agg_obj
  summary(agg_10m)

  results_10m <- run_cs_by_scheme(panel_10m, schemes_all,
                                  id = "uid_10m", xformla = xformla)

  panel_ppml_10m <- panel_ppml %>% filter(road_length_km >= 0.01) %>%
    mutate(ppml_id = as.integer(factor(ppml_id)))
  ppml_10m <- run_ppml(panel_ppml_10m, schemes_all)

  # =========================================================================
  # STEP 9 — SENSITIVITY: EXCLUDE 2024 (ADJUSTED C&S + PPML)
  # =========================================================================

  cat("\n== STEP 9 — EXCLUDE 2024 ==\n\n")

  panel_no24 <- panel_did %>% filter(qtr_int <= 36) %>%
    mutate(uid_no24 = as.integer(factor(uid)))

  fit_no24 <- run_cs(panel_no24, id = "uid_no24", xformla = xformla,
                     label = "Overall excl 2024 (adjusted)")
  agg_no24 <- fit_no24$agg_obj
  summary(agg_no24)

  results_no24 <- run_cs_by_scheme(panel_no24, schemes_all,
                                   id = "uid_no24", xformla = xformla)

  panel_ppml_no24 <- panel_ppml %>% filter(qtr_num < 2024) %>%
    mutate(ppml_id = as.integer(factor(ppml_id)))
  ppml_no24 <- run_ppml(panel_ppml_no24, schemes_all)

  # (COVID exclusion results computed in Step 5b above)
  # STEP 10 — ROBUSTNESS SUMMARY TABLE
  # =========================================================================

  cat("\n== STEP 10 — ROBUSTNESS SUMMARY ==\n\n")

  pre_means_full <- panel_did %>%
    filter(treat_group == 1, qtr_int < g) %>%
    group_by(scheme) %>%
    summarise(pm = mean(outcome_pkm, na.rm = TRUE), .groups = "drop")

  pre_means_10m <- panel_10m %>%
    filter(treat_group == 1, qtr_int < g) %>%
    group_by(scheme) %>%
    summarise(pm = mean(outcome_pkm, na.rm = TRUE), .groups = "drop")

  rob <- tibble(scheme = schemes_all) %>%
    left_join(results_main %>% left_join(pre_means_full, by = "scheme") %>%
      transmute(scheme, cs_adj = att/pm*100, cs_adj_p = pval), by = "scheme") %>%
    left_join(results_unadj %>% left_join(pre_means_full, by = "scheme") %>%
      transmute(scheme, cs_unadj = att/pm*100, cs_unadj_p = pval), by = "scheme") %>%
    left_join(results_10m %>% left_join(pre_means_10m, by = "scheme") %>%
      transmute(scheme, cs_10m = att/pm*100, cs_10m_p = pval), by = "scheme") %>%
    left_join(results_no24 %>% left_join(pre_means_full, by = "scheme") %>%
      transmute(scheme, cs_no24 = att/pm*100, cs_no24_p = pval), by = "scheme") %>%
    left_join(results_nocovid_safe %>% left_join(pre_means_full, by = "scheme") %>%
      transmute(scheme, cs_nocov = att/pm*100, cs_nocov_p = pval), by = "scheme") %>%
    left_join(ppml_fit$result %>% filter(scheme != "Overall") %>%
      transmute(scheme, ppml = pct_change, ppml_p = pval), by = "scheme") %>%
    left_join(ppml_10m$result %>% filter(scheme != "Overall") %>%
      transmute(scheme, ppml_10m = pct_change, ppml_10m_p = pval), by = "scheme") %>%
    left_join(ppml_no24$result %>% filter(scheme != "Overall") %>%
      transmute(scheme, ppml_no24 = pct_change, ppml_no24_p = pval), by = "scheme")

  rob_display <- rob %>%
    transmute(Scheme = scheme,
              `C&S adj` = fmt_cell(cs_adj, cs_adj_p),
              `C&S unadj` = fmt_cell(cs_unadj, cs_unadj_p),
              `C&S adj >=10m` = fmt_cell(cs_10m, cs_10m_p),
              `C&S adj no24` = fmt_cell(cs_no24, cs_no24_p),
              `C&S adj noCOV` = fmt_cell(cs_nocov, cs_nocov_p),
              PPML = fmt_cell(ppml, ppml_p),
              `PPML >=10m` = fmt_cell(ppml_10m, ppml_10m_p),
              `PPML no24` = fmt_cell(ppml_no24, ppml_no24_p))

  cat("--- Per-scheme (% change) ---\n\n")
  print(rob_display, n = 7, width = Inf)

  # Overall row
  ppml_ov     <- ppml_fit$result %>% filter(scheme == "Overall")
  ppml_10m_ov <- ppml_10m$result %>% filter(scheme == "Overall")
  ppml_n24_ov <- ppml_no24$result %>% filter(scheme == "Overall")

  # COVID overall ATT (may be NULL if failed)
  nocov_att <- if (!is.null(fit_nocovid)) agg_nocovid$overall.att else NA_real_
  nocov_se  <- if (!is.null(fit_nocovid)) agg_nocovid$overall.se else NA_real_
  nocov_p   <- if (!is.null(fit_nocovid)) {
    2 * pnorm(-abs(nocov_att / nocov_se))
  } else NA_real_

  overall_tbl <- tibble(
    Spec = c("C&S adj", "C&S unadj", "C&S adj >=10m", "C&S adj no24",
             "C&S adj noCOV",
             "PPML", "PPML >=10m", "PPML no24"),
    Estimate = c(
      sprintf("%.4f", agg_main$overall.att),
      sprintf("%.4f", agg_unadj$overall.att),
      sprintf("%.4f", agg_10m$overall.att),
      sprintf("%.4f", agg_no24$overall.att),
      sprintf("%.4f", nocov_att),
      sprintf("%+.1f%%", ppml_ov$pct_change),
      sprintf("%+.1f%%", ppml_10m_ov$pct_change),
      sprintf("%+.1f%%", ppml_n24_ov$pct_change)),
    p = round(c(
      2*pnorm(-abs(agg_main$overall.att/agg_main$overall.se)),
      2*pnorm(-abs(agg_unadj$overall.att/agg_unadj$overall.se)),
      2*pnorm(-abs(agg_10m$overall.att/agg_10m$overall.se)),
      2*pnorm(-abs(agg_no24$overall.att/agg_no24$overall.se)),
      nocov_p,
      ppml_ov$pval, ppml_10m_ov$pval, ppml_n24_ov$pval), 3)
  )

  cat("\n--- Overall ATT ---\n\n")
  print(overall_tbl, width = Inf)

  write_csv(rob_display, file.path(outdir, "robustness_by_scheme.csv"))
  write_csv(overall_tbl, file.path(outdir, "robustness_overall.csv"))

  cat("\nSaved to:", outdir, "\n")

  invisible(list(
    agg_main = agg_main, agg_unadj = agg_unadj,
    results_main = results_main, results_unadj = results_unadj,
    rob_display = rob_display, overall_tbl = overall_tbl
  ))
}

# =============================================================================
# RUN ALL THREE OUTCOMES
# =============================================================================

outdir_base <- here("output", "pooled")

all_results <- list()

for (oc in outcomes) {
  cat("\n\n")
  cat(paste(rep("=", 72), collapse = ""), "\n")
  cat(paste(rep("=", 72), collapse = ""), "\n")
  cat("  RUNNING:", outcome_labels[oc], "\n")
  cat(paste(rep("=", 72), collapse = ""), "\n")
  cat(paste(rep("=", 72), collapse = ""), "\n\n")

  all_results[[oc]] <- run_full_analysis(
    outcome_name    = oc,
    base_panel      = road_panel_matched,
    covars          = matched_covars,
    xformla         = xformla_adj,
    scheme_timing_df = scheme_timing,
    outdir_base     = outdir_base
  )
}

# =============================================================================
# CROSS-OUTCOME COMPARISON TABLE
# =============================================================================

cat("\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("CROSS-OUTCOME COMPARISON\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

cross_overall <- map_df(outcomes, function(oc) {
  all_results[[oc]]$overall_tbl %>%
    filter(Spec == "C&S adj") %>%
    mutate(Outcome = outcome_labels[oc], .before = 1)
})

cat("--- Primary specification (C&S adjusted) across outcomes ---\n\n")
print(cross_overall, width = Inf)

cross_scheme <- map_df(outcomes, function(oc) {
  all_results[[oc]]$rob_display %>%
    select(Scheme, `C&S adj`) %>%
    mutate(Outcome = oc)
}) %>%
  pivot_wider(names_from = Outcome, values_from = `C&S adj`)

cat("\n--- Per-scheme primary estimates across outcomes ---\n\n")
print(cross_scheme, width = Inf)

write_csv(cross_overall, file.path(outdir_base, "cross_outcome_overall.csv"))
write_csv(cross_scheme, file.path(outdir_base, "cross_outcome_by_scheme.csv"))

cat("\n=== ALL ANALYSES COMPLETE ===\n")
cat("Output root:", outdir_base, "\n")
cat("End of script.\n")