# =============================================================================
# DiD SENSITIVITY ANALYSES — POOLED MATCHING
# =============================================================================
#
# Sensitivity checks for the primary specification in 22-DiD_pooled.R.
# All specifications use the same aligned panel: raw injury counts,
# structural zeros dropped, notyettreated control group, no clustering.
#
# Checks:
#   5  Unadjusted vs adjusted C&S
#  5b  Exclude COVID period (2020 Q1 – 2021 Q2)
#   6  Include structural zeros (all road links)
#   7  PPML TWFE (count model with road-length offset)
#   8  Exclude 2024 / STATS19 specification change
#   9  COVID-interaction TWFE
#  10  Mode-stratified outcomes (Vehicle, ActiveTravel)
#  11  Robustness summary table
#
# Input:  road_panel_matched_pooled.parquet
#         OA_matched_full_pooled.rds
#         roads_caz_props.rds
# Output: output/pooled/sensitivity/
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(fixest)

# =============================================================================
# HELPER FUNCTIONS (aligned with script 22)
# =============================================================================

run_cs <- function(data, yname = "outcome", id = "uid_int",
                   xformla = NULL, label = "", panel = TRUE) {
  if (nchar(label) > 0) cat("  ", label, "\n")
  att <- att_gt(
    yname = yname, tname = "qtr_int", idname = id, gname = "g",
    data = data, control_group = "notyettreated",
    xformla = xformla,
    bstrap = TRUE, anticipation = 0, panel = panel
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

run_cs_by_scheme <- function(data, schemes, yname = "outcome",
                             id = "uid_int", xformla = NULL, panel = TRUE) {
  map_df(schemes, function(s) {
    cat("  Scheme:", s, "\n")
    d <- data %>% filter(scheme == s) %>%
      mutate(!!id := as.integer(factor(.data[[id]])), g = as.numeric(g))
    fit <- tryCatch(run_cs(d, yname = yname, id = id, xformla = xformla,
                           panel = panel),
                    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL })
    if (is.null(fit)) return(tibble(scheme = s, att = NA_real_, se = NA_real_,
                                    ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_))
    bind_cols(tibble(scheme = s), fit$result)
  }) %>%
    mutate(sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                           pval < 0.1 ~ "*", TRUE ~ ""),
           sig_label = paste0(sprintf("%.4f", att), sig))
}

fmt_cell <- function(est, pval) {
  stars <- case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                     pval < 0.1 ~ "*", TRUE ~ "")
  sprintf("%.4f%s", est, stars)
}

# =============================================================================
# LOAD DATA (identical to script 22)
# =============================================================================

road_panel_matched <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(quarter_year = as.yearqtr(quarter_year),
         caz_start_q  = as.yearqtr(caz_start_q))

cat("Rows:", nrow(road_panel_matched),
    "| Roads:", n_distinct(road_panel_matched$identifier),
    "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n")

# --- Majority-rule treatment quarter adjustment ---
road_caz_props <- readRDS(here("data", "processed", "roads_caz_props.rds"))

scheme_start_adj <- road_caz_props %>%
  distinct(scheme, startDt, caz_start_q) %>%
  filter(!is.na(startDt)) %>%
  mutate(
    start_date = dmy(startDt),
    raw_q      = as.yearqtr(start_date),
    q_start    = as.Date(raw_q),
    q_end      = as.Date(raw_q + 0.25) - 1,
    q_mid      = q_start + as.integer(difftime(q_end, q_start, units = "days")) / 2,
    caz_start_q_adj = if_else(start_date > q_mid, raw_q + 0.25, raw_q)
  ) %>%
  select(scheme, caz_start_q_adj)

road_panel_matched <- road_panel_matched %>%
  left_join(scheme_start_adj, by = "scheme") %>%
  mutate(caz_start_q = coalesce(caz_start_q_adj, caz_start_q)) %>%
  select(-caz_start_q_adj)

rm(road_caz_props)

scheme_timing <- road_panel_matched %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q)

# OA-level covariates
matched_full_cov <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))
matched_covars <- matched_full_cov %>%
  mutate(log1p_pop_density = log1p(pmax(pop_density, 0)),
         log1p_road_density_m_km2 = log1p(pmax(road_density_m_km2, 0))) %>%
  select(OA, log1p_pop_density, IMD, log1p_road_density_m_km2) %>%
  distinct(OA, .keep_all = TRUE)
rm(matched_full_cov)

xformla_adj <- ~ log1p_pop_density + IMD + log1p_road_density_m_km2

# =============================================================================
# BUILD MODEL PANEL (identical to script 22)
# =============================================================================

outdir <- here("output", "pooled", "sensitivity")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

min_qtr <- min(as.numeric(road_panel_matched$quarter_year), na.rm = TRUE)

model_panel <- road_panel_matched %>%
  select(panel_id, identifier, OA, quarter_year, scheme,
         caz_start_q, treat_group, length, total_inj_adj_All,
         total_inj_adj_Vehicle, total_inj_adj_ActiveTravel) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  left_join(matched_covars, by = "OA") %>%
  mutate(
    outcome   = total_inj_adj_All,
    road_length_km = pmax(length / 1000, 1e-6),
    group     = if_else(treat_group == 1, "CAZ roads", "Matched controls"),
    post_flag = as.integer(quarter_year >= ref_start),
    uid       = paste0(panel_id, "_", scheme),
    uid_int   = as.integer(factor(uid)),
    qtr_int   = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    ),
    covid_period = as.integer(quarter_year >= as.yearqtr("2020 Q1") &
                              quarter_year <= as.yearqtr("2021 Q2"))
  ) %>%
  filter(!is.na(outcome))

rm(road_panel_matched)

# Keep full panel for Step 6 (include structural zeros)
model_panel_full <- model_panel

# Drop structural zeros (same as script 22)
units_with_inj <- model_panel %>%
  group_by(uid_int) %>%
  summarise(total = sum(outcome, na.rm = TRUE), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(uid_int)

n_before <- n_distinct(model_panel$uid_int)

model_panel <- model_panel %>%
  filter(uid_int %in% units_with_inj) %>%
  mutate(uid_int = as.integer(factor(uid_int)))

cat("\nDropped structural zeros:",
    n_before - n_distinct(model_panel$uid_int), "of", n_before, "units",
    sprintf("(%.1f%%)\n", 100 * (1 - n_distinct(model_panel$uid_int) / n_before)))

schemes_all <- sort(unique(model_panel$scheme))

covid_start_int <- as.integer(round((as.numeric(as.yearqtr("2020 Q1")) - min_qtr) * 4)) + 1L
covid_end_int   <- as.integer(round((as.numeric(as.yearqtr("2021 Q2")) - min_qtr) * 4)) + 1L

cat("Model panel:", nrow(model_panel), "rows |",
    n_distinct(model_panel$uid_int), "units |",
    ncol(model_panel), "columns\n")

# =============================================================================
# STEP 5 — UNADJUSTED VS ADJUSTED C&S
# =============================================================================

cat("\n== STEP 5 — UNADJUSTED VS ADJUSTED ==\n\n")

main_fit <- run_cs(model_panel, xformla = xformla_adj, label = "Overall (adjusted)")
agg_main <- main_fit$agg_obj

unadj_fit <- run_cs(model_panel, xformla = NULL, label = "Overall (unadjusted)")
agg_unadj <- unadj_fit$agg_obj

results_main  <- run_cs_by_scheme(model_panel, schemes_all, xformla = xformla_adj)
results_unadj <- run_cs_by_scheme(model_panel, schemes_all, xformla = NULL)

cat(sprintf("  Unadjusted: %.4f (SE %.4f)\n  Adjusted:   %.4f (SE %.4f)\n",
            agg_unadj$overall.att, agg_unadj$overall.se,
            agg_main$overall.att, agg_main$overall.se))

# =============================================================================
# STEP 5b — EXCLUDE COVID (2020 Q1 – 2021 Q2)
# =============================================================================

cat("\n== STEP 5b — EXCLUDE COVID (2020 Q1 – 2021 Q2) ==\n\n")

panel_nocovid <- model_panel %>%
  filter(!(qtr_int >= covid_start_int & qtr_int <= covid_end_int)) %>%
  mutate(uid_nocov = as.integer(factor(uid)))

cat("Excluded qtr_int", covid_start_int, "to", covid_end_int,
    "| Retained", n_distinct(panel_nocovid$qtr_int), "of",
    n_distinct(model_panel$qtr_int), "quarters\n")

fit_nocovid <- tryCatch({
  att_nc <- att_gt(
    yname = "outcome", tname = "qtr_int", idname = "uid_nocov", gname = "g",
    data = panel_nocovid, control_group = "notyettreated",
    xformla = xformla_adj,
    bstrap = TRUE, anticipation = 0, panel = FALSE
  )
  agg_nc <- aggte(att_nc, type = "simple", na.rm = TRUE)
  list(att_obj = att_nc, agg_obj = agg_nc,
       result = tibble(att = agg_nc$overall.att, se = agg_nc$overall.se,
                       pval = 2 * pnorm(-abs(agg_nc$overall.att / agg_nc$overall.se))))
}, error = function(e) { cat("  Overall FAILED:", e$message, "\n"); NULL })

agg_nocovid <- if (!is.null(fit_nocovid)) fit_nocovid$agg_obj else NULL

if (!is.null(agg_nocovid)) {
  cat("\n=== Overall ATT (excl COVID, adjusted) ===\n")
  summary(agg_nocovid)
}

results_nocovid <- run_cs_by_scheme(panel_nocovid, schemes_all,
                                    id = "uid_nocov", xformla = xformla_adj,
                                    panel = FALSE)

# =============================================================================
# STEP 6 — INCLUDE STRUCTURAL ZEROS (FULL PANEL)
# =============================================================================

cat("\n== STEP 6 — INCLUDE STRUCTURAL ZEROS ==\n\n")

panel_full_reindexed <- model_panel_full %>%
  mutate(uid_full = as.integer(factor(uid)))

cat("Full panel units:", n_distinct(panel_full_reindexed$uid_full),
    "vs primary:", n_distinct(model_panel$uid_int), "\n")

fit_full <- run_cs(panel_full_reindexed, id = "uid_full", xformla = xformla_adj,
                   label = "All roads incl. structural zeros")
agg_full <- fit_full$agg_obj
summary(agg_full)

results_full <- run_cs_by_scheme(panel_full_reindexed, schemes_all,
                                 id = "uid_full", xformla = xformla_adj)

rm(model_panel_full, panel_full_reindexed)

# =============================================================================
# STEP 7 — PPML TWFE (count model with road-length offset)
# =============================================================================

cat("\n== STEP 7 — PPML TWFE ==\n\n")

panel_ppml <- model_panel %>%
  mutate(treated  = as.integer(treat_group),
         post_int = post_flag,
         qtr_num  = as.numeric(quarter_year),
         log_length = log(road_length_km),
         ppml_id  = as.integer(factor(uid_int)))

ppml_overall <- tryCatch({
  m <- fepois(outcome ~ treated:post_int | ppml_id + qtr_num,
              data = panel_ppml, offset = ~log_length)
  b  <- coef(m)["treated:post_int"]
  se <- sqrt(vcov(m)["treated:post_int", "treated:post_int"])
  tibble(scheme = "Overall",
         pct_change = (exp(b)-1)*100, ci_lo_pct = (exp(b-1.96*se)-1)*100,
         ci_hi_pct = (exp(b+1.96*se)-1)*100, pval = 2*pnorm(-abs(b/se)))
}, error = function(e) { cat("  PPML overall FAILED:", e$message, "\n"); NULL })

ppml_by_scheme <- map_df(schemes_all, function(s) {
  cat("  PPML:", s, "\n")
  d <- panel_ppml %>% filter(scheme == s) %>%
    mutate(ppml_id = as.integer(factor(ppml_id)), qtr_num = as.numeric(quarter_year))
  tryCatch({
    m <- fepois(outcome ~ treated:post_int | ppml_id + qtr_num,
                data = d, offset = ~log_length)
    b  <- coef(m)["treated:post_int"]
    se <- sqrt(vcov(m)["treated:post_int", "treated:post_int"])
    tibble(scheme = s, pct_change = (exp(b)-1)*100,
           ci_lo_pct = (exp(b-1.96*se)-1)*100,
           ci_hi_pct = (exp(b+1.96*se)-1)*100,
           pval = 2*pnorm(-abs(b/se)))
  }, error = function(e) {
    cat("    FAILED\n")
    tibble(scheme = s, pct_change = NA_real_, ci_lo_pct = NA_real_,
           ci_hi_pct = NA_real_, pval = NA_real_)
  })
})

ppml_results <- bind_rows(ppml_overall, ppml_by_scheme) %>%
  mutate(sig = case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                         pval < 0.1 ~ "*", TRUE ~ ""),
         sig_label = paste0(sprintf("%+.1f%%", pct_change), sig))
print(ppml_results)

# =============================================================================
# STEP 8 — EXCLUDE 2024
# =============================================================================

cat("\n== STEP 8 — EXCLUDE 2024 ==\n\n")

panel_no24 <- model_panel %>%
  filter(as.numeric(quarter_year) < 2024) %>%
  mutate(uid_no24 = as.integer(factor(uid)))

fit_no24 <- run_cs(panel_no24, id = "uid_no24", xformla = xformla_adj,
                   label = "Overall excl 2024 (adjusted)")
agg_no24 <- fit_no24$agg_obj
summary(agg_no24)

results_no24 <- run_cs_by_scheme(panel_no24, schemes_all,
                                 id = "uid_no24", xformla = xformla_adj)

# =============================================================================
# STEP 9 — COVID-INTERACTION TWFE
# =============================================================================

cat("\n== STEP 9 — COVID-INTERACTION TWFE ==\n\n")

panel_reg <- model_panel %>%
  mutate(treated = as.integer(treat_group), post_int = post_flag,
         qtr_num = as.numeric(quarter_year))

m_twfe_covid <- feols(outcome ~ treated:post_int + treated:covid_period |
                        identifier + qtr_num,
                      data = panel_reg)

m_twfe_nocovid <- feols(outcome ~ treated:post_int | identifier + qtr_num,
                        data = panel_reg)

cat("--- TWFE comparison: with vs without COVID interaction ---\n")
etable(m_twfe_nocovid, m_twfe_covid,
       headers = c("Without COVID", "With COVID interaction"),
       dict = c("treated:post_int" = "Treated x Post",
                 "treated:covid_period" = "Treated x COVID"),
       digits = 5, digits.stats = 3)

# =============================================================================
# STEP 10 — MODE-STRATIFIED OUTCOMES
# =============================================================================

cat("\n== STEP 10 — MODE-STRATIFIED OUTCOMES ==\n\n")

mode_outcomes <- c(Vehicle = "total_inj_adj_Vehicle",
                   ActiveTravel = "total_inj_adj_ActiveTravel")

mode_results <- map_df(names(mode_outcomes), function(mode_name) {
  cat("  Mode:", mode_name, "\n")
  col <- mode_outcomes[mode_name]

  d <- model_panel %>%
    mutate(outcome_mode = .data[[col]]) %>%
    # Drop structural zeros for this outcome
    group_by(uid_int) %>%
    mutate(unit_total = sum(outcome_mode, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(unit_total > 0) %>%
    mutate(uid_mode = as.integer(factor(uid_int)))

  cat("    Units with injuries:", n_distinct(d$uid_mode), "\n")

  fit <- tryCatch(
    run_cs(d, yname = "outcome_mode", id = "uid_mode",
           xformla = xformla_adj, label = paste0("  ", mode_name)),
    error = function(e) { cat("    FAILED:", e$message, "\n"); NULL })

  if (is.null(fit)) return(tibble(mode = mode_name, att = NA_real_,
                                   se = NA_real_, pval = NA_real_))
  fit$result %>% mutate(mode = mode_name, .before = 1)
})

cat("\nMode-stratified overall ATTs:\n")
print(mode_results)

# =============================================================================
# STEP 11 — ROBUSTNESS SUMMARY TABLE
# =============================================================================

cat("\n== STEP 11 — ROBUSTNESS SUMMARY ==\n\n")

# --- Overall ATT comparison ---
nocov_att <- if (!is.null(agg_nocovid)) agg_nocovid$overall.att else NA_real_
nocov_se  <- if (!is.null(agg_nocovid)) agg_nocovid$overall.se else NA_real_
nocov_p   <- if (!is.null(agg_nocovid)) {
  2 * pnorm(-abs(nocov_att / nocov_se))
} else NA_real_

ppml_ov <- ppml_results %>% filter(scheme == "Overall")

overall_tbl <- tibble(
  Specification = c(
    "Primary (C&S adjusted)",
    "C&S unadjusted",
    "C&S excl COVID",
    "C&S incl structural zeros",
    "C&S excl 2024",
    "PPML TWFE",
    "Vehicle only (C&S adj)",
    "Active travel only (C&S adj)"
  ),
  ATT = c(
    sprintf("%.4f", agg_main$overall.att),
    sprintf("%.4f", agg_unadj$overall.att),
    sprintf("%.4f", nocov_att),
    sprintf("%.4f", agg_full$overall.att),
    sprintf("%.4f", agg_no24$overall.att),
    if (!is.null(ppml_ov) && nrow(ppml_ov) > 0)
      sprintf("%+.1f%%", ppml_ov$pct_change) else "FAILED",
    sprintf("%.4f", mode_results$att[mode_results$mode == "Vehicle"]),
    sprintf("%.4f", mode_results$att[mode_results$mode == "ActiveTravel"])
  ),
  SE = c(
    sprintf("%.4f", agg_main$overall.se),
    sprintf("%.4f", agg_unadj$overall.se),
    sprintf("%.4f", nocov_se),
    sprintf("%.4f", agg_full$overall.se),
    sprintf("%.4f", agg_no24$overall.se),
    "",
    sprintf("%.4f", mode_results$se[mode_results$mode == "Vehicle"]),
    sprintf("%.4f", mode_results$se[mode_results$mode == "ActiveTravel"])
  ),
  p = round(c(
    2 * pnorm(-abs(agg_main$overall.att / agg_main$overall.se)),
    2 * pnorm(-abs(agg_unadj$overall.att / agg_unadj$overall.se)),
    nocov_p,
    2 * pnorm(-abs(agg_full$overall.att / agg_full$overall.se)),
    2 * pnorm(-abs(agg_no24$overall.att / agg_no24$overall.se)),
    if (!is.null(ppml_ov) && nrow(ppml_ov) > 0) ppml_ov$pval else NA_real_,
    mode_results$pval[mode_results$mode == "Vehicle"],
    mode_results$pval[mode_results$mode == "ActiveTravel"]
  ), 3)
)

cat("\n--- Overall robustness comparison ---\n\n")
print(overall_tbl, width = Inf)

# --- Per-scheme comparison ---
rob_scheme <- tibble(Scheme = schemes_all) %>%
  left_join(results_main %>% transmute(Scheme = scheme, Primary = fmt_cell(att, pval)),
            by = "Scheme") %>%
  left_join(results_unadj %>% transmute(Scheme = scheme, Unadjusted = fmt_cell(att, pval)),
            by = "Scheme") %>%
  left_join(results_nocovid %>% transmute(Scheme = scheme, `Excl COVID` = fmt_cell(att, pval)),
            by = "Scheme") %>%
  left_join(results_full %>% transmute(Scheme = scheme, `Incl zeros` = fmt_cell(att, pval)),
            by = "Scheme") %>%
  left_join(results_no24 %>% transmute(Scheme = scheme, `Excl 2024` = fmt_cell(att, pval)),
            by = "Scheme") %>%
  left_join(ppml_results %>% filter(scheme != "Overall") %>%
              transmute(Scheme = scheme,
                        PPML = sprintf("%+.1f%%%s", pct_change,
                          case_when(pval < 0.01 ~ "***", pval < 0.05 ~ "**",
                                    pval < 0.1 ~ "*", TRUE ~ ""))),
            by = "Scheme")

cat("\n--- Per-scheme robustness comparison ---\n\n")
print(rob_scheme, width = Inf)

# --- Save ---
write_csv(overall_tbl, file.path(outdir, "robustness_overall.csv"))
write_csv(rob_scheme, file.path(outdir, "robustness_by_scheme.csv"))

cat("\n=== ALL SENSITIVITY ANALYSES COMPLETE ===\n")
cat("Output:", outdir, "\n")
cat("End of script.\n")
