# =============================================================================
# DiD — POOLED MATCHING — PRIMARY ANALYSIS (TOTAL INJURIES)
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(fixest)
library(lubridate)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

run_cs <- function(data, id = "uid_int", xformla = NULL, label = "") {
  if (nchar(label) > 0) cat("  ", label, "\n")
  
  att <- att_gt(
    yname = "outcome",
    tname = "qtr_int",
    idname = id,
    gname = "g",
    data = data,
    control_group = "notyettreated",
    xformla = xformla,
    bstrap = TRUE,
    anticipation = 0,
    panel = TRUE
  )
  
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  
  list(
    att_obj = att,
    agg_obj = agg,
    result = tibble(
      att   = agg$overall.att,
      se    = agg$overall.se,
      ci_lo = agg$overall.att - 1.96 * agg$overall.se,
      ci_hi = agg$overall.att + 1.96 * agg$overall.se,
      pval  = 2 * pnorm(-abs(agg$overall.att / agg$overall.se))
    )
  )
}

run_cs_by_scheme <- function(data, schemes, id = "uid_int", xformla = NULL) {
  map_df(schemes, function(s) {
    cat("  Scheme:", s, "\n")
    
    d <- data %>%
      filter(scheme == s) %>%
      mutate(
        !!id := as.integer(factor(.data[[id]])),
        g = as.numeric(g)
      )
    
    fit <- tryCatch(
      run_cs(d, id = id, xformla = xformla),
      error = function(e) {
        cat("    FAILED:", e$message, "\n")
        NULL
      }
    )
    
    if (is.null(fit)) {
      return(tibble(
        scheme = s,
        att = NA_real_,
        se = NA_real_,
        ci_lo = NA_real_,
        ci_hi = NA_real_,
        pval = NA_real_
      ))
    }
    
    bind_cols(tibble(scheme = s), fit$result)
  }) %>%
    mutate(
      sig = case_when(
        pval < 0.01 ~ "***",
        pval < 0.05 ~ "**",
        pval < 0.1  ~ "*",
        TRUE ~ ""
      ),
      sig_label = paste0(sprintf("%.4f", att), sig)
    )
}

es_from_att <- function(att_obj, trim = c(-12, 12)) {
  agg_dyn <- aggte(att_obj, type = "dynamic", na.rm = TRUE)
  
  tibble(
    event_time = agg_dyn$egt,
    att = agg_dyn$att.egt,
    se = agg_dyn$se.egt
  ) %>%
    mutate(
      ci_lo = att - 1.96 * se,
      ci_hi = att + 1.96 * se
    ) %>%
    filter(event_time >= trim[1], event_time <= trim[2])
}

plot_es <- function(es_df, title, subtitle = "") {
  ggplot(es_df, aes(x = event_time, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Quarters relative to CAZ implementation",
      y = "ATT injury count"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

plot_coefs <- function(res_df, title, subtitle = "") {
  ggplot(res_df, aes(x = att, y = fct_reorder(scheme, att))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(
      aes(xmin = ci_lo, xmax = ci_hi),
      width = 0.25,
      linewidth = 0.7,
      orientation = "y"
    ) +
    geom_point(size = 3, colour = "#2E6FAB") +
    geom_text(aes(label = sig_label), hjust = -0.15, size = 3) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "ATT injury count",
      y = NULL,
      caption = "Error bars: 95% CI  |  * p<0.1  ** p<0.05  *** p<0.01"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
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

cat(
  "Rows:", nrow(road_panel_matched),
  "| Roads:", n_distinct(road_panel_matched$identifier),
  "| Quarters:", n_distinct(road_panel_matched$quarter_year), "\n"
)

# =============================================================================
# MAJORITY-RULE TREATMENT QUARTER ADJUSTMENT
# =============================================================================

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

cat("\nScheme timing majority-rule adjusted:\n")
print(scheme_timing %>% arrange(caz_start_q))

# =============================================================================
# OA-LEVEL COVARIATES
# =============================================================================

matched_full_cov <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))

matched_covars <- matched_full_cov %>%
  mutate(
    log1p_pop_density = log1p(pmax(pop_density, 0)),
    log1p_road_density_m_km2 = log1p(pmax(road_density_m_km2, 0))
  ) %>%
  select(OA, log1p_pop_density, IMD, log1p_road_density_m_km2) %>%
  distinct(OA, .keep_all = TRUE)

rm(matched_full_cov)

outcome_label <- "Total injuries"

# =============================================================================
# BUILD MODEL PANEL
# =============================================================================

outdir <- here("output", "pooled", "All")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

min_qtr <- min(as.numeric(road_panel_matched$quarter_year), na.rm = TRUE)

model_panel <- road_panel_matched %>%
  select(
    panel_id, identifier, OA, quarter_year, scheme,
    caz_start_q, treat_group, total_inj_adj_All
  ) %>%
  left_join(scheme_timing %>% rename(ref_start = caz_start_q), by = "scheme") %>%
  left_join(matched_covars, by = "OA") %>%
  mutate(
    outcome_raw = total_inj_adj_All,
    outcome = total_inj_adj_All,
    
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls"),
    post_flag = as.integer(quarter_year >= ref_start),
    uid = paste0(panel_id, "_", scheme),
    
    uid_int = as.integer(factor(uid)),
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    
    g = case_when(
      treat_group == 1 & !is.na(caz_start_q) ~
        as.numeric(round((as.numeric(caz_start_q) - min_qtr) * 4)) + 1,
      TRUE ~ 0
    ),
    
    season_qtr = factor(as.integer(cycle(quarter_year))),
    
    covid_period = as.integer(
      quarter_year >= as.yearqtr("2020 Q2") &
        quarter_year <= as.yearqtr("2021 Q4")
    )
  ) %>%
  filter(!is.na(outcome)) %>%
  select(
    uid, uid_int, identifier, OA, scheme, quarter_year, qtr_int,
    treat_group, g, post_flag, group,
    outcome_raw, outcome,
    season_qtr, covid_period,
    log1p_pop_density, IMD, log1p_road_density_m_km2
  )

rm(road_panel_matched)

# =============================================================================
# DROP STRUCTURAL ZERO ROAD LINKS
# =============================================================================

units_with_inj <- model_panel %>%
  group_by(uid_int) %>%
  summarise(total = sum(outcome_raw, na.rm = TRUE), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(uid_int)

n_before <- n_distinct(model_panel$uid_int)

model_panel <- model_panel %>%
  filter(uid_int %in% units_with_inj) %>%
  mutate(uid_int = as.integer(factor(uid_int)))

cat(
  "\nDropped structural zeros:",
  n_before - n_distinct(model_panel$uid_int), "of", n_before, "units",
  sprintf("(%.1f%%)\n", 100 * (1 - n_distinct(model_panel$uid_int) / n_before))
)

schemes_all <- sort(unique(model_panel$scheme))

cat(
  "Model panel before seasonal/COVID residualisation:",
  nrow(model_panel), "rows |",
  n_distinct(model_panel$uid_int), "units |",
  ncol(model_panel), "columns\n"
)

# =============================================================================
# SEASONALITY + COVID ADJUSTMENT OUTSIDE att_gt()
# =============================================================================

season_covid_model <- feols(
  outcome_raw ~ season_qtr + covid_period,
  data = model_panel
)

model_panel <- model_panel %>%
  mutate(
    outcome = resid(season_covid_model) + mean(outcome_raw, na.rm = TRUE)
  )

# C&S covariates only — no seasonality/COVID/calendar factors inside xformla
xformla_adj <- ~
  log1p_pop_density +
  IMD +
  log1p_road_density_m_km2

cat(
  "Model panel after seasonal/COVID residualisation:",
  nrow(model_panel), "rows |",
  n_distinct(model_panel$uid_int), "units |",
  ncol(model_panel), "columns\n"
)

# =============================================================================
# STEP 0 — SUMMARY STATISTICS
# =============================================================================

cat("\n== STEP 0 — SUMMARY STATISTICS ==\n\n")

model_panel %>%
  group_by(group) %>%
  summarise(
    n_obs = format(n(), big.mark = ","),
    mean_raw = round(mean(outcome_raw, na.rm = TRUE), 4),
    mean_adjusted = round(mean(outcome, na.rm = TRUE), 4),
    pct_zero_raw = round(100 * mean(outcome_raw == 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# STEP 1 — DESCRIPTIVE DiD
# =============================================================================

cat("\n== STEP 1 — DESCRIPTIVE DiD USING RAW OUTCOME ==\n\n")

ss <- model_panel %>%
  group_by(group, post_flag) %>%
  summarise(
    mean_inj = round(mean(outcome_raw, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  mutate(period = if_else(post_flag == 1, "Post", "Pre"))

print(ss)

dt <- ss %>%
  pivot_wider(id_cols = group, names_from = period, values_from = mean_inj) %>%
  mutate(change = Post - Pre)

did_est <- dt$change[dt$group == "CAZ roads"] -
  dt$change[dt$group == "Matched controls"]

cat("\nDescriptive DiD raw outcome:", round(did_est, 5), "\n")

# =============================================================================
# STEP 2 — TREND PLOT
# =============================================================================

pd <- model_panel %>%
  group_by(group, quarter_year) %>%
  summarise(mean_inj = mean(outcome_raw, na.rm = TRUE), .groups = "drop")

scheme_dates <- scheme_timing %>%
  mutate(x = as.Date(caz_start_q)) %>%
  arrange(x)

ggplot(pd, aes(x = as.Date(quarter_year), y = mean_inj, colour = group)) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    data = scheme_dates,
    aes(xintercept = x),
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  geom_text(
    data = scheme_dates,
    aes(x = x, y = max(pd$mean_inj) * 1.02, label = scheme),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 1,
    vjust = -0.3,
    size = 2.5,
    colour = "grey40"
  ) +
  scale_colour_manual(values = c(
    "CAZ roads" = "#E74C3C",
    "Matched controls" = "#2C3E50"
  )) +
  labs(
    title = "Trend: Total injuries staggered treatment timing",
    x = NULL,
    y = "Mean injury count",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "trend_plot_raw.png"), width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 3 — TWFE WITH ROAD AND QUARTER FIXED EFFECTS
# =============================================================================

cat("\n== STEP 3 — TWFE ==\n\n")

panel_reg <- model_panel %>%
  mutate(
    treated = as.integer(treat_group),
    post_int = post_flag,
    qtr_num = as.numeric(quarter_year)
  )

m_twfe <- feols(
  outcome_raw ~ treated:post_int | identifier + qtr_num,
  data = panel_reg
)

etable(
  m_twfe,
  headers = outcome_label,
  dict = c("treated:post_int" = "Treated x Post"),
  digits = 5,
  digits.stats = 3
)

# =============================================================================
# STEP 4 — MAIN: COVARIATE-ADJUSTED C&S ON SEASON/COVID-ADJUSTED OUTCOME
# =============================================================================

cat("\n== STEP 4 — MAIN: ADJUSTED C&S ==\n\n")

main_fit <- run_cs(
  model_panel,
  xformla = xformla_adj,
  label = "Overall adjusted after seasonality and COVID residualisation"
)

agg_main <- main_fit$agg_obj
summary(agg_main)

es_main <- es_from_att(main_fit$att_obj)

plot_es(
  es_main,
  paste0("Event study: ", outcome_label, " season/COVID adjusted")
)

ggsave(
  file.path(outdir, "event_study_adjusted.png"),
  width = 10,
  height = 7,
  dpi = 300
)

results_main <- run_cs_by_scheme(
  model_panel,
  schemes_all,
  xformla = xformla_adj
)

print(results_main)

plot_coefs(
  results_main,
  paste0("ATT by scheme: ", outcome_label, " season/COVID adjusted")
)

ggsave(
  file.path(outdir, "coef_plot_adjusted.png"),
  width = 10,
  height = 7,
  dpi = 300
)

# =============================================================================
# PER-SCHEME EVENT STUDIES
# =============================================================================

es_by_scheme <- map_df(schemes_all, function(s) {
  d <- model_panel %>%
    filter(scheme == s) %>%
    mutate(
      uid_int = as.integer(factor(uid)),
      g = as.numeric(g)
    )
  
  att <- tryCatch(
    att_gt(
      yname = "outcome",
      tname = "qtr_int",
      idname = "uid_int",
      gname = "g",
      data = d,
      control_group = "notyettreated",
      xformla = xformla_adj,
      bstrap = TRUE,
      anticipation = 0,
      panel = TRUE
    ),
    error = function(e) {
      cat("  Event study failed for", s, ":", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(att)) return(NULL)
  
  tryCatch(
    es_from_att(att) %>% mutate(scheme = s),
    error = function(e) NULL
  )
})

if (nrow(es_by_scheme) > 0) {
  ggplot(es_by_scheme, aes(x = event_time, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
    facet_wrap(~scheme, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Event study by scheme: ", outcome_label),
      x = "Quarters relative to CAZ",
      y = "ATT injury count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(outdir, "event_study_by_scheme.png"),
    width = 12,
    height = 10,
    dpi = 300
  )
}

cat("\n=== PRIMARY ANALYSIS COMPLETE ===\n")
cat("Output:", outdir, "\n")
cat("End of script.\n")





# =============================================================================
# SAVE FINAL MODEL DATA
# =============================================================================

final_model_data <- model_panel %>%
  select(
    uid, uid_int, identifier, OA, scheme,
    quarter_year, qtr_int,
    treat_group, g, post_flag, group,
    outcome_raw, outcome,
    season_qtr, covid_period,
    log1p_pop_density, IMD, log1p_road_density_m_km2
  )

arrow::write_parquet(
  final_model_data,
  here("data", "processed", "final_model_data_pooled_All.parquet")
)

saveRDS(
  final_model_data,
  here("data", "processed", "final_model_data_pooled_All.rds")
)

cat(
  "\nFinal model data saved:",
  nrow(final_model_data), "rows |",
  n_distinct(final_model_data$uid_int), "road-scheme units\n"
)





