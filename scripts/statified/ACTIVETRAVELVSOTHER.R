# =============================================================================
# CAZ DiD — MODEL 1 BY ROAD-USER GROUP
# Active-travel injuries vs vehicle injuries
# COVID-adjusted models only
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(lubridate)
library(fixest)

# =============================================================================
# USER SETTINGS
# =============================================================================

outdir <- here("output", "pooled", "All_clean", "active_vs_vehicle_model1")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

road_panel <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q  = as.yearqtr(caz_start_q)
  )

road_caz_props <- readRDS(
  here("data", "processed", "roads_caz_props.rds")
)

# =============================================================================
# CREATE OUTCOME VARIABLES
# =============================================================================
# Active travel = pedestrian + cyclist injuries
# Vehicle = existing vehicle injury count

road_panel <- road_panel %>%
  mutate(
    total_inj_adj_ActiveTravel_new = case_when(
      is.na(total_inj_adj_Pedestrian) & is.na(total_inj_adj_Cyclist) ~ NA_real_,
      TRUE ~ coalesce(total_inj_adj_Pedestrian, 0) +
        coalesce(total_inj_adj_Cyclist, 0)
    ),
    
    total_inj_adj_Vehicle_new = total_inj_adj_Vehicle
  )

outcome_groups <- c(
  active_travel = "total_inj_adj_ActiveTravel_new",
  vehicle       = "total_inj_adj_Vehicle_new"
)

# Check outcome columns
missing_outcomes <- setdiff(unname(outcome_groups), names(road_panel))

if (length(missing_outcomes) > 0) {
  cat("\nMissing outcome columns:\n")
  print(missing_outcomes)
  stop("Outcome variables not found.")
}

# =============================================================================
# MAJORITY-QUARTER RULE FOR CAZ START DATE
# =============================================================================

scheme_start <- road_caz_props %>%
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

road_panel <- road_panel %>%
  left_join(scheme_start, by = "scheme") %>%
  mutate(
    caz_start_q = coalesce(caz_start_q_adj, caz_start_q)
  ) %>%
  select(-caz_start_q_adj)

scheme_timing <- road_panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q)

cat("\nScheme timing after majority-quarter adjustment:\n")
print(scheme_timing)

rm(road_caz_props, scheme_start)

# =============================================================================
# BUILD MODEL PANEL
# =============================================================================

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

model_panel <- road_panel %>%
  select(
    panel_id,
    identifier,
    OA,
    scheme,
    quarter_year,
    caz_start_q,
    treat_group,
    all_of(unname(outcome_groups))
  ) %>%
  left_join(
    scheme_timing %>% rename(ref_start = caz_start_q),
    by = "scheme"
  ) %>%
  mutate(
    uid     = paste0(panel_id, "_", scheme),
    uid_int = as.integer(factor(uid)),
    
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    
    post_flag = as.integer(treat_group == 1 & quarter_year >= ref_start),
    
    covid_period = factor(
      case_when(
        quarter_year >= as.yearqtr("2020 Q2") &
          quarter_year <= as.yearqtr("2021 Q1") ~ "lockdown",
        
        quarter_year >= as.yearqtr("2021 Q2") &
          quarter_year <= as.yearqtr("2021 Q3") ~ "recovery",
        
        TRUE ~ "pre_covid"
      ),
      levels = c("pre_covid", "lockdown", "recovery")
    ),
    
    group = if_else(treat_group == 1, "CAZ roads", "Matched controls")
  )

schemes_all <- sort(unique(model_panel$scheme))

# =============================================================================
# BASIC SAMPLE CHECKS
# =============================================================================

cat("\nNumber of OAs:", n_distinct(model_panel$OA), "\n")
cat("Number of road links:", n_distinct(model_panel$identifier), "\n")
cat("Number of road-link-by-scheme units:", n_distinct(model_panel$uid_int), "\n")

sample_summary <- model_panel %>%
  group_by(group) %>%
  summarise(
    units = n_distinct(uid_int),
    observations = n(),
    across(
      all_of(unname(outcome_groups)),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        pct_zero = ~ 100 * mean(.x == 0, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

cat("\nSample summary by group:\n")
print(sample_summary)

write_csv(
  sample_summary,
  file.path(outdir, "sample_summary_active_vs_vehicle.csv")
)

# =============================================================================
# BUILD STACKED DATA
# =============================================================================

stacked <- map_dfr(schemes_all, function(sc) {
  
  sc_start <- scheme_timing %>%
    filter(scheme == sc) %>%
    pull(caz_start_q)
  
  if (length(sc_start) == 0 || is.na(sc_start)) return(NULL)
  
  sc_start_int <- as.integer(round((as.numeric(sc_start) - min_qtr) * 4)) + 1L
  
  model_panel %>%
    filter(scheme == sc) %>%
    mutate(
      stack_scheme = sc,
      event_time   = qtr_int - sc_start_int,
      uid_stack    = paste0(uid_int, "_", sc),
      post_stack   = as.integer(treat_group == 1 & event_time >= 0)
    )
}) %>%
  mutate(
    stack_scheme = factor(stack_scheme)
  )

cat("\nStacked rows:", nrow(stacked), "\n")
cat("Stacked units:", n_distinct(stacked$uid_stack), "\n")

# =============================================================================
# HELPER FUNCTION — RUN COVID-ADJUSTED MODEL 1 FOR ONE OUTCOME
# =============================================================================

run_model1_outcome <- function(outcome_name, outcome_col, data) {
  
  d <- data %>%
    filter(!is.na(.data[[outcome_col]]))
  
  cat("\n============================================================\n")
  cat("Outcome:", outcome_name, "\n")
  cat("Column:", outcome_col, "\n")
  cat("Rows:", nrow(d), "\n")
  cat("Units:", n_distinct(d$uid_stack), "\n")
  cat("Mean outcome:", mean(d[[outcome_col]], na.rm = TRUE), "\n")
  cat("% zero:", 100 * mean(d[[outcome_col]] == 0, na.rm = TRUE), "\n")
  cat("============================================================\n")
  
  fml <- as.formula(
    paste0(
      outcome_col,
      " ~ post_stack + i(covid_period, treat_group, ref = 'pre_covid') | ",
      "uid_stack + stack_scheme[qtr_int]"
    )
  )
  
  fit <- feglm(
    fml = fml,
    data = d,
    family = "poisson",
    cluster = ~OA,
    lean = TRUE
  )
  
  cat("\nCOVID-adjusted Model 1 for:", outcome_name, "\n")
  print(
    etable(
      fit,
      dict = c(
        "post_stack" = "CAZ post-treatment"
      )
    )
  )
  
  result <- tibble(
    outcome = outcome_name,
    outcome_col = outcome_col,
    estimate_log_irr = coef(fit)["post_stack"],
    se = se(fit)["post_stack"]
  ) %>%
    mutate(
      ci_lo = estimate_log_irr - 1.96 * se,
      ci_hi = estimate_log_irr + 1.96 * se,
      irr = exp(estimate_log_irr),
      irr_lo = exp(ci_lo),
      irr_hi = exp(ci_hi),
      pct_change = 100 * (irr - 1),
      pct_lo = 100 * (irr_lo - 1),
      pct_hi = 100 * (irr_hi - 1)
    )
  
  list(
    fit = fit,
    result = result
  )
}

# =============================================================================
# RUN ACTIVE-TRAVEL AND VEHICLE MODELS
# =============================================================================

model1_fits <- imap(
  outcome_groups,
  ~ run_model1_outcome(
    outcome_name = .y,
    outcome_col  = .x,
    data = stacked
  )
)

model1_results <- map_dfr(model1_fits, "result") %>%
  mutate(
    outcome_label = case_when(
      outcome == "active_travel" ~ "Active travel",
      outcome == "vehicle" ~ "Vehicle",
      TRUE ~ outcome
    )
  )

cat("\nCombined active-travel vs vehicle results:\n")
print(model1_results)

# =============================================================================
# SIDE-BY-SIDE TABLE
# =============================================================================

model1_table <- model1_results %>%
  transmute(
    outcome = outcome_label,
    log_irr = estimate_log_irr,
    se = se,
    irr = irr,
    irr_lo = irr_lo,
    irr_hi = irr_hi,
    pct_change = pct_change,
    pct_lo = pct_lo,
    pct_hi = pct_hi
  )

cat("\nSide-by-side result table:\n")
print(model1_table)

# =============================================================================
# PLOT SIDE BY SIDE
# =============================================================================

p_active_vs_vehicle <- ggplot(
  model1_results,
  aes(
    x = pct_change,
    y = fct_reorder(outcome_label, pct_change)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = pct_lo, xmax = pct_hi),
    width = 0.18
  ) +
  geom_point(size = 3) +
  labs(
    title = "Model 1: CAZ effects by road-user group",
    subtitle = "COVID-adjusted stacked PPML estimates; active travel vs vehicle injuries",
    x = "% change in injuries",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(p_active_vs_vehicle)

ggsave(
  file.path(outdir, "model1_active_vs_vehicle_covid_adjusted.png"),
  p_active_vs_vehicle,
  width = 8,
  height = 4.5,
  dpi = 300
)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

write_csv(
  model1_results,
  file.path(outdir, "model1_active_vs_vehicle_covid_adjusted_results.csv")
)

write_csv(
  model1_table,
  file.path(outdir, "model1_active_vs_vehicle_covid_adjusted_table.csv")
)

saveRDS(
  model1_fits,
  file.path(outdir, "model1_active_vs_vehicle_covid_adjusted_fits.rds")
)

arrow::write_parquet(
  stacked,
  here("data", "processed", "final_stacked_data_model1_active_vs_vehicle.parquet")
)

cat("\nDone. Outputs saved to:\n")
cat(outdir, "\n")