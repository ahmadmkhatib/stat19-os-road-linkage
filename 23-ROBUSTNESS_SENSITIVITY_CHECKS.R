

# =============================================================================
# CAZ DiD ANALYSIS — ROBUSTNESS, SENSITIVITY AND CHECKS
# =============================================================================
#
# This script builds on the main PPML script.
#
# It contains:
#   1. Callaway & Sant'Anna robustness analysis
#   2. C&S event study excluding the COVID period
#   3. Additional checks based on saved main-analysis data
#
# =============================================================================


# =============================================================================
# 0. PACKAGES AND SETTINGS
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(patchwork)

outdir <- here("output", "pooled", "All_clean", "robustness_checks")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

run_cs <- function(data, xformla = NULL) {
  
  att <- att_gt(
    yname         = "outcome_raw",
    tname         = "qtr_int",
    idname        = "uid_int",
    gname         = "g",
    data          = data,
    control_group = "notyettreated",
    xformla       = xformla,
    bstrap        = TRUE,
    anticipation  = 0,
    panel         = TRUE
  )
  
  agg <- aggte(att, type = "simple", na.rm = TRUE)
  dyn <- aggte(att, type = "dynamic", na.rm = TRUE)
  
  list(att = att, agg = agg, dyn = dyn)
}

extract_cs_event_study <- function(dyn, trim = c(-12, 12)) {
  
  tibble(
    event_time = dyn$egt,
    estimate = dyn$att.egt,
    se = dyn$se.egt
  ) %>%
    mutate(
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se
    ) %>%
    filter(event_time >= trim[1], event_time <= trim[2])
}

plot_event_study <- function(df, title, subtitle, ylab, colour = "#2E6FAB") {
  
  ggplot(df, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = colour) +
    geom_line(linewidth = 0.8, colour = colour) +
    geom_point(size = 1.8, colour = colour) +
    scale_x_continuous(breaks = seq(-12, 12, by = 2)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Quarters relative to CAZ implementation",
      y = ylab
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}


# =============================================================================
# 2. LOAD MAIN-ANALYSIS DATA
# =============================================================================
# These files are created by the main PPML script.

model_panel <- arrow::read_parquet(
  here("data", "processed", "final_model_panel_main.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year)
  )

main_results <- readRDS(
  here("data", "processed", "caz_ppml_main_results.rds")
)


# =============================================================================
# 3. C&S ROBUSTNESS ANALYSIS — FULL SAMPLE
# =============================================================================
# C&S is used as a robustness estimator on the additive injury-count scale.
# Zero-only road-link series are removed because they do not help identify
# this estimator.

cs_data <- model_panel %>%
  group_by(uid_int) %>%
  filter(sum(outcome_raw, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  mutate(uid_int = as.integer(factor(uid_int)))

cs_fit <- run_cs(
  cs_data,
  xformla = ~ log1p_pop_density + IMD + log1p_road_density_m_km2
)

summary(cs_fit$agg)
summary(cs_fit$att)

cs_es <- extract_cs_event_study(cs_fit$dyn)

p_cs <- plot_event_study(
  cs_es,
  title = "Robustness: C&S event study",
  subtitle = "Additive ATT scale; full sample",
  ylab = "ATT, injury count",
  colour = "#2E6FAB"
)

print(p_cs)

ggsave(
  file.path(outdir, "cs_event_study_full_sample.png"),
  p_cs,
  width = 10,
  height = 7,
  dpi = 300
)


# =============================================================================
# 4. C&S CHECK EXCLUDING COVID PERIOD
# =============================================================================
# Removes 2020 Q1–2021 Q4 and reruns C&S to assess whether COVID drives
# pre-treatment diagnostics.

cs_data_nocovid <- cs_data %>%
  filter(
    quarter_year < as.yearqtr("2020 Q1") |
      quarter_year > as.yearqtr("2021 Q4")
  )

cat("Rows before removing COVID:", nrow(cs_data), "\n")
cat("Rows after removing COVID:", nrow(cs_data_nocovid), "\n")

cs_fit_nocovid <- run_cs(
  cs_data_nocovid,
  xformla = ~ log1p_pop_density + IMD + log1p_road_density_m_km2
)

summary(cs_fit_nocovid$agg)
summary(cs_fit_nocovid$att)

cs_es_nocovid <- extract_cs_event_study(cs_fit_nocovid$dyn)

p_cs_nocovid <- plot_event_study(
  cs_es_nocovid,
  title = "Robustness: C&S event study excluding COVID period",
  subtitle = "Excludes 2020 Q1–2021 Q4",
  ylab = "ATT, injury count",
  colour = "#27AE60"
)

print(p_cs_nocovid)

ggsave(
  file.path(outdir, "cs_event_study_excluding_covid.png"),
  p_cs_nocovid,
  width = 10,
  height = 7,
  dpi = 300
)

cs_immediate_pre_nocovid <- cs_es_nocovid %>%
  filter(event_time >= -8, event_time <= -1) %>%
  select(event_time, estimate, ci_lo, ci_hi)

print(cs_immediate_pre_nocovid)


# =============================================================================
# 5. SAVE ROBUSTNESS OUTPUTS
# =============================================================================

cs_summary <- tibble(
  model = "C&S full sample",
  estimate = cs_fit$agg$overall.att,
  se = cs_fit$agg$overall.se,
  ci_lo = estimate - 1.96 * se,
  ci_hi = estimate + 1.96 * se
)

cs_nocovid_summary <- tibble(
  model = "C&S excluding COVID",
  estimate = cs_fit_nocovid$agg$overall.att,
  se = cs_fit_nocovid$agg$overall.se,
  ci_lo = estimate - 1.96 * se,
  ci_hi = estimate + 1.96 * se
)

cs_comparison <- bind_rows(cs_summary, cs_nocovid_summary)

print(cs_comparison)

write_csv(cs_comparison, file.path(outdir, "cs_summary_comparison.csv"))
write_csv(cs_es, file.path(outdir, "cs_event_study_full_sample.csv"))
write_csv(cs_es_nocovid, file.path(outdir, "cs_event_study_excluding_covid.csv"))
write_csv(cs_immediate_pre_nocovid, file.path(outdir, "cs_immediate_pre_period_excluding_covid.csv"))

saveRDS(
  list(
    cs_summary = cs_summary,
    cs_nocovid_summary = cs_nocovid_summary,
    cs_comparison = cs_comparison,
    cs_event_study = cs_es,
    cs_event_study_nocovid = cs_es_nocovid,
    cs_immediate_pre_nocovid = cs_immediate_pre_nocovid,
    ppml_parallel_trends_tests = list(
      pt_8_2 = main_results$pt_8_2,
      pt_6_2 = main_results$pt_6_2,
      pt_4_2 = main_results$pt_4_2
    )
  ),
  here("data", "processed", "caz_robustness_sensitivity_checks.rds")
)

cat("\nROBUSTNESS, SENSITIVITY AND CHECKS COMPLETE.\n")
cat("Outputs saved to:", outdir, "\n")











####



### H0: β−4 = β−3 = β−2 = 0; p = 0.19 
#The treated roads and matched control roads do not show a detectable difference in injury trends before CAZ implementation.
#. The rejection in the broader window is entirely attributable to two noisy quarters 5–6 periods before treatment — 
#which is also where Bath and Bradford have very few post-COVID comparison observations, so sampling noise is amplified.
# Although a joint test of all lead coefficients from six quarters before treatment was statistically significant (p = 0.003), 
#lead coefficients in the year immediately before implementation (−4 to −2 quarters) were jointly insignificant (p = 0.19) 
#and exhibited no systematic pattern. 
#Combined with visual inspection of the event-study plots, 
#the matched control roads  providing evidence consistent with the parallel trends assumption for causal identification..


####    #### what is happening with Bath 

stacked %>%
  filter(stack_scheme == "Bath") %>%
  count(event_time)
stacked %>%
  filter(stack_scheme == "Bath") %>%
  group_by(event_time) %>%
  summarise(
    n_obs = n(),
    injuries = sum(outcome_raw),
    mean_inj = mean(outcome_raw),
    positive_obs = sum(outcome_raw > 0),
    .groups = "drop"
  )

bath_es <- feglm(
  outcome_raw ~ i(event_time_f, treat_group, ref = "-1") +
    covid_lockdown_treated + covid_recovery_treated |
    uid_stack + qtr_int,
  data = stacked %>% filter(stack_scheme == "Bath"),
  family = "poisson",
  cluster = ~OA
)

etable(bath_es)
coef(bath_es)


stacked %>%
  filter(stack_scheme == "Bath")

stacked %>%
  filter(
    stack_scheme == "Bath",
    event_time == 8
  ) %>%
  group_by(treat_group) %>%
  summarise(
    n = n(),
    injuries = sum(outcome_raw),
    mean_inj = mean(outcome_raw),
    positive_obs = sum(outcome_raw > 0),
    .groups = "drop"
  )

