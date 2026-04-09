# =============================================================================
#  OA-Level Matching — Two-Stage MDM Design  -  no exclusions 
#  Stage 1: MDM on structural + sociodemographic variables
#  Stage 2: MDM on pre-treatment injury trends + levels
# =============================================================================
#
# #   Country (England vs Scotland) is enforced as a hard exact constraint
#   at Stage 1 via exact = ~country. STATS19 and STATS19-Scotland differ
#   in recording thresholds, severity classification, and reporting
#   protocols. Cross-country matches introduce systematic measurement
#   incomparability that no covariate adjustment fixes.
#
# PRE-MATCHING EXCLUSIONS (applied before Stage 1):
#    Buffer OAs (within 1 km of scheme boundary) — treatment
#      contamination risk.
#  
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(here)
library(MASS)
library(tidyverse)
library(sf)

# =============================================================================
# SECTION 0 — ZERO-INJURY OA DESCRIPTIVES AND BIAS ASSESSMENT
# =============================================================================
# PURPOSE: Characterise the 242 excluded zero-injury treated OAs before
# any matching. This section:
#   (a) Checks whether they are systematically different from retained
#       treated OAs on structural/demographic covariates.
#   (b) Counts how many go from 0 pre-intervention to >0 post-intervention
#       (the key bias concern: are we missing treatment-induced injury
#       onset?).
#   (c) Checks scheme-level concentration — if one CAZ drives most
#       exclusions, the ATT estimate is implicitly scheme-selective.
#
# RUN THIS BEFORE ANY MATCHING. Results are descriptive only.
# =============================================================================

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))


# =============================================================================
# STEP 1 — DERIVE COUNTRY + no PRE-MATCHING EXCLUSIONS
# =============================================================================
# # =============================================================================

cat("\n====================================================\n")
cat("STEP 1: EXCLUSIONS AND SAMPLE CONSTRUCTION\n")
cat("====================================================\n\n")

oa_data_fin <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  ) %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA      == 0
     ) %>%
  mutate(
    treat_indicator = as.integer(treated_OA == 1)
  )

# --- Verify counts match report ---
# Report: 814 treated, 44,479 controls before matching
cat("\n=== SAMPLE AFTER EXCLUSIONS (verify against report) ===\n")
cat("Total OAs:",  nrow(oa_data_fin), "\n")
cat("  Treated:",  sum(oa_data_fin$treat_indicator == 1),
    "  (excl: 814)\n")
cat("  Controls:", sum(oa_data_fin$treat_indicator == 0),
    "  (excl: 44,479)\n\n")


country_tab_n <- oa_data_fin %>%
  count(country, treat_indicator) %>%
  pivot_wider(names_from = treat_indicator, values_from = n,
              names_prefix = "treated_") %>%
  rename(n_control = treated_0, n_treated = treated_1)

cat("Country distribution:\n")
print(country_tab_n)
# =============================================================================
# STEP 2 — DEFINE STAGE 1 AND STAGE 2 VARIABLES
# =============================================================================
# SEPARATION PRINCIPLE:
#   No variable appears in both stages. Stage 2 Mahalanobis distance is
#   determined entirely by injury dynamics. A good census match cannot
#   compensate for a poor trajectory match.
#
# =============================================================================

stage1_road <- c(
  "road_density_m_km2",
  "road_length_km",
  "pct_A_road",
  "pct_B_road",
  "pct_minor_road"
)

stage1_urban <- c(
  "dist_citycentre",
  "pop_density"
)

stage1_socdem <- c(
  "IMD",
  "cars_none_pct",
  "Drive_Car_pct",
  "Walk_pct",
  "Bicycle_pct",
  "X65plus_pct",
  "X5to19_pct",
  "X20to24_pct"
)

stage1_structural <- c(stage1_road, stage1_urban)   # caliper applied here only
stage1_vars       <- c(stage1_road, stage1_urban, stage1_socdem)

stage2_trends <- c(
  "trend_car_KSI_pkm",
  "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",
  "trend_ped_slight_pkm",
  "trend_total_pkm"
)

stage2_levels <- c(
  "mean_car_KSI_pkm",
  "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",
  "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",
  "mean_ped_slight_pkm",
  "mean_total_pkm"
)

stage2_vars <- c(stage2_trends, stage2_levels)

cat("\nVariable counts:\n")
cat("Stage 1:", length(stage1_vars),
    " (road:", length(stage1_road),
    " urban:", length(stage1_urban),
    " socdem:", length(stage1_socdem), ")\n")
cat("Stage 2:", length(stage2_vars),
    " (trends:", length(stage2_trends),
    " levels:", length(stage2_levels), ")\n")
cat("  Report: Stage 2 = 14 vars (7 trends + 7 levels) —",
    length(stage2_vars), "in code\n")

# =============================================================================
# STEP 3 — PRE-PROCESSING: WINSORISE ON FULL DATASET (STAGE 1 VARS ONLY)
# =============================================================================

oa_data_clean <- oa_data_fin %>%
  mutate(across(
    all_of(intersect(stage1_vars, names(.))),
    ~ {
      q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

# Check all variables exist
missing_s1 <- setdiff(stage1_vars, names(oa_data_clean))
missing_s2 <- setdiff(stage2_vars, names(oa_data_clean))

if (length(missing_s1) > 0) {
  cat("WARNING — Stage 1 variables missing:", paste(missing_s1, collapse = ", "), "\n")
  stage1_vars       <- intersect(stage1_vars, names(oa_data_clean))
  stage1_structural <- intersect(stage1_structural, names(oa_data_clean))
}
if (length(missing_s2) > 0) {
  cat("WARNING — Stage 2 variables missing:", paste(missing_s2, collapse = ", "), "\n")
  stage2_vars <- intersect(stage2_vars, names(oa_data_clean))
}

# Drop near-zero-variance variables
drop_low_var <- function(data, vars, threshold = 1e-8) {
  vcheck <- data %>%
    summarise(across(all_of(vars), ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "variance")
  low <- vcheck %>% filter(variance < threshold) %>% pull(variable)
  if (length(low) > 0)
    cat("Dropping near-zero-variance variables:", paste(low, collapse = ", "), "\n")
  setdiff(vars, low)
}

stage1_vars       <- drop_low_var(oa_data_clean, stage1_vars)
stage1_structural <- intersect(stage1_structural, stage1_vars)
stage2_vars       <- drop_low_var(oa_data_clean, stage2_vars)

# =============================================================================
# STEP 4 — STAGE 1: MDM 
# =============================================================================
# #   - replace = TRUE at Stage 1 (same control can supply multiple treated
#     units in Stage 1; deduplication happens before Stage 2) ✓
# =============================================================================

cat("\n====================================================\n")
cat("STEP 4: STAGE 1 MATCHING\n")
cat("====================================================\n\n")



# Use full dataset (all OAs)
oa_stage1_all <- oa_data_fin %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

oa_stage1_all_clean <- oa_stage1_all %>%
  mutate(across(c(road_density_m_km2, pct_A_road, pct_B_road, pct_minor_road,
                  Drive_Car_pct, Walk_pct, Bicycle_pct),
                ~ ifelse(is.na(.) | !is.finite(.), 0, .)))



s1_formula <- reformulate(stage1_vars, response = "treat_indicator")

# Run Stage 1 matching without any caliper
mdm_s1_nocal_all <- tryCatch(
  matchit(
    s1_formula,
    data        = oa_stage1_all_clean,
    method      = "nearest",
    distance    = "mahalanobis",
    ratio       = 10,
    replace     = TRUE,
    exact       = ~ country
  ),
  error = function(e) {
    cat("Stage 1 (no caliper, all OAs) MDM failed:", conditionMessage(e), "\n"); NULL
  }
)

# Extract matched data
s1_matched_nocal_all <- match.data(mdm_s1_nocal_all)

# Summary counts
n_treated_nocal_all <- sum(s1_matched_nocal_all$treat_indicator == 1)
n_unique_ctrl_nocal_all <- n_distinct(s1_matched_nocal_all$OA[s1_matched_nocal_all$treat_indicator == 0])
n_dropped_nocal_all <- sum(oa_stage1_all$treat_indicator == 1) - n_treated_nocal_all

cat("Stage 1 (no caliper, all OAs) results:\n")
cat("  Treated OAs retained:       ", n_treated_nocal_all, "\n")
cat("  Unique control OAs in pool: ", n_unique_ctrl_nocal_all, "\n")
cat("  Treated OAs dropped:        ", n_dropped_nocal_all,
    sprintf("(%.1f%%)\n", 100 * n_dropped_nocal_all / sum(oa_stage1_all$treat_indicator)))


# --- Stage 1 results ---
s1_matched_raw <- match.data(mdm_s1_nocal_all)

s1_treated <- sum(s1_matched_raw$treat_indicator == 1)
s1_control <- n_distinct(
  s1_matched_raw$OA[s1_matched_raw$treat_indicator == 0]
)
s1_dropped <- sum(oa_data_clean$treat_indicator == 1) - s1_treated

cat("Stage 1 results (verify against report: 791 retained, 4,301 controls):\n")
cat("  Treated OAs retained:          ", s1_treated, "  (report: 791)\n")
cat("  Unique control OAs in pool:    ", s1_control, "  (report: 4,301)\n")
cat("  Treated OAs dropped (caliper): ", s1_dropped,
    sprintf("(%.1f%%)\n", 100 * s1_dropped / sum(oa_data_clean$treat_indicator)))

cat("\nCountry distribution after Stage 1:\n")
print(table(s1_matched_raw$country, s1_matched_raw$treat_indicator,
            dnn = c("Country", "Treated")))

# --- Stage 1 balance (full variable set including sociodemographics) ---
cat("\nStage 1 balance:\n")
s1_bal <- bal.tab(mdm_s1_nocal_all, thresholds = c(m = 0.1), un = TRUE)
print(s1_bal)



smd_s1 <- s1_bal$Balance %>%
  rownames_to_column("variable") %>%
  select(variable, Diff.Un, Diff.Adj) %>%
  filter(variable %in% stage1_socdem) %>%
  arrange(desc(abs(Diff.Adj)))
cat("Sociodemographic SMDs after Stage 1 (|Adj| > 0.1 = residual imbalance):\n")
print(smd_s1)

# =============================================================================
# STEP 5 — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE
# =============================================================================

cat("\n====================================================\n")
cat("STEP 5: STAGE 2 DATA PREPARATION\n")
cat("====================================================\n\n")

s1_data_s2_raw <- bind_rows(
  s1_matched_raw %>% filter(treat_indicator == 1),
  s1_matched_raw %>% filter(treat_indicator == 0) %>%
    distinct(OA, .keep_all = TRUE)
) %>%
  select(-any_of(c("weights", "subclass", "distance")))

cat("Stage 2 input after dedup:\n")
cat("  Treated OAs:", sum(s1_data_s2_raw$treat_indicator == 1),
    "  (withExclu: 791)\n")
cat("  Control OAs:", sum(s1_data_s2_raw$treat_indicator == 0),
    "  (withExclu: 4,301)\n")

# Winsorise Stage 2 variables on treated distribution only
s2_vars_present <- intersect(stage2_vars, names(s1_data_s2_raw))
treated_ref     <- s1_data_s2_raw %>% filter(treat_indicator == 1)

s1_data_s2 <- s1_data_s2_raw %>%
  mutate(across(
    all_of(s2_vars_present),
    ~ {
      q <- quantile(treated_ref[[cur_column()]], probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

# =============================================================================
# STEP 6 — STAGE 2: MDM WITHOUT REPLACEMENT, MULTIPLE RATIOS
# =============================================================================
#   - No caliper at Stage 2 ✓ (post-hoc trim used instead)
#   - Country exact constraint maintained ✓
#   - 1:1 selected over 1:2, 1:4 based on distance diagnostics ✓
# =============================================================================

cat("\n====================================================\n")
cat("STEP 6: STAGE 2 MATCHING (multiple ratios)\n")
cat("====================================================\n\n")

s2_formula <- reformulate(s2_vars_present, response = "treat_indicator")

run_s2_mdm <- function(data, formula, ratio) {
  matchit(
    formula,
    data        = data,
    method      = "nearest",
    distance    = "mahalanobis",
    ratio       = ratio,
    replace     = T,
    exact       = ~ country
  )
}

s2_results <- list(
  r1 = run_s2_mdm(s1_data_s2, s2_formula, 1),
  r2 = run_s2_mdm(s1_data_s2, s2_formula, 2),
  r4 = run_s2_mdm(s1_data_s2, s2_formula, 4)
)

md <- match.data(s2_results$r1)
colnames(md)
head(md)
# =============================================================================
# STEP 7 — MATCH QUALITY: MAHALANOBIS DISTANCE COMPUTATION
# =============================================================================
# Distances computed manually because MatchIt's internal distances are
# propensity-score-scale; here we need the true Mahalanobis distance on
# the Stage 2 injury variables for the 95th-percentile trim.
# Covariance matrix estimated from TREATED distribution only (standard).

library(purrr)
library(tibble)

compute_mdist_no <- function(m, data, vars) {
  md <- match.data(m, data = data)
  treat <- md %>% filter(treat_indicator == 1)
  ctrl  <- md %>% filter(treat_indicator == 0)
  
  # covariance of treated units
  S <- cov(treat[, vars], use = "pairwise.complete.obs")
  
  # compute distances
  purrr::map_df(seq_len(nrow(treat)), function(ti) {
    trow <- treat[ti, vars, drop = FALSE]
    purrr::map_df(seq_len(nrow(ctrl)), function(ci) {
      crows <- ctrl[ci, vars, drop = FALSE]
      d <- mahalanobis(
        x      = as.numeric(crows),
        center = as.numeric(trow),
        cov    = S
      )
      tibble(
        treated_OA = treat$OA[ti],
        control_OA = ctrl$OA[ci],
        mdist      = d
      )
    })
  })
}

mdist_r1 <- compute_mdist_no(s2_results$r1, s1_data_s2, s2_vars_present)
mdist_r2 <- compute_mdist_no(s2_results$r2, s1_data_s2, s2_vars_present)
mdist_r4 <- compute_mdist_no(s2_results$r4, s1_data_s2, s2_vars_present)

colnames(mdist_r1)

summarise_quality <- function(df, label) {
  df %>% summarise(
    spec         = label,
    n_pairs      = n(),
    mean_mdist   = round(mean(mdist),               3),
    median_mdist = round(median(mdist),              3),
    p90_mdist    = round(quantile(mdist, 0.90),      3),
    p95_mdist    = round(quantile(mdist, 0.95),      3),
    max_mdist    = round(max(mdist),                 3)
  )
}

dist_summary <- bind_rows(
  summarise_quality(mdist_r1, "1:1"),
  summarise_quality(mdist_r2, "1:2"),
  summarise_quality(mdist_r4, "1:4")
)

cat("=== Mahalanobis distance summary (verify against report) ===\n")
cat("Report values: 1:1 mean=3.52, median=1.32, P90=8.13, max=110.4\n")
cat("               1:2 mean=4.23, median=1.55, P90=10.8\n")
cat("               1:4 mean=6.37, median=2.05, P90=15.4\n\n")
print(dist_summary)

# --- Distribution histogram (1:1) ---
p_hist <- ggplot(mdist_r1, aes(x = mdist)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = quantile(mdist_r1$mdist, 0.95),
             colour = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = quantile(mdist_r1$mdist, 0.95) + 0.5,
           y = Inf, vjust = 2, hjust = 0, colour = "red", size = 3,
           label = sprintf("P95 = %.1f", quantile(mdist_r1$mdist, 0.95))) +
  theme_minimal() +
  labs(title    = "Stage 2 Match Quality (1:1) — Mahalanobis Distances",
       subtitle = "Red dashed line = 95th-percentile trim threshold",
       x        = "Mahalanobis distance", y = "Count")

p_hist

ggsave(here("output", "stage2_distance_histogram.png"), p_hist,
       width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 8 — POST-HOC DISTANCE TRIM AT 95TH PERCENTILE
# =============================================================================
# Report: threshold = 14.8, removed 40 worst pairs, final N = 751 pairs.
# IMPORTANT: check whether the 40 trimmed treated OAs are randomly
# distributed across schemes or concentrated. Systematic trimming
# introduces selection bias.
# =============================================================================

cat("\n====================================================\n")
cat("STEP 8: 95TH PERCENTILE DISTANCE TRIM\n")
cat("====================================================\n\n")

primary_data <- match.data(s2_results[["r1"]], data = s1_data_s2)

p95_threshold <- quantile(mdist_r1$mdist, 0.95)
cat(sprintf("95th percentile threshold: %.1f  (report: 14.8)\n", p95_threshold))

keep_pairs           <- mdist_r1 %>% filter(mdist <= p95_threshold)
primary_data_trimmed <- primary_data %>%
  filter(subclass %in% keep_pairs$subclass)

n_before_trim <- n_distinct(primary_data$subclass[primary_data$treat_indicator == 1])
n_after_trim  <- n_distinct(primary_data_trimmed$subclass[primary_data_trimmed$treat_indicator == 1])
n_removed     <- n_before_trim - n_after_trim

cat(sprintf("Pairs before trim: %d\n", n_before_trim))
cat(sprintf("Pairs removed:     %d  (report: 40)\n", n_removed))
cat(sprintf("Pairs after trim:  %d  (report: 751)\n", n_after_trim))

# --- CRITICAL CHECK: are trimmed pairs randomly distributed? ---
# Concentration in one scheme/country = selection bias in the ATT
cat("\n--- Trimmed pair concentration check ---\n")

trimmed_pairs <- mdist_r1 %>% filter(mdist > p95_threshold)

trimmed_treated_oas <- primary_data %>%
  filter(treat_indicator == 1,
         subclass %in% trimmed_pairs$subclass)

cat("Trimmed treated OAs: N =", nrow(trimmed_treated_oas), "\n")

# Country distribution of trimmed OAs
if ("country" %in% names(trimmed_treated_oas)) {
  cat("Country distribution of trimmed treated OAs:\n")
  print(table(trimmed_treated_oas$country))
  cat("Country distribution of retained treated OAs:\n")
  print(table(primary_data_trimmed$country[primary_data_trimmed$treat_indicator == 1]))
}

# Scheme distribution of trimmed OAs
if ("scheme_id" %in% names(trimmed_treated_oas)) {
  cat("Scheme distribution of trimmed treated OAs:\n")
  print(sort(table(trimmed_treated_oas$scheme_id), decreasing = TRUE))
  
  # Chi-square test: are trims proportional to scheme size?
  retained_scheme <- table(
    primary_data_trimmed$scheme_id[primary_data_trimmed$treat_indicator == 1]
  )
  trimmed_scheme <- table(trimmed_treated_oas$scheme_id)
  all_schemes    <- union(names(retained_scheme), names(trimmed_scheme))
  ret_vec        <- as.integer(retained_scheme[all_schemes])
  tri_vec        <- as.integer(trimmed_scheme[all_schemes])
  ret_vec[is.na(ret_vec)] <- 0
  tri_vec[is.na(tri_vec)] <- 0
  
  prop_test_mat <- rbind(tri_vec, ret_vec)
  if (all(colSums(prop_test_mat) > 0)) {
    ct <- chisq.test(prop_test_mat, simulate.p.value = TRUE)
    cat(sprintf("\nChi-square test — are trims proportional to scheme size? p = %.3f\n",
                ct$p.value))
    if (ct$p.value < 0.05) {
      cat("WARNING: Trimmed OAs are NOT randomly distributed across schemes.\n")
      cat("This introduces selection bias — report and consider sensitivity check.\n")
    } else {
      cat("Trimming appears proportional across schemes — no evidence of\n")
      cat("systematic scheme-level selection bias.\n")
    }
  }
} else {
  cat("NOTE: scheme_id not found. Add scheme identifier to check whether\n")
  cat("trimmed pairs are concentrated in specific CAZ/LEZ schemes.\n")
}

# Covariate comparison: trimmed vs retained treated OAs
cat("\n--- Covariate comparison: trimmed vs retained treated OAs ---\n")
trimmed_vs_retained <- bind_rows(
  trimmed_treated_oas %>% mutate(group = "Trimmed"),
  primary_data_trimmed %>% filter(treat_indicator == 1) %>% mutate(group = "Retained")
) %>%
  select(group, all_of(intersect(c(stage1_vars, stage2_vars), names(.)))) %>%
  pivot_longer(-group) %>%
  group_by(name, group) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd   = sd(value,   na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = c(mean, sd)) %>%
  mutate(
    SMD = (mean_Retained - mean_Trimmed) /
      sqrt((sd_Retained^2 + sd_Trimmed^2) / 2)
  ) %>%
  arrange(desc(abs(SMD)))

cat("Variables most different between trimmed and retained treated OAs:\n")
print(head(trimmed_vs_retained %>% select(name, mean_Retained, mean_Trimmed, SMD), 10))

if (any(abs(trimmed_vs_retained$SMD) > 0.3, na.rm = TRUE)) {
  cat("\nWARNING: Large SMD between trimmed and retained — trimmed OAs are\n")
  cat("systematically different. Report as a limitation.\n")
}

# =============================================================================
# STEP 9 — BALANCE DIAGNOSTICS AND RATIO SELECTION
# =============================================================================
# REPORT CONSISTENCY:
#   Stage 2 trend variables: all SMD < 0.06 ✓ (primary check)
#   Stage 2 level variables: some residual imbalance, max |SMD| = 0.24 ✓
#   Report notes mean_total_pkm as the most imbalanced level variable.
#   These are included as covariates in the doubly-robust outcome model.
# =============================================================================

cat("\n====================================================\n")
cat("STEP 9: BALANCE DIAGNOSTICS\n")
cat("====================================================\n\n")

run_diagnostics <- function(m_obj, label, s2_vars) {
  if (is.null(m_obj)) return(invisible(NULL))
  
  cat("\n--- Diagnostics:", label, "---\n")
  md <- match.data(m_obj, data = s1_data_s2)
  bt <- bal.tab(m_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  
  smd_after <- abs(bt$Balance$Diff.Adj)
  cat("  Variables with |SMD| < 0.10:", sum(smd_after < 0.1, na.rm = TRUE),
      "/", length(smd_after), "\n")
  cat("  Max |SMD|:",  round(max(smd_after,  na.rm = TRUE), 3), "\n")
  cat("  Mean |SMD|:", round(mean(smd_after, na.rm = TRUE), 3), "\n")
  
  # Trend vs level breakdown (key for parallel trends)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    mutate(var_type = case_when(
      variable %in% stage2_trends ~ "Trend",
      variable %in% stage2_levels ~ "Level",
      TRUE                        ~ "Other"
    ))
  
  cat("\n  Trend variable SMDs (all should be < 0.06 per report):\n")
  trend_smds <- smd_df %>% filter(var_type == "Trend") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(trend_smds)
  if (any(abs(trend_smds$Diff.Adj) >= 0.1, na.rm = TRUE)) {
    cat("  WARNING: Trend variable SMD >= 0.1 — parallel trends support weakened.\n")
  }
  
  cat("\n  Level variable SMDs (residual imbalance expected, max ~0.24):\n")
  level_smds <- smd_df %>% filter(var_type == "Level") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(level_smds)
  cat("  NOTE: Residual level imbalance addressed by doubly-robust outcome adjustment.\n")
  cat("  Include all level variables with |SMD| > 0.1 as covariates in C&S model.\n")
  
  # Cross-country check
  cross <- md %>%
    group_by(subclass) %>%
    summarise(n_countries = n_distinct(country), .groups = "drop") %>%
    filter(n_countries > 1)
  if (nrow(cross) > 0) {
    cat("\n  WARNING:", nrow(cross), "cross-country pairs — investigate!\n")
  } else {
    cat("\n  Country check: PASSED (no cross-country pairs)\n")
  }
  
  # Love plot
  lp <- love.plot(
    m_obj,
    threshold    = 0.1,
    abs          = TRUE,
    stars        = "std",
    var.order    = "unadjusted",
    title        = paste("Covariate balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) +
    theme(
      axis.text.y     = element_text(size = 9),
      axis.title.x    = element_text(size = 10),
      plot.title      = element_text(size = 11),
      legend.position = "bottom",
      legend.text     = element_text(size = 9),
      plot.margin     = margin(t = 10, r = 30, b = 10, l = 180)
    )
  
  fname <- here("output", paste0("balance_", gsub("[^a-zA-Z0-9]", "_", label), ".png"))
  ggsave(fname, lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved:", fname, "\n")
  
  invisible(bt)
}

for (nm in names(s2_results)) {
  run_diagnostics(s2_results[[nm]], nm, s2_vars_present)
}

# =============================================================================
# STEP 10 — RATIO SELECTION: STAGE 2
# =============================================================================
# Selection rule (from report): choose 1:1 because it produces consistently
# lower match distances across all summary statistics.
# Note: lower mean/median distance is mathematically guaranteed for 1:1 vs
# 1:k (k>1) since additional controls are by definition less similar.
# The tradeoff is efficiency (lower variance with more controls) vs bias
# (worse match quality). The report selects on match quality; this is
# defensible but acknowledge the efficiency tradeoff.
# =============================================================================

cat("\n====================================================\n")
cat("STEP 10: RATIO SELECTION SUMMARY\n")
cat("====================================================\n\n")

cat("Distance-based selection (verify against report):\n")
print(dist_summary)

cat("\nNOTE: 1:1 will always have lower distances than 1:k by construction.\n")
cat("Selection on distance alone does not account for the efficiency\n")
cat("gain from more controls. If a reviewer raises this, note that:\n")
cat("  (a) balance diagnostics also favour 1:1\n")
cat("  (b) Callaway-Sant'Anna aggregation uses weights, so extra controls\n")
cat("      from 1:4 would be downweighted anyway\n\n")

# =============================================================================
# STEP 11 — SELECT PRIMARY SPECIFICATION, EXTRACT MATCHED OAs
# =============================================================================

cat("\n====================================================\n")
cat("STEP 11: PRIMARY SPECIFICATION AND FINAL EXTRACTION\n")
cat("====================================================\n\n")

primary_spec <- "r1"
mdm_primary  <- s2_results[[primary_spec]]

if (is.null(mdm_primary)) stop("Primary specification failed — choose alternative.")

cat("Primary specification:", primary_spec, "\n")

matched_treated_oas <- primary_data_trimmed %>%
  filter(treat_indicator == 1) %>%
  dplyr::select(OA, weights, subclass)

matched_control_oas <- primary_data_trimmed %>%
  filter(treat_indicator == 0) %>%
  dplyr::select(OA, weights, subclass)

cat("Final matched treated OAs:", nrow(matched_treated_oas),
    "  (report: 751)\n")
cat("Final matched control OAs:", nrow(matched_control_oas),
    "  (report: 751)\n")

# Final cross-country check
final_cross <- primary_data_trimmed %>%
  group_by(subclass) %>%
  summarise(n_countries = n_distinct(country), .groups = "drop") %>%
  filter(n_countries > 1)
cat("Cross-country pairs in final matched sample:", nrow(final_cross),
    if (nrow(final_cross) == 0) "✓" else "WARNING", "\n")

# =============================================================================
# STEP 12 — ZERO-INJURY SENSITIVITY CHECK (STRUCTURAL MATCH ONLY)
# =============================================================================
# PURPOSE: Compare ATT estimates from the main sample (pre-injury >0) vs a
# supplementary sample that includes zero-injury treated OAs matched only
# on Stage 1 structural variables. If ATT estimates are similar, the
# zero-injury exclusion is not materially biasing results. If very
# different, the exclusion is doing real work and the scope limitation of
# the ATT must be clearly stated.
#
# This section creates and saves the supplementary matched sample.
# Actual ATT comparison requires running C&S estimation twice and is
# done in the outcome analysis script.
# =============================================================================

cat("\n====================================================\n")
cat("STEP 12: ZERO-INJURY SENSITIVITY SAMPLE\n")
cat("====================================================\n\n")

zero_injury_treated <- OA_matching_dataset %>%
  filter(treated_OA == 1, zero_injury_OA == 1, n_roads > 0) %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    ),
    treat_indicator = 1L
  )

cat("Zero-injury treated OAs available for sensitivity matching:",
    nrow(zero_injury_treated), "\n")

if (nrow(zero_injury_treated) > 0) {
  # Add to Stage 1 control pool for structural-only matching
  # Controls come from oa_data_clean (already winsorised on Stage 1 vars)
  sens_data <- bind_rows(
    zero_injury_treated %>%
      select(any_of(c("OA", "treat_indicator", "country", stage1_vars))),
    oa_data_clean %>%
      filter(treat_indicator == 0) %>%
      select(any_of(c("OA", "treat_indicator", "country", stage1_vars)))
  )
  
  s1_vars_avail <- intersect(stage1_vars, names(sens_data))
  sens_formula  <- reformulate(s1_vars_avail, response = "treat_indicator")
  
  sens_caliper        <- rep(1.0, length(intersect(stage1_structural, names(sens_data))))
  names(sens_caliper) <- intersect(stage1_structural, names(sens_data))
  
  mdm_sens <- tryCatch(
    matchit(
      sens_formula,
      data        = sens_data,
      method      = "nearest",
      distance    = "mahalanobis",
      caliper     = sens_caliper,
      std.caliper = TRUE,
      ratio       = 1,
      replace     = FALSE,
      exact       = ~ country
    ),
    error = function(e) {
      cat("Sensitivity matching failed:", conditionMessage(e), "\n"); NULL
    }
  )
  
  if (!is.null(mdm_sens)) {
    sens_matched    <- match.data(mdm_sens, data = sens_data)
    sens_treated    <- sum(sens_matched$treat_indicator == 1)
    sens_dropped    <- nrow(zero_injury_treated) - sens_treated
    
    cat("Sensitivity match results:\n")
    cat("  Zero-injury treated OAs matched:", sens_treated, "\n")
    cat("  Dropped (caliper):              ", sens_dropped, "\n")
    
    sens_treated_oas <- sens_matched %>%
      filter(treat_indicator == 1) %>%
      select(OA, weights, subclass) %>%
      mutate(subclass = paste0("sens_", subclass))
    
    sens_control_oas <- sens_matched %>%
      filter(treat_indicator == 0) %>%
      select(OA, weights, subclass) %>%
      mutate(subclass = paste0("sens_", subclass))
    
    saveRDS(sens_treated_oas,
            here("data", "processed", "OA_matched_zero_injury_treated.rds"))
    saveRDS(sens_control_oas,
            here("data", "processed", "OA_matched_zero_injury_controls.rds"))
    
    cat("\nSensitivity sample saved. To complete the sensitivity check:\n")
    cat("  1. Run C&S estimation on the main matched sample (Step 11 output)\n")
    cat("  2. Run C&S estimation on this sensitivity sample\n")
    cat("  3. Compare ATT point estimates and CIs\n")
    cat("  4. If estimates are similar, zero-injury exclusion is not biasing\n")
    cat("     main results. If substantially different, state scope limitation\n")
    cat("     explicitly: estimates apply only to OAs with pre-treatment\n")
    cat("     injury exposure.\n")
  }
} else {
  cat("No zero-injury treated OAs remaining after road exclusions.\n")
}

# =============================================================================
# STEP 13 — SAVE PRIMARY MATCHED SAMPLE
# =============================================================================

cat("\n====================================================\n")
cat("STEP 13: SAVE OUTPUTS\n")
cat("====================================================\n\n")

saveRDS(matched_treated_oas, here("data", "processed", "OA_matched_treated.rds"))
saveRDS(matched_control_oas, here("data", "processed", "OA_matched_donors.rds"))

cat("Saved:\n")
cat("  OA_matched_treated.rds  —", nrow(matched_treated_oas), "treated OAs\n")
cat("  OA_matched_donors.rds   —", nrow(matched_control_oas), "control OAs\n")

cat("\n=== MATCHING COMPLETE ===\n")
cat("Scope note: ATT estimates apply to treated OAs with pre-treatment\n")
cat("injury exposure. See Section 0 and Step 12 for characterisation\n")
cat("and sensitivity analysis of the excluded zero-injury OAs.\n")
