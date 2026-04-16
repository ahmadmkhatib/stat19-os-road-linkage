# =============================================================================
# DIAGNOSTICS — SKEWNESS, SCOTLAND, ZERO-INJURY OA PERSISTENCE
# =============================================================================
#
# 
#    Q1. Scotland: are the problematic matches actually contributing outcome
#        variation? Produce N-roads-per-OA table/map after removing zero-
#        injury roads across the whole period.
#    Q2. Zero-injury OA structural distinctiveness: does it persist after
#        removing zero-total-period roads? I.e. is this a matching problem
#        or a modelling problem?
#    Q3. Are zero-pre-injury OAs structurally unusual ACROSS the full sample
#        (not just Scotland)?
#    Q4. Interaction / separate model tests for zero vs non-zero pre-injury
#        OAs in the final DiD.
#
# =============================================================================

library(tidyverse)
library(here)
library(arrow)
library(sf)
library(lubridate)
library(fixest)
library(ggrepel)
library(patchwork)   # 

road_panel_model <- open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

road_attributes <- st_read(
  here("data", "processed", "road_attributes_OA.gpkg"),
  quiet = TRUE
) %>%
  st_drop_geometry() %>%
  select(identifier, OA)

road_panel_model <- road_panel_model %>%
  left_join(road_attributes, by = "identifier")

road_panel_model <- road_panel_model %>%
  mutate(
    total_inj_adj_all =
      total_inj_adj_Car.Van +
      total_inj_adj_Pedestrian +
      total_inj_adj_Cyclist +
      total_inj_adj_Other,
    country = case_when(
      substr(OA, 1, 1) == "E" ~ "England",
      substr(OA, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  )

# --- CAZ start dates ---
caz <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
)

caz_dates <- caz %>%
  st_drop_geometry() %>%
  mutate(caz_start_date = dmy(startDt)) %>%
  group_by(scheme) %>%
  summarise(
    caz_start_date = min(caz_start_date, na.rm = TRUE),
    .groups = "drop"
  )

road_panel_model <- road_panel_model %>%
  left_join(caz_dates, by = "scheme")

# =============================================================================
# PART 1 — SKEWNESS ASSESSMENT FOR MATCHING VARIABLES
# =============================================================================
# Winsorisation trims extreme values but does NOT change the shape of the
# distribution between the 1st and 99th percentiles.  For Mahalanobis
# distance matching (MDM), skewness matters in two ways:
#
#   (a) Covariance matrix estimation: the sample covariance is sensitive to
#       skewness because it is a second-moment estimator. Heavy right tails
#       inflate variance estimates and compress the relative contribution of
#       that variable to the Mahalanobis distance.
#
#   (b) "Closeness" is not symmetric for skewed variables: a unit 2 SD above
#       the mean is much further in raw units than a unit 2 SD below it, but
#       MDM treats them symmetrically. This can cause poor matches on the
#       right tail.
#
# RULE OF THUMB applied here:
#   |skewness| > 2 AND max/median > 10 after winsorisation → log-transform.
#   These thresholds are not universal; the key test is whether log-
#   transforming produces meaningfully more symmetric distributions and
#   better post-matching balance. We test both and compare SMDs.
#
# NOTE: log transformation requires a shift for variables with zeros.
#   We use log1p (log(x+1)) for count-like variables and
#   log(x + 0.001) for rate variables that can be exactly zero.

cat("\n====================================================\n")
cat("PART 1: SKEWNESS ASSESSMENT FOR MATCHING VARIABLES\n")
cat("====================================================\n\n")

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Load the winsorised datasets from the matching pipeline if available,
# otherwise recreate them here for diagnosis
stage1_vars <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)

stage2_trends <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm", "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm", "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_total_pkm"
)
all_matching_vars <- c(stage1_vars, stage2_trends, stage2_levels)

# Compute skewness (moment-based: E[(X-mu)^3] / sd^3)
skew <- function(x) {
  x   <- x[!is.na(x)]
  n   <- length(x)
  mu  <- mean(x)
  s   <- sd(x)
  if (s == 0) return(NA_real_)
  sum((x - mu)^3) / (n * s^3)
}

skew_table <- map_df(
  intersect(all_matching_vars, names(OA_matching_dataset)),
  function(v) {
    x_raw <- OA_matching_dataset[[v]]
    # Winsorise at 1/99
    q  <- quantile(x_raw, c(0.01, 0.99), na.rm = TRUE)
    x_w <- pmin(pmax(x_raw, q[1]), q[2])
    tibble(
      variable     = v,
      stage        = case_when(
        v %in% stage1_vars     ~ "Stage 1",
        v %in% stage2_trends   ~ "Stage 2 trend",
        v %in% stage2_levels   ~ "Stage 2 level"
      ),
      skew_raw     = round(skew(x_raw), 2),
      skew_winsor  = round(skew(x_w), 2),
      median_raw   = round(median(x_raw, na.rm = TRUE), 4),
      max_raw      = round(max(x_raw, na.rm = TRUE), 4),
      max_med_ratio = round(max(x_raw, na.rm = TRUE) /
                              (median(x_raw, na.rm = TRUE) + 1e-9), 1),
      flag_log     = abs(skew(x_w)) > 2 & max(x_raw, na.rm = TRUE) /
        (median(x_raw, na.rm = TRUE) + 1e-9) > 10
    )
  }
)

cat("Skewness table (raw vs winsorised; flag_log = TRUE → consider log):\n\n")
print(skew_table %>% arrange(desc(abs(skew_winsor))), n = Inf)

cat("\nVariables flagged for potential log transformation:\n")
print(skew_table %>% filter(flag_log) %>% select(variable, stage, skew_raw, skew_winsor, max_med_ratio))

# --- Visualise before/after for flagged variables ---
flagged_vars <- skew_table %>% filter(flag_log) %>% pull(variable)

if (length(flagged_vars) > 0) {
  
  plots_skew <- map(flagged_vars, function(v) {
    x <- OA_matching_dataset[[v]]
    q <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
    x_w <- pmin(pmax(x, q[1]), q[2])
    # log1p is safe for zero values for rate-like variables
    x_log <- log1p(pmax(x_w, 0))
    df <- bind_rows(
      tibble(val = x_w,   transform = "Winsorised only"),
      tibble(val = x_log, transform = "Winsorised + log1p")
    )
    ggplot(df, aes(x = val)) +
      geom_histogram(bins = 50, fill = "#2980B9", alpha = 0.7) +
      facet_wrap(~ transform, scales = "free") +
      labs(title = v, x = NULL, y = "Count") +
      theme_minimal(base_size = 9)
  })
  
  # Save grid of skewness plots
  patchwork_skew <- wrap_plots(plots_skew, ncol = 2)
  ggsave(here("output", "skewness_flagged_vars.png"),
         patchwork_skew, width = 14, height = 3 * ceiling(length(flagged_vars) / 2),
         dpi = 200)
  cat("\nSkewness distribution plots saved to output/skewness_flagged_vars.png\n")
}

cat("\n--- RECOMMENDATION ---\n")
cat("For MDM specifically, the key question is not whether variables are\n")
cat("skewed in the marginal distribution, but whether the Mahalanobis distance\n")
cat("computed in the raw space is a meaningful measure of 'closeness' for the\n")
cat("matching purpose.\n\n")
cat("For the STAGE 1 structural variables:\n")
cat("  road_density_m_km2, road_length_km, pop_density, dist_citycentre:\n")
cat("  These are right-skewed even after winsorisation (typical for area-level\n")
cat("  spatial data). Log-transforming them before matching is RECOMMENDED\n")
cat("  because it makes 'closeness' more meaningful — a doubling of road density\n")
cat("  is equally meaningful whether you are at 5,000 or 10,000 m/km2, but MDM\n")
cat("  in raw space treats the latter as twice as 'far'.\n\n")
cat("For STAGE 2 trend variables:\n")
cat("  These are log-slopes (already on log scale). Skewness arises from\n")
cat("  structural zeros, not from the scale. Do NOT log-transform again.\n")
cat("  The zero-spike treatment (dropping variables where >60% are structural\n")
cat("  zeros) from FIX 4 is the correct approach here.\n\n")
cat("For STAGE 2 level variables (mean_*_pkm):\n")
cat("  These are injury rates per km — right-skewed due to OAs with high\n")
cat("  rates. Log1p transformation is RECOMMENDED for the same reason as\n")
cat("  Stage 1 density variables.\n\n")
cat("ACTION: Re-run matching pipeline (pipeline_v2.R) with log1p applied\n")
cat("  to the flagged variables AFTER winsorisation in both Stage 1 and\n")
cat("  Stage 2. The recommended call sequence is:\n")
cat("    x_transformed <- log1p(pmax(x_winsorised, 0))\n")
cat("  This handles exact zeros safely. Verify that SMDs improve after\n")
cat("  log transformation using the balance tables.\n\n")


cat("\n====================================================\n")
cat("PART 2: FIXED DIAGNOSTIC CODE\n")
cat("====================================================\n\n")

# --- Convert caz_start_date to decimal year (consistent with quarter_year) ---
road_panel_model <- road_panel_model %>%
  mutate(
    # decimal year: 2019.0 = Jan, 2019.25 = Apr, 2019.5 = Jul, 2019.75 = Oct
    caz_start_decimal = year(caz_start_date) +
      (month(caz_start_date) - 1) / 12,
    # BUG 1 FIX: both sides are now numeric decimal years
    pre_period = quarter_year < caz_start_decimal
  )

cat("Pre/post period distribution after bug fix:\n")
print(table(road_panel_model$pre_period, useNA = "ifany"))
cat("(Should now have both TRUE and FALSE rows)\n\n")

# --- Zero-injury roads across the WHOLE period ---
road_totals <- road_panel_model %>%
  group_by(identifier) %>%
  summarise(
    total_injuries_alltime = sum(total_inj_adj_all, na.rm = TRUE),
    .groups = "drop"
  )

zero_roads_alltime <- road_totals %>% filter(total_injuries_alltime == 0)
cat("Roads with zero injuries across WHOLE period:",
    nrow(zero_roads_alltime), "/", nrow(road_totals),
    sprintf("(%.1f%%)\n\n", 100 * nrow(zero_roads_alltime) / nrow(road_totals)))

# Keep only roads with at least one injury ever
roads_retained <- road_totals %>% filter(total_injuries_alltime > 0)
cat("Roads retained (at least one injury ever):", nrow(roads_retained), "\n\n")

analysis_sample <- road_panel_model %>%
  filter(identifier %in% roads_retained$identifier)

# --- Pre-treatment injuries per OA (using fixed pre_period flag) ---
oa_pre_injury <- analysis_sample %>%
  filter(pre_period == TRUE) %>%
  group_by(OA) %>%
  summarise(
    pre_injuries       = sum(total_inj_adj_all, na.rm = TRUE),
    n_roads_pre        = n_distinct(identifier),
    country            = first(country),
    .groups = "drop"
  ) %>%
  mutate(zero_pre_injury_OA = as.integer(pre_injuries == 0))

cat("Zero-pre-injury OA rate AFTER bug fix AND after removing zero-alltime roads:\n")
cat("  Total OAs in analysis sample:", nrow(oa_pre_injury), "\n")
cat("  OAs with zero pre-period injuries:", sum(oa_pre_injury$zero_pre_injury_OA), "\n")
cat("  Percent:", round(mean(oa_pre_injury$zero_pre_injury_OA) * 100, 1), "%\n\n")

cat("By country:\n")
oa_pre_injury %>%
  group_by(country) %>%
  summarise(
    n_OAs          = n(),
    n_zero_pre     = sum(zero_pre_injury_OA),
    pct_zero_pre   = round(mean(zero_pre_injury_OA) * 100, 1),
    .groups = "drop"
  ) %>% print()

# Merge zero_pre_injury_OA back into analysis_sample
analysis_sample <- analysis_sample %>%
  left_join(oa_pre_injury %>% select(OA, pre_injuries, zero_pre_injury_OA),
            by = "OA")

# --- BUG 2 FIX: post-treatment injuries for zero-pre-injury OAs ---
post_injury_check <- analysis_sample %>%
  filter(zero_pre_injury_OA == 1, pre_period == FALSE) %>%  # BUG 2 FIX: use pre_period
  group_by(OA) %>%
  summarise(
    post_injuries = sum(total_inj_adj_all, na.rm = TRUE),
    n_roads_post  = n_distinct(identifier),
    country       = first(country),
    .groups = "drop"
  )

cat("\n--- Do zero-pre-injury OAs contribute injuries post-treatment? ---\n")
cat("Zero-pre-injury OAs with post-treatment injuries:",
    sum(post_injury_check$post_injuries > 0),
    "/", nrow(post_injury_check),
    sprintf("(%.1f%%)\n", 100 * mean(post_injury_check$post_injuries > 0)))
cat("By country:\n")
post_injury_check %>%
  group_by(country) %>%
  summarise(
    n_OAs          = n(),
    n_with_post    = sum(post_injuries > 0),
    pct_with_post  = round(mean(post_injuries > 0) * 100, 1),
    mean_post_inj  = round(mean(post_injuries), 2),
    .groups = "drop"
  ) %>% print()

# =============================================================================
# PART 3 — COLLEAGUE QUESTIONS
# =============================================================================

cat("\n====================================================\n")
cat("PART 3: COLLEAGUE QUESTIONS\n")
cat("====================================================\n\n")

# ----------------------------------------------------------------------------
# Q1. SCOTLAND DIAGNOSIS: What % of the analysis sample comes from the
#     problematic matched Scottish zero-pre-injury OAs?
#     "If it's very low, no issue."
# ----------------------------------------------------------------------------

cat("--- Q1: Scotland contribution to final analysis sample ---\n\n")

# Roads per OA in the retained (non-zero-alltime) sample
oa_road_counts <- analysis_sample %>%
  group_by(OA, country) %>%
  summarise(
    n_roads_analysis    = n_distinct(identifier),
    total_inj_alltime   = sum(total_inj_adj_all, na.rm = TRUE),
    pre_injuries        = sum(total_inj_adj_all[pre_period == TRUE], na.rm = TRUE),
    post_injuries       = sum(total_inj_adj_all[pre_period == FALSE], na.rm = TRUE),
    zero_pre_injury_OA  = first(zero_pre_injury_OA),
    .groups = "drop"
  )

scotland_summary <- oa_road_counts %>%
  group_by(country) %>%
  summarise(
    n_OAs            = n(),
    n_roads          = sum(n_roads_analysis),
    pct_roads        = round(100 * sum(n_roads_analysis) /
                               sum(oa_road_counts$n_roads_analysis), 1),
    total_injuries   = round(sum(total_inj_alltime), 1),
    pct_injuries     = round(100 * sum(total_inj_alltime) /
                               sum(oa_road_counts$total_inj_alltime), 1),
    .groups = "drop"
  )
cat("Overall country breakdown (after removing zero-alltime-injury roads):\n")
print(scotland_summary)

# Focus on the PROBLEMATIC subset: Scottish zero-pre-injury OAs
# These are the OAs whose matched controls had near-degenerate weights (max 93)
scotland_zero_pre <- oa_road_counts %>%
  filter(country == "Scotland", zero_pre_injury_OA == 1)

total_roads_sample <- sum(oa_road_counts$n_roads_analysis)
total_inj_sample   <- sum(oa_road_counts$total_inj_alltime)

cat("\nProblematic subset: Scottish zero-pre-injury OAs\n")
cat("  N OAs:             ", nrow(scotland_zero_pre), "\n")
cat("  N roads:           ", sum(scotland_zero_pre$n_roads_analysis), "\n")
cat("  % of total roads:  ", round(100 * sum(scotland_zero_pre$n_roads_analysis) /
                                     total_roads_sample, 2), "%\n")
cat("  Total injuries:    ", round(sum(scotland_zero_pre$total_inj_alltime), 1), "\n")
cat("  % of total injuries:", round(100 * sum(scotland_zero_pre$total_inj_alltime) /
                                      total_inj_sample, 2), "%\n")
cat("  Post-treatment injuries:", round(sum(scotland_zero_pre$post_injuries), 1), "\n")
cat("  OAs with any post injuries:",
    sum(scotland_zero_pre$post_injuries > 0), "/", nrow(scotland_zero_pre), "\n\n")

cat("INTERPRETATION:\n")
cat("  If % of total injuries < 2%: problematic matched set contributes\n")
cat("  negligible outcome variation — Scotland sparsity is a matching\n")
cat("  artefact that does not materially affect ATT estimation.\n")
cat("  If % > 5%: exclusion sensitivity analysis is warranted.\n\n")

# Produce a summary table for the paper
oa_road_counts_summary <- oa_road_counts %>%
  mutate(
    oa_type = case_when(
      country == "Scotland" & zero_pre_injury_OA == 1 ~ "Scotland: zero-pre-injury",
      country == "Scotland" & zero_pre_injury_OA == 0 ~ "Scotland: non-zero-pre-injury",
      country == "England"  & zero_pre_injury_OA == 1 ~ "England: zero-pre-injury",
      country == "England"  & zero_pre_injury_OA == 0 ~ "England: non-zero-pre-injury"
    )
  ) %>%
  group_by(oa_type) %>%
  summarise(
    n_OAs              = n(),
    n_roads            = sum(n_roads_analysis),
    pct_total_roads    = round(100 * sum(n_roads_analysis) / total_roads_sample, 2),
    total_inj          = round(sum(total_inj_alltime), 1),
    pct_total_inj      = round(100 * sum(total_inj_alltime) / total_inj_sample, 2),
    mean_roads_per_OA  = round(mean(n_roads_analysis), 1),
    .groups = "drop"
  )
cat("Summary table by OA type:\n")
print(oa_road_counts_summary, n = Inf)

# Roads per OA distribution plot
p_roads_oa <- oa_road_counts %>%
  ggplot(aes(x = n_roads_analysis, fill = country)) +
  geom_histogram(bins = 40, position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("England" = "#2980B9", "Scotland" = "#E74C3C")) +
  labs(
    title = "Distribution of roads per OA in analysis sample",
    subtitle = "After removing roads with zero injuries across whole period",
    x = "Number of roads with ≥1 injury",
    y = "Count of OAs",
    fill = "Country"
  ) +
  theme_minimal()
ggsave(here("output", "roads_per_OA_country.png"),
       p_roads_oa, width = 10, height = 5, dpi = 200)
cat("Roads-per-OA plot saved.\n\n")

# ----------------------------------------------------------------------------
# Q2. ZERO-PRE-INJURY OA STRUCTURAL DISTINCTIVENESS: Does it persist after
#     removing zero-alltime-injury roads?
#     "Are retained roads from zero-pre-injury OAs still structurally distinct?"
# ----------------------------------------------------------------------------

cat("--- Q2: Structural distinctiveness of zero-pre-injury OAs ---\n")
cat("    (after removing zero-alltime-injury roads)\n\n")

OA_matching_dataset_full <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Join zero_pre_injury_OA flag to matching dataset
# (only for OAs that appear in the analysis sample after road exclusions)
oa_flags <- oa_pre_injury %>% select(OA, zero_pre_injury_OA)

match_with_flags <- OA_matching_dataset_full %>%
  inner_join(oa_flags, by = "OA") %>%
  filter(treated_OA == 1)   # focus on treated OAs — that's what affects DiD

stage1_structural <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)

cat("Treated OAs in matching dataset with road-exclusion flag:\n")
print(table(match_with_flags$zero_pre_injury_OA))
cat("\n")

# SMD between zero-pre and non-zero-pre treated OAs
smd_zero_vs_nonzero <- map_df(
  intersect(stage1_structural, names(match_with_flags)),
  function(v) {
    g0 <- match_with_flags[[v]][match_with_flags$zero_pre_injury_OA == 0]
    g1 <- match_with_flags[[v]][match_with_flags$zero_pre_injury_OA == 1]
    pooled_sd <- sqrt((var(g0, na.rm = TRUE) + var(g1, na.rm = TRUE)) / 2)
    tibble(
      variable = v,
      mean_nonzero_pre = round(mean(g0, na.rm = TRUE), 3),
      mean_zero_pre    = round(mean(g1, na.rm = TRUE), 3),
      SMD              = round((mean(g0, na.rm = TRUE) - mean(g1, na.rm = TRUE)) /
                                 (pooled_sd + 1e-9), 3)
    )
  }
) %>% arrange(desc(abs(SMD)))

cat("SMD between zero-pre-injury and non-zero-pre-injury TREATED OAs\n")
cat("(after removing zero-alltime-injury roads from the analysis sample):\n\n")
print(smd_zero_vs_nonzero, n = Inf)

n_large_smd <- sum(abs(smd_zero_vs_nonzero$SMD) > 0.2)
cat("\nVariables with |SMD| > 0.2:", n_large_smd, "\n")
cat("Variables with |SMD| > 0.5:", sum(abs(smd_zero_vs_nonzero$SMD) > 0.5), "\n\n")

cat("INTERPRETATION:\n")
cat("  If structural differences persist at |SMD| > 0.5 for key variables\n")
cat("  (road_length, pop_density, pct_minor_road), the zero-pre-injury\n")
cat("  subgroup remains structurally distinct even in the retained sample.\n")
cat("  This matters for the DiD ONLY if these OAs contribute non-trivial\n")
cat("  post-treatment injury variation (see Q1 above).\n\n")

# Q3: Is this England-wide or mainly Scotland?
cat("--- Q3: Zero-pre-injury OA structural issue — England vs Scotland ---\n\n")

match_with_flags %>%
  filter(treated_OA == 1) %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE ~ "Unknown"
    )
  ) %>%
  group_by(country, zero_pre_injury_OA) %>%
  summarise(
    n_OAs             = n(),
    mean_road_length  = round(mean(road_length_km, na.rm = TRUE), 3),
    mean_pop_density  = round(mean(pop_density, na.rm = TRUE), 1),
    mean_pct_minor    = round(mean(pct_minor_road, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()

cat("\nProportion of treated OAs that are zero-pre-injury, by country:\n")
match_with_flags %>%
  filter(treated_OA == 1) %>%
  mutate(country = case_when(
    substr(LAD24CD, 1, 1) == "E" ~ "England",
    substr(LAD24CD, 1, 1) == "S" ~ "Scotland"
  )) %>%
  group_by(country) %>%
  summarise(
    n_total       = n(),
    n_zero_pre    = sum(zero_pre_injury_OA == 1),
    pct_zero_pre  = round(mean(zero_pre_injury_OA == 1) * 100, 1),
    .groups = "drop"
  ) %>% print()

# ----------------------------------------------------------------------------
# Q4. TWO-STAGE MATCHING REMINDER + INTERACTION TEST CODE
# ----------------------------------------------------------------------------

cat("\n--- Q4: Two-stage matching reminder and interaction test code ---\n\n")

cat("MATCHING DESIGN REMINDER:\n")
cat("  Yes — the pipeline uses two-stage matching:\n")
cat("  Stage 1 matches on road network + sociodemographic structure.\n")
cat("  Stage 2 matches on pre-treatment injury trajectory WITHIN the\n")
cat("  Stage 1 pool. Crucially, Analysis B sets all Stage 2 variables\n")
cat("  to 0 for zero-pre-injury treated OAs, so they are matched to\n")
cat("  OTHER zero-pre-injury controls (near-zero distance).\n")
cat("  This means zero-pre-injury treated OAs are compared to structurally\n")
cat("  similar zero-pre-injury controls — not to injury-exposed controls.\n")
cat("  So structural differences between zero and non-zero subgroups do NOT\n")
cat("  invalidate the comparison within Analysis B, as long as the Stage 1\n")
cat("  structural match is adequate within the zero-pre-injury subgroup.\n\n")

cat("INTERACTION TEST CODE (for final DiD model):\n")
cat("Run this in the DiD script once the panel is assembled:\n\n")

cat("# Test 1: Interaction — does CAZ effect differ by zero vs non-zero pre-injury?\n")
cat("model_interact <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) *\n")
cat("                      zero_pre_injury_OA |\n")
cat("    identifier + quarter_year,\n")
cat("  data   = analysis_sample,\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n")
cat("summary(model_interact)\n\n")

cat("# Test 2: Separate models\n")
cat("model_nonzero <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) |\n")
cat("    identifier + quarter_year,\n")
cat("  data    = filter(analysis_sample, zero_pre_injury_OA == 0),\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n")
cat("model_zero <- feols(\n")
cat("  total_inj_adj_all ~ i(post, treated_group, ref = FALSE) |\n")
cat("    identifier + quarter_year,\n")
cat("  data    = filter(analysis_sample, zero_pre_injury_OA == 1),\n")
cat("  weights = ~ weights,\n")
cat("  cluster = ~ OA\n")
cat(")\n\n")

cat("# CLUSTER NOTE (colleague's question on Scotland):\n")
cat("# Clustering on OA (cluster = ~OA) is already the right approach.\n")
cat("# It accounts for within-OA correlation across roads and time,\n")
cat("# which absorbs the Scotland weight-concentration problem at the\n")
cat("# variance estimation stage. The few high-weight Scottish control OAs\n")
cat("# are already clustered together — their contribution to the variance\n")
cat("# is correctly inflated by the cluster SE estimator, so the t-statistics\n")
cat("# for Scotland-driven effects will be appropriately conservative.\n")
cat("# No additional correction is needed beyond cluster = ~OA.\n\n")

# ----------------------------------------------------------------------------
# SUMMARY TABLE: everything your colleague asked for in one place
# ----------------------------------------------------------------------------

cat("====================================================\n")
cat("SUMMARY FOR COLLEAGUE\n")
cat("====================================================\n\n")

cat("1. % of sample roads from Scottish zero-pre-injury OAs (the problematic\n")
cat("   matched set): see oa_road_counts_summary table above.\n\n")

cat("2. % of total injuries from this group: see pct_total_inj column.\n")
cat("   If < 2%, these matches are a matching artefact and do not drive ATT.\n\n")

cat("3. Do zero-pre-injury OAs remain structurally distinct after removing\n")
cat("   zero-alltime-injury roads? See smd_zero_vs_nonzero table above.\n")
cat("   Key variables to check: road_length_km, pop_density, pct_minor_road.\n\n")

cat("4. Is this mainly Scotland? See Q3 table above (pct_zero_pre by country).\n\n")

cat("5. Two-stage matching design already handles zero-pre-injury OAs by\n")
cat("   comparing them to structurally similar zero-pre-injury controls.\n\n")

cat("6. Clustering on OA (cluster = ~OA) in the final model absorbs the\n")
cat("   Scotland weight-concentration problem at the SE level.\n\n")

cat("7. Interaction and separate-model tests coded above — run once panel\n")
cat("   is assembled. Only needed if zero-pre-injury OAs contribute\n")
cat("   non-trivial post-treatment outcome variation.\n\n")