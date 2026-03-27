# =============================================================================
# CAZ/LEZ Staggered DiD: OA-Level Matching — Two-Stage MDM Design
# Stage 1: MDM (wide caliper) on structural + sociodemographic variables
# Stage 2: MDM (tight caliper) on pre-treatment injury trends + levels
# =============================================================================
#
# DESIGN RATIONALE:
#
#   Stage 1 — Structural restriction (MDM, wide caliper)
#     Removes OAs that are too structurally dissimilar to be plausible
#     comparators, before any injury outcome data are examined.
#     Uses road network, urban form, and sociodemographic variables.
#     Wide caliper (1.0 SD) — the goal is to cut the incomparable tail,
#     not to find tight matches. A control OA far outside the structural
#     distribution of treated OAs is not a valid comparator regardless
#     of its injury trends.
#
#     WHY MDM NOT CEM AT STAGE 1:
#     CEM with 20 variables produces 3^20 possible cells — almost every
#     treated OA ends up alone in its cell and is dropped. MDM with a
#     wide caliper achieves the same structural filtering without the
#     cell-sparsity problem. The caliper still enforces a hard exclusion:
#     any control OA beyond 1.0 SD is dropped entirely, not penalised.
#
#   Stage 2 — Outcome matching (MDM, tight caliper)
#     Within the Stage 1 restricted pool, finds the closest match on
#     pre-treatment injury dynamics — trends and levels by mode and severity.
#     These variables are the direct empirical content of the parallel
#     trends assumption. Nothing from Stage 1 enters Stage 2, so census
#     and structural variables cannot compensate for a trajectory mismatch.
#     Tight caliper (0.25 SD primary) — treated OAs with no sufficiently
#     close trajectory match are dropped rather than badly matched.
#
#     RATIO FLEXIBILITY:
#     Stage 2 is run at multiple ratios (1:1, 1:2, 1:4) and the best
#     is chosen based on balance diagnostics and pre-trend overlay plots.
#     More controls = larger donor pool and more power in C&S estimation,
#     but at the cost of slightly weaker average match quality on trends.
#
#   Country (England vs Scotland) is enforced as a hard exact constraint
#   at BOTH stages via exact = ~ country in MatchIt. This is non-negotiable:
#   STATS19 and STATS19-Scotland differ in recording thresholds, severity
#   classification, and reporting protocols. Cross-country matches introduce
#   systematic measurement incomparability that no covariate adjustment fixes.
#
#   Matching is cross-sectional at OA level — no timing information.
#   Scheme-specific ATTs come from C&S aggregation, not from matching.
#   A pooled donor pool after Stage 1 is an advantage: more structurally
#   comparable OAs available for Stage 2 MDM.
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(here)
library(MASS)
library(tidyverse)

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# =============================================================================
# STEP 0 — DERIVE COUNTRY + PRE-MATCHING EXCLUSIONS
# =============================================================================

oa_data_clean <- OA_matching_dataset %>%
  
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  ) %>%
  
  { . ->> oa_data_with_country } %>%
  filter(country != "Unknown") %>%
  
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA      == 0,
    zero_injury_OA == 0    # sensitivity check: re-run with this removed
  ) %>%
  
  mutate(treat_indicator = as.integer(treated_OA == 1))

# --- Composition check ---
cat("=== SAMPLE AFTER EXCLUSIONS ===\n")
cat("Unknown LAD24CD prefixes:", sum(oa_data_with_country$country == "Unknown"), "\n")
cat("Total OAs:", nrow(oa_data_clean), "\n")
cat("  Treated:", sum(oa_data_clean$treat_indicator), "\n")
cat("  Controls:", sum(oa_data_clean$treat_indicator == 0), "\n\n")

country_tab <- oa_data_clean %>%
  count(country, treat_indicator) %>%
  pivot_wider(names_from = treat_indicator, values_from = n,
              names_prefix = "treated_") %>%
  rename(n_control = treated_0, n_treated = treated_1)
print(country_tab)

# =============================================================================
# STEP 1 — DEFINE STAGE 1 AND STAGE 2 VARIABLES
# =============================================================================
# Variable names mapped directly to OA_matching_census columns.
# Stage 1: structural + sociodemographic — defines comparability of environment
# Stage 2: injury trends + levels — defines parallel trends directly
#
# SEPARATION PRINCIPLE:
#   No variable appears in both stages. This ensures that in Stage 2,
#   the Mahalanobis distance is determined entirely by injury dynamics.
#   A good census match cannot compensate for a poor trajectory match.
# =============================================================================

# --- STAGE 1: Structural + sociodemographic (MDM, wide caliper) ---
# Road network
stage1_road <- c(
  "road_density_m_km2",  # road density — urban intensity proxy
  "road_length_km",      # total road length — exposure scale
  "pct_A_road",          # % A-road segments
  "pct_B_road",          # % B-road segments
  "pct_minor_road"       # % minor road segments
  # speed limit distribution: add when variable available
  # n_roads available but collinear with road_length_km — excluded
)

# Urban form
stage1_urban <- c(
  "dist_citycentre",     # distance to city centre — urban gradient
  "pop_density"          # population density
  # retail density: add when variable available
)

# Sociodemographic — these are confounders of both treatment and injury trends
# They belong in Stage 1 so they cannot compensate for trend mismatches in Stage 2
stage1_socdem <- c(
  "IMD",                 # Index of Multiple Deprivation
  "cars_none_pct",       # % households with no car — active travel proxy
  "Drive_Car_pct",       # % driving to work
  "Walk_pct",            # % walking to work
  "Bicycle_pct",         # % cycling to work
  "X65plus_pct",         # % aged 65+ — high-risk pedestrian group
  "X5to19_pct",          # % aged 5-19 — high-risk pedestrian group
  "X20to24_pct"          # % aged 20-24 — high-risk road user group
  # ethnicity pcts available — add if theory suggests they predict injury trends
)

stage1_vars <- c(stage1_road, stage1_urban, stage1_socdem)

# --- STAGE 2: Pre-treatment injury dynamics (MDM, tight caliper) ---
# Trends: quasi-Poisson GLM slope (log-rate change per quarter per road-km)
stage2_trends <- c(
  "trend_car_KSI_pkm",
  "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",
  "trend_ped_slight_pkm",
  "trend_total_pkm"
)

# Levels: mean quarterly casualties per road-km over pre-treatment period
# Included alongside trends to anchor on absolute scale as well as trajectory.
# Two OAs with the same trend but very different absolute levels may still
# diverge post-treatment in ways that matter for the counterfactual.
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

cat("Stage 1 variables:", length(stage1_vars), "\n")
cat("  Road network:", length(stage1_road), "\n")
cat("  Urban form:", length(stage1_urban), "\n")
cat("  Sociodemographic:", length(stage1_socdem), "\n")
cat("Stage 2 variables:", length(stage2_vars), "\n")
cat("  Injury trends:", length(stage2_trends), "\n")
cat("  Injury levels:", length(stage2_levels), "\n")

# =============================================================================
# STEP 2 — PRE-PROCESSING
# =============================================================================
# Winsorise all matching variables at 1st/99th percentile.
# Prevents outlier OAs from distorting the covariance matrix at both stages.
# Applied once to the full dataset before either matching stage.

all_match_vars <- c(stage1_vars, stage2_vars)

oa_data_clean <- oa_data_clean %>%
  mutate(across(
    all_of(intersect(all_match_vars, names(.))),
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
  stage1_vars <- intersect(stage1_vars, names(oa_data_clean))
}
if (length(missing_s2) > 0) {
  cat("WARNING — Stage 2 variables missing:", paste(missing_s2, collapse = ", "), "\n")
  stage2_vars <- intersect(stage2_vars, names(oa_data_clean))
}

# Drop near-zero-variance variables (cause near-singularity in covariance matrix)
drop_low_var <- function(data, vars, threshold = 1e-8) {
  vcheck <- data %>%
    summarise(across(all_of(vars), ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "variance")
  low <- vcheck %>% filter(variance < threshold) %>% pull(variable)
  if (length(low) > 0)
    cat("Dropping near-zero-variance variables:", paste(low, collapse = ", "), "\n")
  setdiff(vars, low)
}

stage1_vars <- drop_low_var(oa_data_clean, stage1_vars)
stage2_vars <- drop_low_var(oa_data_clean, stage2_vars)

# =============================================================================
# STEP 3 — STAGE 1: MDM WITH WIDE CALIPER ON STRUCTURAL VARIABLES
# =============================================================================
# Purpose: remove OAs too structurally dissimilar to be valid comparators.
# Wide caliper (1.0 SD) — cuts the incomparable tail, not tight matching.
# Country enforced as exact constraint — hard non-negotiable.
#
# Output: oa_stage1_pool — the restricted donor pool for Stage 2.

s1_formula <- reformulate(stage1_vars, response = "treat_indicator")

cat("\n=== STAGE 1: MDM (wide caliper = 1.0 SD) ===\n")

mdm_s1 <- tryCatch(
  matchit(
    s1_formula,
    data        = oa_data_clean,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = 1.0,          # wide — remove incomparable tail only
    std.caliper = TRUE,
    ratio       = 10,           # generous ratio — retain large pool for Stage 2
    replace     = TRUE,         # with replacement — maximises pool size
    exact       = ~ country     # hard constraint — no cross-country matches
  ),
  error = function(e) {
    cat("Stage 1 MDM failed:", conditionMessage(e), "\n"); NULL
  }
)

if (is.null(mdm_s1)) stop("Stage 1 failed — cannot proceed.")

# Extract Stage 1 pool (all OAs in matched set)
s1_data    <- match.data(mdm_s1)
s1_treated <- sum(s1_data$treat_indicator == 1)
s1_control <- sum(s1_data$treat_indicator == 0)
s1_dropped <- sum(oa_data_clean$treat_indicator == 1) - s1_treated

cat("\nStage 1 results:\n")
cat("  Treated OAs retained:", s1_treated, "\n")
cat("  Control OAs in pool:", s1_control, "\n")
cat("  Treated OAs dropped (caliper):", s1_dropped,
    sprintf("(%.1f%%)\n", 100 * s1_dropped / sum(oa_data_clean$treat_indicator)))

cat("\nCountry distribution after Stage 1:\n")
print(table(s1_data$country, s1_data$treat_indicator,
            dnn = c("Country", "Treated")))

# Stage 1 balance — structural variables should be well balanced
cat("\nStage 1 balance (structural variables):\n")
print(bal.tab(mdm_s1, thresholds = c(m = 0.1), un = TRUE))

# Sensitivity: also run Stage 1 at caliper = 0.75 and 1.25
mdm_s1_075 <- tryCatch(
  matchit(s1_formula, data = oa_data_clean, method = "nearest",
          distance = "mahalanobis", caliper = 0.75, std.caliper = TRUE,
          ratio = 10, replace = TRUE, exact = ~ country),
  error = function(e) NULL
)
mdm_s1_125 <- tryCatch(
  matchit(s1_formula, data = oa_data_clean, method = "nearest",
          distance = "mahalanobis", caliper = 1.25, std.caliper = TRUE,
          ratio = 10, replace = TRUE, exact = ~ country),
  error = function(e) NULL
)

cat("\nStage 1 caliper sensitivity — treated OAs retained:\n")
for (obj in list(
  list(label = "0.75 SD", m = mdm_s1_075),
  list(label = "1.00 SD (primary)", m = mdm_s1),
  list(label = "1.25 SD", m = mdm_s1_125)
)) {
  if (!is.null(obj$m)) {
    n <- sum(match.data(obj$m)$treat_indicator == 1)
    cat(sprintf("  Caliper %-22s: %d treated OAs\n", obj$label, n))
  }
}

# =============================================================================
# STEP 4 — STAGE 2: MDM WITH TIGHT CALIPER ON INJURY DYNAMICS
# =============================================================================
# Purpose: within the Stage 1 pool, find the closest match on pre-treatment
# injury trajectories — the direct empirical content of parallel trends.
#
# Run at multiple ratios (1:1, 1:2, 1:4) — choose based on diagnostics:
#   - Pre-trend overlay plot (primary criterion)
#   - SMD balance on Stage 2 variables
#   - Number of treated OAs retained after caliper
#
# Higher ratio = larger donor pool for C&S = more power, but weaker
# average match quality on injury trends. The diagnostics tell you
# where the quality-power tradeoff is acceptable.

# Winsorise Stage 2 variables within Stage 1 pool
s2_vars_present <- intersect(stage2_vars, names(s1_data))

s1_data_s2 <- s1_data %>%
  mutate(across(
    all_of(s2_vars_present),
    ~ {
      q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

s2_formula <- reformulate(s2_vars_present, response = "treat_indicator")

cat("\n=== STAGE 2: MDM (tight caliper, multiple ratios) ===\n")
cat("Matching on", length(s2_vars_present), "injury variables\n")

# Helper: run Stage 2 MDM for a given ratio and caliper
run_s2_mdm <- function(data, formula, ratio, caliper) {
  tryCatch(
    matchit(
      formula,
      data        = data,
      method      = "nearest",
      distance    = "mahalanobis",
      caliper     = caliper,
      std.caliper = TRUE,
      ratio       = ratio,
      replace     = FALSE,       # without replacement in Stage 2
      exact       = ~ country    # country constraint maintained
    ),
    error = function(e) {
      cat(sprintf("  MDM ratio=%d caliper=%.2f failed: %s\n",
                  ratio, caliper, conditionMessage(e)))
      NULL
    }
  )
}

# --- Primary caliper = 0.25 SD, ratios 1:1, 1:2, 1:4 ---
s2_results <- list(
  r1_c025 = run_s2_mdm(s1_data_s2, s2_formula, ratio = 1, caliper = 0.25),
  r2_c025 = run_s2_mdm(s1_data_s2, s2_formula, ratio = 2, caliper = 0.25),
  r4_c025 = run_s2_mdm(s1_data_s2, s2_formula, ratio = 4, caliper = 0.25),
  
  # --- Caliper sensitivity at ratio = 1 ---
  r1_c020 = run_s2_mdm(s1_data_s2, s2_formula, ratio = 1, caliper = 0.20),
  r1_c030 = run_s2_mdm(s1_data_s2, s2_formula, ratio = 1, caliper = 0.30)
)

# --- Ratio/caliper comparison table ---
cat("\nStage 2 results summary:\n")
cat(sprintf("  %-25s  %8s  %8s  %8s\n",
            "Specification", "Treated", "Controls", "Dropped"))
cat(strrep("-", 60), "\n")

for (nm in names(s2_results)) {
  m <- s2_results[[nm]]
  if (!is.null(m)) {
    md    <- match.data(m)
    nt    <- sum(md$treat_indicator == 1)
    nc    <- sum(md$treat_indicator == 0)
    ndrop <- s1_treated - nt
    cat(sprintf("  %-25s  %8d  %8d  %8d (%.1f%%)\n",
                nm, nt, nc, ndrop, 100 * ndrop / s1_treated))
  }
}

# =============================================================================
# STEP 5 — BALANCE DIAGNOSTICS AND RATIO SELECTION
# =============================================================================
# The right ratio is chosen after examining:
#   1. Pre-trend overlay plot — do treated and control trajectories track?
#   2. SMD on Stage 2 injury variables — are trends and levels balanced?
#   3. Exclusion rate — is the caliper dropping too many treated OAs?
#
# Rule: choose the highest ratio at which:
#   - All Stage 2 SMDs remain below 0.10
#   - Pre-trend trajectories visually overlay
#   - Exclusion rate does not exceed ~25% of Stage 1 sample

run_diagnostics <- function(m_obj, label, s1_data, stage2_vars) {
  if (is.null(m_obj)) return(invisible(NULL))
  
  cat("\n--- Diagnostics:", label, "---\n")
  md <- match.data(m_obj)
  
  # SMD balance
  bt <- bal.tab(m_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_after <- abs(bt$Balance$Diff.Adj)
  cat("  Variables with |SMD| < 0.10:",
      sum(smd_after < 0.1, na.rm = TRUE), "/", length(smd_after), "\n")
  cat("  Max |SMD|:", round(max(smd_after, na.rm = TRUE), 3), "\n")
  cat("  Mean |SMD|:", round(mean(smd_after, na.rm = TRUE), 3), "\n")
  
  # Cross-country check
  cross <- md %>%
    group_by(subclass) %>%
    summarise(n_countries = n_distinct(country), .groups = "drop") %>%
    filter(n_countries > 1)
  if (nrow(cross) > 0) {
    cat("  WARNING:", nrow(cross), "cross-country pairs — investigate!\n")
  } else {
    cat("  Country check: PASSED\n")
  }
  
  # Love plot
  lp <- love.plot(
    m_obj,
    threshold    = 0.1,
    abs          = TRUE,
    var.order    = "unadjusted",
    title        = paste("Balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  )
  fname <- paste0("balance_", gsub("[^a-zA-Z0-9]", "_", label), ".png")
  ggsave(fname, lp, width = 10, height = 8)
  cat("  Love plot saved:", fname, "\n")
  
  invisible(bt)
}

# Run diagnostics for each Stage 2 specification
for (nm in names(s2_results)) {
  run_diagnostics(s2_results[[nm]], nm, s1_data_s2, s2_vars_present)
}

# =============================================================================
# STEP 6 — SELECT PRIMARY SPECIFICATION AND EXTRACT MATCHED OAs
# =============================================================================
# Set primary_spec to whichever ratio/caliper performed best in Step 5.
# Default: ratio = 1, caliper = 0.25 (most conservative starting point).
# Change after inspecting diagnostics.

primary_spec <- "r1_c025"   # <-- UPDATE after reviewing Step 5 diagnostics

mdm_primary <- s2_results[[primary_spec]]

if (is.null(mdm_primary)) stop("Primary specification failed — choose alternative.")

cat("\n=== PRIMARY SPECIFICATION:", primary_spec, "===\n")
primary_data <- match.data(mdm_primary)

matched_control_oas <- primary_data %>%
  filter(treat_indicator == 0) %>%
  dplyr::select(OA, weights, subclass)

matched_treated_oas <- primary_data %>%
  filter(treat_indicator == 1) %>%
  dplyr::select(OA, weights, subclass)

cat("Matched treated OAs:", nrow(matched_treated_oas), "\n")
cat("Matched control OAs:", nrow(matched_control_oas), "\n")

saveRDS(matched_control_oas, here("data", "processed", "OA_matched_donors.rds"))
saveRDS(matched_treated_oas, here("data", "processed", "OA_matched_treated.rds"))

# =============================================================================
# STEP 7 — BUILD RESTRICTED ROAD-LINK PANEL
# =============================================================================

road_link_panel_restricted <- road_link_panel %>%
  filter(
    treated_link == 1 |
      (treated_link == 0 &
         OA %in% matched_control_oas$OA &
         buffer_OA == FALSE)
  ) %>%
  left_join(
    oa_data_clean %>% dplyr::select(
      OA, road_density_m_km2, pct_A_road, pct_B_road,
      pct_minor_road, road_length_km, dist_citycentre,
      IMD, cars_none_pct, pop_density
    ),
    by = "OA"
  )

cat("\nRoad links in restricted panel:", nrow(road_link_panel_restricted), "\n")
cat("Unique OAs in panel:", n_distinct(road_link_panel_restricted$OA), "\n")

# =============================================================================
# STEP 8 — CALLAWAY & SANT'ANNA ESTIMATION
# =============================================================================
# Scheme-specific ATTs come from the group aggregation here —
# not from the matching step. The pooled donor pool does not
# preclude individual CAZ effect estimates.

library(did)

outcomes <- c("KSI_adj", "Slight_adj")
modes    <- c("All", "Pedestrian", "Car.Van", "Cyclist", "Other")

run_cs <- function(yvar, data) {
  att_gt(
    yname         = yvar,
    tname         = "quarter_num",
    idname        = "road_link_id",
    gname         = "cohort_g",
    
    xformla       = ~ road_density_m_km2 + pct_A_road + pct_minor_road +
      road_length_km + dist_citycentre + IMD + cars_none_pct,
    
    est_method    = "dr",
    control_group = "notyettreated",
    clustervars   = "OA",
    panel                  = TRUE,
    allow_unbalanced_panel = TRUE,
    data                   = data
  )
}

cs_results <- list()
for (outcome in outcomes) {
  for (mode in modes) {
    col_name <- if (mode == "All") outcome else paste0(outcome, "_", mode)
    cat("Estimating:", col_name, "\n")
    cs_results[[col_name]] <- tryCatch(
      run_cs(col_name, road_link_panel_restricted),
      error = function(e) { cat("  Failed:", conditionMessage(e), "\n"); NULL }
    )
  }
}

# Aggregations
cs_agg <- lapply(cs_results, function(res) {
  if (is.null(res)) return(NULL)
  list(
    simple  = aggte(res, type = "simple"),   # overall ATT
    dynamic = aggte(res, type = "dynamic"),  # event-study
    group   = aggte(res, type = "group")     # by-scheme individual CAZ ATTs
  )
})

# =============================================================================
# STEP 9 — SENSITIVITY CHECKS
# =============================================================================

# 9a. Stage 2 ratio sensitivity
#   Already run in Step 4 — compare diagnostics across r1/r2/r4

# 9b. Stage 1 caliper sensitivity
#   Already run in Step 3 — compare sample sizes at 0.75/1.00/1.25 SD

# 9c. Stage 2 tight caliper sensitivity
#   Already run in Step 4 — compare 0.20/0.25/0.30 SD at ratio = 1

# 9d. Zero-injury OA sensitivity
#   Re-run full pipeline without zero_injury_OA == 0 filter in Step 0

# 9e. Control group sensitivity
#   Re-run C&S with control_group = "nevertreated"

# 9f. Treatment definition sensitivity
#   Primary: treated_50pct; sensitivity: treated_any

# 9g. HonestDiD pre-trend sensitivity (Rambachan & Roth 2023)
# library(HonestDiD)
# Apply to each dynamic aggregation

# =============================================================================
# STEP 10 — PUBLICATION-READY BALANCE REPORT
# =============================================================================

if (!is.null(mdm_primary)) {
  
  bt_final <- bal.tab(mdm_primary, thresholds = c(m = 0.1), un = TRUE)
  
  balance_report <- data.frame(
    variable   = rownames(bt_final$Balance),
    smd_before = round(abs(bt_final$Balance$Diff.Un),  3),
    smd_after  = round(abs(bt_final$Balance$Diff.Adj), 3)
  ) %>%
    mutate(
      balanced = smd_after < 0.1,
      stage    = case_when(
        variable %in% stage1_vars ~ "Stage 1 — structural",
        variable %in% stage2_vars ~ "Stage 2 — injury dynamics",
        TRUE                      ~ "Other"
      )
    ) %>%
    arrange(stage, desc(smd_after))
  
  cat("\n=== FINAL BALANCE REPORT (", primary_spec, ") ===\n", sep = "")
  cat("Stage 2 variables with |SMD| < 0.10:",
      sum(balance_report$balanced[balance_report$stage == "Stage 2 — injury dynamics"],
          na.rm = TRUE), "/",
      sum(balance_report$stage == "Stage 2 — injury dynamics"), "\n")
  cat("Max |SMD| (Stage 2):",
      round(max(balance_report$smd_after[
        balance_report$stage == "Stage 2 — injury dynamics"], na.rm = TRUE), 3), "\n")
  
  print(balance_report)
  
  write.csv(balance_report,
            here("output", "balance_report_twostage.csv"),
            row.names = FALSE)
}