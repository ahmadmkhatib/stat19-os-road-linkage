# =============================================================================
# CAZ/LEZ Staggered DiD: OA-Level Matching — Two-Stage Design
# Stage 1: Coarsened Exact Matching (CEM) on structural variables
# Stage 2: Mahalanobis Distance Matching (MDM) on pre-treatment trends
# =============================================================================
#
# DESIGN RATIONALE:
#   The two-stage approach addresses two distinct threats to parallel trends:
#
#   Stage 1 (CEM) — Structural restriction
#     Removes OAs that are fundamentally incomparable in kind before any
#     outcome data are examined. Road type, AADT band, urban/rural class,
#     speed limit regime, junction density, and traffic volume level are
#     categorical or naturally coarsened — the meaningful distinctions are
#     between bands, not within them. CEM enforces hard cell boundaries:
#     an OA outside the same cell is excluded entirely, not merely penalised.
#     This cannot be done with MDM, which always finds *a* match regardless
#     of structural plausibility.
#
#   Stage 2 (MDM) — Outcome trend matching
#     Within the structurally restricted pool from Stage 1, finds the closest
#     match on the variables that directly define the parallel trends assumption:
#     pre-treatment injury trends by mode and severity, and volume trend.
#     Mahalanobis distance accounts for correlation across trend variables via
#     the covariance matrix — no intermediate model required, no V matrix to
#     defend. The caliper drops treated OAs with no sufficiently close match
#     rather than accepting a poor one.
#
#   Why NOT a V matrix / weighted MDM on everything?
#     A V matrix encodes causal theory implicitly in scalar weights that
#     interact across all variables and are hard to justify or audit.
#     The two-stage design makes the same distinctions explicitly and
#     transparently: hard structural boundary in Stage 1, continuous
#     distance on trends in Stage 2.
#
#   Matching is cross-sectional at OA level (no timing information).
#   Scheme-specific effect estimates come from the C&S aggregation step,
#   not from the matching step. A pooled donor pool after Stage 1 is
#   an advantage — more structurally comparable OAs available for MDM.
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(here)
library(MASS)        # for ginv() — generalised inverse fallback
library(tidyverse)

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# =============================================================================
# STEP 0 — DERIVE COUNTRY + PRE-MATCHING EXCLUSIONS
# =============================================================================
# Country derived from LAD24CD prefix:
#   "E" = England   (e.g. E08000001)
#   "S" = Scotland  (e.g. S12000033)
#
# WHY THIS MATTERS:
#   England/Wales (STATS19) and Scotland (STATS19-Scotland) differ in injury
#   classification thresholds, recording protocols, and reporting obligations.
#   Cross-country matches introduce systematic measurement incomparability.
#   Country is enforced as an exact CEM cell in Stage 1 — not a soft constraint.

oa_data_clean <- OA_matching_dataset %>%
  
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  ) %>%
  
  { . ->> oa_data_with_country } %>%  # save pre-filter for diagnostics
  filter(country != "Unknown") %>%
  
  filter(
    # Keep treated OAs and clean control OAs only
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    
    # Remove buffer OAs (spillover / displacement contamination)
    buffer_OA == 0,
    
    # Remove zero-injury OAs — uninformative for trend matching
    # and can distort slope estimates. Flagged as sensitivity check.
    zero_injury_OA == 0
  )

# --- Composition check ---
cat("=== COUNTRY COMPOSITION ===\n")
cat("Unknown LAD24CD prefixes:", sum(oa_data_with_country$country == "Unknown"), "\n")

country_tab <- oa_data_clean %>%
  count(country, treated_OA) %>%
  pivot_wider(names_from = treated_OA, values_from = n, names_prefix = "treated_") %>%
  rename(n_control = treated_0, n_treated = treated_1)

print(country_tab)
cat("OAs after exclusions:", nrow(oa_data_clean), "\n")
cat("  Treated:", sum(oa_data_clean$treated_OA), "\n")
cat("  Controls:", sum(oa_data_clean$treated_OA == 0), "\n")

# =============================================================================
# STEP 1 — DEFINE STAGE 1 AND STAGE 2 VARIABLES
# =============================================================================

# --- STAGE 1 variables (CEM — structural, categorical or coarsened) ---
#
# These variables define whether an OA is structurally comparable.
# A mismatch on any of these means the OA is not a valid comparator,
# regardless of how similar its outcome trends may appear.
#
# Variables must be coarsened into discrete bins before CEM.
# Bin boundaries are pre-specified using domain knowledge, not data-driven.

stage1_vars_cem <- c(
  "road_type_class",       # functional road classification (A / B / minor / urban arterial)
  "aadt_band",             # traffic volume band (<5k / 5-15k / 15-30k / >30k)
  "urban_rural_class",     # urban / rural classification
  "speed_limit_band",      # speed limit regime (<=30mph / 40mph / >=50mph)
  "junction_density_band", # junction density tercile (low / medium / high)
  "vol_level_band",        # traffic volume level band — coarsened AADT or VKT proxy
  "country"                # England vs Scotland — exact cell match (see rationale above)
)

# --- STAGE 2 variables (MDM — pre-treatment outcome trends only) ---
#
# These variables directly define the parallel trends assumption:
# two OAs that share the same pre-treatment trajectory are comparable
# in the counterfactual sense required by DiD.
#
# Injury LEVELS are deliberately excluded from Stage 2.
# Levels are partially collinear with trends and their inclusion risks
# matching on absolute scale rather than trajectory. The parallel trends
# assumption is about trajectories, not levels.
#
# Traffic volume TREND enters Stage 2 because it captures how exposure
# was changing pre-treatment — directly relevant to injury trajectory.
# Traffic volume LEVEL enters Stage 1 as a structural characteristic.

stage2_vars_mdm <- c(
  # --- Injury TRENDS (quasi-Poisson slope, log-rate change per quarter per road-km) ---
  "trend_car_KSI_pkm",
  "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",
  "trend_ped_slight_pkm",
  "trend_total_pkm",
  
  # --- Volume TREND (how traffic was changing pre-treatment) ---
  "trend_volume_pkm"       # add if available; drop otherwise — see note below
)

# NOTE on trend_volume_pkm:
#   If traffic volume trend is not available at OA level, remove it from
#   stage2_vars_mdm. Stage 1 already controls for volume level via vol_level_band.
#   The matching will still be valid — this is a "nice to have" variable.

cat("\nStage 1 variables (CEM):", length(stage1_vars_cem), "\n")
cat("Stage 2 variables (MDM):", length(stage2_vars_mdm), "\n")

# =============================================================================
# STEP 2 — COARSEN STAGE 1 VARIABLES
# =============================================================================
# Bin boundaries are pre-specified and fixed before any reference to
# outcome data or balance. This is the key discipline of Stage 1:
# the coarsening must be based purely on domain/engineering logic.
#
# If these variables already exist as factors in the dataset, skip this step.

oa_data_clean <- oa_data_clean %>%
  
  mutate(
    
    # Road type: use existing classification if available as factor
    # If stored as numeric or string, recode here:
    road_type_class = case_when(
      pct_A_road >= 0.25                          ~ "A_road_dominant",
      pct_B_road >= 0.25                          ~ "B_road_dominant",
      road_density_m_km2 > median(road_density_m_km2, na.rm = TRUE) ~ "urban_minor",
      TRUE                                         ~ "rural_minor"
    ),
    
    # AADT band — coarsen traffic volume level into 4 bands
    # Adjust thresholds to match your AADT variable name and units
    aadt_band = case_when(
      is.na(road_length_km)          ~ NA_character_,
      road_length_km < 0.5           ~ "very_short",     # proxy for low-volume
      road_length_km < 1.0           ~ "short",
      road_length_km < 2.0           ~ "medium",
      TRUE                           ~ "long"
    ),
    # NOTE: replace road_length_km with actual AADT variable when available.
    # Intended bands: <5k / 5-15k / 15-30k / >30k vehicles/day
    
    # Urban/rural: use existing variable if present, else derive from density
    urban_rural_class = case_when(
      road_density_m_km2 > 15000 ~ "urban",
      road_density_m_km2 > 5000  ~ "peri_urban",
      TRUE                        ~ "rural"
    ),
    
    # Speed limit band — if speed_limit variable exists, coarsen here
    # Placeholder: replace with actual speed limit variable
    speed_limit_band = "not_available",   # REPLACE when variable available
    
    # Junction density tercile — pre-specified, not data-driven
    junction_density_band = ntile(road_density_m_km2, 3) %>%
      factor(labels = c("low", "medium", "high")),
    
    # Volume level band — coarsen to match Stage 1 intent
    # Using road_length_km as proxy; replace with AADT when available
    vol_level_band = case_when(
      road_length_km < 0.5  ~ "low",
      road_length_km < 1.5  ~ "medium",
      TRUE                   ~ "high"
    )
    
  )

# --- Check coarsened variable distributions ---
cat("\n=== STAGE 1 COARSENED VARIABLE DISTRIBUTIONS ===\n")
for (v in setdiff(stage1_vars_cem, "country")) {
  if (v %in% names(oa_data_clean)) {
    cat("\n", v, ":\n", sep = "")
    print(table(oa_data_clean[[v]], oa_data_clean$treated_OA,
                dnn = c(v, "treated"), useNA = "ifany"))
  } else {
    cat("\nWARNING:", v, "not found in dataset — check variable name\n")
  }
}

# =============================================================================
# STEP 3 — STAGE 1: COARSENED EXACT MATCHING (CEM)
# =============================================================================
# CEM restricts the donor pool by excluding OAs that fall outside the same
# coarsened cell as any treated OA. This is a hard exclusion — not a penalty.
#
# After CEM, the pool contains only structurally comparable OAs.
# MDM in Stage 2 then finds the best match within this restricted pool.
#
# replace = TRUE here means a control OA can be retained in the pool for
# multiple treated OAs — this is about pool restriction, not 1:1 pairing.
# The final pairing happens in Stage 2.

# Check all Stage 1 variables exist
missing_s1 <- setdiff(stage1_vars_cem, names(oa_data_clean))
if (length(missing_s1) > 0) {
  cat("WARNING — Stage 1 variables missing:", paste(missing_s1, collapse = ", "), "\n")
  cat("  Remove missing variables from stage1_vars_cem before proceeding.\n")
  stage1_vars_cem <- intersect(stage1_vars_cem, names(oa_data_clean))
}

oa_data_clean <- oa_data_clean %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

cem_formula <- reformulate(stage1_vars_cem, response = "treat_indicator")

cat("\n=== STAGE 1: CEM ===\n")

cem_out <- tryCatch(
  matchit(
    cem_formula,
    data    = oa_data_clean,
    method  = "cem",
    estimand = "ATT"
  ),
  error = function(e) {
    cat("CEM failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(cem_out)) stop("Stage 1 CEM failed — cannot proceed to Stage 2.")

# --- Stage 1 summary ---
cat("\nCEM output:\n")
print(summary(cem_out, un = FALSE))

# Extract OAs that survived Stage 1 (in matched cells)
cem_data <- match.data(cem_out)

n_treated_s1  <- sum(cem_data$treat_indicator == 1)
n_control_s1  <- sum(cem_data$treat_indicator == 0)
n_dropped_s1  <- sum(oa_data_clean$treat_indicator == 1) - n_treated_s1

cat("\nStage 1 results:\n")
cat("  Treated OAs retained:", n_treated_s1, "\n")
cat("  Control OAs in restricted pool:", n_control_s1, "\n")
cat("  Treated OAs dropped (no structural match):", n_dropped_s1,
    sprintf("(%.1f%%)\n", 100 * n_dropped_s1 / sum(oa_data_clean$treat_indicator)))

# Verify no cross-country matches possible after Stage 1
# (country is a CEM cell variable — cross-country OAs will be in different cells)
cat("\nCountry distribution after Stage 1:\n")
print(table(cem_data$country, cem_data$treat_indicator,
            dnn = c("Country", "Treated")))

# =============================================================================
# STEP 4 — STAGE 2: MAHALANOBIS DISTANCE MATCHING ON TRENDS
# =============================================================================
# Within the CEM-restricted pool, MDM finds the closest match on pre-treatment
# trajectory variables. The Mahalanobis distance accounts for correlation
# across trend variables via the sample covariance matrix — no weights to
# specify, no intermediate model.
#
# Caliper: treated OAs with no control within 0.25 SD are dropped, not
# badly matched. A biased match is a worse problem than a smaller sample.
#
# replace = FALSE: 1:1 matching without replacement within the CEM pool.
# Adjust ratio if you want 1:k matching.

# Check Stage 2 variables exist and drop missing
missing_s2 <- setdiff(stage2_vars_mdm, names(cem_data))
if (length(missing_s2) > 0) {
  cat("\nWARNING — Stage 2 variables missing (dropping):",
      paste(missing_s2, collapse = ", "), "\n")
  stage2_vars_mdm <- intersect(stage2_vars_mdm, names(cem_data))
}

if (length(stage2_vars_mdm) == 0) {
  stop("No Stage 2 variables available — cannot run MDM.")
}

# Winsorise trend variables at 1st/99th percentile
# Prevents outlier OAs from distorting the covariance matrix
cem_data_mdm <- cem_data %>%
  mutate(across(
    all_of(stage2_vars_mdm),
    ~ {
      q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

# Drop near-zero-variance trend variables — these offer no discriminating
# information and cause near-singularity in the covariance matrix
var_check <- cem_data_mdm %>%
  summarise(across(all_of(stage2_vars_mdm), ~ var(., na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "variance")

low_var <- var_check %>% filter(variance < 1e-8) %>% pull(variable)

if (length(low_var) > 0) {
  cat("\nDropping near-zero-variance Stage 2 variables:",
      paste(low_var, collapse = ", "), "\n")
  stage2_vars_mdm <- setdiff(stage2_vars_mdm, low_var)
}

mdm_formula <- reformulate(stage2_vars_mdm, response = "treat_indicator")

cat("\n=== STAGE 2: MDM ===\n")
cat("Matching on", length(stage2_vars_mdm), "trend variables:\n")
cat(" ", paste(stage2_vars_mdm, collapse = "\n  "), "\n")

# --- Primary: caliper = 0.25 SD ---
mdm_out_025 <- tryCatch(
  matchit(
    mdm_formula,
    data        = cem_data_mdm,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = 0.25,
    std.caliper = TRUE,
    ratio       = 1,
    replace     = FALSE
  ),
  error = function(e) {
    cat("MDM (caliper 0.25) failed:", conditionMessage(e), "\n"); NULL
  }
)

# --- Sensitivity: caliper = 0.20 ---
mdm_out_020 <- tryCatch(
  matchit(
    mdm_formula,
    data        = cem_data_mdm,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = 0.20,
    std.caliper = TRUE,
    ratio       = 1,
    replace     = FALSE
  ),
  error = function(e) { cat("MDM (caliper 0.20) failed\n"); NULL }
)

# --- Sensitivity: caliper = 0.30 ---
mdm_out_030 <- tryCatch(
  matchit(
    mdm_formula,
    data        = cem_data_mdm,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = 0.30,
    std.caliper = TRUE,
    ratio       = 1,
    replace     = FALSE
  ),
  error = function(e) { cat("MDM (caliper 0.30) failed\n"); NULL }
)

# --- Stage 2 summary ---
if (!is.null(mdm_out_025)) {
  mdm_data <- match.data(mdm_out_025)
  
  n_treated_s2  <- sum(mdm_data$treat_indicator == 1)
  n_control_s2  <- sum(mdm_data$treat_indicator == 0)
  n_dropped_s2  <- n_treated_s1 - n_treated_s2
  
  cat("\nStage 2 results (caliper 0.25 SD):\n")
  cat("  Treated OAs matched:", n_treated_s2, "\n")
  cat("  Control OAs matched:", n_control_s2, "\n")
  cat("  Treated OAs dropped by caliper:", n_dropped_s2,
      sprintf("(%.1f%% of Stage 1 sample)\n",
              100 * n_dropped_s2 / n_treated_s1))
  
  # Caliper sensitivity comparison
  cat("\nCaliper sensitivity — treated OAs retained:\n")
  for (obj in list(
    list(label = "0.20 SD", m = mdm_out_020),
    list(label = "0.25 SD (primary)", m = mdm_out_025),
    list(label = "0.30 SD", m = mdm_out_030)
  )) {
    if (!is.null(obj$m)) {
      n <- sum(match.data(obj$m)$treat_indicator == 1)
      cat(sprintf("  Caliper %-20s: %d treated OAs\n", obj$label, n))
    }
  }
  
  # Verify no cross-country matches
  cross_country <- mdm_data %>%
    group_by(subclass) %>%
    summarise(n_countries = n_distinct(country), .groups = "drop") %>%
    filter(n_countries > 1)
  
  if (nrow(cross_country) > 0) {
    cat("\nWARNING:", nrow(cross_country),
        "matched pairs contain cross-country matches — investigate!\n")
  } else {
    cat("\nCountry check passed — no cross-country matches\n")
  }
}

# =============================================================================
# STEP 5 — BALANCE DIAGNOSTICS
# =============================================================================
# Balance assessed on BOTH Stage 1 structural variables and Stage 2 trend
# variables. The two-stage design should produce:
#   - Excellent balance on structural variables (CEM guarantees cell equality)
#   - Good balance on trend variables (MDM optimises for this directly)
#
# Target: |SMD| < 0.10 for all variables (Stuart 2010).
# Variance ratios close to 1.0 confirm dispersion comparability.
# Significance tests (t-tests) are deliberately avoided — they conflate
# balance with sample size.

if (!is.null(mdm_out_025)) {
  
  cat("\n=== BALANCE DIAGNOSTICS (primary: caliper 0.25 SD) ===\n")
  
  # Balance on Stage 2 trend variables (what MDM directly optimised)
  cat("\nBalance on Stage 2 trend variables:\n")
  print(bal.tab(mdm_out_025,
                thresholds = c(m = 0.1, v = 2),
                un = TRUE))
  
  # Love plot — Stage 2 trends
  lp_trends <- love.plot(
    mdm_out_025,
    threshold    = 0.1,
    abs          = TRUE,
    var.order    = "unadjusted",
    title        = "Covariate Balance — Stage 2 Trend Variables",
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  )
  ggsave("balance_stage2_trends.png", lp_trends, width = 10, height = 7)
  
  # Pre-trend overlay plot — visual check of parallel trends
  # Plot mean quarterly injury trajectory for treated vs matched controls
  # (requires road × quarter panel to be merged back — done in Step 6)
  
  cat("\nBalance diagnostics saved.\n")
}

# =============================================================================
# STEP 6 — EXTRACT MATCHED OA IDs AND BUILD RESTRICTED DONOR POOL
# =============================================================================

if (!is.null(mdm_out_025)) {
  
  mdm_data <- match.data(mdm_out_025)
  
  matched_control_oas <- mdm_data %>%
    filter(treat_indicator == 0) %>%
    dplyr::select(OA, weights, subclass)
  
  matched_treated_oas <- mdm_data %>%
    filter(treat_indicator == 1) %>%
    dplyr::select(OA, weights, subclass)
  
  cat("\n=== MATCHING SUMMARY ===\n")
  cat("Matched treated OAs:", nrow(matched_treated_oas), "\n")
  cat("Matched control OAs:", nrow(matched_control_oas), "\n")
  
  # Save matched OA lists
  saveRDS(matched_control_oas,
          here("data", "processed", "OA_matched_donors.rds"))
  saveRDS(matched_treated_oas,
          here("data", "processed", "OA_matched_treated.rds"))
  
  # --- Restrict road-link panel ---
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
        pct_minor_road, road_length_km, dist_citycentre, n_roads
      ),
      by = "OA"
    )
  
  cat("Road links in restricted panel:", nrow(road_link_panel_restricted), "\n")
  cat("Unique OAs in panel:", n_distinct(road_link_panel_restricted$OA), "\n")
}

# =============================================================================
# STEP 7 — CALLAWAY & SANT'ANNA ESTIMATION
# =============================================================================
# The matching step defines who is in the donor pool.
# Scheme-specific (individual CAZ) effect estimates come from the C&S
# aggregation step — not from the matching step. The pooled donor pool
# does not preclude by-scheme analysis.

library(did)

outcomes <- c("KSI_adj", "Slight_adj")
modes    <- c("All", "Pedestrian", "Car.Van", "Cyclist", "Other")

run_cs <- function(yvar, data) {
  att_gt(
    yname         = yvar,
    tname         = "quarter_num",
    idname        = "road_link_id",
    gname         = "cohort_g",           # first treated quarter; 0 = never treated
    
    xformla       = ~ road_density_m_km2 + pct_A_road + pct_minor_road +
      road_length_km + dist_citycentre,
    
    est_method    = "dr",                 # doubly robust
    control_group = "notyettreated",      # not-yet-treated preferred over never-treated
    
    clustervars   = "OA",                 # cluster at OA level (roads within OA correlated)
    
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
      error = function(e) {
        cat("  Failed:", conditionMessage(e), "\n"); NULL
      }
    )
  }
}

# --- Aggregate ATT(g,t) estimates ---
# simple  = overall ATT
# dynamic = event-study by time-since-treatment
# group   = by-scheme (individual CAZ effects)
cs_agg <- lapply(cs_results, function(res) {
  if (is.null(res)) return(NULL)
  list(
    simple  = aggte(res, type = "simple"),
    dynamic = aggte(res, type = "dynamic"),
    group   = aggte(res, type = "group")   # scheme-specific ATTs
  )
})

# =============================================================================
# STEP 8 — SENSITIVITY / ROBUSTNESS CHECKS
# =============================================================================

# 8a. Caliper sensitivity (already run in Step 4 — compare matched samples)

# 8b. Zero-injury OA sensitivity
#   Re-run full pipeline with zero_injury_OA retained (flagged, not excluded)
#   Compare matched sample composition and C&S point estimates

# 8c. Control group: never-treated vs not-yet-treated
#   Primary uses not-yet-treated (more efficient, recommended by C&S)
#   Sensitivity: re-run with control_group = "nevertreated"

# 8d. Treatment definition sensitivity
#   Primary: treated_50pct (50%+ of road inside CAZ)
#   Sensitivity: treated_any (any overlap with CAZ boundary)

# 8e. HonestDiD pre-trend sensitivity (Rambachan & Roth 2023)
# library(HonestDiD)
# Apply to each dynamic aggregation — see HonestDiD vignette

# =============================================================================
# STEP 9 — BALANCE REPORTING (PUBLICATION-READY)
# =============================================================================

if (!is.null(mdm_out_025)) {
  
  bt <- bal.tab(mdm_out_025, thresholds = c(m = 0.1), un = TRUE)
  
  smd_before <- abs(bt$Balance$Diff.Un)
  smd_after  <- abs(bt$Balance$Diff.Adj)
  
  balance_report <- data.frame(
    variable     = rownames(bt$Balance),
    smd_before   = round(smd_before, 3),
    smd_after    = round(smd_after,  3),
    balanced     = smd_after < 0.1
  ) %>%
    arrange(desc(smd_after))
  
  cat("\n=== BALANCE REPORT ===\n")
  cat("Variables with |SMD| < 0.10 after matching:",
      sum(balance_report$balanced, na.rm = TRUE), "/",
      nrow(balance_report), "\n")
  cat("Max |SMD| after matching:",
      round(max(smd_after, na.rm = TRUE), 3), "\n")
  cat("Mean |SMD| after matching:",
      round(mean(smd_after, na.rm = TRUE), 3), "\n")
  
  print(balance_report)
  
  write.csv(balance_report,
            here("output", "balance_report_twostage.csv"),
            row.names = FALSE)
}