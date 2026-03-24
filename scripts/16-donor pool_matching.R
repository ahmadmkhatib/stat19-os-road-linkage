# =============================================================================
# CAZ/LEZ Staggered DiD: OA-Level Matching — Abadie et al. (2010) 
# Data-Driven Optimal Weight Matrix V
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(WeightIt)
library(quadprog)   # ADD: needed for V matrix optimisation
library(here)
library(MASS)        # for ginv() — generalised matrix inverse
library(tidyverse)

OA_matching_dataset<- readRDS(here("data","processed","OA_matching_census.rds"))
names(OA_matching_dataset)

# =============================================================================
# CAZ/LEZ Staggered DiD: OA-Level Matching for Road-Link Donor Pool Restriction
# Callaway & Sant'Anna (2021) — Doubly Robust Estimator
# Abadie, Diamond & Hainmueller (2010) — Data-Driven Weight Matrix V
# =============================================================================
# VARIABLE NAMES: mapped to actual dataset columns
# STRATEGY:
#   1. Match at OA level to define a comparable donor pool
#   2. All road links within matched control OAs enter C&S estimation
#   3. Matching performed SEPARATELY per treatment cohort
#   4. Weight matrix V estimated data-driven to prioritise parallel trends
#   5. Country derived from LAD24CD — exact matching enforced within loop
# =============================================================================

library(tidyverse)   # includes dplyr, ggplot2, tidyr, magrittr (%>%)
library(MatchIt)
library(cobalt)
library(WeightIt)
library(MASS)        # for ginv() — generalised matrix inverse

# -----------------------------------------------------------------------------
# STEP 0 — DERIVE COUNTRY + PRE-MATCHING EXCLUSIONS
# -----------------------------------------------------------------------------
# Country is derived from the LAD24CD prefix:
#   "E" = England  (e.g. E08000001)
#   "S" = Scotland (e.g. S12000033)
#
# WHY THIS MATTERS:
#   England and Scotland use different road casualty recording systems.
#   STATS19 (England/Wales) and STATS19-Scotland differ in injury thresholds,
#   classification protocols, and reporting obligations.
#   Cross-country matches would introduce systematic measurement incomparability
#   into the control group — exact matching on country prevents this entirely.
#
# NOTE: We do NOT split the code or run two separate pipelines.
#   Country is enforced as an exact matching constraint within the cohort loop.
#   This keeps the pipeline unified, avoids code duplication, and ensures
#   balance diagnostics are reported consistently across all schemes.

oa_data_clean <- OA_matching_dataset %>%
  
  # --- Derive country from LAD24CD prefix ---
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"   # flag unexpected prefixes
    )
  ) %>%
  
  # --- Verify no unknown country codes ---
  # If you see warnings here, check LAD24CD values for unexpected prefixes
  { . ->> oa_data_with_country } %>%   # save pre-filter for diagnostics
  filter(country != "Unknown") %>%
  
  filter(
    # 1. Keep only treated OAs or clean control OAs
    #    treated_OA == 1: inside CAZ/LEZ zone
    #    control_group1_OA or control_group2_OA == 1: eligible controls
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    
    # 2. Remove buffer OAs (spillover / displacement contamination)
    buffer_OA == 0,
    
    # 3. Remove zero-injury OAs — uninformative for outcome matching
    #    and can distort trend estimates
    zero_injury_OA == 0
  )

# --- Country composition check ---
cat("=== COUNTRY COMPOSITION ===\n")
cat("Unknown LAD24CD prefixes flagged:",
    sum(oa_data_with_country$country == "Unknown"), "\n")

country_tab <- oa_data_clean %>%
  count(country, treated_OA) %>%
  tidyr::pivot_wider(names_from = treated_OA,
                     values_from = n,
                     names_prefix = "treated_") %>%
  rename(n_control = treated_0, n_treated = treated_1)

cat("\nOAs by country and treatment status:\n")
print(country_tab)

cat("\nOAs after exclusions:", nrow(oa_data_clean), "\n")
cat("  Treated:", sum(oa_data_clean$treated_OA), "\n")
cat("  Controls:", sum(oa_data_clean$treated_OA == 0), "\n")

# -----------------------------------------------------------------------------
# STEP 0b — DEFINE MATCHING VARIABLES BY TIER
# -----------------------------------------------------------------------------
# Tier 1: Pre-treatment injury outcomes — highest priority for parallel trends
# Tier 2: Road network + structural exposure — proximal confounders
# Note: traffic exposure (AADF/VKT) not in dataset — road_length_km used instead

tier1_vars <- c(
  # --- Injury LEVELS (mean pre-treatment counts per road-km) ---
  "mean_car_KSI_pkm",
  "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",
  "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",
  "mean_ped_slight_pkm",
  "mean_total_pkm",
  
  # --- Injury TRENDS (slope of log-rate change per quarter per road-km) ---
  "trend_car_KSI_pkm",
  "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",
  "trend_ped_slight_pkm",
  "trend_total_pkm"
)

tier2_vars <- c(
  # --- Road network characteristics ---
  "road_density_m_km2",   # road density
  "road_length_km",        # total road length — exposure denominator
  "pct_A_road",            # % A-road
  "pct_B_road",            # % B-road
  "pct_minor_road",        # % minor road
  "n_roads",               # total number of road segments
  "dist_citycentre"        # distance to city centre — urban gradient proxy
)

# Note: sociodemographic vars (IMD, census) not in current dataset
# Add them here if/when merged: imd_score, pct_car_ownership, etc.

all_match_vars <- c(tier1_vars, tier2_vars)

cat("\nMatching variables:\n")
cat("  Tier 1 (outcome levels + trends):", length(tier1_vars), "\n")
cat("  Tier 2 (road network):", length(tier2_vars), "\n")
cat("  Total:", length(all_match_vars), "\n")

# -----------------------------------------------------------------------------
# STEP 0c — ESTIMATE OPTIMAL WEIGHT MATRIX V (Abadie et al. 2010)
# -----------------------------------------------------------------------------
# V is estimated by minimising mean squared prediction error of pre-treatment
# OUTCOME variables (tier 1 only) using all matching variables as predictors.
# This is the Abadie et al. data-driven approach — weights are not manually set.
#
# Implementation:
#   For each outcome variable in tier 1, regress on all predictors.
#   V weight for each predictor = average squared coefficient across outcomes.
#   This reflects how much each predictor contributes to predicting outcomes.
#   Higher V weight = more influence on distance = closer matching required.
# -----------------------------------------------------------------------------

estimate_V_matrix <- function(data, outcome_vars, predictor_vars) {
  
  # Check all variables exist
  missing_vars <- setdiff(c(outcome_vars, predictor_vars), names(data))
  if (length(missing_vars) > 0) {
    warning("Missing variables: ", paste(missing_vars, collapse = ", "))
    outcome_vars   <- intersect(outcome_vars, names(data))
    predictor_vars <- intersect(predictor_vars, names(data))
  }
  
  # Standardise all variables (mean 0, SD 1)
  data_std <- data %>%
    mutate(across(all_of(c(outcome_vars, predictor_vars)),
                  ~ as.numeric(scale(.))))
  
  X <- as.matrix(data_std[, predictor_vars])
  Y <- as.matrix(data_std[, outcome_vars])
  
  n_pred    <- ncol(X)
  v_weights <- numeric(n_pred)
  names(v_weights) <- predictor_vars
  
  n_outcomes_used <- 0
  
  for (y_var in colnames(Y)) {
    y     <- Y[, y_var]
    valid <- complete.cases(cbind(y, X))
    
    if (sum(valid) < 20) {
      cat("    Skipping outcome", y_var, "— too few complete cases\n")
      next
    }
    
    fit   <- lm(y ~ X[valid, ] - 1)  # no intercept — standardised vars
    coefs <- coef(fit)
    coefs[is.na(coefs)] <- 0          # handle multicollinearity
    
    # Weight = squared coefficient normalised by number of outcomes
    v_weights <- v_weights + coefs^2
    n_outcomes_used <- n_outcomes_used + 1
  }
  
  if (n_outcomes_used == 0) stop("No valid outcomes for V estimation")
  
  v_weights <- v_weights / n_outcomes_used
  
  # Normalise so weights sum to number of predictors (preserves scale)
  v_weights <- v_weights / sum(v_weights) * n_pred
  
  # Build diagonal weight matrix
  V <- diag(v_weights)
  dimnames(V) <- list(predictor_vars, predictor_vars)
  
  return(list(V = V, v_weights = v_weights, n_outcomes = n_outcomes_used))
}

# -----------------------------------------------------------------------------
# STEP 0d — BUILD WEIGHTED MAHALANOBIS DISTANCE MATRIX USING V
# -----------------------------------------------------------------------------
# Weighted Mahalanobis distance = (x-y)' (V^{1/2} S^{-1} V^{1/2}) (x-y)
# where S = covariance matrix of predictors, V = diagonal weight matrix
# Variables with higher V weight contribute more to the distance metric.

build_weighted_mahal_vcov <- function(data, match_vars, v_weights) {
  
  X     <- as.matrix(data[, match_vars])
  X_std <- scale(X)
  
  # Pooled covariance matrix
  S <- cov(X_std, use = "pairwise.complete.obs")
  
  # Regularise if nearly singular (add small ridge)
  S <- S + diag(1e-6, nrow(S))
  
  # Inverse — use generalised inverse for safety
  S_inv <- tryCatch(solve(S), error = function(e) {
    cat("    Warning: singular covariance matrix — using generalised inverse\n")
    MASS::ginv(S)
  })
  
  # Scale by V: variables with higher weight get smaller effective distance
  v  <- v_weights[match_vars]
  v[is.na(v)] <- min(v, na.rm = TRUE)  # safety fallback
  
  V_sqrt    <- diag(sqrt(v))
  S_inv_V   <- V_sqrt %*% S_inv %*% V_sqrt
  
  return(S_inv_V)
}

# -----------------------------------------------------------------------------
# STEP 1 — COHORT-STRATIFIED MATCHING LOOP WITH DATA-DRIVEN V
# -----------------------------------------------------------------------------
# Cohort defined by 'scheme' variable — each CAZ/LEZ scheme is a cohort.
# Never-treated OAs (treated_OA == 0) are eligible controls for all cohorts.

schemes <- oa_data_clean %>%
  filter(treated_OA == 1) %>%
  pull(scheme) %>%
  unique() %>%
  sort()

cat("\nSchemes to match:", paste(schemes, collapse = ", "), "\n")

matched_oa_list <- list()

for (g in schemes) {
  
  cat("\n--- Matching scheme:", g, "---\n")
  
  # Cohort dataset: treated OAs for this scheme + all never-treated OAs
  oa_cohort <- oa_data_clean %>%
    filter(
      (treated_OA == 1 & scheme == g) |   # treated: this scheme only
        (treated_OA == 0)                   # controls: never treated
    ) %>%
    mutate(treat_indicator = as.integer(treated_OA == 1 & scheme == g))
  
  # --- Check country composition within this cohort ---
  cat("  Country composition:\n")
  print(table(oa_cohort$country, oa_cohort$treat_indicator,
              dnn = c("Country", "Treated")))
  
  # IMPORTANT: if a scheme is Scotland-only (e.g. Glasgow LEZ),
  # all treated OAs will be Scottish. Control OAs from England will be
  # excluded by the exact = ~ country constraint automatically.
  # Check that enough Scottish controls exist:
  n_scot_controls <- sum(oa_cohort$treat_indicator == 0 &
                           oa_cohort$country == "Scotland")
  n_eng_controls  <- sum(oa_cohort$treat_indicator == 0 &
                           oa_cohort$country == "England")
  cat("  Scottish controls available:", n_scot_controls, "\n")
  cat("  English controls available:", n_eng_controls, "\n")
  
  n_treated  <- sum(oa_cohort$treat_indicator)
  n_controls <- sum(oa_cohort$treat_indicator == 0)
  cat("  Treated:", n_treated, "| Controls available:", n_controls, "\n")
  
  if (n_treated < 5) {
    cat("  Skipping — too few treated OAs (n =", n_treated, ")\n")
    next
  }
  
  # Check all matching variables exist in this cohort
  missing <- setdiff(all_match_vars, names(oa_cohort))
  if (length(missing) > 0) {
    cat("  WARNING — missing variables:", paste(missing, collapse = ", "), "\n")
    match_vars_g <- intersect(all_match_vars, names(oa_cohort))
  } else {
    match_vars_g <- all_match_vars
  }
  
  # ---------------------------------------------------------------------------
  # ESTIMATE V ON CONTROL OAs ONLY
  # Using controls only avoids treatment contamination in V estimation
  # ---------------------------------------------------------------------------
  control_data <- oa_cohort %>%
    filter(treat_indicator == 0) %>%
    dplyr::select(all_of(match_vars_g)) %>%
    filter(complete.cases(pick(everything())))
  
  cat("  Estimating V matrix on", nrow(control_data), "control OAs...\n")
  
  V_result <- tryCatch(
    estimate_V_matrix(
      data           = control_data,
      outcome_vars   = intersect(tier1_vars, match_vars_g),
      predictor_vars = match_vars_g
    ),
    error = function(e) {
      cat("  V estimation failed:", conditionMessage(e),
          "— falling back to standard MDM\n")
      NULL
    }
  )
  
  # Top 5 variables by V weight (for reporting)
  if (!is.null(V_result)) {
    cat("  Top 5 variables by V weight:\n")
    top5 <- sort(V_result$v_weights, decreasing = TRUE)[1:5]
    print(round(top5, 4))
  }
  
  # ---------------------------------------------------------------------------
  # RUN MATCHING
  # ---------------------------------------------------------------------------
  # If V estimation succeeded: use weighted Mahalanobis distance
  # If V estimation failed: fall back to standard Mahalanobis (still valid)
  # ---------------------------------------------------------------------------
  # RUN MATCHING
  # ---------------------------------------------------------------------------
  # If V estimation succeeded: use weighted Mahalanobis distance
  # If V estimation failed: fall back to standard Mahalanobis (still valid)
  
  match_formula <- reformulate(match_vars_g, response = "treat_indicator")
  
  m_out <- tryCatch({
    
    if (!is.null(V_result)) {
      # --- Option A: Weighted Mahalanobis using estimated V ---
      # Build distance matrix manually and pass to MatchIt
      oa_complete <- oa_cohort %>%
        dplyr::select(treat_indicator, all_of(match_vars_g)) %>%
        filter(complete.cases(pick(everything())))
      
      S_inv_V <- build_weighted_mahal_vcov(
        data       = oa_complete,
        match_vars = match_vars_g,
        v_weights  = V_result$v_weights
      )
      
      # MatchIt API note:
      # When distance = "mahalanobis", the formula variables define the
      # Mahalanobis space. mahvars is NOT valid here — instead we pass our
      # custom weighted inverse covariance matrix via the vcov argument.
      # This is the correct way to implement Abadie et al. weighted MDM.
      matchit(
        match_formula,
        data        = oa_cohort,
        method      = "nearest",
        distance    = "mahalanobis",
        vcov        = S_inv_V,       # custom V-weighted inverse covariance matrix
        caliper     = 0.25,
        std.caliper = TRUE,
        ratio       = 4,
        replace     = FALSE,
        
        # EXACT MATCHING ON COUNTRY — non-negotiable
        # England (STATS19) and Scotland (STATS19-Scotland) differ in injury
        # classification thresholds and recording protocols. Cross-country
        # matches would introduce systematic measurement incomparability.
        # This constraint forbids them outright.
        exact = ~ country
        # NOTE: if you later add urban_rural_class, extend to:
        # exact = ~ country + urban_rural_class
      )
      
    } else {
      # --- Option B: Standard Mahalanobis fallback (country exact still enforced) ---
      matchit(
        match_formula,
        data        = oa_cohort,
        method      = "nearest",
        distance    = "mahalanobis",
        caliper     = 0.25,
        std.caliper = TRUE,
        ratio       = 4,
        replace     = FALSE,
        exact       = ~ country   # always enforce — even in fallback
      )
    }
    
  }, error = function(e) {
    cat("  Matching failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(m_out)) next
  
  # Store result with V weights for appendix reporting
  matched_oa_list[[as.character(g)]] <- list(
    matchit_obj = m_out,
    scheme      = g,
    V_weights   = if (!is.null(V_result)) V_result$v_weights else NULL,
    n_treated   = n_treated,
    n_controls  = n_controls
  )
  
  # --- Verify no cross-country matches slipped through ---
  matched_data <- match.data(m_out)
  cross_country <- matched_data %>%
    group_by(subclass) %>%
    summarise(n_countries = n_distinct(country), .groups = "drop") %>%
    filter(n_countries > 1)
  
  if (nrow(cross_country) > 0) {
    cat("  WARNING:", nrow(cross_country),
        "matched pairs contain cross-country matches — investigate!\n")
  } else {
    cat("  Country check passed — no cross-country matches\n")
  }
  
  # ---------------------------------------------------------------------------
  # BALANCE DIAGNOSTICS — non-negotiable for every scheme
  # ---------------------------------------------------------------------------
  cat("\n  Balance summary for scheme", g, ":\n")
  print(bal.tab(m_out, thresholds = c(m = 0.1, v = 2)))
  
  lp <- love.plot(
    m_out,
    threshold    = 0.1,
    abs          = TRUE,
    var.order    = "unadjusted",
    title        = paste0("Covariate Balance — Scheme: ", g),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  )
  ggsave(
    filename = paste0("balance_scheme_", gsub(" ", "_", g), ".png"),
    plot     = lp,
    width    = 10, height = 8
  )
}

# -----------------------------------------------------------------------------
# STEP 2 — EXTRACT MATCHED CONTROL OA IDs (ACROSS ALL SCHEMES)
# -----------------------------------------------------------------------------

matched_control_oas <- lapply(names(matched_oa_list), function(g_chr) {
  m <- matched_oa_list[[g_chr]]$matchit_obj
  
  match.data(m) %>%
    filter(treat_indicator == 0) %>%
    mutate(matched_for_scheme = g_chr) %>%
    dplyr::select(OA, matched_for_scheme, weights, subclass)
}) %>%
  bind_rows() %>%
  # OA may be matched as control for multiple schemes — keep unique
  distinct(OA, .keep_all = TRUE)

matched_treated_oas <- lapply(names(matched_oa_list), function(g_chr) {
  m <- matched_oa_list[[g_chr]]$matchit_obj
  
  match.data(m) %>%
    filter(treat_indicator == 1) %>%
    mutate(matched_for_scheme = g_chr) %>%
    dplyr::select(OA, matched_for_scheme, weights, subclass)
}) %>%
  bind_rows()

cat("\n=== MATCHING SUMMARY ===\n")
cat("Matched control OAs:", nrow(matched_control_oas), "\n")
cat("Matched treated OAs:", nrow(matched_treated_oas), "\n")

# -----------------------------------------------------------------------------
# STEP 2b — REPORT V WEIGHTS ACROSS SCHEMES (for Appendix Table S1.3)
# -----------------------------------------------------------------------------
# Average V weights across schemes — shows which variables drove matching
# Replace manual 4x/3x/2x/1x weights in Table S1.3 with these estimates

V_weights_summary <- lapply(names(matched_oa_list), function(g_chr) {
  vw <- matched_oa_list[[g_chr]]$V_weights
  if (is.null(vw)) return(NULL)
  as.data.frame(t(vw)) %>% mutate(scheme = g_chr)
}) %>%
  bind_rows()

if (nrow(V_weights_summary) > 0) {
  V_weights_avg <- V_weights_summary %>%
    dplyr::select(-scheme) %>%
    summarise(across(everything(), ~ mean(., na.rm = TRUE))) %>%
    tidyr::pivot_longer(
      everything(),
      names_to  = "variable",
      values_to = "mean_V_weight"
    ) %>%
    arrange(desc(mean_V_weight))
  
  cat("\n=== AVERAGE V WEIGHTS ACROSS SCHEMES (Table S1.3) ===\n")
  print(V_weights_avg, n = 30)
  
  # Tag each variable with its tier for table
  V_weights_avg <- V_weights_avg %>%
    mutate(tier = case_when(
      variable %in% tier1_vars ~ "Tier 1: Injury outcomes",
      variable %in% tier2_vars ~ "Tier 2: Road network",
      TRUE                     ~ "Other"
    ))
  
  write.csv(V_weights_avg,
            "V_weights_table_S1_3.csv",
            row.names = FALSE)
  cat("V weights saved to V_weights_table_S1_3.csv\n")
}

# -----------------------------------------------------------------------------
# STEP 3 — BUILD RESTRICTED ROAD-LINK PANEL
# -----------------------------------------------------------------------------

road_link_panel_restricted <- road_link_panel %>%
  filter(
    # Keep all treated road links
    treated_link == 1 |
      # Keep control road links ONLY if OA is in matched donor pool
      (treated_link == 0 &
         OA %in% matched_control_oas$OA &
         buffer_OA == FALSE)
  ) %>%
  # Merge OA-level covariates for C&S covariate adjustment (xformla)
  left_join(
    oa_data_clean %>% dplyr::select(
      OA, road_density_m_km2, pct_A_road, pct_B_road, pct_minor_road,
      road_length_km, dist_citycentre, n_roads
    ),
    by = "OA"
  )

cat("\n=== ROAD LINK PANEL ===\n")
cat("Road links in restricted panel:", nrow(road_link_panel_restricted), "\n")
cat("Unique OAs in panel:", n_distinct(road_link_panel_restricted$OA), "\n")

# -----------------------------------------------------------------------------
# STEP 4 — CALLAWAY & SANT'ANNA ESTIMATION
# -----------------------------------------------------------------------------

library(did)

outcomes <- c("KSI_adj", "Slight_adj")
modes    <- c("All", "Pedestrian", "Car.Van", "Cyclist", "Other")

run_cs <- function(yvar, data) {
  att_gt(
    yname         = yvar,
    tname         = "quarter_num",
    idname        = "road_link_id",
    gname         = "cohort_g",          # first treated quarter; 0 = never treated
    
    xformla       = ~ road_density_m_km2 + pct_A_road + pct_minor_road +
      road_length_km + dist_citycentre,
    
    est_method    = "dr",                # doubly robust
    control_group = "nevertreated",
    
    clustervars   = "OA",                # cluster at OA level
    
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

# Aggregate
cs_agg <- lapply(cs_results, function(res) {
  if (is.null(res)) return(NULL)
  list(
    simple  = aggte(res, type = "simple"),
    dynamic = aggte(res, type = "dynamic"),
    group   = aggte(res, type = "group")
  )
})

# -----------------------------------------------------------------------------
# STEP 5 — SENSITIVITY / ROBUSTNESS CHECKS
# -----------------------------------------------------------------------------

# 5a. Alternative calipers
# Re-run matching loop with caliper = 0.15 (tighter) and 0.35 (looser)

# 5b. Alternative ratio
# Try ratio = 1 (1:1) and ratio = 6 (1:6) — test precision vs quality tradeoff

# 5c. Control group: re-run C&S with control_group = "notyettreated"
# (use control_group2_OA as not-yet-treated pool)

# 5d. HonestDiD — pre-trend sensitivity (Rambachan & Roth 2023)
library(HonestDiD)
# Apply to each dynamic aggregation — see HonestDiD vignette

# 5e. Manual weights as robustness check
# Re-run matching with fixed weights (Tier1=4x, Tier2=3x) instead of V matrix
# If results are stable, this strengthens the data-driven V choice

# -----------------------------------------------------------------------------
# STEP 6 — BALANCE REPORTING (PUBLICATION-READY)
# -----------------------------------------------------------------------------

balance_summary <- sapply(names(matched_oa_list), function(g_chr) {
  bt <- bal.tab(
    matched_oa_list[[g_chr]]$matchit_obj,
    thresholds = c(m = 0.1),
    un = TRUE
  )
  smd_after <- abs(bt$Balance$Diff.Adj)
  c(
    scheme          = g_chr,
    n_treated       = matched_oa_list[[g_chr]]$n_treated,
    pct_balanced    = round(mean(smd_after < 0.1, na.rm = TRUE) * 100, 1),
    max_smd         = round(max(smd_after, na.rm = TRUE), 3),
    mean_smd        = round(mean(smd_after, na.rm = TRUE), 3)
  )
}, simplify = FALSE) %>%
  bind_rows()

cat("\n=== BALANCE SUMMARY ACROSS SCHEMES ===\n")
print(balance_summary)

write.csv(balance_summary, "balance_summary_all_schemes.csv", row.names = FALSE)
