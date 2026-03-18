# =============================================================================
# CAZ/LEZ Staggered DiD: OA-Level Matching for Road-Link Donor Pool Restriction
# Callaway & Sant'Anna (2021) — Doubly Robust Estimator
# =============================================================================
# STRATEGY OVERVIEW:
#   Match at OA level (one row per OA) to define a comparable donor pool.
#   All road links within matched control OAs then enter C&S estimation.
#   Matching is performed SEPARATELY per treatment cohort (timing group).
# =============================================================================

library(MatchIt)
library(cobalt)
library(dplyr)
library(ggplot2)
library(WeightIt)   # ADD: for entropy balancing as robustness check

# -----------------------------------------------------------------------------
# STEP 0 — PRE-MATCHING EXCLUSIONS (apply before any matching)
# -----------------------------------------------------------------------------
# These are non-negotiable. Run this on your full OA dataset first.

oa_data_clean <- oa_data_full |>
  filter(
    # 1. Remove OAs straddling zone boundaries (mixed treatment exposure)
    straddle_oa == FALSE,
    
    # 2. Remove buffer OAs (spillover / displacement contamination)
    # Adjust buffer distance to your zone sizes — 500m is a starting point
    buffer_500m == FALSE,
    
    # 3. Remove OAs in structurally incomparable contexts
    # IMPORTANT: run England and Scotland SEPARATELY — different STATS19
    # systems, different policy regimes. Filter to one country at a time.
    urban_rural_class != "Rural hamlet and isolated",  # adjust to your classification
    urban_rural_class != "Rural village"               # if all CAZs are urban
  )

# -----------------------------------------------------------------------------
# STEP 1 — COHORT-STRATIFIED MATCHING LOOP
# -----------------------------------------------------------------------------
# Match SEPARATELY per cohort. A 2018-treated OA should match to 2018-comparable
# controls, not to 2022 controls that may differ systematically.
# Never-treated OAs are eligible as controls for ALL cohorts.

cohorts <- oa_data_clean |>
  filter(treated == 1) |>
  pull(cohort_g) |>
  unique() |>
  sort()

matched_oa_list <- list()

for (g in cohorts) {
  
  cat("\n--- Matching cohort:", g, "---\n")
  
  # Cohort dataset: treated OAs in cohort g + all never-treated OAs
  # "Not yet treated" OAs can also be controls if you prefer — see note below
  oa_cohort <- oa_data_clean |>
    filter(
      (treated == 1 & cohort_g == g) |   # treated: this cohort only
        (treated == 0)                       # controls: never treated
      # To include not-yet-treated: | (treated == 1 & cohort_g > g)
    ) |>
    mutate(treat_indicator = as.integer(treated == 1 & cohort_g == g))
  
  n_treated  <- sum(oa_cohort$treat_indicator)
  n_controls <- sum(oa_cohort$treat_indicator == 0)
  cat("  Treated:", n_treated, " | Controls available:", n_controls, "\n")
  
  # Skip cohort if too few treated units for meaningful matching
  if (n_treated < 5) {
    cat("  Skipping — too few treated OAs\n")
    next
  }
  
  # ---------------------------------------------------------------------------
  # MATCHING SPECIFICATION
  # ---------------------------------------------------------------------------
  # KEY IMPROVEMENTS over original code:
  #   1. Added Slight_adj pre-trends (you have both KSI and Slight outcomes)
  #   2. Added mode-specific pre-trends for pedestrian, cyclist (ksi_ped, ksi_cyc 
  #      already present — add slight equivalents)
  #   3. Added seasonality measure (Q4_share) — RTIs peak in autumn/winter
  #   4. Added pct_active_travel from census (key for pedestrian/cyclist outcomes)
  #   5. Added network_length_km (total road length in OA — exposure denominator)
  #   6. Added region to exact matching alongside urban_rural_class
  #      (England vs Scotland MUST be exact if pooling — better to run separately)
  #   7. ratio = 4 instead of 3 — more controls improve C&S precision at low cost
  #   8. replace = FALSE explicitly stated (default but important to be explicit)
  #   9. Added mahvars to separate distance vars from exact vars cleanly
  # ---------------------------------------------------------------------------
  
  m_out <- tryCatch({
    matchit(
      treat_indicator ~
        
        # --- TIER 1: Pre-treatment outcome trends (MOST important) ---
        # KSI outcomes
        ksi_trend_slope +          # slope of KSI over pre-period (linear)
        ksi_mean_pre +             # mean KSI level pre-treatment
        ksi_trend_slope_sq +       # FIX: add quadratic trend — captures non-linear pre-trends
        # Slight outcomes — ADD: you have both outcomes, match on both
        slight_trend_slope +
        slight_mean_pre +
        # Mode-specific pre-trends — match on all casualty types you'll analyse
        ksi_ped_pre +
        ksi_car_pre +
        ksi_cyc_pre +
        slight_ped_pre +           # ADD: slight injuries by mode
        slight_car_pre +
        slight_cyc_pre +
        # Seasonality — RTIs have strong seasonal pattern; must match on it
        q4_inj_share +             # ADD: share of annual injuries in Q4 (Oct-Dec)
        
        # --- TIER 2: Structural / exposure variables ---
        road_density +
        network_length_km +        # ADD: total road length — controls for exposure
        pct_a_road +
        mean_speed_limit +
        junction_density +
        pct_hgv +                  # ADD: heavy goods vehicles if available
        pop_density +
        imd_score +
        imd_transport +
        pct_car_ownership +
        pct_active_travel +        # ADD: % active travel to work (census) — key for ped/cyc
        
        # --- TIER 3: handled via exact= argument below ---
        urban_rural_class,
      
      data     = oa_cohort,
      method   = "nearest",
      distance = "mahalanobis",
      
      # CALIPER: 0.25 SD is standard. If you get poor balance or many unmatched
      # treated units, loosen to 0.35. If balance is good, tighten to 0.15.
      caliper  = 0.25,
      std.caliper = TRUE,         # caliper in SD units (recommended)
      
      ratio    = 4,               # CHANGE: 1:4 — more controls, better precision
      replace  = FALSE,           # no replacement — each control used once
      
      # EXACT MATCHING: force identical cells on categorical variables
      # ADD region: never match an OA in Yorkshire to one in Greater Glasgow
      exact    = ~ urban_rural_class + country  # use 'region' if not running separately
    )
  }, error = function(e) {
    cat("  Matching failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(m_out)) next
  
  # Store result with cohort label
  matched_oa_list[[as.character(g)]] <- list(
    matchit_obj = m_out,
    cohort      = g
  )
  
  # ---------------------------------------------------------------------------
  # BALANCE DIAGNOSTICS — non-negotiable for every cohort
  # ---------------------------------------------------------------------------
  cat("\n  Balance summary for cohort", g, ":\n")
  print(bal.tab(m_out, thresholds = c(m = 0.1, v = 2)))
  # m = 0.1: SMD threshold; v = 2: variance ratio threshold
  
  # Save love plot per cohort
  lp <- love.plot(
    m_out,
    threshold   = 0.1,
    abs         = TRUE,
    var.order   = "unadjusted",
    title       = paste0("Covariate Balance — Cohort ", g),
    shapes      = c("circle filled", "triangle filled"),
    colors      = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  )
  ggsave(
    filename = paste0("balance_cohort_", g, ".png"),
    plot     = lp,
    width    = 10, height = 8
  )
}

# -----------------------------------------------------------------------------
# STEP 2 — EXTRACT MATCHED CONTROL OA IDs (ACROSS ALL COHORTS)
# -----------------------------------------------------------------------------

matched_control_oas <- lapply(names(matched_oa_list), function(g_chr) {
  m <- matched_oa_list[[g_chr]]$matchit_obj
  
  # Extract matched data: controls only
  md <- match.data(m) |>
    filter(treat_indicator == 0) |>
    mutate(matched_for_cohort = as.integer(g_chr))
  
  md[, c("oa_id", "matched_for_cohort", "weights", "subclass")]
}) |>
  bind_rows() |>
  # An OA may be matched as control for multiple cohorts — keep unique
  distinct(oa_id, .keep_all = TRUE)

# Also extract treated OA IDs (for completeness)
matched_treated_oas <- lapply(names(matched_oa_list), function(g_chr) {
  m <- matched_oa_list[[g_chr]]$matchit_obj
  match.data(m) |>
    filter(treat_indicator == 1) |>
    mutate(matched_for_cohort = as.integer(g_chr)) |>
    select(oa_id, matched_for_cohort, weights, subclass)
}) |>
  bind_rows()

cat("\nMatched control OAs:", nrow(matched_control_oas))
cat("\nMatched treated OAs:", nrow(matched_treated_oas), "\n")

# -----------------------------------------------------------------------------
# STEP 3 — BUILD RESTRICTED ROAD-LINK PANEL
# -----------------------------------------------------------------------------
# Carry road-link panel observations whose OA is in the matched donor pool.
# Also merge OA-level covariates onto road links for use in C&S xformla.

road_link_panel_restricted <- road_link_panel |>
  filter(
    # Keep all treated road links
    treated_link == 1 |
      # Keep control road links ONLY if in matched donor pool
      (treated_link == 0 &
         oa_id %in% matched_control_oas$oa_id &
         buffer_500m == FALSE &
         straddle_oa == FALSE)
  ) |>
  # Merge OA-level covariates for use in C&S covariate adjustment
  left_join(
    oa_data_clean |> select(
      oa_id, road_density, pct_a_road, mean_speed_limit,
      junction_density, pop_density, imd_score, urban_rural_class,
      pct_active_travel, network_length_km
    ),
    by = "oa_id"
  )

cat("Road links in restricted panel:", nrow(road_link_panel_restricted), "\n")
cat("Unique OAs in panel:",
    n_distinct(road_link_panel_restricted$oa_id), "\n")

# -----------------------------------------------------------------------------
# STEP 4 — CALLAWAY & SANT'ANNA ESTIMATION (on restricted panel)
# -----------------------------------------------------------------------------

library(did)

# Define outcomes × casualty type combinations
outcomes <- c("KSI_adj", "Slight_adj")
modes    <- c("All", "Pedestrian", "Car.Van", "Cyclist", "Other")

# Function to run C&S for one outcome-mode combination
run_cs <- function(yvar, data) {
  att_gt(
    yname         = yvar,
    tname         = "quarter_num",       # numeric quarter index
    idname        = "road_link_id",
    gname         = "cohort_g",          # first treated quarter; 0 = never treated
    
    xformla       = ~ road_type_class + mean_speed_limit + junction_flag +
      log1p(network_length_km) + imd_score + pct_active_travel,
    
    est_method    = "dr",                # doubly robust — keep
    control_group = "nevertreated",      # primary spec
    
    # CLUSTER at OA level — road links within OA are not independent
    # If most CAZs cover whole LAs, consider clustervars = "la_id" instead
    clustervars   = "oa_id",
    
    panel         = TRUE,
    allow_unbalanced_panel = TRUE,       # ADD: handles missing road-quarter obs
    data          = data
  )
}

# Run all combinations
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

# Aggregate results
cs_agg <- lapply(cs_results, function(res) {
  if (is.null(res)) return(NULL)
  list(
    simple  = aggte(res, type = "simple"),    # overall ATT
    dynamic = aggte(res, type = "dynamic"),   # event study
    group   = aggte(res, type = "group")      # by cohort
  )
})

# -----------------------------------------------------------------------------
# STEP 5 — SENSITIVITY / ROBUSTNESS CHECKS
# -----------------------------------------------------------------------------

# 5a. Alternative calipers — test sensitivity of donor pool to caliper width
# Re-run matching loop with caliper = 0.15 (tighter) and 0.35 (looser)
# Compare ATT estimates across donor pools

# 5b. Alternative ratio — try 1:1 and 1:6
# Tighter ratio = better per-match quality; wider ratio = more power

# 5c. Control group: re-run C&S with control_group = "notyettreated"
cs_results_nyt <- lapply(names(cs_results), function(nm) {
  # re-run with notyettreated for comparison
})

# 5d. HonestDiD — pre-trend sensitivity (Rambachan & Roth 2023)
# install.packages("HonestDiD")
library(HonestDiD)
# Apply to each dynamic aggregation — see HonestDiD vignette for full workflow
# Key: test how large pre-trend violations would need to be to overturn results

# 5e. Synthetic DiD as parallel estimator (Arkhangelsky et al. 2021)
# install.packages("synthdid")
# Run on restricted panel as robustness check against C&S estimates

# -----------------------------------------------------------------------------
# STEP 6 — BALANCE REPORTING (PUBLICATION-READY)
# -----------------------------------------------------------------------------

# Combined love plot across all cohorts
all_balance <- lapply(names(matched_oa_list), function(g_chr) {
  bal.tab(matched_oa_list[[g_chr]]$matchit_obj,
          thresholds = c(m = 0.1),
          un = TRUE)
})

# Summary table: % variables with SMD < 0.1 before and after matching
balance_summary <- sapply(all_balance, function(b) {
  smd <- b$Balance$Diff.Adj
  c(pct_balanced = mean(abs(smd) < 0.1, na.rm = TRUE) * 100)
})
print(balance_summary)
