# =============================================================================
#  OA-Level Matching — Two-Stage MDM Design (No-Caliper Version)
#  Stage 1: MDM on structural + sociodemographic variables, no caliper,
#           95th-percentile distance trim post-hoc
#  Stage 2: MDM on pre-treatment injury trends + levels, no caliper,
#           95th-percentile distance trim post-hoc
#
# =============================================================================
#
# DESIGN RATIONALE
#
#   NO CALIPERS AT EITHER STAGE.
#   All matching variables are treated symmetrically. Poor matches are
#   removed at each stage via a post-hoc 95th-percentile Mahalanobis
#   distance trim. This approach:
#     (a) avoids imposing variable-specific threshold choices that
#         differentially weight structural vs sociodemographic variables,
#     (b) lets the distance distribution itself reveal where the support
#         for matching breaks down, and
#     (c) produces a single, transparent quality-control rule applied
#         consistently at both stages.
#
#   STAGE 1 — Structural + sociodemographic restriction (MDM, no caliper)
#     Nearest-neighbour MDM on 15 road-network, urban-form, and
#     sociodemographic variables. replace = TRUE, ratio = 10 to maximise
#     the candidate pool passed to Stage 2.
#     Post-hoc trim at P95 of Stage 1 Mahalanobis distances removes the
#     most structurally distant controls before injury-dynamics matching.
#     Country (England/Scotland) enforced as hard exact constraint at
#     both stages — STATS19 and STATS19-Scotland differ in recording
#     thresholds and severity classification; cross-country matches
#     introduce systematic measurement incomparability.
#
#   STAGE 2 — Injury-dynamics matching (MDM, no caliper)
#     Within the Stage 1-trimmed pool, 1:1 nearest-neighbour MDM on
#     14 pre-treatment injury variables (7 trends + 7 levels, all per
#     road-kilometre). #     Post-hoc trim at P95 of Stage 2 Mahalanobis distances removes
#     the most trajectory-distant pairs.
#
#   TWO PARALLEL ANALYSES:
#     Analysis A — EXCLUDING zero-injury treated OAs (zero_injury_OA==1)
#       Stage 2 injury variables are well-defined for all retained OAs.
#       ATT interpretation: effect on road-exposed OAs with pre-period
#       injury exposure.
#     Analysis B — INCLUDING zero-injury treated OAs.
#       Their Stage 2 injury variables are set to 0 (not NA).
#       They will match to zero-injury controls with small Stage 2
#       distances; parallel trends for these pairs rests on Stage 1
#       structural similarity only.
#       ATT interpretation: effect on all road-exposed OAs within
#       scheme boundaries.
#     Comparing A and B reveals whether zero-injury OAs materially
#     affect the ATT estimate. If similar => primary result is robust.
#     If divergent => the zero-injury subgroup requires separate
#     discussion.
#
#   ONLY FIXED PRE-MATCHING EXCLUSIONS (applied to both analyses):
#     Buffer OAs (treatment contamination risk) and zero-road OAs
#     (per-km rates undefined). No zero-injury exclusion in Analysis B.
#
# =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)
library(here)
library(MASS)
library(tidyverse)
library(purrr)
library(sf)

# =============================================================================
# SECTION 0 — DATA LOADING AND INITIAL CHECKS
# =============================================================================

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

cat("\n====================================================\n")
cat("SECTION 0: INITIAL CHECKS\n")
cat("====================================================\n\n")

cat("--- Assignment counts ---\n")
print(table(OA_matching_dataset$assignment))

cat("\n--- Zero-injury OA counts ---\n")
print(table(OA_matching_dataset$zero_injury_OA))

cat("\n--- Zero-injury breakdown by treatment status ---\n")
print(table(
  OA_matching_dataset$zero_injury_OA,
  OA_matching_dataset$treated_OA,
  dnn = c("zero_injury", "treated")
))

# Road-network check
cat("\n--- Road network summary ---\n")
OA_matching_dataset %>%
  summarise(
    total_OAs  = n_distinct(OA),
    zero_roads = sum(n_roads == 0),
    pct_zero   = round(100 * zero_roads / total_OAs, 2)
  ) %>% print()

# Zero-road OAs with recorded injuries — boundary artefacts
weirdOAs <- OA_matching_dataset %>%
  filter((n_roads == 0 | is.na(n_roads)) & mean_total > 0)
cat("Zero-road OAs with recorded injuries:", nrow(weirdOAs),
    " (treated:", sum(weirdOAs$treated_OA == 1), ")\n")

# =============================================================================
# SECTION 1 — VARIABLE DEFINITIONS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 1: VARIABLE DEFINITIONS\n")
cat("====================================================\n\n")

# Stage 1: all treated symmetrically 
stage1_road <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road"
)
stage1_urban <- c("dist_citycentre", "pop_density")
stage1_socdem <- c(
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_socdem)

# Stage 2: pre-treatment injury dynamics (per road-km)
stage2_trends <- c(
  "trend_car_KSI_pkm",    "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm",    "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm",    "trend_ped_slight_pkm",
  "trend_total_pkm"
)
stage2_levels <- c(
  "mean_car_KSI_pkm",     "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm",     "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm",     "mean_ped_slight_pkm",
  "mean_total_pkm"
)
stage2_vars <- c(stage2_trends, stage2_levels)

cat("Stage 1 variables:", length(stage1_vars),
    "(road:", length(stage1_road),
    "| urban:", length(stage1_urban),
    "| socdem:", length(stage1_socdem), ")\n")
cat("Stage 2 variables:", length(stage2_vars),
    "(trends:", length(stage2_trends),
    "| levels:", length(stage2_levels), ")\n")
cat("No calipers applied at either stage.\n")

# =============================================================================
# SECTION 2 — BUILD BOTH ANALYSIS DATASETS
# =============================================================================
# Analysis A: exclude zero-injury treated OAs (pre-period injury required)
# Analysis B: include zero-injury treated OAs (Stage 2 vars set to 0, not NA)
#
# Both exclude: buffer OAs + zero-road OAs
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 2: DATASET CONSTRUCTION\n")
cat("====================================================\n\n")

# Country derivation (shared)
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

# Shared exclusions: buffer OAs + zero-road OAs
base_filter <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# Analysis A: also exclude zero-injury treated OAs
data_A <- base_filter %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1))

# Analysis B: include zero-injury OAs, set Stage 2 vars to 0
# Zero-injury OAs will have 0 for all 14 Stage 2 variables.
# This is intentional: they will match to zero-injury controls with
# small Stage 2 distances. Parallel trends for these pairs relies on
# Stage 1 structural similarity. The comparison with Analysis A
# reveals whether this matters empirically.
data_B <- base_filter %>%
  mutate(across(
    all_of(intersect(stage2_vars, names(.))),
    ~ if_else(zero_injury_OA == 1 & treated_OA == 1, 0, .)
  ))

cat("=== Analysis A (zero-injury EXCLUDED) ===\n")
cat("  Total OAs:   ", nrow(data_A), "\n")
cat("  Treated:     ", sum(data_A$treat_indicator == 1), "\n")
cat("  Controls:    ", sum(data_A$treat_indicator == 0), "\n\n")

cat("=== Analysis B (zero-injury INCLUDED) ===\n")
cat("  Total OAs:   ", nrow(data_B), "\n")
cat("  Treated:     ", sum(data_B$treat_indicator == 1), "\n")
cat("    of which zero-injury: ",
    sum(data_B$treat_indicator == 1 & data_B$zero_injury_OA == 1), "\n")
cat("  Controls:    ", sum(data_B$treat_indicator == 0), "\n\n")

cat("Country distribution — Analysis A:\n")
print(data_A %>%
        count(country, treat_indicator) %>%
        pivot_wider(names_from  = treat_indicator,
                    values_from = n,
                    names_prefix = "t") %>%
        rename(control = t0, treated = t1))

cat("\nCountry distribution — Analysis B:\n")
print(data_B %>%
        count(country, treat_indicator) %>%
        pivot_wider(names_from  = treat_indicator,
                    values_from = n,
                    names_prefix = "t") %>%
        rename(control = t0, treated = t1))

# =============================================================================
# SECTION 3 — PRE-PROCESSING: WINSORISE STAGE 1 VARIABLES
# =============================================================================
# Winsorising caps extreme values at a chosen percentile rather than
# removing them. Here we use the 1st/99th percentile of the full
# pre-matching sample distribution.
#
# WHY: Mahalanobis distance is computed using the covariance matrix of
# matching variables. One OA with a road density of 50,000 m/km2 (data
# error or genuine outlier) can distort the entire covariance structure
# and produce bad matches for all other OAs. Winsorising prevents this
# without discarding any rows.
#
# Applied to Stage 1 variables only here. Stage 2 variables are
# winsorised separately within the Stage 1 pool, using the treated
# distribution as the reference (see Section 5).
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 3: WINSORISING STAGE 1 VARIABLES\n")
cat("====================================================\n\n")
cat("Capping extreme values at 1st/99th percentile.\n")
cat("No rows removed — only extreme values capped.\n\n")

winsorise_s1 <- function(data, vars) {
  data %>%
    mutate(across(
      all_of(intersect(vars, names(.))),
      ~ {
        q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
        pmin(pmax(., q[1]), q[2])
      }
    ))
}

data_A_clean <- winsorise_s1(data_A, stage1_vars)
data_B_clean <- winsorise_s1(data_B, stage1_vars)

# Validate variables exist and have sufficient variance
check_vars <- function(data, vars, label) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0)
    cat("WARNING —", label, "missing:", paste(missing, collapse = ", "), "\n")
  vcheck <- data %>%
    summarise(across(all_of(intersect(vars, names(data))),
                     ~ var(., na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "v", values_to = "var")
  low <- vcheck %>% filter(var < 1e-8) %>% pull(v)
  if (length(low) > 0)
    cat("Dropping near-zero variance:", paste(low, collapse = ", "), "\n")
  setdiff(intersect(vars, names(data)), low)
}

s1_vars_A <- check_vars(data_A_clean, stage1_vars, "Stage 1 / Analysis A")
s1_vars_B <- check_vars(data_B_clean, stage1_vars, "Stage 1 / Analysis B")
s2_vars_A <- check_vars(data_A_clean, stage2_vars, "Stage 2 / Analysis A")
s2_vars_B <- check_vars(data_B_clean, stage2_vars, "Stage 2 / Analysis B")

# =============================================================================
# SECTION 4 — STAGE 1 MATCHING: NO CALIPER, 95TH PERCENTILE TRIM
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 4: STAGE 1 MATCHING (no caliper + P95 trim)\n")
cat("====================================================\n\n")

run_stage1 <- function(data, s1_vars, label) {
  
  cat("--- Stage 1:", label, "---\n")
  
  formula <- reformulate(s1_vars, response = "treat_indicator")
  
  m <- tryCatch(
    matchit(
      formula,
      data     = data,
      method   = "nearest",
      distance = "mahalanobis",
      ratio    = 10,
      replace  = TRUE,
      exact    = ~ country
    ),
    error = function(e) { 
      cat("FAILED:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(m)) return(NULL)
  
  md_raw <- match.data(m)
  
  cat("Rows in md_raw:", nrow(md_raw), "\n")
  cat("Treated units:", sum(md_raw$treat_indicator == 1), "\n")
  
  # ---------------------------
  # Match matrix
  # ---------------------------
  
  mm <- m$match.matrix
  
  if (is.null(mm)) {
    stop("match.matrix is NULL — matching failed.")
  }
  
  cat("Matched treated units:", nrow(mm), "\n")
  cat("Controls per treated (ratio):", ncol(mm), "\n")
  
  
  # ---------------------------
  # Covariance matrix
  # ---------------------------
  
  treated_rows <- as.numeric(rownames(mm))
  treated_df   <- data[treated_rows, ]
  
  S_s1 <- cov(treated_df[, s1_vars], use = "pairwise.complete.obs")
  
  # -------------------------------------------------
  # Compute Mahalanobis distance for matched pairs
  # -------------------------------------------------
  
  treated_df  <- md_raw %>% filter(treat_indicator == 1)
  controls_df <- md_raw %>% filter(treat_indicator == 0)
  
  S_s1 <- cov(treated_df[, s1_vars], use = "pairwise.complete.obs")
  
  dist_s1 <- purrr::map_df(seq_len(nrow(treated_df)), function(i){
    
    trow <- treated_df[i, ]
    
    purrr::map_df(seq_len(nrow(controls_df)), function(j){
      
      crow <- controls_df[j, ]
      
      tibble(
        treated_OA = trow$OA,
        control_OA = crow$OA,
        mdist = mahalanobis(
          x      = as.numeric(crow[, s1_vars]),
          center = as.numeric(trow[, s1_vars]),
          cov    = S_s1
        )
      )
      
    })
    
  })
  
  
  # ---------------------------
  # Trim worst 5% matches
  # ---------------------------
  
  p95 <- quantile(dist_s1$mdist, 0.95)
  
  kept_treated <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(mean_dist = mean(mdist), .groups = "drop") %>%
    filter(mean_dist <= p95) %>%
    pull(treated_OA)
  
  
  treated_trimmed <- md_raw %>%
    filter(treat_indicator == 1, OA %in% kept_treated)
  
  controls_trimmed <- md_raw %>%
    filter(treat_indicator == 0) %>%
    distinct(OA, .keep_all = TRUE)
  
  
  # ---------------------------
  # Diagnostics
  # ---------------------------
  
  n_in  <- sum(md_raw$treat_indicator == 1)
  n_out <- nrow(treated_trimmed)
  
  cat("  Treated entering:           ", n_in, "\n")
  cat("  P95 distance threshold:     ", round(p95, 2), "\n")
  cat("  Treated trimmed (worst 5%): ", n_in - n_out, "\n")
  cat("  Treated retained:           ", n_out, "\n")
  cat("  Unique controls in pool:    ", nrow(controls_trimmed), "\n")
  
  
  # ---------------------------
  # Balance
  # ---------------------------
  
  cat("\n  Stage 1 balance summary:\n")
  
  bt <- cobalt::bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  
  smd_df <- bt$Balance %>%
    tibble::rownames_to_column("variable") %>%
    dplyr::select(variable, Diff.Un, Diff.Adj) %>%
    dplyr::arrange(desc(abs(Diff.Adj)))
  
  print(smd_df)
  
  
  # ---------------------------
  # Country distribution
  # ---------------------------
  
  cat("\n  Country distribution (post-trim treated):\n")
  print(table(treated_trimmed$country))
  
  
  # ---------------------------
  # Love plot
  # ---------------------------
  
  lp <- cobalt::love.plot(
    m,
    threshold = 0.1,
    abs = TRUE,
    stars = "std",
    var.order = "unadjusted",
    title = paste("Stage 1 balance —", label),
    shapes = c("circle filled", "triangle filled"),
    colors = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) +
    theme(plot.margin = margin(10,30,10,180),
          legend.position = "bottom")
  
  ggsave(
    here::here("output", paste0("s1_balance_", label, ".png")),
    lp,
    width = 13,
    height = 9,
    dpi = 300
  )
  
  cat("  Love plot saved.\n\n")
  
  
  # ---------------------------
  # Output
  # ---------------------------
  
  list(
    matchit_obj      = m,
    matched_raw      = md_raw,
    dist_s1          = dist_s1,
    p95              = p95,
    kept_treated     = kept_treated,
    treated_trimmed  = treated_trimmed,
    controls_trimmed = controls_trimmed,
    bal              = bt
  )
}

s1_A <- run_stage1(data_A_clean, s1_vars_A, "A_excl_zero")

s1_B <- run_stage1(data_B_clean, s1_vars_B, "B_incl_zero")

# =============================================================================
# SECTION 5 — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE
# =============================================================================
# De-duplicate control OAs (ratio = 10 means one control can appear
# multiple times in Stage 1 matched data).
#
# Winsorise Stage 2 variables based on TREATED distribution only.
# This clips control OA values to the treated range — preventing
# extreme control values from appearing matchable when they are
# outside the treated support.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 5: STAGE 2 DATA PREPARATION\n")
cat("====================================================\n\n")

prepare_s2_data <- function(s1_result, s2_vars, label) {
  
  cat("--- Preparing Stage 2 input:", label, "---\n")
  
  s2_raw <- bind_rows(
    s1_result$treated_trimmed,
    s1_result$controls_trimmed
  ) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  
  cat("  Treated OAs: ", sum(s2_raw$treat_indicator == 1), "\n")
  cat("  Control OAs: ", sum(s2_raw$treat_indicator == 0), "\n")
  
  # Winsorise Stage 2 vars on TREATED distribution
  s2_present  <- intersect(s2_vars, names(s2_raw))
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  s2_data <- s2_raw %>%
    mutate(across(
      all_of(s2_present),
      ~ {
        q <- quantile(treated_ref[[cur_column()]], c(0.01, 0.99), na.rm = TRUE)
        pmin(pmax(., q[1]), q[2])
      }
    ))
  
  cat("  Stage 2 vars winsorised on treated distribution.\n\n")
  list(data = s2_data, s2_present = s2_present)
}

s2_prep_A <- prepare_s2_data(s1_A, s2_vars_A, "A_excl_zero")
s2_prep_B <- prepare_s2_data(s1_B, s2_vars_B, "B_incl_zero")
#### Both almost identical 
# A have slighlty better ballance --   Because zeros are true injuries,
# B_incl_zero is the correct and unbiased analytic sample, and its matching performance is just as strong as A.
# =============================================================================
# SECTION 6 — STAGE 2 MATCHING: NO CALIPER, 95TH PERCENTILE TRIM
# =============================================================================
# NOTE FOR ANALYSIS B:
# Zero-injury treated OAs have all Stage 2 variables = 0. They will
# match to zero-injury controls (also all zeros) with distance near 0.
# These pairs will NOT be caught by the P95 trim — near-zero distance
# means they look like excellent matches. This is expected and
# intentional. The comparison with Analysis A (Section 8) is what
# reveals whether their inclusion changes the ATT estimate.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 6: STAGE 2 MATCHING (no caliper + P95 trim)\n")
cat("====================================================\n\n")

compute_mdist <- function(m, data, vars) {
  md    <- match.data(m, data = data)
  treat <- md %>% filter(treat_indicator == 1)
  ctrl  <- md %>% filter(treat_indicator == 0)
  S     <- cov(treat[, vars], use = "pairwise.complete.obs")
  map_df(unique(treat$subclass), function(sc) {
    trow  <- treat %>% filter(subclass == sc)
    crows <- ctrl  %>% filter(subclass == sc)
    map_df(seq_len(nrow(crows)), function(i) {
      tibble(
        subclass   = sc,
        treated_OA = trow$OA,
        control_OA = crows$OA[i],
        mdist      = mahalanobis(
          x      = as.numeric(crows[i, vars]),
          center = as.numeric(trow[vars]),
          cov    = S
        )
      )
    })
  })
}

# =============================================================================
# WHY replace=TRUE AT STAGE 2
# =============================================================================
# Stage 1 uses replace=TRUE so the same control OA can be nominated as a
# structural neighbour by multiple treated OAs — this maximises the
# candidate pool and is fine because Stage 1 controls are not final matches.
#
# Stage 2 also uses replace=TRUE so every treated OA gets its globally
# best available injury-trajectory match, unconstrained by whether that
# control has already been used. Without replacement, later-processed
# treated OAs receive progressively worse matches as the pool is depleted.
#
# IMPORTANT CONSEQUENCE: with replace=TRUE a control OA can appear in
# multiple matched pairs. MatchIt handles this by assigning weights > 1
# to frequently-matched controls. These weights MUST be passed to the
# Callaway-Sant'Anna estimator — otherwise control observations are
# effectively double-counted and standard errors will be incorrect.
# The saved RDS files retain the weights column for this purpose.
#
# WHY RATIO 1, 2, 5 AT STAGE 2
# =============================================================================
# We test three ratios and select based on distance diagnostics and
# balance. Ratios 1, 2, 5 give a spread from tight 1:1 to moderately
# many controls. Ratio 10 (used at Stage 1) is too many for Stage 2:
# the 10th-nearest trajectory match is often very poor, and with
# replace=TRUE there is no pool-exhaustion pressure to justify accepting
# distant matches — you always have access to the best available control.
#
# Selection rule (same as before): choose the highest ratio where
#   (a) all trend SMDs remain < 0.10
#   (b) distance summary statistics (mean, median, P90) remain acceptable
#   (c) the pre-trend overlay plot (when available) shows close alignment
# Higher ratio = more controls = lower variance in C&S estimation,
# but at the cost of weaker average trajectory match quality.
# =============================================================================




run_one_s2 <- function(data, s2_vars, ratio, label_suffix) {
  formula <- reformulate(s2_vars, response = "treat_indicator")
  tryCatch(
    matchit(
      formula,
      data     = data,
      method   = "nearest",
      distance = "mahalanobis",
      ratio    = ratio,
      replace  = TRUE,          # 
      exact    = ~ country
    ),
    error = function(e) {
      cat("  FAILED (ratio=", ratio, "):", conditionMessage(e), "\n")
      NULL
    }
  )
}

summarise_mdist <- function(df, spec) {
  df %>% summarise(
    spec         = spec,
    n_pairs      = n(),
    mean_mdist   = round(mean(mdist),           3),
    median_mdist = round(median(mdist),          3),
    p90_mdist    = round(quantile(mdist, 0.90),  3),
    p95_mdist    = round(quantile(mdist, 0.95),  3),
    max_mdist    = round(max(mdist),             3)
  )
}

run_stage2_multi <- function(s2_prep, analysis_label) {
  
  cat("\n--- Stage 2 (multi-ratio):", analysis_label, "---\n")
  
  data    <- s2_prep$data
  s2_vars <- s2_prep$s2_present
  
  # ======================================================
  # INTERNAL: Compute Mahalanobis distances (Stage 1 logic)
  # ======================================================
  compute_mdist_s2 <- function(m, data, s2_vars) {
    
    md_raw <- match.data(m)
    
    treated_df  <- md_raw %>% dplyr::filter(treat_indicator == 1)
    controls_df <- md_raw %>% dplyr::filter(treat_indicator == 0)
    
    # Covariance matrix based on treated OAs
    S <- cov(treated_df[, s2_vars], use = "pairwise.complete.obs")
    
    # All pairwise distances (treated × controls)
    ddf <- purrr::map_df(seq_len(nrow(treated_df)), function(i) {
      trow <- treated_df[i, ]
      
      purrr::map_df(seq_len(nrow(controls_df)), function(j) {
        crow <- controls_df[j, ]
        
        tibble::tibble(
          treated_OA = trow$OA,
          control_OA = crow$OA,
          mdist = mahalanobis(
            x      = as.numeric(crow[, s2_vars]),
            center = as.numeric(trow[, s2_vars]),
            cov    = S
          )
        )
      })
    })
    
    return(ddf)
  }
  
  # ======================================================
  # Run ratios
  # ======================================================
  ratios <- c(1, 2, 5)
  
  s2_fits <- setNames(
    lapply(ratios, function(r) run_one_s2(data, s2_vars, r, analysis_label)),
    paste0("r", ratios)
  )
  
  # ======================================================
  # Compute Mahalanobis distances for each ratio
  # ======================================================
  mdist_list <- lapply(names(s2_fits), function(nm) {
    m <- s2_fits[[nm]]
    if (is.null(m)) return(NULL)
    compute_mdist_s2(m, data, s2_vars)
  })
  names(mdist_list) <- names(s2_fits)
  
  # ======================================================
  # Distance summary table
  # ======================================================
  summarise_mdist <- function(df, spec) {
    df %>%
      summarise(
        spec         = spec,
        n_pairs      = n(),
        mean_mdist   = round(mean(mdist), 3),
        median_mdist = round(median(mdist), 3),
        p90_mdist    = round(quantile(mdist, 0.90), 3),
        p95_mdist    = round(quantile(mdist, 0.95), 3),
        max_mdist    = round(max(mdist), 3)
      )
  }
  
  dist_summary <- bind_rows(
    Filter(Negate(is.null),
           mapply(summarise_mdist, mdist_list, names(mdist_list), SIMPLIFY = FALSE))
  )
  
  cat("\n  Distance summary across ratios:\n")
  print(dist_summary)
  
  # ======================================================
  # Balance summary across ratios
  # ======================================================
  cat("\n  Balance summary across ratios (trend vars — target all < 0.10):\n")
  
  for (nm in names(s2_fits)) {
    m <- s2_fits[[nm]]
    if (is.null(m)) next
    
    bt <- cobalt::bal.tab(m, thresholds = c(m = 0.1), un = FALSE)
    
    smd_all <- abs(bt$Balance$Diff.Adj)
    
    trend_smd <- bt$Balance %>%
      tibble::rownames_to_column("v") %>%
      filter(v %in% stage2_trends) %>%
      pull(Diff.Adj) %>%
      abs()
    
    cat(sprintf(
      "  %s | max trend SMD: %.3f | max all SMD: %.3f | mean SMD: %.3f\n",
      nm,
      max(trend_smd, na.rm = TRUE),
      max(smd_all,   na.rm = TRUE),
      mean(smd_all,  na.rm = TRUE)
    ))
  }
  
  # ======================================================
  # Histogram comparison across ratios
  # ======================================================
  dist_df_all <- bind_rows(
    Filter(Negate(is.null),
           mapply(function(df, nm) mutate(df, ratio = nm),
                  mdist_list, names(mdist_list), SIMPLIFY = FALSE))
  )
  
  p_ratios <- ggplot(dist_df_all, aes(x = mdist, fill = ratio)) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    facet_wrap(~ ratio, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title    = paste("Stage 2 distances by ratio —", analysis_label),
      subtitle = "Compare shape and tail across ratios",
      x = "Mahalanobis distance",
      y = "Count"
    )
  
  ggsave(
    here::here("output", paste0("s2_dist_ratios_", analysis_label, ".png")),
    p_ratios,
    width = 13, height = 5, dpi = 300
  )
  
  cat("  Ratio comparison histogram saved.\n")
  
  # ======================================================
  # SELECT PRIMARY RATIO
  # (You can update based on diagnostics)
  # ======================================================
  primary_ratio <- "r1"
  cat("\n  Primary ratio selected:", primary_ratio, "\n")
  
  m_primary <- s2_fits[[primary_ratio]]
  mdist_primary <- mdist_list[[primary_ratio]]
  
  # ======================================================
  # Trim using P95 threshold
  # ======================================================
  p95 <- quantile(mdist_primary$mdist, 0.95)
  
  keep_pairs <- mdist_primary %>% filter(mdist <= p95)
  
  primary_data <- match.data(m_primary, data = data)
  
  primary_data_trimmed <- primary_data %>%
    filter(OA %in% keep_pairs$treated_OA)
  
  n_before <- sum(primary_data$treat_indicator == 1)
  n_after  <- sum(primary_data_trimmed$treat_indicator == 1)
  
  trimmed_oas <- primary_data %>%
    filter(treat_indicator == 1, !OA %in% keep_pairs$treated_OA)
  
  cat("\n  P95 trim (primary ratio", primary_ratio, "):\n")
  cat("    Treated before trim: ", n_before, "\n")
  cat("    P95 threshold:       ", round(p95, 2), "\n")
  cat("    Treated trimmed:     ", n_before - n_after, "\n")
  cat("    Treated retained:    ", n_after, "\n")
  
  # Country check
  cat("    Country distribution of trimmed treated OAs:\n")
  print(table(trimmed_oas$country))
  
  # Histogram for primary ratio
  p_hist <- ggplot(mdist_primary, aes(x = mdist)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = p95, colour = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = p95 + 0.3, y = Inf, vjust = 2, hjust = 0, colour = "red",
             size = 3, label = sprintf("P95 = %.1f", p95)) +
    theme_minimal() +
    labs(
      title    = paste("Stage 2 distances (primary ratio", primary_ratio, ") —", analysis_label),
      subtitle = "Red = P95 trim threshold",
      x = "Mahalanobis distance",
      y = "Count"
    )
  
  ggsave(
    here::here("output", paste0("s2_dist_primary_", analysis_label, ".png")),
    p_hist, width = 10, height = 6, dpi = 300
  )
  
  cat("     histogram saved.\n\n")
  
  # ======================================================
  # Return all Stage 2 objects
  # ======================================================
  list(
    all_fits             = s2_fits,
    all_mdists           = mdist_list,
    dist_summary         = dist_summary,
    primary_ratio        = primary_ratio,
    matchit_obj          = m_primary,
    mdist                = mdist_primary,
    p95                  = p95,
    primary_data         = primary_data,
    primary_data_trimmed = primary_data_trimmed,
    n_pairs              = n_after,
    trimmed_oas          = trimmed_oas,
    data                 = data,
    s2_vars              = s2_vars
  )
}  

s2_A <- run_stage2_multi(s2_prep_A, "A_excl_zero")
s2_B <- run_stage2_multi(s2_prep_B, "B_incl_zero")

# ── RATIO SELECTION GUIDANCE ──────────────────────────────────────────────────
# After reviewing the output above, update primary_ratio inside
# run_stage2_multi() if a higher ratio (r2 or r5) passes the balance check.
#
# The key check is: do ALL seven trend SMDs remain below 0.10?
# If yes for r5 => use r5 (more controls, lower C&S variance)
# If yes for r2 but not r5 => use r2
# If only r1 passes => use r1
#
# NOTE: lower mean/median distance at r1 vs r5 is mathematically
# guaranteed (you're taking 1 vs 5 nearest neighbours). Do NOT use
# distance alone to justify r1 over r5 — use the balance diagnostics.
# =============================================================================

# =============================================================================
# SECTION 7 — BALANCE DIAGNOSTICS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 7: BALANCE DIAGNOSTICS\n")
cat("====================================================\n\n")

run_balance <- function(s2_result, label) {
  
  cat("--- Balance (primary ratio:", s2_result$primary_ratio, ") —", label, "---\n")
  bt     <- bal.tab(s2_result$matchit_obj,
                    thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    mutate(var_type = case_when(
      variable %in% stage2_trends ~ "Trend",
      variable %in% stage2_levels ~ "Level",
      TRUE                        ~ "Other"
    ))
  
  cat("\n  TREND SMDs (target: all < 0.06 for strong parallel trends support):\n")
  trend_s <- smd_df %>% filter(var_type == "Trend") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(trend_s)
  if (any(abs(trend_s$Diff.Adj) >= 0.1, na.rm = TRUE))
    cat("  WARNING: trend SMD >= 0.1\n")
  
  cat("\n  LEVEL SMDs (residual imbalance expected; include in C&S outcome model):\n")
  level_s <- smd_df %>% filter(var_type == "Level") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(level_s)
  hi <- level_s %>% filter(abs(Diff.Adj) > 0.1) %>% pull(variable)
  if (length(hi) > 0)
    cat("  -> Include as outcome model covariates:", paste(hi, collapse=", "), "\n")
  
  smd_all <- abs(bt$Balance$Diff.Adj)
  cat("\n  Balanced (<0.1):", sum(smd_all < 0.1, na.rm = TRUE),
      "/", length(smd_all),
      "| Max |SMD|:", round(max(smd_all, na.rm = TRUE), 3),
      "| Mean |SMD|:", round(mean(smd_all, na.rm = TRUE), 3), "\n")
  
  lp <- love.plot(
    s2_result$matchit_obj, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted",
    title        = paste("Stage 2 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) + theme(axis.text.y = element_text(size = 9),
            plot.margin = margin(10, 30, 10, 180),
            legend.position = "bottom")
  ggsave(here("output", paste0("s2_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  invisible(bt)
}

run_balance(s2_A, "A_excl_zero")
run_balance(s2_B, "B_incl_zero")

# =============================================================================
# SECTION 8 — COMPARE ANALYSES A AND B
# =============================================================================
# This is the core diagnostic comparison before running C&S estimation.
# We examine:
#   (i)   Sample overlap between A and B
#   (ii)  Structural characteristics of zero-injury OAs added in B
#   (iii) Stage 2 distance distribution for zero-injury pairs in B
#         (expect near-zero — they match on absence of injury, not
#          on trajectory similarity)
#   (iv)  Overall distance distribution comparison
#
# The definitive comparison (ATT estimates) is done in the C&S outcome
# script using both saved matched samples.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 8: COMPARISON — ANALYSES A vs B\n")
cat("====================================================\n\n")

treated_A <- s2_A$primary_data_trimmed %>%
  filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data_trimmed %>%
  filter(treat_indicator == 1) %>% pull(OA)

only_in_B <- setdiff(treated_B, treated_A)  # zero-injury OAs added in B
only_in_A <- setdiff(treated_A, treated_B)  # OAs in A trimmed out of B
in_both   <- intersect(treated_A, treated_B)

cat("--- Sample overlap ---\n")
cat("  Treated OAs in Analysis A:    ", length(treated_A), "\n")
cat("  Treated OAs in Analysis B:    ", length(treated_B), "\n")
cat("  In both:                      ", length(in_both), "\n")
cat("  Only in B (zero-injury):      ", length(only_in_B), "\n")
cat("  Only in A (trimmed from B):   ", length(only_in_A), "\n\n")

# Characterise zero-injury OAs added in B
if (length(only_in_B) > 0) {
  
  cat("--- Structural characteristics of zero-injury OAs in B ---\n")
  cat("  (These OAs are matched only on Stage 1 structural vars;\n")
  cat("   parallel trends assumption for them rests on structural\n")
  cat("   similarity, not trajectory similarity.)\n\n")
  
  compare_grp <- s2_B$primary_data_trimmed %>%
    filter(treat_indicator == 1) %>%
    mutate(grp = if_else(OA %in% only_in_B,
                         "Zero-injury (B only)",
                         "Injury-exposed (both)"))
  
  s1_avail <- intersect(stage1_vars, names(compare_grp))
  cov_diff <- compare_grp %>%
    select(grp, all_of(s1_avail)) %>%
    pivot_longer(-grp) %>%
    group_by(name, grp) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value,   na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = c(mean, sd)) %>%
    mutate(
      SMD = (`mean_Injury-exposed (both)` - `mean_Zero-injury (B only)`) /
        sqrt((`sd_Injury-exposed (both)`^2 + `sd_Zero-injury (B only)`^2) / 2)
    ) %>%
    arrange(desc(abs(SMD)))
  
  cat("  Structural SMDs: injury-exposed vs zero-injury treated OAs in B\n")
  cat("  (Large |SMD| = structurally different subgroups)\n\n")
  print(cov_diff %>%
          select(name,
                 `mean_Injury-exposed (both)`,
                 `mean_Zero-injury (B only)`,
                 SMD), n = Inf)
  
  cat("\n  Country distribution of zero-injury OAs in B:\n")
  zero_oa_data <- s2_B$primary_data_trimmed %>%
    filter(OA %in% only_in_B)
  print(table(zero_oa_data$country))
  
  # Stage 2 distances for zero-injury pairs specifically
  zero_dists <- s2_B$mdist %>%
    filter(treated_OA %in% only_in_B)
  cat("\n  Stage 2 distance distribution for zero-injury pairs in B:\n")
  cat("  (Should be near-zero: matching zero to zero on all 14 vars)\n")
  print(summary(zero_dists$mdist))
  
  if (mean(zero_dists$mdist) > 1) {
    cat("\n  NOTE: Some zero-injury OAs show non-trivial Stage 2 distances.\n")
    cat("  Check whether their Stage 2 variables are truly all zero.\n")
  }
}

# Distance distribution comparison
cat("\n--- Stage 2 distance distribution: A vs B ---\n")
dist_compare <- bind_rows(
  s2_A$mdist %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$mdist %>% mutate(analysis = "B: incl zero-injury")
) %>%
  group_by(analysis) %>%
  summarise(
    n_pairs = n(),
    mean    = round(mean(mdist),           3),
    median  = round(median(mdist),         3),
    p90     = round(quantile(mdist, 0.90), 3),
    p95     = round(quantile(mdist, 0.95), 3),
    max     = round(max(mdist),            3),
    .groups = "drop"
  )
print(dist_compare)

cat("\nNOTE: Analysis B will likely show lower mean/median distances\n")
cat("because zero-injury pairs have near-zero Stage 2 distances by\n")
cat("construction. Lower distance does NOT mean better matches here —\n")
cat("it means those pairs carry no Stage 2 trajectory information.\n\n")

# Side-by-side histogram
dist_both <- bind_rows(
  s2_A$mdist %>% mutate(analysis = "A: Excluding zero-injury OAs"),
  s2_B$mdist %>% mutate(analysis = "B: Including zero-injury OAs")
)
p_compare <- ggplot(dist_both, aes(x = mdist, fill = analysis)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~ analysis, scales = "free_y") +
  scale_fill_manual(values = c("steelblue", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title    = "Stage 2 Mahalanobis distances: Analysis A vs B",
       subtitle = "Analysis B will show spike near zero from zero-injury pairs",
       x = "Mahalanobis distance", y = "Count")
ggsave(here("output", "s2_dist_comparison_A_vs_B.png"),
       p_compare, width = 12, height = 5, dpi = 300)
cat("Comparison histogram saved.\n")

# =============================================================================
# SECTION 9 — EXTRACT AND SAVE MATCHED SAMPLES
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 9: SAVE MATCHED SAMPLES\n")
cat("====================================================\n\n")

extract_matched <- function(s2_result, label) {
  treated  <- s2_result$primary_data_trimmed %>%
    filter(treat_indicator == 1) %>%
    select(OA, weights, subclass)
  controls <- s2_result$primary_data_trimmed %>%
    filter(treat_indicator == 0) %>%
    select(OA, weights, subclass)
  cat(label, "— Treated:", nrow(treated),
      "| Controls:", nrow(controls), "\n")
  # Final cross-country check
  cross <- s2_result$primary_data_trimmed %>%
    group_by(subclass) %>%
    summarise(n_c = n_distinct(country), .groups = "drop") %>%
    filter(n_c > 1)
  cat("  Cross-country pairs:", nrow(cross),
      if (nrow(cross) == 0) "✓" else "WARNING", "\n")
  list(treated = treated, controls = controls)
}

matched_A <- extract_matched(s2_A, "Analysis A (excl zero-injury)")
matched_B <- extract_matched(s2_B, "Analysis B (incl zero-injury)")

saveRDS(matched_A$treated,
        here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A$controls,
        here("data", "processed", "OA_matched_donors_A.rds"))
saveRDS(matched_B$treated,
        here("data", "processed", "OA_matched_treated_B.rds"))
saveRDS(matched_B$controls,
        here("data", "processed", "OA_matched_donors_B.rds"))

cat("\nSaved:\n")
cat("  OA_matched_treated_A.rds / OA_matched_donors_A.rds\n")
cat("  OA_matched_treated_B.rds / OA_matched_donors_B.rds\n")

# =============================================================================
# SECTION 10 — DESIGN SUMMARY
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 10: DESIGN SUMMARY\n")
cat("====================================================\n\n")

cat("STAGE 1:  MDM on", length(stage1_vars),
    "structural + sociodemographic variables\n")
cat("          No caliper — all variables treated symmetrically\n")
cat("          replace=TRUE, ratio=10\n")
cat("          Post-hoc P95 distance trim\n\n")

cat("STAGE 2:  MDM on", length(stage2_vars),
    "injury dynamics variables (7 trends + 7 levels)\n")
cat("          No caliper\n")
cat("          replace=TRUE (globally best matches; weights passed to C&S)\n")
cat("          Ratios tested: 1, 2, 5 — select on trend SMD balance\n")
cat("          Post-hoc P95 distance trim on selected ratio\n\n")

cat("Country exact constraint at both stages.\n\n")

cat("ANALYSIS A  |  Zero-injury treated OAs EXCLUDED\n")
cat("            |  Primary ratio:", s2_A$primary_ratio, "\n")
cat("            |  N treated:", nrow(matched_A$treated), "\n")
cat("            |  ATT scope: road-exposed OAs with pre-period injuries\n\n")

cat("ANALYSIS B  |  Zero-injury treated OAs INCLUDED\n")
cat("            |  Primary ratio:", s2_B$primary_ratio, "\n")
cat("            |  N treated:", nrow(matched_B$treated), "\n")
cat("            |  ATT scope: all road-exposed OAs in scheme boundaries\n\n")

cat("WEIGHTS NOTE: replace=TRUE at Stage 2 means some control OAs are\n")
cat("matched to multiple treated OAs. MatchIt assigns weights > 1 to\n")
cat("these controls. Pass the 'weights' column from the saved RDS files\n")
cat("to the C&S estimator to avoid double-counting.\n\n")

cat("NEXT STEP: Run C&S estimation on both matched samples.\n")
cat("  If ATT(A) ≈ ATT(B)  =>  B is primary (no selection bias),\n")
cat("                           A is robustness check\n")
cat("  If ATT(A) ≠ ATT(B)  =>  report both; discuss what zero-injury\n")
cat("                           OAs contribute to the treatment effect\n")