# =============================================================================
# OA-LEVEL TWO-STAGE MAHALANOBIS DISTANCE MATCHING
# =============================================================================
#
# PURPOSE:
#   Construct matched comparison groups for a Difference-in-Differences (DiD)
#   #
# DESIGN LOGIC — TWO-STAGE MDM:
#   The two-stage approach separates two distinct matching goals:
#
#   STAGE 1 — Structural plausibility (15 variables)
#     Ensures controls are structurally credible substitutes for treated OAs
#     in terms of road network, urban form, and sociodemographic composition.
#     ratio=10 deliberately casts a wide net — this is not the final matched
#     sample, just a restricted candidate pool for Stage 2.
#     replace=TRUE: every treated OA draws from the full control pool.
#     exact=~country: Accounting for  different contexts of England and Scotland 
#
#   STAGE 2 — Injury-dynamics matching (14 variables: 7 trends + 7 levels)
#     Within the Stage 1 pool, selects controls with the most similar
#     pre-treatment injury trajectories. Trend variables are the primary
#     target because parallel trends is the identifying assumption in C&S DiD.
#     #  replace=TRUE: globally best trajectory match for every treated OA.
#     ratio= 1, 2, 5 are examoined for the best trend SMD balance
#
# TWO  ANALYSES:
#   Analysis A — EXCLUDING zero-injury treated OAs 
#     All Stage 2 variables are well-defined and observed.
#     Parallel trends rests on empirical pre-period injury dynamics.
#     ATT scope: road-exposed OAs with observed injury history.
#     Limitation: excludes 191 treated OAs — scope bias if zero-injury
#     areas respond differently to the intervension
#
#   Analysis B — INCLUDING zero-injury treated OAs 
#     Stage 2 variables set to 0 for zero-injury treated OAs.
#     Broader ATT scope but weaker identification for the zero-injury subset:
#     (a) no empirical parallel trends evidence for zero-injury OAs;
#     (b) Scottish zero-injury OAs have near-degenerate matching 
#
# FIXED EXCLUSIONS (both analyses):
#   Buffer OAs: contamination risk — may have been indirectly affected.
#   Zero-road OAs: injury rates per km undefined; 74 have recorded injuries
#                  suggesting data linkage errors — excluded for safety.
#
# OUTPUTS (Section 9):
#   OA_matched_treated_A.rds  — treated OA IDs + weights, Analysis A
#   OA_matched_donors_A.rds   — control OA IDs + weights, Analysis A
#   OA_matched_treated_B.rds  — treated OA IDs + weights, Analysis B (capped)
#   OA_matched_donors_B.rds   — control OA IDs + weights, Analysis B (capped)
#   OA_matched_full_A.rds     — full matched ready dataset, Analysis A
#   OA_matched_full_B.rds     — full matched ready dataset, Analysis B
#
# 
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
library(ggrepel)   


OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

print(table(OA_matching_dataset$assignment))


cat("\n--- Zero-injury OA counts ---\n")
print(table(OA_matching_dataset$zero_injury_OA))
# 0 = has observed injuries in pre-period; 1 = no injuries recorded

cat("\n--- Zero-injury breakdown by treatment status ---\n")
print(table(
  OA_matching_dataset$zero_injury_OA,
  OA_matching_dataset$treated_OA,
  dnn = c("zero_injury", "treated")
))
cat("\n--- Road network summary ---\n")
OA_matching_dataset %>%
  summarise(
    total_OAs  = n_distinct(OA),
    zero_roads = sum(n_roads == 0 | is.na(n_roads)),
    pct_zero   = round(100 * zero_roads / total_OAs, 2)
  ) %>% print()

# Flag zero-road OAs that somehow have recorded injuries — data anomalies
weirdOAs <- OA_matching_dataset %>%
  filter((n_roads == 0 | is.na(n_roads)) & mean_total > 0)
cat("Zero-road OAs with recorded injuries:", nrow(weirdOAs),
    " (treated:", sum(weirdOAs$treated_OA == 1), ")\n")
# These are excluded in Section 2 — injury rates per km are undefined
# without road length, and the injury records suggest data linkage errors.

# =============================================================================
# VARIABLE DEFINITIONS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 1: VARIABLE DEFINITIONS\n")
cat("====================================================\n\n")

# --- Stage 1: Structural + sociodemographic (15 variables) ---
# Road variables: capture exposure and network type
stage1_road <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road"
)
# Urban form: distance to city centre and population density
# Both affect baseline injury risk and travel mode mix
stage1_urban <- c("dist_citycentre", "pop_density")

# Sociodemographic: affect travel behaviour and injury exposure
# IMD = deprivation; car ownership; mode share; age groups
stage1_socdem <- c(
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_socdem)

# --- Stage 2: Pre-treatment injury dynamics (14 variables) ---
# TRENDS (7): slopes of log-injury rates per road-km over pre-period
# These are the primary matching target — parallel trends in C&S DiD
# requires treated and controls to have moved together before treatment.
# Matching on trends directly maximises parallel trends plausibility.
stage2_trends <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm",
  "trend_cyc_KSI_pkm", "trend_cyc_slight_pkm",
  "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_total_pkm"
)
# LEVELS (7): mean injury rates per road-km over pre-period
# Secondary target — residual level imbalance is expected and acceptable;
# remaining imbalance enters the C&S outcome model as covariates.
stage2_levels <- c(
  "mean_car_KSI_pkm", "mean_car_slight_pkm",
  "mean_cyc_KSI_pkm", "mean_cyc_slight_pkm",
  "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
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
cat("No calipers at either stage — all variables treated symmetrically.\n")
# No caliper rationale: calipers on individual variables impose implicit
# relative weights. MDM already accounts for variable scale and correlation
# via the covariance matrix. Dropping units via caliper risks exclusion of
# valid matches and reduces ATT representativeness.

# =============================================================================
# SECTION 2 — BUILD BOTH ANALYSIS DATASETS
# =============================================================================

# Derive country from LAD code prefix (E = England, S = Scotland)
# Used for exact matching constraint in Stage 1
OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  )

# --- Shared exclusions applied to both analyses ---
# buffer_OA == 0: buffer OAs excluded — risk of indirect contamination
#  n_roads > 0: zero-road OAs excluded — injury rates per km undefined.
#   74 zero-road OAs have recorded injuries
base_filter <- OA_matching_dataset %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA == 0,
    n_roads   >  0
  ) %>%
  mutate(treat_indicator = as.integer(treated_OA == 1))

# --- Analysis A: exclude zero-injury treated OAs ---
# Rationale: Stage 2 matching on injury trends requires observed pre-period
# injury dynamics. Zero-injury OAs have no injury history to match on,
# so parallel trends cannot be empirically supported for this subgroup.
# ATT scope is narrowed to injury-exposed OAs — a scope limitation,

data_A <- base_filter %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1))

# --- Analysis B: include zero-injury OAs ---
# Rationale: provides a broader ATT estimate covering all treated OAs.
# Zero-injury treated OAs have their Stage 2 variables set to 0 so they
# match to zero-injury controls with near-zero distances.
# Parallel trends for this subgroup rests on Stage 1 structural similarity

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
cat("    of which zero-injury:",
    sum(data_B$treat_indicator == 1 & data_B$zero_injury_OA == 1), "\n")
cat("  Controls:    ", sum(data_B$treat_indicator == 0), "\n\n")

cat("Country distribution — Analysis A:\n")
print(data_A %>% count(country, treat_indicator) %>%
        pivot_wider(names_from = treat_indicator, values_from = n,
                    names_prefix = "t") %>% rename(control = t0, treated = t1))

cat("\nCountry distribution — Analysis B:\n")
print(data_B %>% count(country, treat_indicator) %>%
        pivot_wider(names_from = treat_indicator, values_from = n,
                    names_prefix = "t") %>% rename(control = t0, treated = t1))

# =============================================================================
# SECTION 3 — WINSORISE STAGE 1 VARIABLES (pre matching) to remove outlier OAs
# =============================================================================
# Winsorising caps extreme values at 1st/99th percentile computed from the
# FULL dataset (treated + controls) — no rows are removed.
#
# Rationale: Mahalanobis distance uses the covariance matrix of matching
# variables. A single extreme outlier can inflate a diagonal element and
# effectively zero-out that variable's contribution to the distance metric,
# biasing match quality for all OAs. Winsorising prevents this without
# dropping observations or imposing ad hoc transformations.
#
# Stage 2 variables are winsorised separately in Section 5, anchored to
# the treated distribution — so controls are winsorised to the treated
# range, not the other way around.
# =============================================================================

### how extreme are the values compared tp the 99th %

map_df(stage1_vars, function(v){
  tibble(
    variable = v,
    p99 = quantile(data_A[[v]],0.99,na.rm=TRUE),
    max = max(data_A[[v]],na.rm=TRUE)
  )
})

##  there are extreme outliers which can distrot the matching especially road and pop density 
###  Without winsorisation, these OAs would inflate the covariance matrix dramatically.
### these are just 1% 
cat("\n====================================================\n")
cat(" WINSORISING STAGE 1 VARIABLES\n")
cat("====================================================\n\n")

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

# Check for near-zero variance variables — these would destabilise the
# covariance matrix inversion in Mahalanobis distance computation.
# Any such variables are dropped from the matching formula.
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

s1_vars_A <- check_vars(data_A_clean, stage1_vars, "Stage 1 / A")
s1_vars_B <- check_vars(data_B_clean, stage1_vars, "Stage 1 / B")
s2_vars_A <- check_vars(data_A_clean, stage2_vars, "Stage 2 / A")
s2_vars_B <- check_vars(data_B_clean, stage2_vars, "Stage 2 / B")

# =============================================================================
# — STAGE 1 MATCHING
# =============================================================================
# nearest-neighbour MDM on structural + sociodemographic variables.
# ratio=10: large pool passed to Stage 2 — not the final sample.
# replace=TRUE: every treated OA draws from the full control pool,
#   preventing early-matched treated OAs from exhausting good controls.
# exact=~country: enforces England/Scotland separation 
# No caliper: all 15 variables contribute equally via covariance weighting.
# =============================================================================

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
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  mm <- m$match.matrix
  cat("  Treated in:", nrow(mm), "\n")
  
  # Compute Mahalanobis distances manually for diagnostics
  # (MatchIt does not store distances when distance="mahalanobis")
  # Covariance matrix computed from treated rows only to avoid
  # control-pool variance dominating the metric
  S_s1 <- cov(data[as.integer(rownames(mm)), s1_vars],
              use = "pairwise.complete.obs")
  
  dist_s1 <- map_df(seq_len(nrow(mm)), function(i) {
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    c_indices  <- mm[i, ]
    c_indices  <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    map_df(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      tibble(
        treated_OA = treated_id,
        control_OA = crow[["OA"]],
        mdist      = mahalanobis(
          x      = as.numeric(crow[s1_vars]),
          center = as.numeric(trow[s1_vars]),
          cov    = S_s1
        )
      )
    })
  })
  
  # Extract matched treated and unique matched controls
  treated_matched     <- data[as.integer(rownames(mm)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  control_row_indices <- unique(as.integer(mm[!is.na(mm)]))
  controls_matched    <- data[control_row_indices, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated retained:", nrow(treated_matched), "\n")
  cat("  Unique controls in pool:", nrow(controls_matched), "\n")
  
  # Balance diagnostics — target: all SMDs < 0.1 after matching
  cat("\n  Stage 1 balance:\n")
  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  print(smd_df)
  
  cat("\n  Country distribution (treated):\n")
  print(table(treated_matched$country))
  
  # Verify exact constraint was respected — no cross-country pairs
  cross_check <- dist_s1 %>%
    left_join(data %>% select(OA, country) %>% rename(control_country = country),
              by = c("control_OA" = "OA")) %>%
    left_join(data %>% select(OA, country) %>% rename(treated_country = country),
              by = c("treated_OA" = "OA")) %>%
    filter(treated_country != control_country)
  cat("  Cross-country pairs:", nrow(cross_check),
      if (nrow(cross_check) == 0) "✓" else "WARNING: cross-country pairs found!", "\n")
  
  lp <- love.plot(
    m, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted",
    title        = paste("Stage 1 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) + theme(plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  
  ggsave(here("output", paste0("s1_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  list(
    matchit_obj = m,
    dist_s1     = dist_s1,
    treated     = treated_matched,
    controls    = controls_matched,
    bal         = bt
  )
}

s1_A <- run_stage1(data_A_clean, s1_vars_A, "A_excl_zero")
s1_B <- run_stage1(data_B_clean, s1_vars_B, "B_incl_zero")

# Print full balance tables for reporting
bal.tab(s1_A$matchit_obj, un = TRUE)
bal.tab(s1_B$matchit_obj, un = TRUE)



# Stage 1 matching restricts the control pool to OAs that are structurally comparable to treated OAs 
# The aim is not covariate balance but removal of extreme non-comparable areas prior to outcome-based matching in stage 2


# =============================================================================
# — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE
# =============================================================================
# Stage 2 data = Stage 1 treated + Stage 1 unique controls.
# Stage 2 variables are winsorised anchored to the TREATED distribution:
#   quantiles computed from treated OAs only, then applied to controls.
# Rationale: prevents extreme control values from pulling winsorisation
# thresholds away from the treated range. Controls are winsorised to the
# support of the treated distribution.
# =============================================================================


prepare_s2_data <- function(s1_result, s2_vars) {
  
  s2_raw <- bind_rows(
    s1_result$treated,
    s1_result$controls
  ) %>%
    select(-any_of(c("weights", "subclass", "distance")))
  # Remove Stage 1 matching artefacts — Stage 2 will create its own
  
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  s2_raw %>%
    mutate(across(
      all_of(intersect(s2_vars, names(.))),
      ~ {
        # Winsorise bounds derived from treated distribution only
        q <- quantile(treated_ref[[cur_column()]], c(0.01, 0.99), na.rm = TRUE)
        pmin(pmax(., q[1]), q[2])
      }
    ))
}

s2_data_A <- prepare_s2_data(s1_A, s2_vars_A)
s2_data_B <- prepare_s2_data(s1_B, s2_vars_B)

# =============================================================================
#  — STAGE 2 MATCHING
# =============================================================================
# MDM on pre-treatment injury dynamics within the Stage 1 candidate pool.
# replace=TRUE: globally optimal trajectory match for each treated OA.
# 
#
# Mahalanobis distances computed manually from match matrix because
# MatchIt does not store distances when distance="mahalanobis".
# Each treated OA is assigned the MEAN distance across its matched controls.
# Each control OA is assigned the MEAN distance across all treated OAs
#
# WEIGHTS (replace=TRUE):
#   Controls matched to multiple treated OAs receive weight > 1.
#   These weights must be passed to att_gt(weightsname = "weights").
#   # =============================================================================

run_stage2 <- function(data, s2_vars, ratio, label) {
  
  cat("\n--- Stage 2:", label, "ratio", ratio, "---\n")
  
  formula <- reformulate(s2_vars, response = "treat_indicator")
  
  m <- matchit(
    formula,
    data     = data,
    method   = "nearest",
    distance = "mahalanobis",
    ratio    = ratio,
    replace  = TRUE
    # No exact constraint at Stage 2 — country separation already enforced
    # at Stage 1; Stage 2 pool is already country-homogeneous per treated OA
  )
  
  mm <- m$match.matrix
  
  # Covariance matrix from treated rows — same rationale as Stage 1
  S_s2 <- cov(data[as.integer(rownames(mm)), s2_vars],
              use = "pairwise.complete.obs")
  
  # Compute pairwise Mahalanobis distances from match matrix
  dist_s2 <- map_df(seq_len(nrow(mm)), function(i) {
    
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    
    c_indices  <- mm[i, ]
    c_indices  <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    
    dists <- map_dbl(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      mahalanobis(
        x      = as.numeric(crow[s2_vars]),
        center = as.numeric(trow[s2_vars]),
        cov    = S_s2
      )
    })
    
    tibble(
      OA            = treated_id,
      control_OA    = data[as.integer(c_indices), "OA", drop = FALSE] %>% pull(OA),
      mdist_control = dists,
      mdist_treated = mean(dists)  # treated OA summary: mean dist to its matches
    )
  })
  
  # Treated OA distance = mean across its matched controls
  treated_dists <- dist_s2 %>%
    distinct(OA, mdist_treated) %>%
    rename(mdist = mdist_treated)
  
  # Control OA distance = mean across all treated OAs it was matched to
  control_dists <- dist_s2 %>%
    group_by(OA = control_OA) %>%
    summarise(mdist = mean(mdist_control), .groups = "drop")
  
  all_dists <- bind_rows(treated_dists, control_dists) %>%
    group_by(OA) %>%
    summarise(mdist = mean(mdist), .groups = "drop")
  
  matched_data <- match.data(m) %>%
    left_join(all_dists, by = "OA")
  
  cat("Treated matched:  ", sum(matched_data$treat_indicator == 1), "\n")
  cat("Controls matched: ", sum(matched_data$treat_indicator == 0), "\n")
  cat("NAs in mdist:     ", sum(is.na(matched_data$mdist)), "\n")
  # NAs should be 0 — if not, distance computation has a bug
  cat("Distance summary (treated only):\n")
  print(summary(matched_data$mdist[matched_data$treat_indicator == 1]))
  
  list(
    matchit_obj   = m,
    primary_ratio = ratio,
    primary_data  = matched_data
  )
}

# Run all ratio variants — primary selected on balance in Section 7
s2_A_r1 <- run_stage2(s2_data_A, s2_vars_A, 1, "A")
s2_A_r2 <- run_stage2(s2_data_A, s2_vars_A, 2, "A")
s2_A_r5 <- run_stage2(s2_data_A, s2_vars_A, 5, "A")

s2_B_r1 <- run_stage2(s2_data_B, s2_vars_B, 1, "B")
s2_B_r2 <- run_stage2(s2_data_B, s2_vars_B, 2, "B")
s2_B_r5 <- run_stage2(s2_data_B, s2_vars_B, 5, "B")

# =============================================================================
#  BALANCE DIAGNOSTICS
# =============================================================================


run_balance <- function(s2_result, label) {
  
  cat("--- Balance (ratio:", s2_result$primary_ratio, ") —", label, "---\n")
  bt     <- bal.tab(s2_result$matchit_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    mutate(var_type = case_when(
      variable %in% stage2_trends ~ "Trend",
      variable %in% stage2_levels ~ "Level",
      TRUE                        ~ "Other"
    ))
  
  cat("\n  TREND SMDs (target: all < 0.06):\n")
  trend_s <- smd_df %>% filter(var_type == "Trend") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(trend_s)
  if (any(abs(trend_s$Diff.Adj) >= 0.1, na.rm = TRUE))
    cat("  WARNING: trend SMD >= 0.1 — parallel trends support weakened\n")
  
  cat("\n  LEVEL SMDs (residual imbalance expected — include in C&S outcome model):\n")
  level_s <- smd_df %>% filter(var_type == "Level") %>%
    select(variable, Diff.Adj) %>% arrange(desc(abs(Diff.Adj)))
  print(level_s)
  hi <- level_s %>% filter(abs(Diff.Adj) > 0.1) %>% pull(variable)
  if (length(hi) > 0)
    cat("  -> Include as C&S outcome model covariates:", paste(hi, collapse = ", "), "\n")
  
  smd_all <- abs(bt$Balance$Diff.Adj)
  cat("\n  Balanced (<0.1):", sum(smd_all < 0.1, na.rm = TRUE), "/", length(smd_all),
      "| Max |SMD|:", round(max(smd_all, na.rm = TRUE), 3),
      "| Mean |SMD|:", round(mean(smd_all, na.rm = TRUE), 3), "\n")
  
  lp <- love.plot(
    s2_result$matchit_obj, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted", title = paste("Stage 2 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) + theme(axis.text.y = element_text(size = 9),
            plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  ggsave(here("output", paste0("s2_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  invisible(bt)
}

# ratio=5 is the primary specification 
##   — Best trend SMD balance across ratios 1, 2, 5
#   — Level balance also improves with higher ratio
#   — Efficiency gain is real and does not come at the cost of worse balance

s2_A <- s2_A_r5
s2_B <- s2_B_r5

run_balance(s2_A, "A_excl_zero")
run_balance(s2_B, "B_incl_zero")

# =============================================================================
#  — COMPARE ANALYSES A AND B
# =============================================================================
# Documents: (a) sample overlap; (b) structural differences of zero-injury
# OAs added in B; (c) distance distributions across analyses.
#
# If ATT(A) ≈ ATT(B) after DiD: A is confirmed as primary; B as robustness.
# If ATT(A) ≠ ATT(B): heterogeneous treatment effects by injury exposure —
#   report both and discuss what the difference reveals substantively.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 8: COMPARISON — ANALYSES A vs B\n")
cat("====================================================\n\n")

treated_A <- s2_A$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)

only_in_B <- setdiff(treated_B, treated_A)   # zero-injury OAs unique to B
only_in_A <- setdiff(treated_A, treated_B)   # should be 0 by construction
in_both   <- intersect(treated_A, treated_B)

cat("--- Sample overlap ---\n")
cat("  Treated in A:              ", length(treated_A), "\n")
cat("  Treated in B:              ", length(treated_B), "\n")
cat("  In both:                   ", length(in_both), "\n")
cat("  Only in B (zero-injury):   ", length(only_in_B), "\n")
cat("  Only in A (trimmed from B):", length(only_in_A), "\n\n")

if (length(only_in_A) > 0)
  cat("WARNING: OAs in A but not B — unexpected, check exclusion logic\n")

if (length(only_in_B) > 0) {
  
  cat("--- Structural characteristics of zero-injury OAs added in B ---\n")
  # Documents WHY Analysis B is weaker — zero-injury OAs are structurally
  # different from injury-exposed OAs (shorter roads, denser, minor roads only)
  compare_grp <- s2_B$primary_data %>%
    filter(treat_indicator == 1) %>%
    mutate(grp = if_else(OA %in% only_in_B,
                         "Zero-injury (B only)", "Injury-exposed (both)"))
  
  s1_avail <- intersect(stage1_vars, names(compare_grp))
  cov_diff <- compare_grp %>%
    select(grp, all_of(s1_avail)) %>% pivot_longer(-grp) %>%
    group_by(name, grp) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value,   na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = grp, values_from = c(mean, sd)) %>%
    mutate(SMD = (`mean_Injury-exposed (both)` - `mean_Zero-injury (B only)`) /
             sqrt((`sd_Injury-exposed (both)`^2 + `sd_Zero-injury (B only)`^2) / 2)) %>%
    arrange(desc(abs(SMD)))
  
  print(cov_diff %>%
          select(name, `mean_Injury-exposed (both)`, `mean_Zero-injury (B only)`, SMD),
        n = Inf)
  
  cat("\n  Country distribution of zero-injury OAs in B:\n")
  print(table(s2_B$primary_data %>% filter(OA %in% only_in_B) %>% pull(country)))
  # Expected: heavily Scottish — drives the weight collapse in B
  
  cat("\n  Stage 2 distances for zero-injury treated OAs:\n")
  print(summary(s2_B$primary_data %>%
                  filter(OA %in% only_in_B, treat_indicator == 1) %>%
                  pull(mdist)))
  # Expected: near-zero (they match to other zero-injury controls trivially)
}

# Distance distribution comparison — documents match quality across analyses
cat("\n--- Stage 2 distance distribution: A vs B (treated OAs only) ---\n")
dist_compare <- bind_rows(
  s2_A$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$primary_data %>% filter(treat_indicator == 1) %>%
    select(OA, mdist) %>% mutate(analysis = "B: incl zero-injury")
) %>%
  group_by(analysis) %>%
  summarise(
    n      = n(),
    mean   = round(mean(mdist,           na.rm = TRUE), 3),
    median = round(median(mdist,         na.rm = TRUE), 3),
    p90    = round(quantile(mdist, 0.90, na.rm = TRUE), 3),
    p95    = round(quantile(mdist, 0.95, na.rm = TRUE), 3),
    max    = round(max(mdist,            na.rm = TRUE), 3),
    .groups = "drop"
  )
print(dist_compare)

# Identify worst-matched treated OAs in A for reporting / sensitivity check
cat("\n--- Worst-matched treated OAs in Analysis A (top 20 by distance) ---\n")
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  arrange(desc(mdist)) %>%
  select(OA, mdist, country, starts_with("trend_")) %>%
  head(20) %>%
  print()

# Distance tail summary — informs caliper sensitivity decision
cat("\n--- Distance tail: Analysis A treated OAs ---\n")
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  summarise(
    n_over_20  = sum(mdist > 20),
    n_over_30  = sum(mdist > 30),
    n_over_50  = sum(mdist > 50),
    n_over_100 = sum(mdist > 100)
  ) %>% print()
# Decision: keep full sample as primary (tail is small: ~11 OAs > 30).
# A caliper-trimmed sensitivity analysis (drop mdist > 30) is available
# but not applied to primary — trimming 1.4% of treated sample risks
# appearing to select results rather than improve identification.

# =============================================================================
# SECTION 9 — WEIGHT DIAGNOSTICS + APPLY CAP TO ANALYSIS B
# =============================================================================
# replace=TRUE produces weights > 1 for controls matched to multiple treated OAs.
# These weights must be passed to att_gt(weightsname = "weights").
#
# Analysis A: max weight = 12, efficiency = 0.57 — acceptable, no cap needed.
#
# Analysis B: max weight = 93 (two Scottish OAs), efficiency = 0.048 — CRITICAL.
#   Root cause: 114 Scottish zero-injury treated OAs matched to a tiny pool
#   of Scottish zero-injury controls. The same 2-5 controls are recycled ~93
#   times each. The estimate for Scottish zero-injury OAs would rest on
#   essentially 2 data points.
#   Weight cap of 5 applied: raises efficiency to 0.672; trend SMDs unchanged
#   (0.038 with cap vs 0.045 without). Cap is applied directly to the
#   primary_data weights column before saving — this is what gets passed to DiD.
#   Note: capping reduces but does not eliminate the Scotland sparsity problem.
#   This is a fundamental data limitation for zero-injury Scottish OAs and
#   is a key reason Analysis B is the robustness check, not the primary.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 9: WEIGHT DIAGNOSTICS\n")
cat("====================================================\n\n")

weight_diagnostics <- function(s2_result, analysis_label) {
  
  controls <- s2_result$primary_data %>% filter(treat_indicator == 0)
  
  cat("\n--- Weight diagnostics:", analysis_label, "---\n")
  
  cat("\nWeight distribution (controls):\n")
  print(quantile(controls$weights,
                 probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)))
  
  cat("\nEffective N:", round(sum(controls$weights)^2 / sum(controls$weights^2), 1),
      "vs nominal N:", nrow(controls), "\n")
  cat("Efficiency ratio:",
      round((sum(controls$weights)^2 / sum(controls$weights^2)) / nrow(controls), 3), "\n")
  
  cat("\nBy country:\n")
  controls %>%
    group_by(country) %>%
    summarise(
      n      = n(),
      mean_w = round(mean(weights), 2),
      max_w  = round(max(weights), 2),
      eff_n  = round(sum(weights)^2 / sum(weights^2), 1),
      .groups = "drop"
    ) %>% print()
  
  # Cap sensitivity: show effect of different caps on efficiency and balance
  caps <- c(5, 10, 20, Inf)
  treated_rows <- s2_result$primary_data %>% filter(treat_indicator == 1)
  
  cap_results <- map_df(caps, function(cap) {
    capped    <- controls %>% mutate(w_cap = pmin(weights, cap))
    eff_n_cap <- sum(capped$w_cap)^2 / sum(capped$w_cap^2)
    trend_smds <- map_dbl(stage2_trends, function(v) {
      t_mean    <- mean(treated_rows[[v]], na.rm = TRUE)
      c_mean    <- weighted.mean(controls[[v]], w = capped$w_cap, na.rm = TRUE)
      pooled_sd <- sqrt((var(treated_rows[[v]], na.rm = TRUE) +
                           var(controls[[v]],   na.rm = TRUE)) / 2)
      if (pooled_sd == 0) return(0)
      abs(t_mean - c_mean) / pooled_sd
    })
    tibble(
      weight_cap    = if_else(is.infinite(cap), "No cap", paste0("Cap ", cap)),
      n_controls    = nrow(controls),
      effective_n   = round(eff_n_cap, 1),
      efficiency    = round(eff_n_cap / nrow(controls), 3),
      max_trend_smd = round(max(trend_smds), 4),
      pct_at_cap    = round(100 * mean(controls$weights >= cap), 2)
    )
  })
  
  cat("\nWeight cap sensitivity:\n")
  print(cap_results)
  invisible(cap_results)
}

wd_A <- weight_diagnostics(s2_A, "A_excl_zero")
wd_B <- weight_diagnostics(s2_B, "B_incl_zero")

# Apply weight cap of 5 to Analysis B
# Cap chosen at 5: maximises efficiency improvement (0.672) while
# keeping trend SMDs low (0.038). Cap of 10 gives less improvement (0.585)
# for no meaningful balance benefit. Analysis A needs no cap (max = 12,
# cap sensitivity shows negligible change — no cap preserves all information).
s2_B$primary_data <- s2_B$primary_data %>%
  mutate(weights = pmin(weights, 5))

cat("\n--- Analysis B weights after cap at 5 ---\n")
s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  summarise(
    max_weight  = max(weights),
    effective_n = round(sum(weights)^2 / sum(weights^2), 1),
    efficiency  = round((sum(weights)^2 / sum(weights^2)) / n(), 3)
  ) %>% print()

# =============================================================================
# SECTION 10 — EXTRACT AND SAVE MATCHED DATASETS
# =============================================================================
# Two output formats per analysis:
#
# (a) OA_matched_treated/donors: OA + weights only — minimal files for
#     joining to panel outcome data in the DiD script.
#
# (b) OA_matched_full: complete matched dataset including all matching
#     variables, distances, and weights — used for balance reporting,
#     sensitivity analyses, and covariate selection in C&S outcome model.
#
# IMPORTANT for DiD:
#   Pass weights to att_gt(..., weightsname = "weights")
#   Covariates flagged by run_balance() should enter xformla in att_gt()
#   Analysis A: primary specification
#   Analysis B: robustness check (weights capped at 5)
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 10: SAVE MATCHED DATASETS\n")
cat("====================================================\n\n")

# --- Analysis A ---
matched_A_treated  <- s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights)

matched_A_controls <- s2_A$primary_data %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

cat("Analysis A — Treated:", nrow(matched_A_treated),
    "| Controls:", nrow(matched_A_controls), "\n")
cat("  Weight range (controls): [",
    round(min(matched_A_controls$weights), 3), ",",
    round(max(matched_A_controls$weights), 3), "]\n")

# --- Analysis B (weights already capped at 5) ---
matched_B_treated  <- s2_B$primary_data %>%
  filter(treat_indicator == 1) %>%
  select(OA, weights)

matched_B_controls <- s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  select(OA, weights)

cat("Analysis B — Treated:", nrow(matched_B_treated),
    "| Controls:", nrow(matched_B_controls), "\n")
cat("  Weight range (controls): [",
    round(min(matched_B_controls$weights), 3), ",",
    round(max(matched_B_controls$weights), 3), "]\n")

# Final integrity checks before saving
stopifnot(
  "A: treated weights should all be 1" =
    all(matched_A_treated$weights == 1),
  "B: treated weights should all be 1" =
    all(matched_B_treated$weights == 1),
  "A: no NA weights in controls" =
    !anyNA(matched_A_controls$weights),
  "B: no NA weights in controls" =
    !anyNA(matched_B_controls$weights),
  "A: no duplicate treated OAs" =
    !anyDuplicated(matched_A_treated$OA),
  "B: no duplicate treated OAs" =
    !anyDuplicated(matched_B_treated$OA),
  "B: max control weight <= 5" =
    max(matched_B_controls$weights) <= 5
)
cat("All integrity checks passed.\n\n")

# Save minimal OA + weight files (for joining to panel data in DiD script)
saveRDS(matched_A_treated,  here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A_controls, here("data", "processed", "OA_matched_donors_A.rds"))
saveRDS(matched_B_treated,  here("data", "processed", "OA_matched_treated_B.rds"))
saveRDS(matched_B_controls, here("data", "processed", "OA_matched_donors_B.rds"))

# Save full matched datasets (for balance reporting and covariate selection)
saveRDS(s2_A$primary_data, here("data", "processed", "OA_matched_full_A.rds"))
saveRDS(s2_B$primary_data, here("data", "processed", "OA_matched_full_B.rds"))

cat("Saved:\n")
cat("  OA_matched_treated_A.rds  —", nrow(matched_A_treated), "treated OAs\n")
cat("  OA_matched_donors_A.rds   —", nrow(matched_A_controls), "control OAs\n")
cat("  OA_matched_treated_B.rds  —", nrow(matched_B_treated), "treated OAs\n")
cat("  OA_matched_donors_B.rds   —", nrow(matched_B_controls), "control OAs\n")
cat("  OA_matched_full_A.rds     — full matched dataset, Analysis A\n")
cat("  OA_matched_full_B.rds     — full matched dataset, Analysis B\n")

# =============================================================================
# SECTION 11 — DESIGN SUMMARY
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 11: DESIGN SUMMARY\n")
cat("====================================================\n\n")

cat("STAGE 1  | MDM on", length(stage1_vars),
    "vars (road + urban + socdem)\n")
cat("         | No caliper — symmetric treatment of all variables\n")
cat("         | replace=TRUE, ratio=10, exact=~country\n")
cat("         | Purpose: structural plausibility filter\n\n")

cat("STAGE 2  | MDM on", length(stage2_vars),
    "vars (7 trends + 7 levels, per road-km)\n")
cat("         | No caliper\n")
cat("         | replace=TRUE, ratio=5 (selected on trend SMD balance)\n")
cat("         | Purpose: parallel trends support via trajectory matching\n\n")

cat("ANALYSIS A | Zero-injury treated OAs EXCLUDED (PRIMARY)\n")
cat("           | N treated:", nrow(matched_A_treated),
    "| N controls:", nrow(matched_A_controls), "\n")
cat("           | Max control weight:", round(max(matched_A_controls$weights), 2),
    "| No weight cap applied\n\n")

cat("ANALYSIS B | Zero-injury treated OAs INCLUDED (ROBUSTNESS CHECK)\n")
cat("           | N treated:", nrow(matched_B_treated),
    "| N controls:", nrow(matched_B_controls), "\n")
cat("           | Control weights capped at 5 (Scottish zero-injury sparsity)\n\n")

cat("NEXT STEPS FOR DiD:\n")
cat("  1. Join matched OA lists to panel outcome data\n")
cat("  2. Run att_gt(..., weightsname = 'weights') on both datasets\n")
cat("  3. Include level covariates flagged in Section 7 in xformla\n")
cat("  4. If ATT(A) ≈ ATT(B): A is primary, B is robustness check\n")
cat("  5. If ATT(A) ≠ ATT(B): report both; discuss heterogeneous effects\n")
cat("     by pre-treatment injury exposure level\n")






#NEXT STEP:
  #   Pass matched datasets to att_gt(..., weightsname = "weights") in did package.
  #   If ATT(A) ≈ ATT(B): A confirmed as primary, B reported as robustness check.
  #   If ATT(A) ≠ ATT(B): report both; difference reveals heterogeneous effects
  #                        between injury-exposed and zero-injury road areas.