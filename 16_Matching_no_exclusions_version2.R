# =============================================================================
#  designe: first establish that comparators are plausible structural substitutes, then select among them on pre-treatment outcome dynamics
#  OA-Level Matching — Two-Stage MDM Design (No-Caliper)
#  Stage 1: MDM on structural + sociodemographic variables
#  Stage 2: MDM on pre-treatment injury trends + levels, no caliper,
#           
### 
## there is 2 analysis 
# Analysis A  without 0 injury OAs: population is a subset of all treated areas.
# Analysis B  effect across all treated areas, the estimate is broader in scope but less credible in execution.
# =============================================================================
#
# DESIGN RATIONALE
#
# #   All matching variables are treated symmetrically. P
#
#   STAGE 1 — Structural + sociodemographic restriction (MDM, no caliper)
#     Nearest-neighbour MDM on 15 road-network, urban-form, and
#     sociodemographic variables. replace=TRUE, ratio=10 to maximise
#     the candidate pool passed to Stage 2.
#     
#
#   STAGE 2 — Injury-dynamics matching (MDM, no caliper)
#     Within the Stage 1 pool, nearest-neighbour MDM on 14
#     pre-treatment injury variables (7 trends + 7 levels, per road-km).
#     replace=TRUE: globally best trajectory match for every treated OA.
#     Ratios 1, 2, 5 tested; primary selected on trend SMD balance.
#     
#
#   TWO PARALLEL ANALYSES:
#     Analysis A — EXCLUDING zero-injury treated OAs
#       Stage 2 variables well-defined for all retained OAs.
#       ATT scope: road-exposed OAs with pre-period injury exposure.
#     Analysis B — INCLUDING zero-injury treated OAs
#       Stage 2 variables set to 0 for zero-injury treated OAs.
#       They match to zero-injury controls with near-zero Stage 2
#       distances; parallel trends rests on Stage 1 structural similarity.
#       ATT scope: all road-exposed OAs within scheme boundaries.
#     Comparing A and B empirically tests whether zero-injury exclusion
#     biases the ATT estimate.
#
#   FIXED PRE-MATCHING EXCLUSIONS (both analyses):
#     Buffer OAs (contamination risk) and zero-road OAs (undefined rates).
#
#   WEIGHTS NOTE:
#     replace=TRUE at Stage 2 assigns weight > 1 to controls matched to
#     multiple treated OAs. Pass the weights column to att_gt() in the C&S estimation:
#       att_gt(..., weightsname = "weights")
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

# =============================================================================
# 
OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))


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

cat("\n--- Road network summary ---\n")
OA_matching_dataset %>%
  summarise(
    total_OAs  = n_distinct(OA),
    zero_roads = sum(n_roads == 0 | is.na(n_roads)),
    pct_zero   = round(100 * zero_roads / total_OAs, 2)
  ) %>% print()

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

stage1_road <- c(
  "road_density_m_km2", "road_length_km",
  "pct_A_road", "pct_B_road", "pct_minor_road"
)
stage1_urban  <- c("dist_citycentre", "pop_density")
stage1_socdem <- c(
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)
stage1_vars <- c(stage1_road, stage1_urban, stage1_socdem)

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
stage2_vars <- c(stage2_trends, stage2_levels)

cat("Stage 1 variables:", length(stage1_vars),
    "(road:", length(stage1_road),
    "| urban:", length(stage1_urban),
    "| socdem:", length(stage1_socdem), ")\n")
cat("Stage 2 variables:", length(stage2_vars),
    "(trends:", length(stage2_trends),
    "| levels:", length(stage2_levels), ")\n")
cat("No calipers at either stage.\n")

# =============================================================================
# SECTION 2 — BUILD BOTH ANALYSIS DATASETS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 2: DATASET CONSTRUCTION\n")
cat("====================================================\n\n")

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

# Analysis A: exclude zero-injury treated OAs
data_A <- base_filter %>%
  filter(!(treated_OA == 1 & zero_injury_OA == 1))

# Analysis B: include zero-injury OAs; set their Stage 2 vars to 0
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
# SECTION 3 — WINSORISE STAGE 1 VARIABLES
# =============================================================================
# Caps extreme values at 1st/99th percentile — no rows removed.
# Prevents outliers from distorting the Mahalanobis covariance matrix.
# Stage 2 variables winsorised separately in Section 5 on treated dist.
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 3: WINSORISING\n")
cat("====================================================\n\n")

winsorise_s1 <- function(data, vars) {
  data %>%
    mutate(across(
      all_of(intersect(vars, names(.))),
      ~ { q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2]) }
    ))
}

data_A_clean <- winsorise_s1(data_A, stage1_vars)
data_B_clean <- winsorise_s1(data_B, stage1_vars)

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
# SECTION 4 — STAGE 1 MATCHING: NO CALIPER 
# =============================================================================
#  Distance computed on matched pairs only (via subclass),
#  #   Subclass-based = 800 × 10 = 8K distances 
##
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 4: STAGE 1 MATCHING (no caliper)\n")
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
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  mm <- m$match.matrix
  
  cat("  Treated in:", nrow(mm), "\n")
  
  # Distance computation for matched pairs
  S_s1 <- cov(data[as.integer(rownames(mm)), s1_vars],
              use = "pairwise.complete.obs")
  
  dist_s1 <- map_df(seq_len(nrow(mm)), function(i) {
    
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    
    c_indices <- mm[i, ]
    c_indices <- c_indices[!is.na(c_indices)]
    
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
  
  # Keep ALL treated and matched controls (no trimming)
  treated_matched <- data[as.integer(rownames(mm)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  
  control_row_indices <- unique(as.integer(mm[!is.na(mm)]))
  
  controls_matched <- data[control_row_indices, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  cat("  Treated retained:", nrow(treated_matched), "\n")
  cat("  Unique controls in pool:", nrow(controls_matched), "\n")
  
  # Balance diagnostics
  cat("\n  Stage 1 balance:\n")
  
  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  
  print(smd_df)
  
  cat("\n  Country distribution (treated):\n")
  print(table(treated_matched$country))
  
  # Cross-country check
  cross_check <- dist_s1 %>%
    left_join(
      data %>% select(OA, country) %>% rename(control_country = country),
      by = c("control_OA" = "OA")
    ) %>%
    left_join(
      data %>% select(OA, country) %>% rename(treated_country = country),
      by = c("treated_OA" = "OA")
    ) %>%
    filter(treated_country != control_country)
  
  cat("  Cross-country pairs:", nrow(cross_check),
      if (nrow(cross_check) == 0) "✓" else "WARNING", "\n")
  
  lp <- love.plot(
    m, threshold = 0.1, abs = TRUE, stars = "std",
    var.order    = "unadjusted",
    title        = paste("Stage 1 balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before", "After")
  ) + theme(plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  
  ggsave(
    here("output", paste0("s1_balance_", label, ".png")),
    lp,
    width = 13,
    height = 9,
    dpi = 300
  )
  
  cat("  Love plot saved.\n\n")
  
  list(
    matchit_obj      = m,
    dist_s1          = dist_s1,
    treated  = treated_matched,
    controls = controls_matched,
    bal              = bt
  )
}
s1_A <- run_stage1(data_A_clean, s1_vars_A, "A_excl_zero")
s1_B <- run_stage1(data_B_clean, s1_vars_B, "B_incl_zero")



bal.tab(s1_A$matchit_obj, un = TRUE)
love.plot(s1_A$matchit_obj)

bal.tab(s1_B$matchit_obj, un = TRUE)
love.plot(s1_B$matchit_obj)




# =============================================================================
# SECTION 5 — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 5: STAGE 2 DATA PREPARATION\n")
cat("====================================================\n\n")

prepare_s2_data <- function(s1_result,s2_vars){
  
  s2_raw <- bind_rows(
    s1_result$treated,
    s1_result$controls
  ) %>%
    select(-any_of(c("weights","subclass","distance")))
  
  treated_ref <- s2_raw %>%
    filter(treat_indicator==1)
  
  s2_raw %>%
    mutate(across(
      all_of(intersect(s2_vars,names(.))),
      ~{
        q <- quantile(treated_ref[[cur_column()]],
                      c(0.01,0.99),
                      na.rm=TRUE)
        
        pmin(pmax(.,q[1]),q[2])
      }
    ))
  
}

s2_data_A <- prepare_s2_data(s1_A, s2_vars_A)
s2_data_B <- prepare_s2_data(s1_B, s2_vars_B)



# =============================================================================
# SECTION 6 — STAGE 2 MATCHING: NO CALIPER 
# =============================================================================
#
# NOTE FOR ANALYSIS B:
#   Zero-injury treated OAs have Stage 2 variables = 0. They match to
#   zero-injury controls with near-zero distances. The A vs B comparison in Section 8
#   reveals whether their inclusion changes the ATT estimate.
#
# WEIGHTS:
#   replace=TRUE assigns weight > 1 to controls matched to multiple treated
#   OAs. Pass weights column to att_gt(weightsname = "weights").
# =============================================================================


#SECTION 6: STAGE 2 MATCHING

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
  )
  
  mm <- m$match.matrix
  
  # Compute Mahalanobis distances from the  matrix
  S_s2 <- cov(data[as.integer(rownames(mm)), s2_vars],
              use = "pairwise.complete.obs")
  
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
      OA             = treated_id,         # treated OA gets mean dist to its controls
      control_OA     = data[as.integer(c_indices), "OA", drop = FALSE] %>% pull(OA),
      mdist_control  = dists,
      mdist_treated  = mean(dists)         # summarise treated OA by mean dist to its matches
    )
  })
  
  # Treated OA distances = mean distance across their matched controls
  treated_dists <- dist_s2 %>%
    distinct(OA, mdist_treated) %>%
    rename(mdist = mdist_treated)
  
  # Control OA distances = mean distance across all treated OAs they were matched to
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
  cat("Distance summary (treated only):\n")
  print(summary(matched_data$mdist[matched_data$treat_indicator == 1]))
  
  list(
    matchit_obj   = m,
    primary_ratio = ratio,
    primary_data  = matched_data
  )
}
s2_A_r1 <- run_stage2(s2_data_A,s2_vars_A,1,"A")
s2_A_r2 <- run_stage2(s2_data_A,s2_vars_A,2,"A")
s2_A_r5 <- run_stage2(s2_data_A,s2_vars_A,5,"A")

s2_B_r1 <- run_stage2(s2_data_B,s2_vars_B,1,"B")
s2_B_r2 <- run_stage2(s2_data_B,s2_vars_B,2,"B")
s2_B_r5 <- run_stage2(s2_data_B,s2_vars_B,5,"B")




# r5 is the best option 
#### r5 has the best trend SMDs of all four specifications (0.034) — the strongest empirical support for parallel trends
# Level balance is better at R5, meaning less work for the doubly-robust correction
#efficiency gained  from r5 is real and doesnt comes at the cost of worse balance and more covariates needed in the outcome model



#
# =============================================================================
# — BALANCE DIAGNOSTICS
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
    var.order = "unadjusted", title = paste("Stage 2 balance —", label),
    shapes = c("circle filled", "triangle filled"),
    colors = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) + theme(axis.text.y = element_text(size = 9),
            plot.margin = margin(10, 30, 10, 180), legend.position = "bottom")
  ggsave(here("output", paste0("s2_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  invisible(bt)
}



s2_A <- s2_A_r5   
s2_B <- s2_B_r5

run_balance(s2_A, "A_excl_zero")
run_balance(s2_B, "B_incl_zero")

## not that different in the balance ..... we should run both and compare the ATT estimates. 
# If they're similar, the zero-injury OAs aren't driving results and can report A as primary with B as a robustness check. 
# If they diverge, that's a meaningful finding about whether the scheme affected previously-uninjured roads differently,
#which is worth discussing substantively rather than resolving by choosing one specification.

# =============================================================================
# — COMPARE ANALYSES A AND B
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 8: COMPARISON — ANALYSES A vs B\n")
cat("====================================================\n\n")

treated_A <- s2_A$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data %>% filter(treat_indicator == 1) %>% pull(OA)

only_in_B <- setdiff(treated_B, treated_A)
only_in_A <- setdiff(treated_A, treated_B)
in_both   <- intersect(treated_A, treated_B)

cat("--- Sample overlap ---\n")
cat("  Treated in A:              ", length(treated_A), "\n")
cat("  Treated in B:              ", length(treated_B), "\n")
cat("  In both:                   ", length(in_both), "\n")
cat("  Only in B (zero-injury):   ", length(only_in_B), "\n")
cat("  Only in A (trimmed from B):", length(only_in_A), "\n\n")

if (length(only_in_B) > 0) {
  
  cat("--- Structural characteristics of zero-injury OAs added in B ---\n")
  compare_grp <- s2_B$primary_data %>%
    filter(treat_indicator == 1) %>%
    mutate(grp = if_else(OA %in% only_in_B, "Zero-injury (B only)", "Injury-exposed (both)"))
  
  s1_avail <- intersect(stage1_vars, names(compare_grp))
  cov_diff <- compare_grp %>%
    select(grp, all_of(s1_avail)) %>% pivot_longer(-grp) %>%
    group_by(name, grp) %>%
    summarise(mean = mean(value, na.rm=T), sd = sd(value, na.rm=T), .groups="drop") %>%
    pivot_wider(names_from = grp, values_from = c(mean, sd)) %>%
    mutate(SMD = (`mean_Injury-exposed (both)` - `mean_Zero-injury (B only)`) /
             sqrt((`sd_Injury-exposed (both)`^2 + `sd_Zero-injury (B only)`^2) / 2)) %>%
    arrange(desc(abs(SMD)))
  
  print(cov_diff %>%
          select(name, `mean_Injury-exposed (both)`, `mean_Zero-injury (B only)`, SMD),
        n = Inf)
  
  cat("\n  Country distribution of zero-injury OAs in B:\n")
  print(table(s2_B$primary_data %>% filter(OA %in% only_in_B) %>% pull(country)))
  
  zero_dists <- s2_B$primary_data %>%
    filter(OA %in% only_in_B)
  
  cat("\n  Stage 2 distances for zero-injury pairs (expect near-zero):\n")
  print(summary(zero_dists$mdist))
}

# Distance comparison table
cat("\n--- Stage 2 distance distribution: A vs B ---\n")

dist_compare <- bind_rows(
  s2_A$primary_data %>% select(OA, mdist) %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$primary_data %>% select(OA, mdist) %>% mutate(analysis = "B: incl zero-injury")
) %>%
  group_by(analysis) %>%
  summarise(
    n      = n(),
    mean   = round(mean(mdist,             na.rm = TRUE), 3),
    median = round(median(mdist,           na.rm = TRUE), 3),
    p90    = round(quantile(mdist, 0.90,   na.rm = TRUE), 3),
    p95    = round(quantile(mdist, 0.95,   na.rm = TRUE), 3),
    max    = round(max(mdist,              na.rm = TRUE), 3),
    .groups = "drop"
  )

print(dist_compare)


# A Mahalanobis distance of 115/123  means some treated OAs are being matched to controls that are very far away#
# in injury-dynamic space, which undermines the parallel trends assumption for those specific units.

# who are those OAs 
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  arrange(desc(mdist)) %>%
  select(OA, mdist, starts_with("trend_"), starts_with("mean_")) %>%
  head(30)

# These are OAs with extreme pre-treatment injury dynamics that sit far from the control distribution.

# See the distribution of treated OA distances
s2_A$primary_data %>%
  filter(treat_indicator == 1) %>%
  summarise(
    n_over_20  = sum(mdist > 20),
    n_over_30  = sum(mdist > 30),
    n_over_50  = sum(mdist > 50),
    n_over_100 = sum(mdist > 100)
  )

# 30 as caliper threshold — will  drop 11 OAs (1.4% of treated sample), 
# this eliminate the genuinely problematic matches,then all retained pairs have reasonable parallel trends support.

# Identify treated OAs to drop
drop_OAs <- s2_A$primary_data %>%
  filter(treat_indicator == 1, mdist > 30) %>%
  pull(OA)
drop_OAs


# Are they concentrated in particular schemes or years?
s2_A$primary_data %>%
  filter(OA %in% drop_OAs) %>%
  select(OA, mdist, country) %>%
  arrange(desc(mdist))





# =============================================================================
# SECTION 9 — EXTRACT AND SAVE MATCHED SAMPLES
# =============================================================================

extract_matched <- function(s2_result, label) {
  treated  <- s2_result$primary_data %>%
    filter(treat_indicator == 1) %>%
    select(OA, weights)
  controls <- s2_result$primary_data %>%
    filter(treat_indicator == 0) %>%
    select(OA, weights)
  cat(label, "— Treated:", nrow(treated),
      "| Controls:", nrow(controls), "\n")
  list(treated = treated, controls = controls)
}

matched_A <- extract_matched(s2_A, "Analysis A (excl zero-injury)")
matched_B <- extract_matched(s2_B, "Analysis B (incl zero-injury)")




#####################################################################################

# checking wieght distribution A
w_checkA <- s2_A$primary_data %>%
  filter(treat_indicator == 0) %>%
  summarise(
    n_controls     = n(),
    mean_weight    = mean(weights),
    max_weight     = max(weights),
    pct_weight_gt5 = mean(weights > 5) * 100,
    effective_n    = sum(weights)^2 / sum(weights^2)
  )
print(w_checkA)

# checking wieght distribution B
w_checkB <- s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  summarise(
    n_controls     = n(),
    mean_weight    = mean(weights),
    max_weight     = max(weights),
    pct_weight_gt5 = mean(weights > 5) * 100,
    effective_n    = sum(weights)^2 / sum(weights^2)
  )
print(w_checkB)


top_controls <- s2_B$primary_data %>%
  filter(treat_indicator == 0) %>%
  arrange(desc(weights)) %>%
  select(OA, weights, country, scheme) %>%
  head(20)
print(top_controls)


table(s2_B$primary_data$weights > 10)
hist(s2_B$primary_data$weights, breaks = 50)
### only 10  have weights above 10 and they dominated the controls 

# =============================================================================
# WEIGHT DIAGNOSTICS + CAP SENSITIVITY
# =============================================================================

weight_diagnostics <- function(s2_result, analysis_label) {
  
  controls <- s2_result$primary_data %>%
    filter(treat_indicator == 0)
  
  cat("\n--- Weight diagnostics:", analysis_label, "---\n")
  
  # Full distribution
  cat("\nWeight distribution (controls):\n")
  print(quantile(controls$weights,
                 probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)))
  
  cat("\nEffective N:", round(sum(controls$weights)^2 / sum(controls$weights^2), 1),
      "vs nominal N:", nrow(controls), "\n")
  cat("Efficiency ratio:",
      round((sum(controls$weights)^2 / sum(controls$weights^2)) / nrow(controls), 3), "\n")
  
  # By country
  cat("\nBy country:\n")
  controls %>%
    group_by(country) %>%
    summarise(
      n          = n(),
      mean_w     = round(mean(weights), 2),
      max_w      = round(max(weights), 2),
      eff_n      = round(sum(weights)^2 / sum(weights^2), 1),
      .groups = "drop"
    ) %>% print()
  
  # Weight cap sensitivity: cap at 10, 20, Inf (no cap)
  caps <- c(5, 10, 20, Inf)
  
  cap_results <- map_df(caps, function(cap) {
    
    capped <- controls %>%
      mutate(w_cap = pmin(weights, cap))
    
    # Recalculate effective N
    eff_n_cap <- sum(capped$w_cap)^2 / sum(capped$w_cap^2)
    
    # Recalculate trend SMDs with capped weights (weighted mean)
    treated_rows <- s2_result$primary_data %>%
      filter(treat_indicator == 1)
    
    trend_smds <- map_dbl(stage2_trends, function(v) {
      t_mean <- mean(treated_rows[[v]], na.rm = TRUE)
      c_mean <- weighted.mean(controls[[v]], w = capped$w_cap, na.rm = TRUE)
      pooled_sd <- sqrt((var(treated_rows[[v]], na.rm=T) +
                           var(controls[[v]], na.rm=T)) / 2)
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



# The zero-injury OAs (B only) are meaningfully different from the injury-exposed OAs on several structural dimensions:
# Much shorter roads (0.24 km vs 1.10 km, SMD 0.86) — these are small, low-traffic OAs
# Higher population density (17,681 vs 10,985, SMD -0.81) — more urban/residential
# Almost entirely minor roads (96% vs 83%, SMD -0.73) — very little A/B road exposure
# Predominantly Scotland (114 of 191 = 60%) — geographically skewed





#save 

#


saveRDS(matched_A$treated,  here("data", "processed", "OA_matched_treated_A.rds"))
saveRDS(matched_A$controls, here("data", "processed", "OA_matched_donors_A.rds"))



saveRDS(matched_B$treated,  here("data", "processed", "OA_matched_treated_B.rds"))
saveRDS(matched_B$controls, here("data", "processed", "OA_matched_donors_B.rds"))







# =============================================================================
# SECTION 10 — DESIGN SUMMARY
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 10: DESIGN SUMMARY\n")
cat("====================================================\n\n")

cat("STAGE 1  | MDM on", length(stage1_vars), "vars (road + urban + socdem)\n")
cat("         | No caliper — all variables treated symmetrically\n")
cat("         | replace=TRUE, ratio=10\n")

cat("STAGE 2  | MDM on", length(stage2_vars), "vars (7 trends + 7 levels, per km)\n")
cat("         | No caliper\n")
cat("         | replace=TRUE, ratios 1/2/5 tested\n")

cat("Country exact constraint at both stages.\n\n")

cat("ANALYSIS A | Zero-injury treated OAs EXCLUDED\n")
cat("           | Ratio:", s2_A$primary_ratio,
    "| N treated:", nrow(matched_A$treated), "\n\n")

cat("ANALYSIS B | Zero-injury treated OAs INCLUDED (Stage 2 vars = 0)\n")
cat("           | Ratio:", s2_B$primary_ratio,
    "| N treated:", nrow(matched_B$treated), "\n\n")

cat("WEIGHTS: Pass weights column to att_gt(weightsname = 'weights')\n\n")

cat("NEXT STEP: Run C&S on both samples.\n")
cat("  ATT(A) ≈ ATT(B) => A is primary, B is robustness check\n")
cat("  ATT(A) ≠ ATT(B) => report both; discuss zero-injury subgroup\n")


########## final matching data #


#####    If ATT(A) ≈ ATT(B): The zero-injury OAs don't change the result. 
##    B is the primary analysis — it makes no upfront exclusion, includes all road-exposed treated OAs,#
# and the broader scope is more defensible. Report A as a robustness check in a footnote or appendix.


### If ATT(A) ≠ ATT(B) meaningfully: The zero-injury OAs are doing something different. 
# This is actually an interesting finding — it means the effect for low-exposure roads differs from high-exposure roads.
# In this case report both and discuss what the difference tells you about heterogeneity in the treatment effect.
# The one scenario where A would be definitively preferred over B is if the Stage 2 balance for B is substantially worse than for A 
#— specifically if trend SMDs in B exceed 0.10 while in A they don't. 
# That would mean the zero-injury OAs are degrading the parallel trends support. 







