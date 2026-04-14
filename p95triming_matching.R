# =============================================================================
#  designe: first establish that comparators are plausible structural substitutes, then select among them on pre-treatment outcome dynamics
#  OA-Level Matching — Two-Stage MDM Design (No-Caliper)
#  Stage 1: MDM on structural + sociodemographic variables, no caliper,
#           95th-percentile distance trim post-hoc
#  Stage 2: MDM on pre-treatment injury trends + levels, no caliper,
#           95th-percentile distance trim post-hoc
### 
# =============================================================================
#
# DESIGN RATIONALE
#
# #   All matching variables are treated symmetrically. Poor matches are
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
#     sociodemographic variables. replace=TRUE, ratio=10 to maximise
#     the candidate pool passed to Stage 2.
#     Post-hoc trim at P95 of Stage 1 Mahalanobis distances (computed
#     on matched pairs only via subclass — NOT all-pairs) removes the
#     most structurally distant controls before injury-dynamics matching.
#     Trim criterion: minimum distance per treated OA (best available
#     structural match). Treated OAs whose closest structural comparator
#     exceeds P95 are removed.
#     Country (England/Scotland) enforced as hard exact constraint at
#     both stages — STATS19 and STATS19-Scotland differ in recording
#     thresholds and severity classification.
#
#   STAGE 2 — Injury-dynamics matching (MDM, no caliper)
#     Within the Stage 1-trimmed pool, nearest-neighbour MDM on 14
#     pre-treatment injury variables (7 trends + 7 levels, per road-km).
#     replace=TRUE: globally best trajectory match for every treated OA.
#     Ratios 1, 2, 5 tested; primary selected on trend SMD balance.
#     Post-hoc P95 trim removes the most trajectory-distant pairs,
#     identified by subclass (keeps both treated AND matched controls).
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
# SECTION 4 — STAGE 1 MATCHING: NO CALIPER + P95 TRIM
# =============================================================================
#  Distance computed on matched pairs only (via subclass),
#  #   Subclass-based = 800 × 10 = 8K distances 
## P95 trim uses min(mdist) per treated OA (best-match distance),
#    Using minimum correctly identifies OAs with no close
#   structural comparator. 
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
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) return(NULL)
  
  mm <- m$match.matrix  # rows = treated (row indices into data), cells = control indices
  
  cat("  Treated in:", nrow(mm), "\n")
  
  # -- Distances via match.matrix (matched pairs only, fast) -----------------
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
  
  # -- P95 trim on MIN distance per treated OA --------------------------------
  best_dist <- dist_s1 %>%
    group_by(treated_OA) %>%
    summarise(best_dist = min(mdist), .groups = "drop")
  
  p95          <- quantile(best_dist$best_dist, 0.95)
  kept_treated <- best_dist %>% filter(best_dist <= p95) %>% pull(treated_OA)
  
  # -- Extract trimmed treated and their controls from match.matrix -----------
  # NO subclass lookup needed — replace=TRUE doesn't produce subclass in md_raw.
  # We work directly from match.matrix row indices.
  kept_rows    <- which(data[as.integer(rownames(mm)), "OA"] %in% kept_treated)
  kept_mm      <- mm[kept_rows, , drop = FALSE]
  
  treated_trimmed <- data[as.integer(rownames(kept_mm)), , drop = FALSE] %>%
    mutate(treat_indicator = 1L)
  
  # All unique control row indices for kept treated OAs
  control_row_indices <- unique(as.integer(kept_mm[!is.na(kept_mm)]))
  controls_trimmed    <- data[control_row_indices, , drop = FALSE] %>%
    mutate(treat_indicator = 0L)
  
  n_in  <- nrow(mm)
  n_out <- nrow(treated_trimmed)
  
  cat("  P95 threshold (best-match):", round(p95, 2), "\n")
  cat("  Treated trimmed:", n_in - n_out, "\n")
  cat("  Treated retained:", n_out, "\n")
  cat("  Unique controls in pool:", nrow(controls_trimmed), "\n")
  
  # Balance (on full matched sample before trim — MatchIt object is pre-trim)
  cat("\n  Stage 1 balance:\n")
  bt     <- bal.tab(m, thresholds = c(m = 0.1), un = TRUE)
  smd_df <- bt$Balance %>%
    rownames_to_column("variable") %>%
    select(variable, Diff.Un, Diff.Adj) %>%
    arrange(desc(abs(Diff.Adj)))
  print(smd_df)
  
  cat("\n  Country distribution (post-trim treated):\n")
  print(table(treated_trimmed$country))
  
  # Cross-country check: no treated-control pair should cross country
  cross_check <- treated_trimmed %>%
    select(OA, country) %>%
    rename(treated_country = country) %>%
    left_join(
      dist_s1 %>% filter(treated_OA %in% kept_treated) %>%
        left_join(data %>% select(OA, country) %>% rename(control_country = country),
                  by = c("control_OA" = "OA")),
      by = c("OA" = "treated_OA")
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
  ggsave(here("output", paste0("s1_balance_", label, ".png")),
         lp, width = 13, height = 9, dpi = 300)
  cat("  Love plot saved.\n\n")
  
  list(
    matchit_obj      = m,
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
  
  cat("  Treated OAs:", sum(s2_raw$treat_indicator == 1), "\n")
  cat("  Control OAs:", sum(s2_raw$treat_indicator == 0), "\n")
  
  # Winsorise Stage 2 vars on TREATED distribution only
  s2_present  <- intersect(s2_vars, names(s2_raw))
  treated_ref <- s2_raw %>% filter(treat_indicator == 1)
  
  s2_data <- s2_raw %>%
    mutate(across(
      all_of(s2_present),
      ~ { q <- quantile(treated_ref[[cur_column()]], c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2]) }
    ))
  
  cat("  Stage 2 vars winsorised on treated distribution.\n\n")
  list(data = s2_data, s2_present = s2_present)
}

s2_prep_A <- prepare_s2_data(s1_A, s2_vars_A, "A_excl_zero")
s2_prep_B <- prepare_s2_data(s1_B, s2_vars_B, "B_incl_zero")

# =============================================================================
# SECTION 6 — STAGE 2 MATCHING: NO CALIPER + P95 TRIM
# =============================================================================
#P95 trim filters on subclass (keeps treated AND their matched
#   controls). #
# NOTE FOR ANALYSIS B:
#   Zero-injury treated OAs have Stage 2 variables = 0. They match to
#   zero-injury controls with near-zero distances and will NOT be caught
#   by the P95 trim. The A vs B comparison in Section 8
#   reveals whether their inclusion changes the ATT estimate.
#
# WEIGHTS:
#   replace=TRUE assigns weight > 1 to controls matched to multiple treated
#   OAs. Pass weights column to att_gt(weightsname = "weights").
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 6: STAGE 2 MATCHING (no caliper + P95 trim)\n")
cat("====================================================\n\n")

# Compute Mahalanobis distances for matched pairs (subclass-based)
# Replace compute_mdist with this match.matrix-based version
compute_mdist <- function(m, data, vars) {
  mm <- m$match.matrix
  S  <- cov(data[as.integer(rownames(mm)), vars],
            use = "pairwise.complete.obs")
  
  map_df(seq_len(nrow(mm)), function(i) {
    t_idx      <- as.integer(rownames(mm)[i])
    trow       <- data[t_idx, , drop = FALSE]
    treated_id <- trow[["OA"]]
    
    c_indices <- mm[i, ]
    c_indices <- c_indices[!is.na(c_indices)]
    if (length(c_indices) == 0) return(tibble())
    
    map_df(seq_along(c_indices), function(j) {
      crow <- data[as.integer(c_indices[j]), , drop = FALSE]
      tibble(
        subclass   = rownames(mm)[i],   # treated row index as subclass ID
        treated_OA = treated_id,
        control_OA = crow[["OA"]],
        mdist      = mahalanobis(
          x      = as.numeric(crow[vars]),
          center = as.numeric(trow[vars]),
          cov    = S
        )
      )
    })
  })
}

# Replace summarise_mdist with a safe version that checks for empty input
summarise_mdist <- function(df, spec) {
  if (is.null(df) || nrow(df) == 0 || !"mdist" %in% names(df)) {
    return(tibble(spec = spec, n_pairs = 0L, mean = NA, median = NA,
                  p90 = NA, p95 = NA, max = NA))
  }
  df %>% summarise(
    spec   = spec,
    n_pairs = n(),
    mean   = round(mean(mdist),           3),
    median = round(median(mdist),         3),
    p90    = round(quantile(mdist, 0.90), 3),
    p95    = round(quantile(mdist, 0.95), 3),
    max    = round(max(mdist),            3)
  )
}
run_one_s2 <- function(data, s2_vars, ratio) {
  formula <- reformulate(s2_vars, response = "treat_indicator")
  tryCatch(
    matchit(formula, data = data, method = "nearest",
            distance = "mahalanobis", ratio = ratio,
            replace = TRUE, exact = ~ country),
    error = function(e) {
      cat("  FAILED (ratio =", ratio, "):", conditionMessage(e), "\n"); NULL
    }
  )
}

run_stage2_multi <- function(s2_prep, analysis_label) {
  
  cat("\n--- Stage 2 (multi-ratio):", analysis_label, "---\n")
  data    <- s2_prep$data
  s2_vars <- s2_prep$s2_present
  
  # Run ratios 1, 2, 5
  ratios  <- c(1, 2, 5)
  s2_fits <- setNames(
    lapply(ratios, function(r) run_one_s2(data, s2_vars, r)),
    paste0("r", ratios)
  )
  
  # Distances for each ratio (subclass-based — fast)
  mdist_list <- lapply(names(s2_fits), function(nm) {
    m <- s2_fits[[nm]]
    if (is.null(m)) return(NULL)
    compute_mdist(m, data, s2_vars)
  })
  names(mdist_list) <- names(s2_fits)
  
  # Distance summary
  summarise_mdist <- function(df, spec) {
    df %>% summarise(
      spec = spec, n_pairs = n(),
      mean   = round(mean(mdist),           3),
      median = round(median(mdist),         3),
      p90    = round(quantile(mdist, 0.90), 3),
      p95    = round(quantile(mdist, 0.95), 3),
      max    = round(max(mdist),            3)
    )
  }
  
  dist_summary <- bind_rows(Filter(Negate(is.null),
                                   mapply(summarise_mdist, mdist_list, names(mdist_list), SIMPLIFY = FALSE)))
  
  cat("\n  Distance summary across ratios:\n")
  print(dist_summary)
  
  # Balance summary (trend SMDs — selection criterion)
  cat("\n  Balance: trend SMDs per ratio (select highest ratio with all < 0.10):\n")
  for (nm in names(s2_fits)) {
    m <- s2_fits[[nm]]
    if (is.null(m)) next
    bt        <- bal.tab(m, thresholds = c(m = 0.1), un = FALSE)
    smd_all   <- abs(bt$Balance$Diff.Adj)
    trend_smd <- bt$Balance %>% rownames_to_column("v") %>%
      filter(v %in% stage2_trends) %>% pull(Diff.Adj) %>% abs()
    cat(sprintf("  %s | max trend SMD: %.3f | max all SMD: %.3f | mean SMD: %.3f\n",
                nm, max(trend_smd, na.rm=T), max(smd_all, na.rm=T), mean(smd_all, na.rm=T)))
  }
  
  # Histogram across ratios
  dist_df_all <- bind_rows(Filter(Negate(is.null),
                                  mapply(function(df, nm) mutate(df, ratio = nm),
                                         mdist_list, names(mdist_list), SIMPLIFY = FALSE)))
  
  p_ratios <- ggplot(dist_df_all, aes(x = mdist, fill = ratio)) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    facet_wrap(~ ratio, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() + theme(legend.position = "none") +
    labs(title    = paste("Stage 2 distances by ratio —", analysis_label),
         subtitle = "Select highest ratio with acceptable trend SMD balance",
         x = "Mahalanobis distance", y = "Count")
  ggsave(here("output", paste0("s2_dist_ratios_", analysis_label, ".png")),
         p_ratios, width = 13, height = 5, dpi = 300)
  cat("  Ratio comparison histogram saved.\n")
  
  # ── SELECT PRIMARY RATIO ──────────────────────────────────────────────────
  primary_ratio <- "r5"   # 
  m_primary     <- s2_fits[[primary_ratio]]
  mdist_primary <- mdist_list[[primary_ratio]]
  
  # ── P95 TRIM 
  p95 <- quantile(mdist_primary$mdist, 0.95)
  
  # kept_treated: OA identifiers of treated units whose best match <= P95
  kept_treated_ids <- mdist_primary %>%
    filter(mdist <= p95) %>%
    pull(treated_OA) %>%
    unique()
  
  # primary_data from match.data() has no subclass with replace=TRUE.
  # Filter treated rows by OA, then pull their matched control OAs
  # from mdist_primary (which was built from match.matrix).
  kept_control_ids <- mdist_primary %>%
    filter(treated_OA %in% kept_treated_ids) %>%
    pull(control_OA) %>%
    unique()
  
  primary_data         <- match.data(m_primary, data = data)
  primary_data_trimmed <- primary_data %>%
    filter(
      (treat_indicator == 1 & OA %in% kept_treated_ids) |
        (treat_indicator == 0 & OA %in% kept_control_ids)
    )
  
  n_before    <- sum(primary_data$treat_indicator == 1)
  n_after     <- sum(primary_data_trimmed$treat_indicator == 1)
  trimmed_oas <- primary_data %>%
    filter(treat_indicator == 1, !OA %in% kept_treated_ids)
  
  cat("\n  P95 trim (primary ratio", primary_ratio, "):\n")
  cat("    Treated before trim:", n_before, "\n")
  cat("    P95 threshold:      ", round(p95, 2), "\n")
  cat("    Treated trimmed:    ", n_before - n_after, "\n")
  cat("    Treated retained:   ", n_after, "\n")
  cat("    Controls retained:  ",
      sum(primary_data_trimmed$treat_indicator == 0), "\n")
  
  cat("    Country distribution of trimmed treated OAs:\n")
  print(table(trimmed_oas$country))
  
  # Cross-country check via mdist_primary (has both OA identifiers)
  cross_check <- mdist_primary %>%
    filter(treated_OA %in% kept_treated_ids) %>%
    left_join(data %>% select(OA, country) %>% rename(t_country = country),
              by = c("treated_OA" = "OA")) %>%
    left_join(data %>% select(OA, country) %>% rename(c_country = country),
              by = c("control_OA" = "OA")) %>%
    filter(t_country != c_country)
  cat("    Cross-country pairs (post-trim):", nrow(cross_check),
      if (nrow(cross_check) == 0) "✓" else "WARNING", "\n")
  # Primary histogram
  p_hist <- ggplot(mdist_primary, aes(x = mdist)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = p95, colour = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = p95 + 0.3, y = Inf, vjust = 2, hjust = 0,
             colour = "red", size = 3, label = sprintf("P95 = %.1f", p95)) +
    theme_minimal() +
    labs(title    = paste("Stage 2 distances (ratio", primary_ratio, ") —", analysis_label),
         subtitle = "Red = P95 trim threshold",
         x = "Mahalanobis distance", y = "Count")
  ggsave(here("output", paste0("s2_dist_primary_", analysis_label, ".png")),
         p_hist, width = 10, height = 6, dpi = 300)
  cat("    Histogram saved.\n\n")
  
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


### Trend balance is  better at r5 than r1

# =============================================================================
# SECTION 7 — BALANCE DIAGNOSTICS
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 7: BALANCE DIAGNOSTICS\n")
cat("====================================================\n\n")

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

run_balance(s2_A, "A_excl_zero")
run_balance(s2_B, "B_incl_zero")



####### compre R1 and R5 

# Temporary r5 result objects for balance inspection
s2_A_r5 <- s2_A
s2_A_r5$matchit_obj <- s2_A$all_fits[["r5"]]
s2_A_r5$primary_ratio <- "r5"

s2_B_r5 <- s2_B
s2_B_r5$matchit_obj <- s2_B$all_fits[["r5"]]
s2_B_r5$primary_ratio <- "r5"

run_balance(s2_A_r5, "A_excl_zero_r5")
run_balance(s2_B_r5, "B_incl_zero_r5")

###    and the winner is B  R5 on every dimension that matters
####B r5 has the best trend SMDs of all four specifications (0.034) — the strongest empirical support for parallel trends
# Level balance is better at R5, meaning less work for the doubly-robust correction
# The efficiency gain from r5 is real and does not  comes at the cost of worse level balance and more covariates needed in the outcome model

# =============================================================================
# SECTION 8 — COMPARE ANALYSES A AND B
# =============================================================================

cat("\n====================================================\n")
cat("SECTION 8: COMPARISON — ANALYSES A vs B\n")
cat("====================================================\n\n")

treated_A <- s2_A$primary_data_trimmed %>% filter(treat_indicator == 1) %>% pull(OA)
treated_B <- s2_B$primary_data_trimmed %>% filter(treat_indicator == 1) %>% pull(OA)

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
  compare_grp <- s2_B$primary_data_trimmed %>%
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
  print(table(s2_B$primary_data_trimmed %>% filter(OA %in% only_in_B) %>% pull(country)))
  
  zero_dists <- s2_B$mdist %>% filter(treated_OA %in% only_in_B)
  cat("\n  Stage 2 distances for zero-injury pairs (expect near-zero):\n")
  print(summary(zero_dists$mdist))
}

# Distance comparison table
cat("\n--- Stage 2 distance distribution: A vs B ---\n")
dist_compare <- bind_rows(
  s2_A$mdist %>% mutate(analysis = "A: excl zero-injury"),
  s2_B$mdist %>% mutate(analysis = "B: incl zero-injury")
) %>% group_by(analysis) %>%
  summarise(n = n(), mean = round(mean(mdist), 3), median = round(median(mdist), 3),
            p90 = round(quantile(mdist, 0.90), 3), p95 = round(quantile(mdist, 0.95), 3),
            max = round(max(mdist), 3), .groups = "drop")
print(dist_compare)
cat("\nNOTE: B will show lower mean/median — zero-injury pairs have near-zero\n")
cat("distances by construction. Lower distance ≠ better match here.\n\n")

dist_both <- bind_rows(
  s2_A$mdist %>% mutate(analysis = "A: Excluding zero-injury OAs"),
  s2_B$mdist %>% mutate(analysis = "B: Including zero-injury OAs")
)
p_compare <- ggplot(dist_both, aes(x = mdist, fill = analysis)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~ analysis, scales = "free_y") +
  scale_fill_manual(values = c("steelblue", "darkorange")) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "Stage 2 distances: A vs B",
       subtitle = "Spike near zero in B = zero-injury pairs",
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
    select(OA, weights)
  controls <- s2_result$primary_data_trimmed %>%
    filter(treat_indicator == 0) %>%
    select(OA, weights)
  cat(label, "— Treated:", nrow(treated),
      "| Controls:", nrow(controls), "\n")
  list(treated = treated, controls = controls)
}

matched_A <- extract_matched(s2_A, "Analysis A (excl zero-injury)")
matched_B <- extract_matched(s2_B, "Analysis B (incl zero-injury)")




#####################################################################################

# checking wieght distribution 
w_check <- s2_B$primary_data_trimmed %>%
  filter(treat_indicator == 0) %>%
  summarise(
    n_controls     = n(),
    mean_weight    = mean(weights),
    max_weight     = max(weights),
    pct_weight_gt5 = mean(weights > 5) * 100,
    effective_n    = sum(weights)^2 / sum(weights^2)
  )
print(w_check)


top_controls <- s2_B$primary_data_trimmed %>%
  filter(treat_indicator == 0) %>%
  arrange(desc(weights)) %>%
  select(OA, weights, country, scheme) %>%
  head(20)
print(top_controls)




# =============================================================================
# WEIGHT DIAGNOSTICS + CAP SENSITIVITY
# =============================================================================

weight_diagnostics <- function(s2_result, analysis_label) {
  
  controls <- s2_result$primary_data_trimmed %>%
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
    treated_rows <- s2_result$primary_data_trimmed %>%
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
cat("         | P95 trim on min(matched-pair distance) per treated OA\n\n")

cat("STAGE 2  | MDM on", length(stage2_vars), "vars (7 trends + 7 levels, per km)\n")
cat("         | No caliper\n")
cat("         | replace=TRUE, ratios 1/2/5 tested\n")
cat("         | P95 trim on subclass (retains treated + matched controls)\n\n")

cat("Country exact constraint at both stages.\n\n")

cat("ANALYSIS A | Zero-injury treated OAs EXCLUDED\n")
cat("           | Ratio:", s2_A$primary_ratio,
    "| N treated:", nrow(matched_A$treated), "\n\n")

cat("ANALYSIS B | Zero-injury treated OAs INCLUDED (Stage 2 vars = 0)\n")
cat("           | Ratio:", s2_B$primary_ratio,
    "| N treated:", nrow(matched_B$treated), "\n\n")

cat("WEIGHTS: Pass weights column to att_gt(weightsname = 'weights')\n\n")

cat("NEXT STEP: Run C&S on both samples.\n")
cat("  ATT(A) ≈ ATT(B) => B is primary, A is robustness check\n")
cat("  ATT(A) ≠ ATT(B) => report both; discuss zero-injury subgroup\n")


########## final matching data #

str(s2_B_r5 )


#####    If ATT(A) ≈ ATT(B): The zero-injury OAs don't change the result. 
##    B is the primary analysis — it makes no upfront exclusion, includes all road-exposed treated OAs,#
# and the broader scope is more defensible. Report A as a robustness check in a footnote or appendix.


### If ATT(A) ≠ ATT(B) meaningfully: The zero-injury OAs are doing something different. 
# This is actually an interesting finding — it means the effect for low-exposure roads differs from high-exposure roads.
# In this case report both and discuss what the difference tells you about heterogeneity in the treatment effect.
# The one scenario where A would be definitively preferred over B is if the Stage 2 balance for B is substantially worse than for A 
#— specifically if trend SMDs in B exceed 0.10 while in A they don't. 
# That would mean the zero-injury OAs are degrading the parallel trends support. 




# =============================================================================
# TRIM SENSITIVITY: Compare P90 / P95 / P99 / No-trim across both analyses
# =============================================================================

run_trim_sensitivity <- function(s2_result, analysis_label) {
  
  cat("\n--- Trim sensitivity:", analysis_label, "---\n")
  
  mdist_primary <- s2_result$mdist
  primary_data  <- s2_result$primary_data
  data          <- s2_result$data
  
  thresholds <- c(
    "No trim" = 1.00,
    "P99"     = 0.99,
    "P95"     = 0.95,
    "P90"     = 0.90
  )
  
  results <- map_df(names(thresholds), function(trim_name) {
    
    pct <- thresholds[trim_name]
    
    if (pct == 1.00) {
      # No trim: use all matched treated OAs
      kept_treated_ids <- unique(mdist_primary$treated_OA)
    } else {
      p_thresh <- quantile(mdist_primary$mdist, pct)
      kept_treated_ids <- mdist_primary %>%
        filter(mdist <= p_thresh) %>%
        pull(treated_OA) %>%
        unique()
    }
    
    kept_control_ids <- mdist_primary %>%
      filter(treated_OA %in% kept_treated_ids) %>%
      pull(control_OA) %>%
      unique()
    
    trimmed_data <- primary_data %>%
      filter(
        (treat_indicator == 1 & OA %in% kept_treated_ids) |
          (treat_indicator == 0 & OA %in% kept_control_ids)
      )
    
    n_treated  <- sum(trimmed_data$treat_indicator == 1)
    n_controls <- sum(trimmed_data$treat_indicator == 0)
    
    # Trend SMD balance on trimmed sample
    # Use the matchit object but reweight — approximate via direct SMD calc
    treat_rows   <- trimmed_data %>% filter(treat_indicator == 1)
    control_rows <- trimmed_data %>% filter(treat_indicator == 0)
    
    trend_smds <- map_df(stage2_trends, function(v) {
      t_vals <- treat_rows[[v]]
      c_vals <- control_rows[[v]]
      pooled_sd <- sqrt((var(t_vals, na.rm=T) + var(c_vals, na.rm=T)) / 2)
      tibble(
        variable = v,
        smd = if_else(pooled_sd > 0,
                      (mean(t_vals, na.rm=T) - mean(c_vals, na.rm=T)) / pooled_sd,
                      0)
      )
    })
    
    level_smds <- map_df(stage2_levels, function(v) {
      t_vals <- treat_rows[[v]]
      c_vals <- control_rows[[v]]
      pooled_sd <- sqrt((var(t_vals, na.rm=T) + var(c_vals, na.rm=T)) / 2)
      tibble(
        variable = v,
        smd = if_else(pooled_sd > 0,
                      (mean(t_vals, na.rm=T) - mean(c_vals, na.rm=T)) / pooled_sd,
                      0)
      )
    })
    
    tibble(
      trim        = trim_name,
      n_treated   = n_treated,
      n_controls  = n_controls,
      n_trimmed   = sum(primary_data$treat_indicator == 1) - n_treated,
      max_trend_smd  = max(abs(trend_smds$smd), na.rm = TRUE),
      mean_trend_smd = mean(abs(trend_smds$smd), na.rm = TRUE),
      max_level_smd  = max(abs(level_smds$smd), na.rm = TRUE),
      mean_level_smd = mean(abs(level_smds$smd), na.rm = TRUE)
    )
  })
  
  cat("\n")
  print(results %>% mutate(across(where(is.double), ~ round(., 3))))
  
  # Plot: balance vs sample size trade-off
  p <- results %>%
    mutate(trim = factor(trim, levels = names(thresholds))) %>%
    pivot_longer(c(max_trend_smd, max_level_smd),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           max_trend_smd = "Max trend |SMD|",
                           max_level_smd = "Max level |SMD|")) %>%
    ggplot(aes(x = n_treated, y = value, colour = metric, label = trim)) +
    geom_point(size = 3) +
    geom_line(aes(group = metric)) +
    ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
    geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("Max trend |SMD|" = "#2ECC71",
                                   "Max level |SMD|" = "#E74C3C")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      title    = paste("Trim threshold sensitivity —", analysis_label),
      subtitle = "Dashed line = 0.10 SMD threshold | Trade-off: tighter trim = better balance but smaller N",
      x        = "N treated retained",
      y        = "Max |SMD|",
      colour   = NULL
    )
  
  ggsave(here("output", paste0("trim_sensitivity_", analysis_label, ".png")),
         p, width = 10, height = 6, dpi = 300)
  cat("  Sensitivity plot saved.\n")
  
  invisible(results)
}

trim_sens_A <- run_trim_sensitivity(s2_A, "A_excl_zero")
trim_sens_B <- run_trim_sensitivity(s2_B, "B_incl_zero")

### i dont need trimming 
## see code version 2 for the updated code without trimming


