# =============================================================================
#  OA-Level Matching — Two-Stage MDM Design
#  Stage 1: MDM (wide caliper) on structural + sociodemographic variables
#  Stage 2: MDM (tight caliper) on pre-treatment injury trends + levels
# =============================================================================
#
# DESIGN RATIONALE:
#
#   Stage 1 — Structural restriction (MDM)
#     Removes OAs that are too structurally dissimilar to be plausible
#     comparators, before any injury outcome data are examined.
#     Uses road network, urban form, and sociodemographic variables.
#     Wide caliper (1.0 SD) — the goal is to cut the incomparable tail,
#     not to find tight matches. A control OA far outside the structural
#     distribution of treated OAs is not a valid comparator regardless
#     of its injury trends.
#
#     
#   Stage 2 — Outcome matching (MDM)
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
#
# =============================================================================

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Road-link panel for C&S estimation — must contain:
#   road_link_id, OA, quarter_num, cohort_g, treated_link, buffer_OA,
#   outcome columns (KSI_adj, Slight_adj etc.)

# =============================================================================
#   CHECKS  
# =============================================================================

table(OA_matching_dataset$assignment)
table(OA_matching_dataset$zero_injury_OA)

OA_matching_dataset %>%
  summarise(
    total_OAs  = n_distinct(OA),
    zero_roads = sum(n_roads == 0 | is.na(n_roads)),
    pct_zero   = 100 * zero_roads / total_OAs
  ) %>%
  print()

# OAs with no roads but recorded injuries — boundary artefacts
weirdOAs <- OA_matching_dataset %>%
  filter((n_roads == 0 | is.na(n_roads)) & mean_total > 0)
cat("Zero-road OAs with injuries:", nrow(weirdOAs), "\n")
cat("  Of which treated:", sum(weirdOAs$treated_OA == 1), "\n")

# =============================================================================
# STEP 1 — DERIVE COUNTRY + PRE-MATCHING EXCLUSIONS
# =============================================================================

oa_data_final <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Unknown"
    )
  ) %>%
  filter(
    (treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1),
    buffer_OA      == 0,
    zero_injury_OA == 0,
    n_roads        >  0    # exclude zero-road OAs: _pkm = 0 distorts Stage 2 distance
  ) %>%
  mutate(
    treat_indicator = as.integer(treated_OA == 1)
  )

cat("\n=== SAMPLE AFTER EXCLUSIONS ===\n")
cat("Total OAs:",    nrow(oa_data_final), "\n")
cat("  Treated:",    sum(oa_data_final$treat_indicator == 1), "\n")
cat("  Controls:",   sum(oa_data_final$treat_indicator == 0), "\n\n")

country_tab <- oa_data_final %>%
  count(country, treat_indicator) %>%
  pivot_wider(names_from = treat_indicator, values_from = n,
              names_prefix = "treated_") %>%
  rename(n_control = treated_0, n_treated = treated_1)
print(country_tab)

# Note: Scotland is ~25% of treated but ~9% of controls — the country exact
# constraint is therefore load-bearing for Scottish CAZ schemes.

# =============================================================================
# STEP 2 — DEFINE STAGE 1 AND STAGE 2 VARIABLES
# =============================================================================
# SEPARATION PRINCIPLE:
#   No variable appears in both stages. In Stage 2, Mahalanobis distance
#   is determined entirely by injury dynamics. A good census match cannot
#   compensate for a poor trajectory match.

# --- STAGE 1: Structural + sociodemographic (MDM, wide caliper) ---
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

stage1_vars <- c(stage1_road, stage1_urban, stage1_socdem)

# --- STAGE 2: Pre-treatment injury dynamics (MDM, tight caliper) ---
# Per-km rates used throughout — consistent scale across OAs of different sizes.
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

cat("Stage 1 variables:", length(stage1_vars), "\n")
cat("  Road network:",     length(stage1_road),   "\n")
cat("  Urban form:",       length(stage1_urban),  "\n")
cat("  Sociodemographic:", length(stage1_socdem), "\n")
cat("Stage 2 variables:", length(stage2_vars), "\n")
cat("  Injury trends:",   length(stage2_trends), "\n")
cat("  Injury levels:",   length(stage2_levels), "\n")

# =============================================================================
# STEP 3 — PRE-PROCESSING: WINSORISE ON FULL DATASET (STAGE 1 VARS ONLY)
# =============================================================================
# Stage 1 variables winsorised at 1st/99th percentile of full pre-matching
# sample. Stage 2 variables are winsorised separately within the Stage 1 pool
# and based on the treated distribution only (see Step 5).

all_match_vars <- c(stage1_vars, stage2_vars)

oa_data_clean <- oa_data_final %>%
  mutate(across(
    all_of(intersect(stage1_vars, names(.))),
    ~ {
      q <- quantile(., probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

# Check all variables exist — warn and drop missing
missing_s1 <- setdiff(stage1_vars, names(oa_data_clean))
missing_s2 <- setdiff(stage2_vars, names(oa_data_clean))

if (length(missing_s1) > 0) {
  cat("WARNING — Stage 1 variables missing:", paste(missing_s1, collapse = ", "), "\n")
  stage1_vars <- intersect(stage1_vars, names(oa_data_clean))
}
names(oa_data_clean)
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
# STEP 4 — STAGE 1: MDM W ON STRUCTURAL VARIABLES
# =============================================================================
# Purpose: remove OAs too structurally dissimilar to be valid comparators.
## Country enforced as hard exact constraint.
# replace = TRUE with ratio = 10 maximises the pool passed to Stage 2.
#


###standardisation 


s1_formula <- reformulate(stage1_vars, response = "treat_indicator")


# structural variables (road + urban form only)
stage1_structural <- c(stage1_road, stage1_urban)

# apply 1.0 SD caliper ONLY to structural vars
caliper_list <- rep(1.0, length(stage1_structural))
names(caliper_list) <- stage1_structural

mdm_s1 <- tryCatch(
  matchit(
    s1_formula,
    data        = oa_data_clean,
    method      = "nearest",
    distance    = "mahalanobis",
    caliper     = caliper_list,   
    std.caliper = TRUE,
    ratio       = 10,
    replace     = TRUE,
    exact       = ~ country
  ),
  error = function(e) {
    cat("Stage 1 MDM failed:", conditionMessage(e), "\n"); NULL
  }
)

# --- Extract Stage 1 results ---
s1_matched_raw <- match.data(mdm_s1)

s1_treated <- sum(s1_matched_raw$treat_indicator == 1)
s1_control <- n_distinct(
  s1_matched_raw$OA[s1_matched_raw$treat_indicator == 0]   # unique controls
)
s1_dropped <- sum(oa_data_clean$treat_indicator == 1) - s1_treated

cat("\nStage 1 results:\n")
cat("  Treated OAs retained:          ", s1_treated, "\n")
cat("  Unique control OAs in pool:    ", s1_control, "\n")
cat("  Treated OAs dropped (caliper): ", s1_dropped,
    sprintf("(%.1f%%)\n", 100 * s1_dropped / sum(oa_data_clean$treat_indicator)))

cat("\nCountry distribution after Stage 1:\n")
print(table(s1_matched_raw$country, s1_matched_raw$treat_indicator,
            dnn = c("Country", "Treated")))

# Stage 1 balance
cat("\nStage 1 balance (structural variables):\n")
print(bal.tab(mdm_s1, thresholds = c(m = 0.1), un = TRUE))



dropped_treated <- oa_data_clean %>%
  filter(treat_indicator == 1 & !(OA %in% s1_matched_raw$OA))

nrow(dropped_treated)
dropped_treated$OA


oa_data_clean$drop_status <- ifelse(oa_data_clean$OA %in% dropped_treated$OA,
                                    "Dropped", "Retained")

pca_s1 <- prcomp(oa_data_clean[, stage1_vars], scale. = TRUE)
scores <- as.data.frame(pca_s1$x)

scores$drop_status <- ifelse(
  oa_data_clean$OA %in% dropped_treated$OA,
  "Dropped", "Retained"
)

df_plot <- data.frame(
  OA = oa_data_clean$OA,
  treat_indicator = oa_data_clean$treat_indicator,
  drop_status = case_when(
    oa_data_clean$treat_indicator == 1 & oa_data_clean$OA %in% dropped_treated$OA ~ "Treated: Dropped",
    oa_data_clean$treat_indicator == 1                                             ~ "Treated: Retained",
    TRUE                                                                           ~ "Control"
  )
)


scores <- as.data.frame(pca_s1$x)
scores <- cbind(scores, df_plot)   # safe merge

ggplot(scores, aes(PC1, PC2, colour = drop_status)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Treated: Dropped"  = "red",
    "Treated: Retained" = "blue",
    "Control"           = "grey70"
  )) +
  theme_minimal() +
  labs(title = "PCA of Structural Variables (Stage 1)",
       colour = "OA Type")



# =============================================================================
# STEP 5 — PREPARE STAGE 2 INPUT: DEDUP + WINSORISE
# =============================================================================
# 

# Deduplicate: one row per OA
# Step 5 — after dedup, before winsorising
s1_data_s2_raw <- bind_rows(
  s1_matched_raw %>% filter(treat_indicator == 1),
  s1_matched_raw %>% filter(treat_indicator == 0) %>%
    distinct(OA, .keep_all = TRUE)
) %>%
  select(-any_of(c("weights", "subclass", "distance")))  # drop Stage 1 matchit columns


cat("\n=== STAGE 2 INPUT AFTER DEDUP ===\n")
cat("Treated OAs:", sum(s1_data_s2_raw$treat_indicator == 1), "\n")
cat("Control OAs:", sum(s1_data_s2_raw$treat_indicator == 0), "\n")

# Winsorise Stage 2 variables based on TREATED distribution only
s2_vars_present <- intersect(stage2_vars, names(s1_data_s2_raw))

treated_ref <- s1_data_s2_raw %>% filter(treat_indicator == 1)

s1_data_s2 <- s1_data_s2_raw %>%
  mutate(across(
    all_of(s2_vars_present),
    ~ {
      q <- quantile(treated_ref[[cur_column()]], probs = c(0.01, 0.99), na.rm = TRUE)
      pmin(pmax(., q[1]), q[2])
    }
  ))

# =============================================================================
# STEP 6 — STAGE 2: MDM WITH TIGHT CALIPER ON INJURY DYNAMICS
# =============================================================================
# FIX: replace=FALSE throughout — each control OA used at most once per
# specification to avoid double-counting in C&S estimation.
# Country exact constraint maintained.

s2_formula <- reformulate(s2_vars_present, response = "treat_indicator")

cat("\n=== STAGE 2: MDM (tight caliper, multiple ratios) ===\n")
cat("Matching on", length(s2_vars_present), "injury variables\n")


run_s2_mdm_nocal <- function(data, formula, ratio) {
  matchit(
    formula,
    data        = data,
    method      = "nearest",
    distance    = "mahalanobis",
    ratio       = ratio,
    replace     = FALSE,
    exact       = ~ country,
    model       = T
  )
}


# Re-run Stage 2 matching on the clean data
s2_results <- list(
  r1 = run_s2_mdm_nocal(s1_data_s2, s2_formula, 1),
  r2 = run_s2_mdm_nocal(s1_data_s2, s2_formula, 2),
  r4 = run_s2_mdm_nocal(s1_data_s2, s2_formula, 4)
)

# extract primary
mdm_primary <- s2_results[["r1"]]
primary_data <- match.data(mdm_primary, data = s1_data_s2)

summary(s2_results$r1)

##### matching qualtity index 

compute_mdist <- function(m, data, vars) {
  
  md <- match.data(m, data = data, weights = "match_wt")   
  
  treat <- md %>% filter(treat_indicator == 1)
  ctrl  <- md %>% filter(treat_indicator == 0)
  
  S <- cov(treat[, vars], use = "pairwise.complete.obs")
  
  dist_list <- purrr::map_df(unique(treat$subclass), function(sc) {
    trow  <- treat %>% filter(subclass == sc)
    crows <- ctrl  %>% filter(subclass == sc)
    
    purrr::map_df(seq_len(nrow(crows)), function(i) {
      d <- mahalanobis(
        x      = as.numeric(crows[i, vars]),
        center = as.numeric(trow[vars]),
        cov    = S
      )
      tibble(
        subclass   = sc,
        treated_OA = trow$OA,
        control_OA = crows$OA[i],
        mdist      = d
      )
    })
  })
  
  return(dist_list)
}


mdist_r1 <- compute_mdist(s2_results$r1, s1_data_s2, s2_vars_present)
mdist_r2 <- compute_mdist(s2_results$r2, s1_data_s2, s2_vars_present)
mdist_r4 <- compute_mdist(s2_results$r4, s1_data_s2, s2_vars_present)



summarise_quality <- function(df) {
  df %>% summarise(
    mean_mdist = mean(mdist),
    median_mdist = median(mdist),
    p90_mdist = quantile(mdist, 0.90),
    max_mdist = max(mdist)
  )
}

summarise_quality(mdist_r1)
summarise_quality(mdist_r2)
summarise_quality(mdist_r4)

### winner: 1:1 matching ?



ggplot(mdist_r1, aes(x = mdist)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Stage 2 Match Quality (1:1) — Mahalanobis Distances",
       x = "Mahalanobis distance", y = "Count")


cat("\nStage 2 results summary:\n")
cat(sprintf("  %-25s  %8s  %8s  %8s\n", "Specification", "Treated", "Controls", "Dropped"))
cat(strrep("-", 60), "\n")

### there are some extreme distances --- ie bad matches 
# see the distribution
quantile(mdist_r1$mdist, probs = c(0.90, 0.95, 0.99))

# trim at 95th percentile
keep_pairs <- mdist_r1 %>% filter(mdist <= quantile(mdist, 0.95))

primary_data_trimmed <- primary_data %>%
  filter(subclass %in% keep_pairs$subclass)
# This removed 40 of the worst matches keeping 750+ clean ones.

for (nm in names(s2_results)) {
  m <- s2_results[[nm]]
  
  if (!inherits(m, "matchit")) {
    cat(sprintf("%-20s FAILED (not a matchit object)\n", nm))
    next
  }
  
  #  reconstruct matched data using match.data
  md <- tryCatch(
    match.data(m, data = s1_data_s2, weights = "match_weight"),
    error = function(e) NULL
  )
  
  if (is.null(md)) {
    cat(sprintf("%-20s FAILED (match.data could not reconstruct)\n", nm))
    next
  }
  
  nt <- sum(md$treat_indicator == 1)
  nc <- sum(md$treat_indicator == 0)
  ndrop <- s1_treated - nt
  
  cat(sprintf("%-20s %8d %8d %8d (%.1f%%)\n",
              nm, nt, nc, ndrop, 100 * ndrop / s1_treated))
}

###   1 to 1 is the best option 
# =============================================================================
# STEP 7 — BALANCE DIAGNOSTICS AND RATIO SELECTION
# =============================================================================
# Selection rule: choose the highest ratio at which:
#   - All Stage 2 SMDs remain below 0.10
#   - Pre-trend trajectories visually overlay (on _pkm scale — see plot below)
#   - Exclusion rate does not exceed ~25% of Stage 1 sample

run_diagnostics <- function(m_obj, label, stage2_vars) {
  if (is.null(m_obj)) return(invisible(NULL))
  
  cat("\n--- Diagnostics:", label, "---\n")
  md <- match.data(m_obj, data = s1_data_s2, weights = "match_weight")
  
  # SMD balance
  bt        <- bal.tab(m_obj, thresholds = c(m = 0.1, v = 2), un = TRUE)
  smd_after <- abs(bt$Balance$Diff.Adj)
  cat("  Variables with |SMD| < 0.10:", sum(smd_after < 0.1, na.rm = TRUE),
      "/", length(smd_after), "\n")
  cat("  Max |SMD|:",  round(max(smd_after,  na.rm = TRUE), 3), "\n")
  cat("  Mean |SMD|:", round(mean(smd_after, na.rm = TRUE), 3), "\n")
  
  # Cross-country check — must be zero
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
    stars        = "std",
    var.order    = "unadjusted",
    title        = paste("Covariate balance —", label),
    shapes       = c("circle filled", "triangle filled"),
    colors       = c("#E74C3C", "#2ECC71"),
    sample.names = c("Before matching", "After matching")
  ) +
    theme(
      axis.text.y   = element_text(size = 9),
      axis.title.x  = element_text(size = 10),
      plot.title    = element_text(size = 11),
      legend.position = "bottom",
      legend.text   = element_text(size = 9),
      plot.margin   = margin(t = 10, r = 30, b = 10, l = 180)
    )
  
  fname <- here("output", paste0("balance_", gsub("[^a-zA-Z0-9]", "_", label), ".png"))
  ggsave(fname, lp, width = 13, height =9, dpi = 300)
  cat("  Love plot saved:", fname, "\n")
  
  invisible(bt)
}

for (nm in names(s2_results)) {
  run_diagnostics(s2_results[[nm]], nm, s2_vars_present)
}

# --- Pre-trend overlay plot (on _pkm scale — consistent with Stage 2 vars) ---
# FIX: previous version plotted raw counts per OA; now plots per road-km rates
# to match the scale on which Stage 2 matching was done.

plot_pretrend <- function(m_obj, label, inj_pre, road_lengths) {
  if (is.null(m_obj)) return(invisible(NULL))
  
  md        <- match.data(m_obj, data = s1_data_s2)
  oas_in    <- md$OA
  
  trend_pkm <- inj_pre %>%
    filter(OA %in% oas_in) %>%
    left_join(road_lengths %>% select(OA, road_length_km), by = "OA") %>%
    mutate(
      total_pkm = if_else(
        !is.na(road_length_km) & road_length_km > 0,
        total_injuries / road_length_km,
        NA_real_
      ),
      group = case_when(
        treated_OA        == 1 ~ "Treated",
        control_group2_OA == 1 ~ "Control",
        TRUE                   ~ NA_character_
      )
    ) %>%
    filter(!is.na(group), !is.na(total_pkm)) %>%
    group_by(group, quarter_year) %>%
    summarise(mean_pkm = mean(total_pkm, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(trend_pkm, aes(x = quarter_year, y = mean_pkm, colour = group)) +
    geom_line() +
    geom_point(size = 1) +
    labs(
      title    = paste("Pre-treatment parallel trends —", label),
      subtitle = "Mean total casualties per road-km (matched OAs only)",
      x        = "Quarter",
      y        = "Mean casualties per road-km",
      colour   = NULL
    ) +
    theme_minimal()
  
  fname <- here("output", paste0("pretrend_", gsub("[^a-zA-Z0-9]", "_", label), ".png"))
  ggsave(fname, p, width = 10, height = 6)
  cat("  Pre-trend plot saved:", fname, "\n")
  invisible(p)
}

# Run pre-trend plots for each specification
# Requires inj_pre and road_lengths objects from the matching dataset script
# inj_pre       <- readRDS(here("data", "processed", "inj_pre.rds"))
# road_lengths  <- readRDS(here("data", "processed", "road_lengths.rds"))
# Uncomment above and run plots once those objects are available:
# for (nm in names(s2_results)) {
#   plot_pretrend(s2_results[[nm]], nm, inj_pre, road_lengths)
# }






# =============================================================================
# STEP 8 — SELECT PRIMARY SPECIFICATION AND EXTRACT MATCHED OAs
# =============================================================================
# Update primary_spec after reviewing Step 7 diagnostics.
# Default: r1_c025 (most conservative).

primary_spec <- "r1"   # <-- UPDATE after reviewing Step 7 diagnostics

mdm_primary <- s2_results[[primary_spec]]

if (is.null(mdm_primary)) stop("Primary specification failed — choose alternative.")

cat("\n=== PRIMARY SPECIFICATION:", primary_spec, "===\n")
primary_data <- match.data(mdm_primary, data = s1_data_s2)

matched_treated_oas <- primary_data_trimmed %>%
  filter(treat_indicator == 1) %>%
  dplyr::select(OA, weights, subclass)

matched_control_oas <- primary_data_trimmed %>%
  filter(treat_indicator == 0) %>%
  dplyr::select(OA, weights, subclass)

cat("Matched treated OAs:", nrow(matched_treated_oas), "\n")
cat("Matched control OAs:", nrow(matched_control_oas), "\n")

# Cross-country final check — must be zero
final_cross <- primary_data %>%
  group_by(subclass) %>%
  summarise(n_countries = n_distinct(country), .groups = "drop") %>%
  filter(n_countries > 1)
cat("Cross-country pairs in primary specification:", nrow(final_cross),
    if (nrow(final_cross) == 0) "✓" else "WARNING", "\n")





saveRDS(matched_control_oas, here("data", "processed", "OA_matched_donors.rds"))
saveRDS(matched_treated_oas, here("data", "processed", "OA_matched_treated.rds"))

