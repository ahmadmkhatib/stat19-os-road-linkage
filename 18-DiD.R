

# =============================================================================
# PRELIMINARY DiD — STAGGERED ATT ESTIMATION
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(did)
library(dplyr)
library(here)
# =============================================================================
# 
# =============================================================================

matched_full <- readRDS(here("data", "processed", "OA_matched_full_A.rds"))

matched_treated_OAs  <- matched_full %>% filter(treat_indicator == 1) %>% pull(OA)
matched_control_OAs  <- matched_full %>% filter(treat_indicator == 0) %>% pull(OA)
all_matched_OAs      <- matched_full %>% pull(OA)

cat("Matched treated OAs: ",  length(matched_treated_OAs), "\n")
cat("Matched control OAs: ",  length(matched_control_OAs), "\n")
cat("Total matched OAs:   ",  length(all_matched_OAs),     "\n")

# =============================================================================
# STEP 2 — Load road-to-OA lookup
# =============================================================================
# road_attributes already has OA attached from script 9 (road_attributes_OA.gpkg)

road_attributes <- st_read(here("data", "processed", "road_attributes_OA.gpkg")) %>%
  st_drop_geometry() %>%
  select(identifier, OA)                  # adjust OA column name if different

cat("Roads with OA linkage:", nrow(road_attributes), "\n")
cat("Missing OA:           ", sum(is.na(road_attributes$OA)), "\n")

# =============================================================================
# STEP 3 — Filter roads to matched OAs only
# =============================================================================

matched_roads <- road_attributes %>%
  filter(OA %in% all_matched_OAs)

cat("Roads in matched OAs:", n_distinct(matched_roads$identifier), "\n")
cat("OAs represented:     ", n_distinct(matched_roads$OA), "\n")

# Check coverage — any matched OAs with no roads?
missing_OA_roads <- setdiff(all_matched_OAs, matched_roads$OA)
cat("Matched OAs with no roads:", length(missing_OA_roads), "\n")

# =============================================================================
# STEP 4 — Load full panel and filter to matched roads
# =============================================================================

road_panel_full <- arrow::open_dataset(
  here("data", "processed", "road_panel_dataset")
) %>% collect()

road_panel_matched <- road_panel_full %>%
  inner_join(matched_roads, by = "identifier") %>%    # attaches OA column
  left_join(
    matched_full %>% select(OA, weights, treat_indicator, baseline_injury_stratum),
    by = "OA"
  )

cat("\nMatched panel dimensions:\n")
cat("  Rows:              ", nrow(road_panel_matched), "\n")
cat("  Unique roads:      ", n_distinct(road_panel_matched$identifier), "\n")
cat("  Unique OAs:        ", n_distinct(road_panel_matched$OA), "\n")
cat("  Treated roads:     ", n_distinct(road_panel_matched$identifier[road_panel_matched$treat_indicator == 1]), "\n")
cat("  Control roads:     ", n_distinct(road_panel_matched$identifier[road_panel_matched$treat_indicator == 0]), "\n")

# =============================================================================
# STEP 5 — Add road length for rate outcomes (injuries per road-km)
# =============================================================================
# You need this for pkm outcomes in DiD

road_lengths <- road_attributes %>%                   # or wherever length lives
  select(identifier, road_length_km)                  # adjust column name

road_panel_matched <- road_panel_matched %>%
  left_join(road_lengths, by = "identifier") %>%
  mutate(
    road_length_km = replace_na(road_length_km, 0),
    
    # Rate outcomes — per road-km (add small constant to avoid /0)
    total_pkm  = total_inj_adj_All  / (road_length_km + 1e-6),
    KSI_pkm    = KSI_adj_All        / (road_length_km + 1e-6),
    slight_pkm = Slight_adj_All     / (road_length_km + 1e-6)
  )

# =============================================================================
# STEP 6 — Treatment timing: convert to integer for att_gt
# =============================================================================
# att_gt needs gname = 0 for never-treated, integer period for treated

# Get all unique quarters as ordered factor → integer
all_quarters <- sort(unique(road_panel_matched$quarter_year))
quarter_int  <- tibble(
  quarter_year  = all_quarters,
  quarter_int   = seq_along(all_quarters)
)

panel_matched <- road_panel_matched %>%
  left_join(quarter_int, by = "quarter_year") %>%
  left_join(
    road_panel_matched %>%
      filter(treat_indicator == 1) %>%
      distinct(OA, caz_start_q) %>%
      left_join(quarter_int %>% rename(caz_start_q = quarter_year,
                                       first_treat_int = quarter_int),
                by = "caz_start_q"),
    by = "OA"
  ) %>%
  mutate(
    # Never-treated controls get 0 (required by att_gt)
    first_treat_int = replace_na(first_treat_int, 0)
  )

cat("\nTreatment timing check:\n")
panel_matched %>%
  filter(treat_indicator == 1) %>%
  distinct(scheme, caz_start_q, first_treat_int) %>%
  arrange(first_treat_int) %>%
  print()

cat("Never-treated (0):", 
    n_distinct(panel_matched$identifier[road_panel_matched$first_treat_int == 0]), 
    "roads\n")

# =============================================================================
# STEP 7 — Final integrity checks
# =============================================================================

stopifnot(
  "No NA in quarter_int"     = !anyNA(panel_matched$quarter_int),
  "No NA in first_treat_int" = !anyNA(panel_matched$first_treat_int),
  "No NA in weights"         = !anyNA(panel_matched$weights),
  "No NA in OA"              = !anyNA(panel_matched$OA)
)

cat("\nOutcome summary:\n")
panel_matched %>%
  summarise(
    mean_total  = round(mean(total_pkm,  na.rm = TRUE), 4),
    mean_KSI    = round(mean(KSI_pkm,    na.rm = TRUE), 4),
    mean_slight = round(mean(slight_pkm, na.rm = TRUE), 4),
    pct_zero    = round(mean(total_inj_adj_All == 0) * 100, 1)
  ) %>% print()

# =============================================================================
# STEP 8 — Save matched panel
# =============================================================================

saveRDS(panel_matched, here("data", "processed", "panel_matched.rds"))


# Check structure
cat("Panel rows:", nrow(panel_matched), "\n")
cat("Treated OAs:", n_distinct(panel_matched$OA[panel_matched$treat_indicator == 1]), "\n")
cat("Control OAs:", n_distinct(panel_matched$OA[panel_matched$treat_indicator == 0]), "\n")
cat("Time periods:", n_distinct(panel_matched$quarter), "\n")

# =============================================================================
# PRIMARY ATT — TOTAL INJURIES PER ROAD-KM
# =============================================================================
# xformla from balance output — log1p-transformed vars + s1 imbalanced vars

xformla <- as.formula(
  paste("~", paste(outcome_covs$analysis_A, collapse = " + "))
)
cat("Covariates in xformla:\n", paste(outcome_covs$analysis_A, collapse = "\n "), "\n\n")

# att_gt: group = first treated quarter, time = calendar quarter
out_total <- att_gt(
  yname         = "injuries_pkm",        # adjust to your outcome column
  gname         = "first_treated_quarter", # adjust to your treatment timing column
  idname        = "OA",
  tname         = "quarter",
  xformla       = xformla,
  weightsname   = "weights",
  data          = panel_matched,
  est_method    = "reg",
  clustervars   = "OA",
  control_group = "notyettreated"        # preferred for staggered adoption
)

summary(out_total)

# Aggregate to simple event-study and overall ATT
es_total  <- aggte(out_total, type = "dynamic",   na.rm = TRUE)
att_total <- aggte(out_total, type = "simple",     na.rm = TRUE)

cat("\n--- Overall ATT (total injuries/km) ---\n")
summary(att_total)

cat("\n--- Event study ---\n")
summary(es_total)

# Quick event-study plot
ggdid(es_total,
      title = "Event study — total injuries per road-km",
      xlab  = "Quarters relative to LEZ implementation",
      ylab  = "ATT (injuries/km)")

ggsave(here("output", "event_study_total.png"), width = 10, height = 6, dpi = 300)

# =============================================================================
# SECONDARY — BY SEVERITY (KSI vs SLIGHT)
# =============================================================================

run_att <- function(outcome_var, label) {
  cat("\n--- ATT:", label, "---\n")
  out <- att_gt(
    yname         = outcome_var,
    gname         = "first_treated_quarter",
    idname        = "OA",
    tname         = "quarter",
    xformla       = xformla,
    weightsname   = "weights",
    data          = panel_matched,
    est_method    = "reg",
    clustervars   = "OA",
    control_group = "notyettreated"
  )
  att   <- aggte(out, type = "simple",  na.rm = TRUE)
  es    <- aggte(out, type = "dynamic", na.rm = TRUE)
  cat("ATT estimate:", round(att$overall.att, 4),
      "| SE:", round(att$overall.se, 4),
      "| p:", round(2 * pnorm(-abs(att$overall.att / att$overall.se)), 3), "\n")
  list(att_gt = out, att = att, es = es, label = label)
}

results <- list(
  total      = run_att("injuries_pkm",       "Total injuries/km"),
  KSI        = run_att("KSI_pkm",            "KSI/km"),
  slight     = run_att("slight_pkm",         "Slight injuries/km"),
  pedestrian = run_att("ped_injuries_pkm",   "Pedestrian injuries/km"),
  cyclist    = run_att("cyc_injuries_pkm",   "Cyclist injuries/km")
)

# =============================================================================
# SUMMARY TABLE FOR ABSTRACT/PRESENTATION
# =============================================================================

summary_table <- map_df(results, function(r) {
  att <- r$att
  tibble(
    outcome  = r$label,
    ATT      = round(att$overall.att, 4),
    SE       = round(att$overall.se,  4),
    CI_lower = round(att$overall.att - 1.96 * att$overall.se, 4),
    CI_upper = round(att$overall.att + 1.96 * att$overall.se, 4),
    p_value  = round(2 * pnorm(-abs(att$overall.att / att$overall.se)), 3)
  )
})

cat("\n========================================\n")
cat("PRELIMINARY DiD RESULTS SUMMARY\n")
cat("========================================\n")
print(summary_table)