# ============================================================
# Merge OA Matching Dataset with Census Characteristics


library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(naniar)

### census OA data 
OA_char_raw <- read.csv(here("data","processed","outputArea_raw.csv"))
OA_char_percent <- read.csv(here("data","processed","outputArea_percent.csv"))

miss_var_summary(OA_char_raw)
miss_var_summary(OA_char_percent)


OA_char_percent <-  OA_char_percent %>% filter(!is.na(Bicycle))

### business  OA data 
OA__EW_businesses <- read.csv(here("data","processed","OA_EW_businesses.csv"))
OA__scot_businesses<- read.csv(here("data","processed","OA_Sco_businesses.csv"))


OA_EW_businesses_clean <- OA__EW_businesses %>%
  rename(OA = OA21CD)
OA_Sco_businesses_clean <- OA__scot_businesses %>%
  rename(OA = OA22) 

# Combine datasets
OA_businesses <- bind_rows(
  OA_EW_businesses_clean,
  OA_Sco_businesses_clean
)


OA_matching_dataset <- readRDS(
  here("data","processed","OA_matching_dataset.rds")
)

# 
# OA shapefile and compute area
# ------------------------------------------------------------

oa_sub <- st_read(
  here("data","processed","shp_files","OA_subset.shp"),
  quiet = TRUE
) %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  mutate(
    area_m2 = as.numeric(st_area(.)),
    area_km2 = area_m2 / 1e6
  )

# --------------- 
# checks
# -------------- 

cat("OA_matching_dataset rows:", nrow(OA_matching_dataset), "\n")
cat("OA_char_raw rows:", nrow(OA_char_raw), "\n")
cat("OA_char_percent rows:", nrow(OA_char_percent), "\n")

cat("OA_char_raw dups:",
    OA_char_raw %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")

cat("OA_char_percent dups:",
    OA_char_percent %>% count(OA) %>% filter(n > 1) %>% nrow(), "\n")
cat("OA__EW_businesses dups:",
    OA__EW_businesses %>% count(OA21CD) %>% filter(n > 1) %>% nrow(), "\n")

OA_businesses %>% count(OA) %>% filter(n > 1)%>% nrow()


OA_businesses %>% 
  count(OA, sort = TRUE)

OA_businesses %>% 
  filter(OA %in% (OA_businesses %>% count(OA) %>% filter(n > 1) %>% pull(OA))) %>%
  head(20)

OA_businesses <- OA_businesses %>% distinct()

OA_businesses %>% 
  count(OA) %>% 
  filter(n > 1)

### still 4 duplicates 
OA_businesses <- OA_businesses %>%
  group_by(OA) %>%
  summarise(
    business_retail_per_km2 = sum(business_retail_per_km2, na.rm = TRUE),
    business_accommodation_food_per_km2 = sum(business_accommodation_food_per_km2, na.rm = TRUE),
    .groups = "drop"
  )



vars_to_rename <- setdiff(
  names(OA_char_raw),
  c("OA","country","Total","IMD")
)

OA_char_raw_renamed <- OA_char_raw %>%
  rename_with(~ paste0(.x,"_n"), all_of(vars_to_rename))

names(OA_char_percent)
OA_char_pct_renamed <- OA_char_percent %>%
  rename_with(~ paste0(.x,"_pct"), all_of(vars_to_rename)) %>%
  select(-country,-Total,-IMD)

# Merge census tables
# ------------------------

OA_census <- OA_char_raw_renamed %>%
  left_join(OA_char_pct_renamed, by="OA")

# ------------------------------------------------------------
# Merge census onto OA matching dataset
# ------------------------------------------------------------

OA_matching_census <- OA_matching_dataset %>%
  left_join(OA_census, by="OA")


#'### ------- add busnesses 
OA_matching_census <- OA_matching_dataset %>%
left_join(OA_census, by = "OA") %>%
  left_join(OA_businesses, by = "OA")


# ------------------------------------------------------------
# Add OA area + population density
# ------------------------------------------------------------

OA_matching_census <- OA_matching_census %>%
  left_join(
    oa_sub %>%
      st_drop_geometry() %>%
      select(OA, area_km2),
    by = "OA"
  ) %>%
  mutate(
    pop_density = Total / area_km2,
    log_pop_density = log1p(pop_density)
  )



summary(OA_matching_census$pop_density)
# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

stopifnot(nrow(OA_matching_census) == nrow(OA_matching_dataset))

stopifnot(
  OA_matching_census %>%
    count(OA) %>%
    filter(n > 1) %>%
    nrow() == 0
)



naniar::miss_var_summary(OA_matching_census)

#### remoive the OA with nas 
OA_matching_census <- OA_matching_census %>%
  drop_na(
    White_pct, Mixed_pct, Asian_pct, Black_pct, Other_ethnicity_pct,
    X4under_pct, X5to9_pct, X10to14_pct, X15to19_pct
     )

OA_matching_census <- OA_matching_census %>%
  mutate(scheme = replace_na(scheme, "Control"))

# ============================================================
# PRE-MATCHING QA
# All checks must pass before saving. A hard stop is raised
# at the end if any fail.
# ============================================================

qa_pass <- TRUE
qa_fail <- function(msg) { cat("  FAIL \u2717", msg, "\n"); qa_pass <<- FALSE }
qa_ok   <- function(msg) { cat("  PASS \u2713", msg, "\n") }

cat("\n================================================\n")
cat("PRE-MATCHING QA\n")
cat("================================================\n")

# ── No duplicate OAs ───────────────────────────────────────────────────────
n_rows <- nrow(OA_matching_census)
n_oa   <- n_distinct(OA_matching_census$OA)
if (n_rows == n_oa) {
  qa_ok(sprintf("No duplicate OAs (%d rows)", n_rows))
} else {
  qa_fail(sprintf("%d rows but only %d distinct OA codes — duplicates present", n_rows, n_oa))
}

# ── Country  ─────────────────────────────────────────────────

country_tbl <- OA_matching_census |>
  mutate(.country = case_when(
    substr(LAD24CD, 1, 1) == "E" ~ "England",
    substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
    substr(LAD24CD, 1, 1) == "W" ~ "Wales",
    TRUE                         ~ "Unknown"
  )) |>
  count(.country)
print(country_tbl)
n_unknown <- country_tbl |> filter(.country == "Unknown") |> pull(n) |> sum()
if (n_unknown == 0) qa_ok("No unrecognised LAD24CD prefixes") else
  qa_fail(sprintf("%d OAs have unrecognised LAD24CD prefix (not E/S/W)", n_unknown))

# ── Assignment  counts by country ─────────────────────────────────────

OA_matching_census |>
  count(treated_OA, control_group1_OA, control_group2_OA, buffer_OA) |>
  arrange(treated_OA, control_group1_OA, control_group2_OA) |>
  print()

n_treated <- sum(OA_matching_census$treated_OA        == 1, na.rm = TRUE)
n_ctrl1   <- sum(OA_matching_census$control_group1_OA == 1, na.rm = TRUE)
n_ctrl2   <- sum(OA_matching_census$control_group2_OA == 1, na.rm = TRUE)
if (n_treated == 0) qa_fail("No treated OAs found")       else qa_ok(sprintf("Treated OAs:       %d", n_treated))
if (n_ctrl1   == 0) qa_fail("No control group 1 OAs")     else qa_ok(sprintf("Control group 1:   %d", n_ctrl1))
if (n_ctrl2   == 0) qa_fail("No control group 2 OAs")     else qa_ok(sprintf("Control group 2:   %d", n_ctrl2))

# injury profiles by group × country ─────────────────
# 

inj_profile <- OA_matching_census |>
  mutate(
    .country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                         ~ "Other"
    ),
    .group = case_when(
      treated_OA        == 1 ~ "Treated",
      control_group1_OA == 1 ~ "Control_1",
      control_group2_OA == 1 ~ "Control_2",
      buffer_OA         == 1 ~ "Buffer",
      TRUE                   ~ "Other"
    )
  ) |>
  group_by(.country, .group) |>
  summarise(
    n                    = n(),
    pct_zero_mean_total  = round(100 * mean(mean_total     == 0, na.rm = TRUE), 1),
    median_mean_total    = round(median(mean_total,              na.rm = TRUE), 4),
    pct_zero_trend_total = round(100 * mean(trend_total_pkm == 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(.country, .group)
print(inj_profile, n = Inf)

scot_c1 <- inj_profile |> filter(.country == "Scotland", .group == "Control_1")
if (nrow(scot_c1) == 0) {
  cat("  NOTE: No Scottish Control_1 OAs found.\n")
} else if (scot_c1$pct_zero_mean_total < 100) {
  qa_ok(sprintf(
    "Scottish Control_1 OAs: %.1f%% zero mean injuries (was 100%% before script-12 fix)",
    scot_c1$pct_zero_mean_total
  ))
} else {
  qa_fail("Scottish Control_1 OAs STILL ALL ZEROS — re-run script 12 with the fix applied")
}

eng_c1 <- inj_profile |> filter(.country == "England", .group == "Control_1")
if (nrow(eng_c1) > 0 && eng_c1$pct_zero_mean_total == 100)
  qa_fail("English Control_1 OAs are ALL ZEROS — script 12 fix may not be applied")

# Stage 1 variable completeness (treated OAs) ───────────────────────────
cat("\n[5] Stage 1 variable completeness (treated OAs)\n")

stage1_expected <- c(
  "road_density_m_km2", "road_length_km", "pct_A_road", "pct_B_road", "pct_minor_road",
  "dist_citycentre", "pop_density",
  "IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
  "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct"
)
missing_cols_s1 <- setdiff(stage1_expected, names(OA_matching_census))
if (length(missing_cols_s1) > 0)
  qa_fail(sprintf("Stage 1 columns absent from dataset: %s",
                  paste(missing_cols_s1, collapse = ", ")))

na_s1 <- OA_matching_census |>
  filter(treated_OA == 1) |>
  summarise(across(all_of(intersect(stage1_expected, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_s1) == 0) {
  qa_ok("No missing Stage 1 values in treated OAs")
} else {
  cat("  WARNING — Stage 1 NAs in treated OAs:\n"); print(na_s1)
  qa_fail("Missing Stage 1 values in treated OAs")
}

#  Stage 2 variable completeness (matching pool, roads > 0) ───────────────
cat("\n[6] Stage 2 variable completeness (matching pool, n_roads > 0)\n")

stage2_expected <- c(
  "trend_car_KSI_pkm", "trend_car_slight_pkm", "trend_cyc_KSI_pkm",
  "trend_cyc_slight_pkm", "trend_ped_KSI_pkm", "trend_ped_slight_pkm",
  "trend_other_KSI_pkm", "trend_other_slight_pkm", "trend_total_pkm",
  "mean_car_KSI_pkm", "mean_car_slight_pkm", "mean_cyc_KSI_pkm",
  "mean_cyc_slight_pkm", "mean_ped_KSI_pkm", "mean_ped_slight_pkm",
  "mean_other_KSI_pkm", "mean_other_slight_pkm", "mean_total_pkm"
)
na_s2 <- OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1,
         buffer_OA == 0, n_roads > 0) |>
  summarise(across(all_of(intersect(stage2_expected, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_s2) == 0) {
  qa_ok("No missing Stage 2 values in matching pool (OAs with roads)")
} else {
  cat("  INVESTIGATE — Stage 2 NAs in matching pool (pkm NAs unexpected for roaded OAs):\n")
  print(na_s2)
  qa_fail("Missing Stage 2 values for roaded OAs in matching pool")
}

# Road variable sanity ───────────────────────────────────────────────────


road_sanity <- OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1 | control_group2_OA == 1) |>
  summarise(
    n_negative_road_length = sum(road_length_km    < 0, na.rm = TRUE),
    n_negative_density     = sum(road_density_m_km2 < 0, na.rm = TRUE),
    n_pct_A_over_100       = sum(pct_A_road         > 100, na.rm = TRUE),
    n_zero_roads           = sum(n_roads == 0 | is.na(n_roads))
  )
print(road_sanity)
if (road_sanity$n_negative_road_length == 0 &&
    road_sanity$n_negative_density     == 0 &&
    road_sanity$n_pct_A_over_100       == 0) {
  qa_ok("Road variables have no impossible values")
} else {
  qa_fail("Impossible values in road variables — check OA_roads_clean pipeline")
}

# ── Per-km injury rate plausibility (treated, roads > 0) ──────────────────
cat("\n[8] Per-km injury rate plausibility (treated OAs with roads)\n")

pkm_sanity <- OA_matching_census |>
  filter(treated_OA == 1, n_roads > 0) |>
  summarise(
    n                    = n(),
    n_na_mean_total_pkm  = sum(is.na(mean_total_pkm)),
    n_negative_pkm       = sum(mean_total_pkm < 0,              na.rm = TRUE),
    median_mean_total    = round(median(mean_total_pkm,         na.rm = TRUE), 4),
    p99_mean_total       = round(quantile(mean_total_pkm, 0.99, na.rm = TRUE), 4),
    max_mean_total       = round(max(mean_total_pkm,           na.rm = TRUE), 4)
  )
print(pkm_sanity)
if (pkm_sanity$n_negative_pkm == 0)
  qa_ok("No negative per-km injury rates in treated OAs") else
  qa_fail("Negative per-km injury rates in treated OAs")
if (pkm_sanity$n_na_mean_total_pkm == 0)
  qa_ok("No NA mean_total_pkm for treated OAs with roads") else
  qa_fail(sprintf(
    "%d treated OAs with roads have NA mean_total_pkm", pkm_sanity$n_na_mean_total_pkm
  ))

# ──Zero-injury flag consistency ──────────────────────────────────────────
cat("\n[9] zero_injury_OA flag\n")

OA_matching_census |>
  filter(treated_OA == 1 | control_group1_OA == 1) |>
  count(treated_OA, control_group1_OA, zero_injury_OA) |>
  print()

n_zero_treated <- sum(OA_matching_census$treated_OA   == 1 &
                      OA_matching_census$zero_injury_OA == 1, na.rm = TRUE)
cat(sprintf("  Zero-injury treated OAs (excluded from Analysis A): %d / %d (%.1f%%)\n",
            n_zero_treated, n_treated, 100 * n_zero_treated / max(n_treated, 1)))
qa_ok("zero_injury_OA flag present and consistent")

# ── Census variable completeness (treated OAs) ───────────────────────────
cat("\n[10] Census variable completeness (treated OAs)\n")

census_vars <- c("IMD", "cars_none_pct", "Drive_Car_pct", "Walk_pct",
                 "Bicycle_pct", "X65plus_pct", "X5to19_pct", "X20to24_pct",
                 "pop_density")
na_census <- OA_matching_census |>
  filter(treated_OA == 1) |>
  summarise(across(all_of(intersect(census_vars, names(OA_matching_census))),
                   ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") |>
  filter(n_NA > 0)
if (nrow(na_census) == 0) qa_ok("No missing census variables in treated OAs") else {
  cat("  Missing census values in treated OAs:\n"); print(na_census)
}

# ── Final result ──────────────────────────────────────────────────────────────
cat("\n================================================\n")
if (qa_pass) {
  cat("ALL QA CHECKS PASSED — safe to save and proceed to matching\n")
} else {
  cat("ONE OR MORE CHECKS FAILED — fix issues above before matching\n")
}
cat("================================================\n\n")

stopifnot("Pre-matching QA failed — see output above" = qa_pass)

# ------------------------------------------------------------

saveRDS(
  OA_matching_census,
  here("data","processed","OA_matching_census.rds")
)

# OA_matching_census<- readRDS(here("data","processed","OA_matching_census.rds"))
names(OA_matching_census)


# ---------------
#  spatial version
# ---------------------------------

OA_matching_census_sf <- oa_sub %>%
  select(OA, geometry) %>%
  left_join(OA_matching_census, by="OA") %>%
  st_as_sf()

st_write(
  OA_matching_census_sf,
  here("data","processed","shp_files","OA_matching_census.gpkg"),
  delete_dsn = TRUE
)


