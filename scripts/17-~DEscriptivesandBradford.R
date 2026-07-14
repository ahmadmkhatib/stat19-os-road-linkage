# =============================================================================
# BRADFORD / WEST YORKSHIRE REPORTING-SHOCK DEEP DIVE
# =============================================================================
# Purpose:
#   Check whether Bradford's post-2021 casualty divergence looks like:
#     1. a West Yorkshire force-wide reporting/recording shift,
#     2. a Bradford-wide local shock,
#     3. a treated-OA-specific pattern,
#     4. a collision-volume change, casualty-composition change, or severity mix.
#
# Inputs:
#   data/processed/injuries_final.rds
#     - best raw STATS19-derived source for LAD, police_force, collision counts.
#   data/processed/injuries_with_oa.rds
#     - injury records assigned to Output Areas; needed for treated-vs-other-OA
#       Bradford checks because injuries_final.rds has no OA column.
#   data/processed/OA_matched_full_pooled.rds
#     - final pooled matched OA sample, including Bradford treated/control OAs.
#   data/processed/OA_matching_pairs_pooled.rds
#     - exact treated-control OA pairs, useful for Bradford matched-control LADs.
#
# Outputs:
#   output/diagnostics/reporting_deep_dive/*.png
#   output/diagnostics/reporting_deep_dive/*.csv
# =============================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(zoo)
library(here)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "reporting_deep_dive")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

WEST_YORKSHIRE_FORCE_CODE <- 13
BRADFORD_LAD <- "E08000032"

wy_lads <- c(
  "E08000032", # Bradford
  "E08000035", # Leeds
  "E08000034", # Kirklees
  "E08000036", # Wakefield
  "E08000033"  # Calderdale
)

# =============================================================================
# HELPERS
# =============================================================================

drop_geom_if_needed <- function(x) {
  if (inherits(x, "sf")) st_drop_geometry(x) else x
}

as_quarter <- function(x) {
  if (inherits(x, "yearqtr")) return(x)
  if (inherits(x, "Date")) return(as.yearqtr(x))
  as.yearqtr(as.Date(x))
}

index_to_2019 <- function(data, group_cols, value_col, index_col = "index") {
  value_sym <- rlang::ensym(value_col)
  index_sym <- rlang::ensym(index_col)
  
  data %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(
      base_2019 = mean(
        (!!value_sym)[
          quarter_year >= as.yearqtr("2019 Q1") &
            quarter_year <= as.yearqtr("2019 Q4")
        ],
        na.rm = TRUE
      ),
      !!index_sym := 100 * (!!value_sym) / base_2019
    ) %>%
    ungroup()
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(file.path(outdir, filename), plot,
         width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", filename)
}

find_first_col <- function(data, candidates, label) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) {
    stop(
      "Could not find ", label, " column. Available columns are:\n",
      paste(names(data), collapse = ", ")
    )
  }
  found[1]
}

theme_check <- function() {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================

injuries_raw <- readRDS(here("data", "processed", "injuries_final.rds")) %>%
  drop_geom_if_needed() %>%
  mutate(
    date = as.Date(date),
    quarter_year = as_quarter(date),
    is_west_yorks = police_force == WEST_YORKSHIRE_FORCE_CODE
  )

injuries_oa <- readRDS(here("data", "processed", "injuries_with_oa.rds")) %>%
  drop_geom_if_needed() %>%
  mutate(
    date = as.Date(date),
    quarter_year = as_quarter(date)
  )

matched_pooled <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds")) %>%
  drop_geom_if_needed()

pairs_pooled <- readRDS(here("data", "processed", "OA_matching_pairs_pooled.rds")) %>%
  drop_geom_if_needed()

cat("\n=== Loaded data ===\n")
cat("injuries_raw rows:    ", nrow(injuries_raw), "\n")
cat("injuries_oa rows:     ", nrow(injuries_oa), "\n")
cat("matched_pooled rows:  ", nrow(matched_pooled), "\n")
cat("pairs_pooled rows:    ", nrow(pairs_pooled), "\n\n")

oa_col_in_injuries <- find_first_col(
  injuries_oa,
  c("OA", "OA11CD", "OA21CD", "oa_code", "geo_code"),
  "OA"
)

cat("Using OA column in injuries_oa:", oa_col_in_injuries, "\n\n")

# =============================================================================
# 1. WEST YORKSHIRE VS ALL OTHER FORCES: CASUALTIES AND COLLISIONS
# =============================================================================

force_group_quarter <- injuries_raw %>%
  filter(!is.na(police_force)) %>%
  group_by(quarter_year, is_west_yorks) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  ) %>%
  mutate(group = if_else(is_west_yorks, "West Yorkshire", "All other forces"))

force_group_indexed <- force_group_quarter %>%
  index_to_2019(c("group"), n_casualties, "casualty_index") %>%
  index_to_2019(c("group"), n_collisions, "collision_index") %>%
  index_to_2019(c("group"), casualties_per_collision, "cpc_index")

write_csv(force_group_indexed, file.path(outdir, "01_west_yorks_vs_other_forces_indexed.csv"))

p_wy_casualty <- ggplot(force_group_indexed, aes(quarter_year, casualty_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  labs(
    title = "Raw casualty count indexed to 2019 average = 100",
    subtitle = "West Yorkshire vs. all other police forces",
    x = "Quarter", y = "Casualty count index", colour = NULL
  ) +
  theme_check()

print(p_wy_casualty)
save_plot(p_wy_casualty, "01_west_yorks_raw_casualty_index.png")

p_wy_collision <- ggplot(force_group_indexed, aes(quarter_year, collision_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  labs(
    title = "Raw collision count indexed to 2019 average = 100",
    subtitle = "West Yorkshire vs. all other police forces",
    x = "Quarter", y = "Collision count index", colour = NULL
  ) +
  theme_check()

print(p_wy_collision)
save_plot(p_wy_collision, "02_west_yorks_raw_collision_index.png")

p_wy_cpc <- ggplot(force_group_indexed, aes(quarter_year, cpc_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  labs(
    title = "Casualties per collision indexed to 2019 average = 100",
    subtitle = "If this jumps without collisions jumping, composition/severity may have changed",
    x = "Quarter", y = "Casualties per collision index", colour = NULL
  ) +
  theme_check()

print(p_wy_cpc)
save_plot(p_wy_cpc, "03_west_yorks_casualties_per_collision_index.png")

# Formal count break tests. The printed interaction term is the useful part.
force_quarter_counts <- injuries_raw %>%
  filter(!is.na(police_force)) %>%
  group_by(police_force, quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    .groups = "drop"
  ) %>%
  mutate(
    is_west_yorks = police_force == WEST_YORKSHIRE_FORCE_CODE,
    post_break = quarter_year >= as.yearqtr("2021 Q1")
  )

wy_casualty_break_test <- glm(
  n_casualties ~ factor(police_force) + factor(quarter_year) + is_west_yorks:post_break,
  data = force_quarter_counts,
  family = quasipoisson()
)

wy_collision_break_test <- glm(
  n_collisions ~ factor(police_force) + factor(quarter_year) + is_west_yorks:post_break,
  data = force_quarter_counts,
  family = quasipoisson()
)

wy_break_results <- bind_rows(
  summary(wy_casualty_break_test)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(str_detect(term, "is_west_yorks")) %>%
    mutate(outcome = "casualties"),
  summary(wy_collision_break_test)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(str_detect(term, "is_west_yorks")) %>%
    mutate(outcome = "collisions")
)

print(wy_break_results)
write_csv(wy_break_results, file.path(outdir, "02_west_yorks_break_tests.csv"))

# =============================================================================
# 2. BRADFORD VS ACTUAL MATCHED-CONTROL LADS
# =============================================================================

bradford_control_lads <- matched_pooled %>%
  filter(scheme == "Bradford", treat_indicator == 0) %>%
  distinct(LAD24CD, LAD24NM) %>%
  pull(LAD24CD)

cat("\nBradford matched-control LAD count:", length(bradford_control_lads), "\n")

bradford_vs_controls <- injuries_raw %>%
  filter(LAD24CD == BRADFORD_LAD | LAD24CD %in% bradford_control_lads) %>%
  mutate(group = if_else(LAD24CD == BRADFORD_LAD, "Bradford", "Bradford's matched-control LADs")) %>%
  group_by(quarter_year, group) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    n_lads = n_distinct(LAD24CD),
    casualties_per_lad = n_casualties / n_lads,
    collisions_per_lad = n_collisions / n_lads,
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  ) %>%
  index_to_2019(c("group"), casualties_per_lad, "casualty_index") %>%
  index_to_2019(c("group"), collisions_per_lad, "collision_index") %>%
  index_to_2019(c("group"), casualties_per_collision, "cpc_index")

write_csv(bradford_vs_controls, file.path(outdir, "03_bradford_vs_matched_control_lads_indexed.csv"))

p_bradford_cas <- ggplot(bradford_vs_controls, aes(quarter_year, casualty_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  labs(
    title = "Raw casualty count per LAD indexed to 2019 average = 100",
    subtitle = "Bradford vs. its actual matched-control LADs",
    x = "Quarter", y = "Casualty count index", colour = NULL
  ) +
  theme_check()

print(p_bradford_cas)
save_plot(p_bradford_cas, "04_bradford_vs_control_lads_casualty_index.png")

p_bradford_col <- ggplot(bradford_vs_controls, aes(quarter_year, collision_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  labs(
    title = "Raw collision count per LAD indexed to 2019 average = 100",
    subtitle = "Bradford vs. its actual matched-control LADs",
    x = "Quarter", y = "Collision count index", colour = NULL
  ) +
  theme_check()

print(p_bradford_col)
save_plot(p_bradford_col, "05_bradford_vs_control_lads_collision_index.png")

# =============================================================================
# 3. WEST YORKSHIRE LADS SIDE BY SIDE
# =============================================================================

wy_lad_quarter <- injuries_raw %>%
  filter(LAD24CD %in% wy_lads) %>%
  group_by(LAD24CD, LAD24NM, quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    .groups = "drop"
  ) %>%
  index_to_2019(c("LAD24CD", "LAD24NM"), n_casualties, "casualty_index") %>%
  index_to_2019(c("LAD24CD", "LAD24NM"), n_collisions, "collision_index")

write_csv(wy_lad_quarter, file.path(outdir, "04_west_yorkshire_lads_indexed.csv"))

p_wy_lads_cas <- ggplot(wy_lad_quarter, aes(quarter_year, casualty_index, colour = LAD24NM)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.1) +
  labs(
    title = "West Yorkshire LAD casualty counts indexed to 2019 = 100",
    x = "Quarter", y = "Casualty index", colour = NULL
  ) +
  theme_check()

print(p_wy_lads_cas)
save_plot(p_wy_lads_cas, "06_west_yorkshire_lads_casualty_index.png")

p_wy_lads_col <- ggplot(wy_lad_quarter, aes(quarter_year, collision_index, colour = LAD24NM)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.1) +
  labs(
    title = "West Yorkshire LAD collision counts indexed to 2019 = 100",
    x = "Quarter", y = "Collision index", colour = NULL
  ) +
  theme_check()

print(p_wy_lads_col)
save_plot(p_wy_lads_col, "07_west_yorkshire_lads_collision_index.png")

# =============================================================================
# 4. BRADFORD TREATED OAS VS OTHER BRADFORD OAS
# =============================================================================
# This section uses injuries_with_oa.rds because injuries_final.rds has no OA.

bradford_treated_oas <- matched_pooled %>%
  filter(scheme == "Bradford", treat_indicator == 1) %>%
  distinct(OA)

bradford_matched_control_oas <- matched_pooled %>%
  filter(scheme == "Bradford", treat_indicator == 0) %>%
  distinct(OA)

bradford_oa_quarter <- injuries_oa %>%
  filter(LAD24CD == BRADFORD_LAD) %>%
  mutate(
    OA_join = .data[[oa_col_in_injuries]],
    group = case_when(
      OA_join %in% bradford_treated_oas$OA ~ "Bradford treated OAs",
      OA_join %in% bradford_matched_control_oas$OA ~ "Bradford matched-control OAs",
      TRUE ~ "Other Bradford OAs"
    )
  ) %>%
  group_by(quarter_year, group) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  ) %>%
  filter(!is.na(group)) %>%
  index_to_2019(c("group"), n_casualties, "casualty_index") %>%
  index_to_2019(c("group"), n_collisions, "collision_index") %>%
  index_to_2019(c("group"), casualties_per_collision, "cpc_index")

write_csv(bradford_oa_quarter, file.path(outdir, "05_bradford_oa_groups_indexed.csv"))

p_bradford_oa_cas <- ggplot(bradford_oa_quarter, aes(quarter_year, casualty_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.1) +
  labs(
    title = "Bradford casualties by OA group indexed to 2019 = 100",
    subtitle = "Tests whether the post-2021 rise is treated-area-specific or Bradford-wide",
    x = "Quarter", y = "Casualty index", colour = NULL
  ) +
  theme_check()

print(p_bradford_oa_cas)
save_plot(p_bradford_oa_cas, "08_bradford_oa_groups_casualty_index.png")

p_bradford_oa_col <- ggplot(bradford_oa_quarter, aes(quarter_year, collision_index, colour = group)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.1) +
  labs(
    title = "Bradford collisions by OA group indexed to 2019 = 100",
    x = "Quarter", y = "Collision index", colour = NULL
  ) +
  theme_check()

print(p_bradford_oa_col)
save_plot(p_bradford_oa_col, "09_bradford_oa_groups_collision_index.png")

# =============================================================================
# 5. BRADFORD VS CONTROL LADS BY SEVERITY
# =============================================================================

severity_col <- intersect(
  c("casualty_severity", "casualty_severity_label", "severity"),
  names(injuries_raw)
)[1]

if (!is.na(severity_col)) {
  bradford_severity <- injuries_raw %>%
    filter(LAD24CD == BRADFORD_LAD | LAD24CD %in% bradford_control_lads) %>%
    mutate(
      group = if_else(LAD24CD == BRADFORD_LAD, "Bradford", "Bradford's matched-control LADs"),
      severity = as.character(.data[[severity_col]])
    ) %>%
    group_by(quarter_year, group, severity) %>%
    summarise(
      n_casualties = n(),
      n_lads = n_distinct(LAD24CD),
      casualties_per_lad = n_casualties / n_lads,
      .groups = "drop"
    ) %>%
    index_to_2019(c("group", "severity"), casualties_per_lad, "casualty_index")
  
  write_csv(bradford_severity, file.path(outdir, "06_bradford_vs_controls_by_severity.csv"))
  
  p_severity <- ggplot(bradford_severity, aes(quarter_year, casualty_index, colour = group)) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.0) +
    facet_wrap(~ severity, scales = "free_y") +
    labs(
      title = "Bradford vs matched-control LADs by casualty severity",
      x = "Quarter", y = "Casualty index", colour = NULL
    ) +
    theme_check()
  
  print(p_severity)
  save_plot(p_severity, "10_bradford_vs_controls_by_severity.png", width = 12, height = 7)
} else {
  cat("No casualty severity column found. Skipping severity split.\n")
}

# =============================================================================
# 6. BRADFORD VS CONTROL LADS BY CASUALTY TYPE / ROAD USER
# =============================================================================

type_col <- intersect(
  c("casualty_type1", "casualty_type", "casualty_class", "road_user_type"),
  names(injuries_raw)
)[1]

if (!is.na(type_col)) {
  top_types <- injuries_raw %>%
    filter(LAD24CD == BRADFORD_LAD | LAD24CD %in% bradford_control_lads) %>%
    count(type = as.character(.data[[type_col]]), sort = TRUE) %>%
    slice_head(n = 8) %>%
    pull(type)
  
  bradford_type <- injuries_raw %>%
    filter(LAD24CD == BRADFORD_LAD | LAD24CD %in% bradford_control_lads) %>%
    mutate(
      group = if_else(LAD24CD == BRADFORD_LAD, "Bradford", "Bradford's matched-control LADs"),
      type = as.character(.data[[type_col]])
    ) %>%
    filter(type %in% top_types) %>%
    group_by(quarter_year, group, type) %>%
    summarise(
      n_casualties = n(),
      n_lads = n_distinct(LAD24CD),
      casualties_per_lad = n_casualties / n_lads,
      .groups = "drop"
    ) %>%
    index_to_2019(c("group", "type"), casualties_per_lad, "casualty_index")
  
  write_csv(bradford_type, file.path(outdir, "07_bradford_vs_controls_by_casualty_type.csv"))
  
  p_type <- ggplot(bradford_type, aes(quarter_year, casualty_index, colour = group)) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.0) +
    facet_wrap(~ type, scales = "free_y") +
    labs(
      title = "Bradford vs matched-control LADs by casualty type",
      subtitle = "Top casualty types only",
      x = "Quarter", y = "Casualty index", colour = NULL
    ) +
    theme_check()
  
  print(p_type)
  save_plot(p_type, "11_bradford_vs_controls_by_casualty_type.png", width = 13, height = 8)
} else {
  cat("No casualty type column found. Skipping casualty-type split.\n")
}

# =============================================================================
# 7. SIMPLE POST-2021 SUMMARY TABLES
# =============================================================================

pre_post_summary <- function(data, group_cols, value_cols) {
  data %>%
    mutate(period = if_else(quarter_year >= as.yearqtr("2021 Q1"), "Post-2021", "Pre-2021")) %>%
    group_by(across(all_of(c(group_cols, "period")))) %>%
    summarise(across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(
      names_from = period,
      values_from = all_of(value_cols),
      names_glue = "{.value}_{period}"
    )
}

summary_bradford_controls <- pre_post_summary(
  bradford_vs_controls,
  group_cols = c("group"),
  value_cols = c("casualty_index", "collision_index", "cpc_index")
)

summary_wy_lads <- pre_post_summary(
  wy_lad_quarter,
  group_cols = c("LAD24CD", "LAD24NM"),
  value_cols = c("casualty_index", "collision_index")
)

summary_bradford_oa <- pre_post_summary(
  bradford_oa_quarter,
  group_cols = c("group"),
  value_cols = c("casualty_index", "collision_index", "cpc_index")
)

write_csv(summary_bradford_controls, file.path(outdir, "08_summary_bradford_vs_controls_pre_post.csv"))
write_csv(summary_wy_lads, file.path(outdir, "09_summary_west_yorkshire_lads_pre_post.csv"))
write_csv(summary_bradford_oa, file.path(outdir, "10_summary_bradford_oa_groups_pre_post.csv"))

cat("\n=== Bradford vs matched-control LADs: pre/post indexed summary ===\n")
print(summary_bradford_controls)

cat("\n=== West Yorkshire LADs: pre/post indexed summary ===\n")
print(summary_wy_lads)

cat("\n=== Bradford OA groups: pre/post indexed summary ===\n")
print(summary_bradford_oa)


# =============================================================================
# CAZ EVENT DATA DESCRIPTIVE WIDE TABLES
# =============================================================================
# Creates tables like:
#   CAZ 1 / CAZ 2 / ... / pooled CAZ / national / CAZ control groups
#   across time 1, time 2, time 3, ...
#
# Outputs:
#   caz_event_descriptive_casualties_wide.csv
#   caz_event_descriptive_collisions_wide.csv
#   caz_event_descriptive_casualties_per_collision_wide.csv
# =============================================================================

# ---- Build CAZ event dataset from your existing objects ----
# injuries_oa = event/casualty data over time
# matched_pooled = tells us scheme + treated/control OAs

oa_col_in_matched <- find_first_col(
  matched_pooled,
  c("OA", "OA11CD", "OA21CD", "oa_code", "geo_code"),
  "OA in matched_pooled"
)

caz_lookup <- matched_pooled %>%
  transmute(
    OA_join = .data[[oa_col_in_matched]],
    scheme = as.character(scheme),
    treat_indicator = as.integer(treat_indicator)
  ) %>%
  distinct()

caz_event_data <- injuries_oa %>%
  mutate(
    OA_join = .data[[oa_col_in_injuries]]
  ) %>%
  inner_join(caz_lookup, by = "OA_join") %>%
  group_by(scheme, treat_indicator, quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  )

# ---- National comparator from all injury data ----
national_event_data <- injuries_raw %>%
  group_by(quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  ) %>%
  mutate(
    scheme = "national",
    treat_indicator = NA_integer_
  )

# ---- Pooled CAZ row: all treated CAZ areas combined ----
pooled_caz_event_data <- caz_event_data %>%
  filter(treat_indicator == 1) %>%
  group_by(quarter_year) %>%
  summarise(
    n_casualties = sum(n_casualties, na.rm = TRUE),
    n_collisions = sum(n_collisions, na.rm = TRUE),
    casualties_per_collision = n_casualties / n_collisions,
    .groups = "drop"
  ) %>%
  mutate(
    scheme = "pooled CAZ",
    treat_indicator = NA_integer_
  )

# ---- Function to make the wide table ----
make_caz_descriptive_wide <- function(value_col, output_filename) {
  
  value_sym <- rlang::ensym(value_col)
  
  treated_control_rows <- caz_event_data %>%
    mutate(
      row_label = case_when(
        treat_indicator == 1 ~ scheme,
        treat_indicator == 0 ~ paste0(scheme, " control grp"),
        TRUE ~ scheme
      )
    ) %>%
    select(row_label, quarter_year, value = !!value_sym)
  
  pooled_row <- pooled_caz_event_data %>%
    transmute(
      row_label = "pooled CAZ",
      quarter_year,
      value = !!value_sym
    )
  
  national_row <- national_event_data %>%
    transmute(
      row_label = "national",
      quarter_year,
      value = !!value_sym
    )
  
  combined_long <- bind_rows(
    treated_control_rows,
    pooled_row,
    national_row
  ) %>%
    arrange(quarter_year) %>%
    mutate(
      time_label = paste0(
        "time ",
        dense_rank(quarter_year)
      )
    )
  
  row_order <- combined_long %>%
    distinct(row_label) %>%
    mutate(
      order_group = case_when(
        row_label == "pooled CAZ" ~ 2,
        row_label == "national" ~ 3,
        str_detect(row_label, "control grp") ~ 4,
        TRUE ~ 1
      )
    ) %>%
    arrange(order_group, row_label)
  
  wide_table <- combined_long %>%
    select(row_label, time_label, value) %>%
    pivot_wider(
      names_from = time_label,
      values_from = value
    ) %>%
    left_join(row_order, by = "row_label") %>%
    arrange(order_group, row_label) %>%
    select(-order_group)
  
  write_csv(
    wide_table,
    file.path(outdir, output_filename)
  )
  
  return(wide_table)
}

# ---- Create descriptive tables ----
caz_casualties_wide <- make_caz_descriptive_wide(
  n_casualties,
  "caz_event_descriptive_casualties_wide.csv"
)

caz_collisions_wide <- make_caz_descriptive_wide(
  n_collisions,
  "caz_event_descriptive_collisions_wide.csv"
)

caz_cpc_wide <- make_caz_descriptive_wide(
  casualties_per_collision,
  "caz_event_descriptive_casualties_per_collision_wide.csv"
)

cat("\n=== CAZ event descriptive table: casualties ===\n")
print(caz_casualties_wide)

cat("\n=== CAZ event descriptive table: collisions ===\n")
print(caz_collisions_wide)

cat("\n=== CAZ event descriptive table: casualties per collision ===\n")
print(caz_cpc_wide)










# =============================================================================
# BRADFORD VS WEST YORKSHIRE EXCLUDING BRADFORD + ENGLISH REGION REFERENCES
# =============================================================================
# Input:
#   data/processed/injuries_final.rds
#
# Output:
#   output/diagnostics/reporting_deep_dive/
#     bradford_vs_wy_excl_bradford_and_regions_casualty_index.png
#     bradford_vs_wy_excl_bradford_and_regions_casualty_index.csv
#
# Notes:
#   - Bradford is plotted as one line.
#   - West Yorkshire excluding Bradford combines Leeds, Kirklees, Calderdale,
#     and Wakefield into one line.
#   - Other English regions are included as reference lines using the STATS19
#     police_force code, because injuries_final.rds may not contain region names.
#   - Series are indexed to their own 2019 quarterly average = 100.
# =============================================================================

library(sf)
library(tidyverse)
library(lubridate)
library(zoo)
library(here)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "reporting_deep_dive")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

BRADFORD_LAD <- "E08000032"

wy_lads <- c(
  "E08000032", # Bradford
  "E08000035", # Leeds
  "E08000034", # Kirklees
  "E08000036", # Wakefield
  "E08000033"  # Calderdale
)

wy_excl_bradford_lads <- setdiff(wy_lads, BRADFORD_LAD)

drop_geom_if_needed <- function(x) {
  if (inherits(x, "sf")) st_drop_geometry(x) else x
}

as_quarter <- function(x) {
  if (inherits(x, "yearqtr")) return(x)
  if (inherits(x, "Date")) return(as.yearqtr(x))
  as.yearqtr(as.Date(x))
}

find_first_col <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) NA_character_ else found[1]
}

police_force_to_region <- function(police_force) {
  dplyr::case_when(
    police_force %in% c(1, 48) ~ "London",
    police_force %in% c(3, 4, 5, 6, 7) ~ "North West",
    police_force %in% c(10, 11, 17) ~ "North East",
    police_force %in% c(12, 13, 14, 16) ~ "Yorkshire and The Humber",
    police_force %in% c(20, 21, 22, 23) ~ "West Midlands",
    police_force %in% c(30, 31, 32, 33, 34) ~ "East Midlands",
    police_force %in% c(35, 36, 37, 40, 41, 42) ~ "East of England",
    police_force %in% c(43, 44, 45, 46, 47) ~ "South East",
    police_force %in% c(50, 52, 53, 54, 55) ~ "South West",
    TRUE ~ NA_character_
  )
}

index_to_2019 <- function(data, group_cols, value_col, index_col = "index") {
  value_sym <- rlang::ensym(value_col)
  index_sym <- rlang::ensym(index_col)
  
  data %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(
      base_2019 = mean(
        (!!value_sym)[
          quarter_year >= as.yearqtr("2019 Q1") &
            quarter_year <= as.yearqtr("2019 Q4")
        ],
        na.rm = TRUE
      ),
      !!index_sym := 100 * (!!value_sym) / base_2019
    ) %>%
    ungroup()
}

theme_check <- function() {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.title = element_blank()
    )
}

injuries_raw <- readRDS(here("data", "processed", "injuries_final.rds")) %>%
  drop_geom_if_needed() %>%
  mutate(
    date = as.Date(date),
    quarter_year = as_quarter(date)
  )

lad_col <- find_first_col(
  injuries_raw,
  c("LAD24CD", "local_authority_ons_district", "local_authority_district")
)

if (is.na(lad_col)) {
  stop(
    "No LAD/local authority code column found. Available columns:\n",
    paste(names(injuries_raw), collapse = ", ")
  )
}

injuries_raw <- injuries_raw %>%
  mutate(
    LAD_CODE = as.character(.data[[lad_col]]),
    region = police_force_to_region(police_force)
  )

message("Using LAD column: ", lad_col)
message("Using police_force-derived English regions for reference lines.")

bradford_quarter <- injuries_raw %>%
  filter(LAD_CODE == BRADFORD_LAD) %>%
  group_by(quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    .groups = "drop"
  ) %>%
  mutate(group = "Bradford")

wy_excl_bradford_quarter <- injuries_raw %>%
  filter(LAD_CODE %in% wy_excl_bradford_lads) %>%
  group_by(quarter_year) %>%
  summarise(
    n_casualties = n(),
    n_collisions = n_distinct(collision_index),
    n_lads = n_distinct(LAD_CODE),
    casualties_per_lad = n_casualties / n_lads,
    collisions_per_lad = n_collisions / n_lads,
    .groups = "drop"
  ) %>%
  transmute(
    quarter_year,
    group = "West Yorkshire excl. Bradford",
    n_casualties = casualties_per_lad,
    n_collisions = collisions_per_lad
  )

region_quarter <- injuries_raw %>%
  filter(
    !is.na(region),
    !LAD_CODE %in% wy_lads
  ) %>%
  group_by(quarter_year, region) %>%
  summarise(
    n_casualties_total = n(),
    n_collisions_total = n_distinct(collision_index),
    n_lads = n_distinct(LAD_CODE),
    .groups = "drop"
  ) %>%
  transmute(
    quarter_year,
    group = region,
    n_casualties = n_casualties_total / n_lads,
    n_collisions = n_collisions_total / n_lads
  )

plot_data <- bind_rows(
  bradford_quarter,
  wy_excl_bradford_quarter,
  region_quarter
) %>%
  filter(!is.na(n_casualties), !is.na(quarter_year)) %>%
  index_to_2019(c("group"), n_casualties, "casualty_index") %>%
  mutate(
    group_type = case_when(
      group == "Bradford" ~ "Bradford",
      group == "West Yorkshire excl. Bradford" ~ "West Yorkshire excluding Bradford",
      TRUE ~ "Other English regions"
    ),
    group = factor(
      group,
      levels = c(
        "Bradford",
        "West Yorkshire excl. Bradford",
        sort(unique(group[!group %in% c("Bradford", "West Yorkshire excl. Bradford")]))
      )
    )
  )

write_csv(
  plot_data,
  file.path(outdir, "bradford_vs_wy_excl_bradford_and_regions_casualty_index.csv")
)

plot_colours <- c(
  "Bradford" = "#C0392B",
  "West Yorkshire excl. Bradford" = "#006D77"
)

region_groups <- setdiff(levels(plot_data$group), names(plot_colours))
region_colours <- rep("#B0B0B0", length(region_groups))
names(region_colours) <- region_groups
plot_colours <- c(plot_colours, region_colours)

p <- ggplot(
  plot_data,
  aes(
    x = quarter_year,
    y = casualty_index,
    colour = group,
    linewidth = group_type,
    alpha = group_type
  )
) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(
    xintercept = as.numeric(as.yearqtr("2021 Q1")),
    linetype = "dotted",
    colour = "grey50"
  ) +
  geom_line() +
  geom_point(aes(size = group_type), show.legend = FALSE) +
  scale_colour_manual(values = plot_colours) +
  scale_linewidth_manual(
    values = c(
      "Bradford" = 1.15,
      "West Yorkshire excluding Bradford" = 1.05,
      "Other English regions" = 0.55
    )
  ) +
  scale_alpha_manual(
    values = c(
      "Bradford" = 1,
      "West Yorkshire excluding Bradford" = 1,
      "Other English regions" = 0.45
    )
  ) +
  scale_size_manual(
    values = c(
      "Bradford" = 1.8,
      "West Yorkshire excluding Bradford" = 1.6,
      "Other English regions" = 0.8
    )
  ) +
  labs(
    title = "Bradford compared with West Yorkshire excluding Bradford and English regions",
    subtitle = "Casualty counts per LAD, indexed to each area's 2019 quarterly average = 100",
    x = "Quarter",
    y = "Casualty index"
  ) +
  guides(linewidth = "none", alpha = "none") +
  theme_check()

print(p)

ggsave(
  file.path(outdir, "bradford_vs_wy_excl_bradford_and_regions_casualty_index.png"),
  p,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

message("Saved plot and CSV to: ", outdir)





# Does the "raw index" divergence show up in Leeds and Kirklees too
# (same force, no CAZ), or is it Bradford-specific even within its own force?

west_yorks_district_check <- injuries_raw %>%
  filter(LAD24CD %in% c("E08000032", "E08000035", "E08000034")) %>%  # Bradford, Leeds, Kirklees
  mutate(district = case_when(
    LAD24CD == "E08000032" ~ "Bradford",
    LAD24CD == "E08000035" ~ "Leeds",
    LAD24CD == "E08000034" ~ "Kirklees"
  )) %>%
  group_by(quarter_year, district) %>%
  summarise(n_casualties = n(), .groups = "drop")

baseline_wy_districts <- west_yorks_district_check %>%
  filter(quarter_year >= as.yearqtr("2019 Q1"), quarter_year <= as.yearqtr("2019 Q4")) %>%
  group_by(district) %>%
  summarise(base = mean(n_casualties), .groups = "drop")

west_yorks_district_indexed <- west_yorks_district_check %>%
  left_join(baseline_wy_districts, by = "district") %>%
  mutate(index = 100 * n_casualties / base)

ggplot(west_yorks_district_indexed, aes(x = quarter_year, y = index, colour = district)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = as.numeric(as.yearqtr("2021 Q1")), linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  labs(title = "Casualty count index, West Yorkshire districts",
       subtitle = "Bradford vs. Leeds vs. Kirklees (all same police force)",
       x = "Quarter", y = "Index (2019 = 100)", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")

#### The evidence now points fairly clearly toward a genuine West Yorkshire-wide reporting/recording shift, not a Bradford-specific artifact and not a genuine CAZ effect.




