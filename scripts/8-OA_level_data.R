#============================================================
# OA-Level Data 
#============================================================
# Purpose:
# Construct OA dataset defining treatment, buffer, and control
# areas.
## Classification uses majority-area (>=50%) spatial rules.
### output :  OA_level_from_polygons.rds
###           OA_subset.shp
#============================================================

library(tidyverse)
library(sf)
library(here)
library(osmdata)


select <- dplyr::select
filter <- dplyr::filter

lads_sub <- readRDS(here("data","processed","LADs_sub.rds"))

oa <- st_read(here("data","processed","shp_files","OAs_comb.shp")) %>%
  st_transform(27700) %>%
  st_make_valid()

caz <- st_read(
  here("data","processed","shp_files","CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_transform(st_crs(27700)) %>%
  st_make_valid()

st_crs(caz)
#============================================================
# Keep OAs within study LADs
#============================================================

##  keeping the LAD with the largest overlap
oa_sub <- st_join(
  oa,
  lads_sub %>% dplyr::select(LAD24CD, LAD24NM),
  join = st_intersects,
  largest = TRUE          # keeps only the LAD with greatest overlap area
) %>%
  filter(!is.na(LAD24CD))



st_crs(oa_sub)$epsg
st_crs(caz)$epsg
#============================================================
# Calculate OA area
#============================================================

oa_area <- oa_sub %>%
  mutate(total_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  dplyr:: select(OA, total_area)

# ============================================================
# OA proportion inside CAZ
# ============================================================

stopifnot(st_crs(oa_sub) == st_crs(caz))

# Dissolve multiple polygon pieces belonging to the same scheme
caz_by_scheme <- caz %>%
  group_by(scheme) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

oa_caz_intersection <- st_intersection(
  oa_sub %>% dplyr::select(OA),
  caz_by_scheme %>% dplyr::select(scheme)
) %>%
  mutate(
    int_area = as.numeric(st_area(geometry))
  ) %>%
  st_drop_geometry()

# Sum all intersected pieces before calculating the OA proportion
oa_caz_prop <- oa_caz_intersection %>%
  left_join(
    oa_area %>%
      mutate(total_area = as.numeric(total_area)),
    by = "OA"
  ) %>%
  group_by(OA, scheme) %>%
  summarise(
    int_area = sum(int_area, na.rm = TRUE),
    total_area = first(total_area),
    .groups = "drop"
  ) %>%
  mutate(
    prop_caz = int_area / total_area,
    pct_inside_caz = 100 * prop_caz
  )

# In the unlikely event that an OA intersects multiple CAZ schemes,
# associate it with the scheme containing the largest area share.
oa_caz_best <- oa_caz_prop %>%
  group_by(OA) %>%
  slice_max(prop_caz, n = 1, with_ties = FALSE) %>%
  ungroup()

# Primary majority-area treatment rule
treated_OAs <- oa_caz_best %>%
  filter(prop_caz >= 0.50) %>%
  select(OA, scheme, prop_caz, pct_inside_caz)
#============================================================
# Create 1km CAZ buffer
#============================================================
caz_buffer <- st_buffer(caz,1000) %>%
  st_difference(caz)

oa_buffer_intersection <- st_intersection(
  oa_sub %>%dplyr::  select(OA),
  caz_buffer
) %>%
  mutate(int_area = st_area(geometry)) %>%
  st_drop_geometry()

oa_buffer_prop <- oa_buffer_intersection %>%
  left_join(oa_area, by="OA") %>%
  mutate(prop_buffer = as.numeric(int_area / total_area))

buffer_OAs <- oa_buffer_prop %>%
  filter(prop_buffer >= 0.5) %>%
  pull(OA)

#============================================================
# Identify CAZ LADs
#============================================================

treated_LADs <- oa_sub %>%
  filter(OA %in% treated_OAs$OA) %>%
  distinct(LAD24CD) %>%
  pull(LAD24CD)

#============================================================
# Distance measures
# Three alternative urbanity proxies, each measuring distance
# from each OA centroid to a different reference point:
#   (a) LAD geometric centroid     — original
#   (b) LAD population-weighted centroid — Census 2021
#   (c) OS Built-Up Area centroid  — Ordnance Survey Open Data
#============================================================

oa_centroids <- st_centroid(oa_sub)

# --- (a) Distance to LAD geometric centroid (original) ---
lad_geom_centroids <- oa_sub |>
  group_by(LAD24CD) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid()

dist_citycentre <- st_distance(
  st_geometry(oa_centroids),
  lad_geom_centroids[match(oa_centroids$LAD24CD,
                           lad_geom_centroids$LAD24CD), ]$geometry,
  by_element = TRUE
) |>
  as.numeric()

# --- (b) Distance to LAD population-weighted centroid ---
# Uses Census 2021 OA population to weight OA centroids
oa_census_pop <- readRDS(here("data", "processed", "OA_matching_census.rds")) |>
  dplyr::select(OA, Total)

oa_coords <- oa_centroids |>
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  left_join(oa_census_pop, by = "OA")

lad_pwc <- oa_coords |>
  filter(!is.na(Total), Total > 0) |>
  group_by(LAD24CD) |>
  summarise(
    pwc_x = sum(x * Total) / sum(Total),
    pwc_y = sum(y * Total) / sum(Total),
    .groups = "drop"
  )

lad_pwc_sf <- st_as_sf(lad_pwc, coords = c("pwc_x", "pwc_y"), crs = 27700)

dist_pop_centre <- st_distance(
  st_geometry(oa_centroids),
  lad_pwc_sf[match(oa_centroids$LAD24CD,
                   lad_pwc_sf$LAD24CD), ] |>
    st_geometry(),
  by_element = TRUE
) |>
  as.numeric()

# --- (c) Distance to OS Built-Up Area centroid ---
# Source: OS Open Built Up Areas (Ordnance Survey, OGL v3)
# Centroids of the named BUA polygon matching each LAD's main settlement.
####   from https://osdatahub.os.uk/data/downloads/open/BuiltUpAreas
os_bua_centres <- read_csv(
  here("data", "processed", "os_bua_city_centres_BNG.csv"),
  show_col_types = FALSE    
)

os_bua_sf <- st_as_sf(os_bua_centres, coords = c("easting", "northing"), crs = 27700)

dist_BUA_centroid <- st_distance(
  st_geometry(oa_centroids),
  os_bua_sf[match(oa_centroids$LAD24CD, os_bua_sf$LAD24CD), ] |>
    st_geometry(),
  by_element = TRUE
) |>
  as.numeric()

# --- (d) Compare the three measures ---
cat("\n=== Distance measure comparison ===\n")
comparison <- tibble(
  OA        = oa_centroids$OA,
  LAD24CD   = oa_centroids$LAD24CD,
  geom      = dist_citycentre,
  pwc       = dist_pop_centre,
  os_bua    = dist_BUA_centroid
)

cat("Correlations:\n")
cat("  geom  vs pwc:    ", round(cor(comparison$geom, comparison$pwc),    3), "\n")
cat("  geom  vs os_bua: ", round(cor(comparison$geom, comparison$os_bua), 3), "\n")
cat("  pwc   vs os_bua: ", round(cor(comparison$pwc,  comparison$os_bua), 3), "\n\n")

comparison |>
  summarise(
    across(c(geom, pwc, os_bua), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")
  ) |>
  print()

#============================================================
# OA classification
####
#dups
oa_sub %>% st_drop_geometry() %>% count(OA) %>% filter(n > 1) %>% nrow()





OA_analysis <- oa_sub %>%
  st_drop_geometry() %>%
  dplyr::select(OA, LAD24CD, LAD24NM) %>%
  mutate(
        treated_OA = if_else(OA %in% treated_OAs$OA,1,0),
    
    buffer_OA = if_else(
      OA %in% buffer_OAs &
        treated_OA == 0,
      1,0
    ),
    
    control_group1_OA = if_else(
      LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    ),
    
    control_group2_OA = if_else(
      !LAD24CD %in% treated_LADs &
        treated_OA == 0 &
        buffer_OA == 0,
      1,0
    )
  )

OA_analysis$dist_citycentre <- dist_citycentre
OA_analysis$dist_pop_centre <- dist_pop_centre
OA_analysis$dist_BUA_centroid  <- dist_BUA_centroid

nrow(OA_analysis)
OA_analysis %>% distinct(OA) %>% nrow()

# Are all oa_sub OAs actually in OA_analysis?
anti_join(
  oa_sub %>% st_drop_geometry() %>% dplyr::select(OA),
  OA_analysis, by = "OA"
) %>% nrow()


OA_analysis |>
  dplyr::select(dist_citycentre, dist_pop_centre, dist_BUA_centroid) |>
  summary()

table(OA_analysis$treated_OA)

OA_analysis |>
  filter(treated_OA == 1) |>
  dplyr::select(dist_citycentre, dist_pop_centre, dist_BUA_centroid) |>
  summary()

OA_analysis <- OA_analysis %>%
  mutate(
    assignment = case_when(
      treated_OA == 1 ~ "treated",
      buffer_OA == 1 ~ "buffer",
      control_group1_OA == 1 ~ "control_same_city",
      control_group2_OA == 1 ~ "control_other_city",
      TRUE ~ NA_character_
    )
  )
nrow(treated_OAs)

OA_analysis <- OA_analysis %>%
  left_join(
    treated_OAs,
    by = "OA"
  )


OA_analysis |>
  group_by(assignment) |>
  summarise(
    median_geom   = median(dist_citycentre),
    median_pwc    = median(dist_pop_centre),
    median_os_bua = median(dist_BUA_centroid)
  )

table(OA_analysis$assignment)

#distri
OA_analysis %>%
  count(assignment) %>%
  mutate(pct = n / sum(n))

glimpse(OA_analysis)






# ============================================================
# Diagnostics for OAs near the 50% treatment threshold
# ============================================================

oa_margin_diagnostic <- oa_caz_best %>%
  mutate(
    treated_50 = prop_caz >= 0.50,
    
    # "Marginal" is defined as within 10 percentage points
    # on either side of the primary 50% threshold.
    marginal_40_60 =
      prop_caz >= 0.40 & prop_caz < 0.60,
    
    marginal_below =
      prop_caz >= 0.40 & prop_caz < 0.50,
    
    marginal_above =
      prop_caz >= 0.50 & prop_caz < 0.60,
    
    overlap_band = case_when(
      prop_caz < 0.40 ~ "Below 40%",
      prop_caz < 0.50 ~ "40% to <50%",
      prop_caz < 0.60 ~ "50% to <60%",
      prop_caz < 1.00 ~ "60% to <100%",
      TRUE ~ "100% inside"
    )
  )

# Overall numbers and percentages
overall_margin_summary <- oa_margin_diagnostic %>%
  summarise(
    intersecting_OAs = n(),
    treated_OAs = sum(treated_50),
    marginal_OAs_40_60 = sum(marginal_40_60),
    marginal_below_50 = sum(marginal_below),
    marginal_above_50 = sum(marginal_above),
    pct_marginal_among_intersecting =
      100 * marginal_OAs_40_60 / intersecting_OAs,
    pct_marginal_among_treated =
      100 * marginal_above_50 / treated_OAs
  )

cat("\n--- OAs around the 50% CAZ threshold ---\n")
print(overall_margin_summary)


scheme_margin_summary <- oa_margin_diagnostic %>%
  group_by(scheme) %>%
  summarise(
    intersecting_OAs = n(),
    treated_OAs = sum(treated_50),
    marginal_40_60 = sum(marginal_40_60),
    marginal_below_50 = sum(marginal_below),
    marginal_above_50 = sum(marginal_above),
    pct_marginal_among_intersecting =
      100 * marginal_40_60 / intersecting_OAs,
    pct_marginal_among_treated =
      if_else(
        treated_OAs > 0,
        100 * marginal_above_50 / treated_OAs,
        NA_real_
      ),
    .groups = "drop"
  )

cat("\n--- Marginal OAs by scheme ---\n")
print(scheme_margin_summary, n = Inf)



threshold_sensitivity <- map_dfr(
  c(40, 45, 50, 55, 60),
  function(threshold_pct) {
    
    cutoff <- threshold_pct / 100
    
    oa_caz_best %>%
      summarise(
        threshold_pct = threshold_pct,
        treated_OAs = sum(prop_caz >= cutoff),
        OAs_switching_vs_50 = sum(
          (prop_caz >= cutoff) != (prop_caz >= 0.50)
        )
      )
  }
)

cat("\n--- Sensitivity to alternative area thresholds ---\n")
print(threshold_sensitivity)

# Create one clean overlap record per OA
oa_overlap_lookup <- oa_caz_best %>%
  st_drop_geometry() %>%
  dplyr::select(OA, prop_caz, pct_inside_caz) %>%
  distinct(OA, .keep_all = TRUE)

# Remove overlap variables created by any previous run
OA_analysis <- OA_analysis %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "prop_caz",
        "pct_inside_caz",
        "prop_caz.x",
        "prop_caz.y",
        "pct_inside_caz.x",
        "pct_inside_caz.y"
      )
    )
  ) %>%
  left_join(
    oa_overlap_lookup,
    by = "OA"
  ) %>%
  mutate(
    prop_caz = tidyr::replace_na(prop_caz, 0),
    pct_inside_caz = tidyr::replace_na(pct_inside_caz, 0)
  )

OA_analysis %>%
  summarise(
    disagreement = sum(
      treated_OA != as.integer(prop_caz >= 0.50)
    )
  ) %>%
  print()

OA_analysis %>%
  dplyr::select(
    OA,
    assignment,
    scheme,
    prop_caz,
    pct_inside_caz
  ) %>%
  dplyr::filter(prop_caz > 0) %>%
  dplyr::arrange(dplyr::desc(prop_caz)) %>%
  tibble::as_tibble() %>%
  print(n = 30)



saveRDS(
  OA_analysis,
  here("data","processed","OA_level_from_polygons.rds")
)

english_schemes <- c(
  "Bath",
  "Birmingham",
  "Bradford",
  "Bristol",
  "Newcastle",
  "Portsmouth",
  "Sheffield"
)

english_oa_margin <- oa_margin_diagnostic %>%
  dplyr::filter(scheme %in% english_schemes) %>%
  mutate(
    distance_from_50 = abs(prop_caz - 0.50)
  )

english_margin_summary <- english_oa_margin %>%
  summarise(
    intersecting_OAs = n(),
    
    treated_OAs =
      sum(prop_caz >= 0.50),
    
    fully_inside_OAs =
      sum(prop_caz >= 0.999999),
    
    partially_inside_treated_OAs =
      sum(prop_caz >= 0.50 & prop_caz < 0.999999),
    
    marginal_OAs_40_60 =
      sum(prop_caz >= 0.40 & prop_caz < 0.60),
    
    marginal_below_50 =
      sum(prop_caz >= 0.40 & prop_caz < 0.50),
    
    marginal_above_50 =
      sum(prop_caz >= 0.50 & prop_caz < 0.60),
    
    pct_marginal_among_intersecting =
      100 * marginal_OAs_40_60 / intersecting_OAs,
    
    pct_marginal_among_treated =
      100 * marginal_above_50 / treated_OAs
  )

print(english_margin_summary)


st_write(
  oa_sub,
  here("data","processed","shp_files","OA_subset.shp"),
  delete_dsn=TRUE
)



