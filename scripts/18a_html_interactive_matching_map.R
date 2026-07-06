# =============================================================================
# EXPORT INTERACTIVE MATCHING MAP
# Recent pooled matching data, England only
# Includes diagnostics for large/control-area outlier OAs and long-distance pairs
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(leaflet)
  library(htmlwidgets)
  library(htmltools)
  library(here)
})

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
count  <- dplyr::count

drop_geometry_if_present <- function(x) {
  if (inherits(x, "sf")) {
    st_drop_geometry(x)
  } else {
    x
  }
}

cat("============================================================\n")
cat(" EXPORTING INTERACTIVE MATCHING MAP: RECENT POOLED, ENGLAND ONLY\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD RECENT MATCHING DATA
# =============================================================================

cat("[1/6] Loading recent pooled matching data...\n")

treated_pooled  <- readRDS(here("data", "processed", "OA_matched_treated_pooled.rds"))
controls_pooled <- readRDS(here("data", "processed", "OA_matched_donors_pooled.rds"))
pairs_pooled    <- readRDS(here("data", "processed", "OA_matching_pairs_pooled.rds"))
full_data       <- readRDS(here("data", "processed", "OA_matching_census.rds"))

pair_treated_col <- case_when(
  "treated_OA" %in% names(pairs_pooled) ~ "treated_OA",
  "treated_oa" %in% names(pairs_pooled) ~ "treated_oa",
  "treated"    %in% names(pairs_pooled) ~ "treated",
  "OA_treated" %in% names(pairs_pooled) ~ "OA_treated",
  TRUE ~ NA_character_
)

pair_control_col <- case_when(
  "control_OA" %in% names(pairs_pooled) ~ "control_OA",
  "control_oa" %in% names(pairs_pooled) ~ "control_oa",
  "control"    %in% names(pairs_pooled) ~ "control",
  "OA_control" %in% names(pairs_pooled) ~ "OA_control",
  TRUE ~ NA_character_
)

if (is.na(pair_treated_col) || is.na(pair_control_col)) {
  stop(
    "Could not identify treated/control OA columns in OA_matching_pairs_pooled.rds. ",
    "Found columns: ", paste(names(pairs_pooled), collapse = ", ")
  )
}

pairs_pooled <- pairs_pooled %>%
  transmute(
    scheme = as.character(scheme),
    treated_OA = as.character(.data[[pair_treated_col]]),
    control_OA = as.character(.data[[pair_control_col]])
  )

treated_pooled <- treated_pooled %>% mutate(OA = as.character(OA))
controls_pooled <- controls_pooled %>% mutate(OA = as.character(OA))
full_data <- full_data %>% mutate(OA = as.character(OA))
full_data_tab <- drop_geometry_if_present(full_data)

if ("LAD24CD" %in% names(full_data_tab)) {
  full_data_tab <- full_data_tab %>% filter(substr(LAD24CD, 1, 1) == "E")
  full_data <- full_data %>% filter(OA %in% full_data_tab$OA)
}

treated_pooled <- treated_pooled %>% filter(OA %in% full_data_tab$OA)
controls_pooled <- controls_pooled %>% filter(OA %in% full_data_tab$OA)

pairs_pooled <- pairs_pooled %>%
  filter(treated_OA %in% treated_pooled$OA, control_OA %in% controls_pooled$OA)

cat("  Treated OAs:", n_distinct(treated_pooled$OA), "\n")
cat("  Control OAs:", n_distinct(controls_pooled$OA), "\n")
cat("  Matched pairs:", nrow(pairs_pooled), "\n\n")

# =============================================================================
# LOAD SPATIAL FILES
# =============================================================================

cat("[2/6] Loading spatial files...\n")

oa_raw <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
)

oa_id <- intersect(c("OA21CD", "OA11CD", "geo_code", "OA"), names(oa_raw))[1]

if (is.na(oa_id)) {
  stop("No OA identifier found in OA_subset.shp.")
}

oa_sf_full <- oa_raw %>%
  rename(OA = all_of(oa_id)) %>%
  mutate(OA = as.character(OA)) %>%
  filter(OA %in% full_data_tab$OA)

# Use British National Grid for area/distance diagnostics, then WGS84 for leaflet.
oa_area_diag <- oa_sf_full %>%
  st_transform(27700) %>%
  mutate(map_geom_area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  select(OA, map_geom_area_km2)

oa_sf <- oa_sf_full %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  st_simplify(dTolerance = 0.0001, preserveTopology = TRUE)

lad_raw <- st_read(
  here("data", "processed", "shp_files", "LADs_filtered.shp"),
  quiet = TRUE
)

lad_id <- intersect(c("LAD21CD", "LAD22CD", "LAD24CD", "lad_code"), names(lad_raw))[1]
lad_nm <- intersect(c("LAD21NM", "LAD22NM", "LAD24NM", "lad_name"), names(lad_raw))[1]

if (is.na(lad_id) || is.na(lad_nm)) {
  stop("Could not identify LAD code/name columns in LADs_filtered.shp.")
}

lad_sf <- lad_raw %>%
  rename(
    LAD_CODE = all_of(lad_id),
    LAD_NAME = all_of(lad_nm)
  ) %>%
  filter(substr(LAD_CODE, 1, 1) == "E") %>%
  st_transform(4326)

# =============================================================================
# LOOKUP TABLES
# =============================================================================

cat("[3/6] Building lookup tables...\n")

oa_meta <- full_data %>%
  drop_geometry_if_present() %>%
  select(
    OA,
    any_of(c(
      "LAD24CD", "LAD24NM", "scheme",
      "area_km2", "n_roads", "road_length_km",
      "road_density_m_km2", "pop_density",
      "mean_total_pkm", "trend_total_pkm"
    ))
  ) %>%
  distinct() %>%
  left_join(oa_area_diag, by = "OA")

if (!"LAD24NM" %in% names(oa_meta) && "LAD24CD" %in% names(oa_meta)) {
  lad_names_df <- lad_sf %>%
    st_drop_geometry() %>%
    select(LAD_CODE, LAD_NAME) %>%
    distinct()
  
  oa_meta <- oa_meta %>%
    left_join(lad_names_df, by = c("LAD24CD" = "LAD_CODE")) %>%
    rename(LAD24NM = LAD_NAME)
}

if (!"LAD24NM" %in% names(oa_meta)) {
  oa_meta <- oa_meta %>% mutate(LAD24NM = NA_character_)
}

treated_scheme_lookup <- pairs_pooled %>%
  distinct(treated_OA, scheme) %>%
  group_by(treated_OA) %>%
  summarise(scheme = paste(sort(unique(scheme)), collapse = "; "), .groups = "drop")

control_scheme_lookup <- pairs_pooled %>%
  distinct(control_OA, scheme) %>%
  group_by(control_OA) %>%
  summarise(scheme = paste(sort(unique(scheme)), collapse = "; "), .groups = "drop")

ctrl_lookup <- pairs_pooled %>%
  group_by(control_OA) %>%
  summarise(
    matched_to = paste(sort(unique(treated_OA)), collapse = ", "),
    matched_schemes = paste(sort(unique(scheme)), collapse = "; "),
    n_pair_rows = n(),
    n_treated_matched = n_distinct(treated_OA),
    .groups = "drop"
  )

# =============================================================================
# PAIR DISTANCE DIAGNOSTICS
# =============================================================================

cat("[4/6] Building area and distance diagnostics...\n")

oa_centroids_sf <- oa_sf_full %>%
  st_transform(27700) %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  select(OA, geometry)

treated_centroids <- tibble(
  treated_OA = oa_centroids_sf$OA,
  treated_geometry = st_geometry(oa_centroids_sf)
)

control_centroids <- tibble(
  control_OA = oa_centroids_sf$OA,
  control_geometry = st_geometry(oa_centroids_sf)
)

pair_distance_diag <- pairs_pooled %>%
  left_join(treated_centroids, by = "treated_OA") %>%
  left_join(control_centroids, by = "control_OA") %>%
  mutate(distance_km = as.numeric(st_distance(treated_geometry, control_geometry, by_element = TRUE)) / 1000) %>%
  select(scheme, treated_OA, control_OA, distance_km)

control_diag <- pairs_pooled %>%
  left_join(
    oa_meta %>%
      select(any_of(c(
        "OA", "LAD24CD", "LAD24NM", "area_km2", "map_geom_area_km2",
        "n_roads", "road_length_km", "road_density_m_km2", "pop_density",
        "mean_total_pkm", "trend_total_pkm"
      ))),
    by = c("control_OA" = "OA")
  ) %>%
  left_join(pair_distance_diag, by = c("scheme", "treated_OA", "control_OA")) %>%
  group_by(control_OA) %>%
  summarise(
    matched_schemes = paste(sort(unique(scheme)), collapse = "; "),
    matched_to = paste(sort(unique(treated_OA)), collapse = ", "),
    n_pair_rows = n(),
    n_treated_matched = n_distinct(treated_OA),
    max_distance_km = max(distance_km, na.rm = TRUE),
    mean_distance_km = mean(distance_km, na.rm = TRUE),
    LAD24CD = first(LAD24CD),
    LAD24NM = first(LAD24NM),
    area_km2 = first(area_km2),
    map_geom_area_km2 = first(map_geom_area_km2),
    n_roads = first(n_roads),
    road_length_km = first(road_length_km),
    road_density_m_km2 = first(road_density_m_km2),
    pop_density = first(pop_density),
    mean_total_pkm = first(mean_total_pkm),
    trend_total_pkm = first(trend_total_pkm),
    .groups = "drop"
  )

area_cutoff_99 <- quantile(control_diag$map_geom_area_km2, 0.99, na.rm = TRUE)
distance_cutoff_99 <- quantile(pair_distance_diag$distance_km, 0.99, na.rm = TRUE)

control_diag <- control_diag %>%
  mutate(
    extreme_area_control = replace_na(map_geom_area_km2 >= area_cutoff_99, FALSE),
    extreme_distance_control = replace_na(max_distance_km >= distance_cutoff_99, FALSE)
  )

control_popup_lookup <- control_diag %>%
  transmute(
    control_OA,
    diagnostic_popup = paste0(
      "<br><b>Pair rows:</b> ", n_pair_rows,
      "<br><b>Unique treated OAs:</b> ", n_treated_matched,
      "<br><b>Max centroid distance:</b> ", round(max_distance_km, 1), " km",
      "<br><b>Map geometry area:</b> ", round(map_geom_area_km2, 2), " km2",
      if_else(extreme_area_control, "<br><b>Area flag:</b> top 1% matched controls", ""),
      if_else(extreme_distance_control, "<br><b>Distance flag:</b> top 1% matched controls", "")
    )
  )

cat("  Area top 1% cutoff:", round(area_cutoff_99, 2), "km2\n")
cat("  Distance top 1% cutoff:", round(distance_cutoff_99, 1), "km\n")
cat("  E00036468 diagnostic:\n")
print(control_diag %>% filter(control_OA == "E00036468"), width = 220)

# =============================================================================
# CREATE SPATIAL OBJECTS
# =============================================================================

cat("[5/6] Preparing spatial layers...\n")

treated_sf <- oa_sf %>%
  filter(OA %in% treated_pooled$OA) %>%
  left_join(oa_meta, by = "OA") %>%
  left_join(treated_scheme_lookup, by = c("OA" = "treated_OA"), suffix = c("_meta", "_matched")) %>%
  mutate(
    scheme = coalesce(scheme_matched, scheme_meta),
    scheme = as.character(scheme)
  ) %>%
  select(-any_of(c("scheme_meta", "scheme_matched")))

control_sf <- oa_sf %>%
  filter(OA %in% controls_pooled$OA) %>%
  left_join(oa_meta, by = "OA") %>%
  left_join(control_scheme_lookup, by = c("OA" = "control_OA"), suffix = c("_meta", "_matched")) %>%
  mutate(
    scheme = coalesce(scheme_matched, scheme_meta),
    scheme = as.character(scheme)
  ) %>%
  select(-any_of(c("scheme_meta", "scheme_matched"))) %>%
  left_join(ctrl_lookup, by = c("OA" = "control_OA")) %>%
  left_join(control_popup_lookup, by = c("OA" = "control_OA"))

treated_lad_popup <- rep("", nrow(treated_sf))
if ("LAD24NM" %in% names(treated_sf)) {
  treated_lad_popup <- paste0(
    "<br><b>LAD:</b> ",
    htmlEscape(coalesce(treated_sf$LAD24NM, ""))
  )
}

control_lad_popup <- rep("", nrow(control_sf))
if ("LAD24NM" %in% names(control_sf)) {
  control_lad_popup <- paste0(
    "<br><b>LAD:</b> ",
    htmlEscape(coalesce(control_sf$LAD24NM, ""))
  )
}

treated_sf$popup <- paste0(
  "<b>Treated OA:</b> ", htmlEscape(treated_sf$OA),
  "<br><b>Scheme:</b> ", htmlEscape(treated_sf$scheme),
  treated_lad_popup
)

control_sf$popup <- paste0(
  "<b>Control OA:</b> ", htmlEscape(control_sf$OA),
  "<br><b>Matched scheme(s):</b> ", htmlEscape(coalesce(control_sf$matched_schemes, "")),
  "<br><b>Matched to treated OA(s):</b> ", htmlEscape(coalesce(control_sf$matched_to, "")),
  control_lad_popup,
  coalesce(control_sf$diagnostic_popup, "")
)

# =============================================================================
# BUILD MAP
# =============================================================================

cat("[6/6] Building interactive map...\n")

m <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(
    data = lad_sf,
    color = "#777777",
    weight = 0.7,
    opacity = 0.5,
    group = "LAD boundaries"
  )

schemes <- pairs_pooled %>%
  distinct(scheme) %>%
  arrange(scheme) %>%
  pull(scheme)

valid_layers <- c()

for (sc in schemes) {
  cat("Processing scheme:", sc, "\n")
  
  tr_oas <- pairs_pooled %>%
    filter(scheme == sc) %>%
    pull(treated_OA) %>%
    unique()
  
  ctrl_oas <- pairs_pooled %>%
    filter(scheme == sc) %>%
    pull(control_OA) %>%
    unique()
  
  tr_sf <- treated_sf %>% filter(OA %in% tr_oas)
  ctrl_sf <- control_sf %>% filter(OA %in% ctrl_oas)
  
  if (nrow(tr_sf) == 0 && nrow(ctrl_sf) == 0) next
  
  layer_name <- paste0("Scheme: ", sc)
  valid_layers <- c(valid_layers, layer_name)
  
  if (nrow(ctrl_sf) > 0) {
    m <- m %>%
      addPolygons(
        data = ctrl_sf,
        fillColor = ifelse(ctrl_sf$OA %in% control_diag$control_OA[control_diag$extreme_area_control],
                           "#6A3D9A", "#1F78B4"),
        fillOpacity = 0.72,
        color = "#0B3C5D",
        weight = 1.1,
        opacity = 1,
        popup = ctrl_sf$popup,
        label = ctrl_sf$OA,
        group = layer_name,
        highlightOptions = highlightOptions(
          weight = 2.5,
          fillOpacity = 0.95,
          bringToFront = TRUE
        )
      )
  }
  
  if (nrow(tr_sf) > 0) {
    m <- m %>%
      addPolygons(
        data = tr_sf,
        fillColor = "#D85A30",
        fillOpacity = 0.9,
        color = "#8B2010",
        weight = 2,
        popup = tr_sf$popup,
        label = tr_sf$OA,
        group = layer_name,
        highlightOptions = highlightOptions(
          weight = 3,
          fillOpacity = 1,
          bringToFront = TRUE
        )
      )
  }
}

m <- m %>%
  addLayersControl(
    overlayGroups = c("LAD boundaries", valid_layers),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("#D85A30", "#1F78B4", "#6A3D9A"),
    labels = c("Treated OA", "Matched control OA", "Matched control, top 1% area"),
    title = "Recent matching groups",
    opacity = 0.9
  ) %>%
  setView(lng = -2.0, lat = 53.2, zoom = 6)

# =============================================================================
# SAVE
# =============================================================================

output_file <- here("outputs", "interactive_matching_map_recent_england.html")
summary_file <- here("outputs", "interactive_matching_map_recent_england_summary.csv")
control_diag_file <- here("outputs", "interactive_matching_control_diagnostics.csv")
pair_distance_file <- here("outputs", "interactive_matching_pair_distances.csv")

dir.create(here("outputs"), showWarnings = FALSE)

saveWidget(
  widget = m,
  file = output_file,
  selfcontained = TRUE
)

map_summary <- pairs_pooled %>%
  group_by(scheme) %>%
  summarise(
    treated_oas = n_distinct(treated_OA),
    control_oas = n_distinct(control_OA),
    matched_pairs = n(),
    .groups = "drop"
  ) %>%
  arrange(scheme)

write_csv(map_summary, summary_file)
write_csv(control_diag %>% arrange(desc(map_geom_area_km2)), control_diag_file)
write_csv(pair_distance_diag %>% arrange(desc(distance_km)), pair_distance_file)

cat("\nMAP SAVED:\n", output_file, "\n")
cat("SUMMARY SAVED:\n", summary_file, "\n")
cat("CONTROL DIAGNOSTICS SAVED:\n", control_diag_file, "\n")





cat("PAIR DISTANCES SAVED:\n", pair_distance_file, "\n")





# =============================================================================
# DIAGNOSTIC: FLAG WEIRD / EXTREME MATCHED CONTROL OAs
# =============================================================================

weird_area_cutoff_99 <- quantile(
  control_diag$map_geom_area_km2,
  probs = 0.99,
  na.rm = TRUE
)

weird_distance_cutoff_99 <- quantile(
  control_diag$max_distance_km,
  probs = 0.99,
  na.rm = TRUE
)

weird_reuse_cutoff_99 <- quantile(
  control_diag$n_pair_rows,
  probs = 0.99,
  na.rm = TRUE
)

control_diag <- control_diag %>%
  mutate(
    weird_area_99 = map_geom_area_km2 >= weird_area_cutoff_99,
    weird_distance_99 = max_distance_km >= weird_distance_cutoff_99,
    weird_reuse_99 = n_pair_rows >= weird_reuse_cutoff_99,
    weird_any_99 = weird_area_99 | weird_distance_99 | weird_reuse_99,
    weird_reason = case_when(
      weird_area_99 & weird_distance_99 & weird_reuse_99 ~
        "top 1% area + distance + reuse",
      weird_area_99 & weird_distance_99 ~
        "top 1% area + distance",
      weird_area_99 & weird_reuse_99 ~
        "top 1% area + reuse",
      weird_distance_99 & weird_reuse_99 ~
        "top 1% distance + reuse",
      weird_area_99 ~
        "top 1% area",
      weird_distance_99 ~
        "top 1% distance",
      weird_reuse_99 ~
        "top 1% reuse",
      TRUE ~ ""
    )
  )

weird_controls_99 <- control_diag %>%
  filter(weird_any_99) %>%
  arrange(desc(weird_area_99), desc(map_geom_area_km2), desc(max_distance_km)) %>%
  select(any_of(c(
    "control_OA", "matched_schemes", "matched_to",
    "LAD24CD", "LAD24NM",
    "weird_reason",
    "map_geom_area_km2", "area_km2",
    "max_distance_km", "mean_distance_km",
    "n_pair_rows", "n_treated_matched",
    "n_roads", "road_length_km",
    "road_density_m_km2", "pop_density",
    "mean_total_pkm", "trend_total_pkm"
  )))

cat("\n=== WEIRD MATCHED CONTROL OA DIAGNOSTICS ===\n")
cat("Top 1% area cutoff:", round(weird_area_cutoff_99, 3), "km2\n")
cat("Top 1% distance cutoff:", round(weird_distance_cutoff_99, 1), "km\n")
cat("Top 1% reuse cutoff:", weird_reuse_cutoff_99, "pair rows\n")
cat("Number of weird controls:", nrow(weird_controls_99), "\n\n")

print(weird_controls_99, n = 100, width = 220)


