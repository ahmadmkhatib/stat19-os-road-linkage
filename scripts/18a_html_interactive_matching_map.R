# =============================================================================
# EXPORT INTERACTIVE MATCHING MAP AS SHAREABLE HTML
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(leaflet)
  library(htmlwidgets)
  library(htmltools)
  library(here)
})

cat("============================================================\n")
cat(" EXPORTING INTERACTIVE MATCHING MAP\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("[1/5] Loading data...\n")

treated_A  <- readRDS(here("data", "processed", "OA_matched_treated_A.rds"))
controls_A <- readRDS(here("data", "processed", "OA_matched_donors_A.rds"))
pairs_A    <- readRDS(here("data", "processed", "OA_matching_pairs_A.rds"))
full_data  <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# =============================================================================
# LOAD SPATIAL FILES
# =============================================================================

cat("[2/5] Loading spatial files...\n")

oa_path  <- here("data", "processed", "shp_files", "OA_subset.shp")
lad_path <- here("data", "processed", "shp_files", "LADs_filtered.shp")

oa_raw <- st_read(oa_path, quiet = TRUE)

oa_id <- intersect(
  c("OA21CD", "OA11CD", "geo_code", "OA"),
  names(oa_raw)
)[1]

oa_sf <- oa_raw |>
  rename(OA = all_of(oa_id))

analysis_oas <- union(treated_A$OA, controls_A$OA)

oa_sf <- oa_sf |>
  filter(OA %in% analysis_oas) |>
  st_transform(4326) |>
  st_simplify(dTolerance = 0.0001, preserveTopology = TRUE)

lad_raw <- st_read(lad_path, quiet = TRUE)

lad_id <- intersect(
  c("LAD21CD","LAD22CD","LAD24CD","lad_code"),
  names(lad_raw)
)[1]

lad_nm <- intersect(
  c("LAD21NM","LAD22NM","LAD24NM","lad_name"),
  names(lad_raw)
)[1]

lad_sf <- lad_raw |>
  rename(
    LAD_CODE = all_of(lad_id),
    LAD_NAME = all_of(lad_nm)
  ) |>
  st_transform(4326)

# =============================================================================
# LOOKUP TABLES
# =============================================================================

cat("[3/5] Building lookup tables...\n")

oa_meta <- full_data |>
  select(OA, LAD24CD, any_of("LAD24NM")) |>
  distinct()

if (!"LAD24NM" %in% names(oa_meta)) {
  
  lad_names_df <- lad_sf |>
    st_drop_geometry() |>
    select(LAD_CODE, LAD_NAME) |>
    distinct()
  
  oa_meta <- oa_meta |>
    left_join(
      lad_names_df,
      by = c("LAD24CD" = "LAD_CODE")
    ) |>
    rename(LAD24NM = LAD_NAME)
}

# =============================================================================
# CREATE TREATED + CONTROL LAYERS
# =============================================================================

cat("[4/5] Preparing map layers...\n")

# ---- treated ----

treated_sf <- oa_sf |>
  filter(OA %in% treated_A$OA) |>
  left_join(oa_meta, by = "OA")

treated_sf$popup <- paste0(
  "<b>Treated OA:</b> ", treated_sf$OA,
  "<br><b>LAD:</b> ", treated_sf$LAD24NM
)

# ---- controls ----

control_sf <- oa_sf |>
  filter(OA %in% controls_A$OA) |>
  left_join(oa_meta, by = "OA")

# lookup: which treated OA(s) each control belongs to
ctrl_lookup <- pairs_A |>
  group_by(control_OA) |>
  summarise(
    matched_to = paste(unique(treated_OA), collapse = ", "),
    .groups = "drop"
  )

control_sf <- control_sf |>
  left_join(
    ctrl_lookup,
    by = c("OA" = "control_OA")
  )

control_sf$popup <- paste0(
  "<b>Control OA:</b> ", control_sf$OA,
  "<br><b>LAD:</b> ", control_sf$LAD24NM,
  "<br><b>Matched to:</b> ", control_sf$matched_to
)

# =============================================================================
# BUILD LEAFLET MAP
# =============================================================================

cat("[5/5] Building interactive map...\n")

m <- leaflet(
  options = leafletOptions(preferCanvas = TRUE)
) |>
  
  addProviderTiles(
    providers$CartoDB.Positron
  ) |>
  
  # LAD boundaries
  addPolylines(
    data = lad_sf,
    color = "#777777",
    weight = 0.7,
    opacity = 0.5,
    group = "LAD boundaries"
  ) |>
  
  # Controls
  addPolygons(
    data = control_sf,
    fillColor = "#2E6FAB",
    fillOpacity = 0.55,
    color = "#1A4F8A",
    weight = 0.8,
    popup = ~popup,
    label = ~OA,
    group = "Matched controls",
    highlightOptions = highlightOptions(
      weight = 2,
      fillOpacity = 0.9,
      bringToFront = TRUE
    )
  ) |>
  
  # Treated
  addPolygons(
    data = treated_sf,
    fillColor = "#D85A30",
    fillOpacity = 0.8,
    color = "#8B2010",
    weight = 1.5,
    popup = ~popup,
    label = ~OA,
    group = "Treated OAs",
    highlightOptions = highlightOptions(
      weight = 3,
      fillOpacity = 1,
      bringToFront = TRUE
    )
  ) |>
  
  addLayersControl(
    overlayGroups = c(
      "Treated OAs",
      "Matched controls",
      "LAD boundaries"
    ),
    options = layersControlOptions(
      collapsed = FALSE
    )
  ) |>
  
  addLegend(
    position = "bottomright",
    colors = c("#D85A30", "#2E6FAB"),
    labels = c("Treated OA", "Matched control"),
    opacity = 0.9,
    title = "Matching groups"
  ) |>
  
  setView(
    lng = -2.0,
    lat = 53.5,
    zoom = 6
  )

# =============================================================================
# SAVE HTML
# =============================================================================

output_file <- here(
  "outputs",
  "interactive_matching_map.html"
)

dir.create(here("outputs"), showWarnings = FALSE)

saveWidget(
  widget = m,
  file = output_file,
  selfcontained = TRUE
)

