# =============================================================================
# EXPORT INTERACTIVE MATCHING MAP (SCHEME LEVEL - CORRECT MATCH LOGIC)
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

oa_raw <- st_read(
  here("data", "processed", "shp_files", "OA_subset.shp"),
  quiet = TRUE
)

oa_id <- intersect(
  c("OA21CD", "OA11CD", "geo_code", "OA"),
  names(oa_raw)
)[1]

oa_sf <- oa_raw |>
  rename(OA = all_of(oa_id)) |>
  st_transform(4326) |>
  st_simplify(dTolerance = 0.0001, preserveTopology = TRUE)

lad_raw <- st_read(
  here("data", "processed", "shp_files", "LADs_filtered.shp"),
  quiet = TRUE
)

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
  select(OA, LAD24CD, any_of("LAD24NM"), scheme) |>
  distinct()

# fallback if LAD name missing
if (!"LAD24NM" %in% names(oa_meta)) {
  
  lad_names_df <- lad_sf |>
    st_drop_geometry() |>
    select(LAD_CODE, LAD_NAME) |>
    distinct()
  
  oa_meta <- oa_meta |>
    left_join(lad_names_df, by = c("LAD24CD" = "LAD_CODE")) |>
    rename(LAD24NM = LAD_NAME)
}

# =============================================================================
# CREATE SPATIAL OBJECTS
# =============================================================================

cat("[4/5] Preparing spatial layers...\n")

treated_sf <- oa_sf |>
  filter(OA %in% treated_A$OA) |>
  left_join(oa_meta, by = "OA")

control_sf <- oa_sf |>
  filter(OA %in% controls_A$OA) |>
  left_join(oa_meta, by = "OA")

# ensure clean types
treated_sf$OA <- as.character(treated_sf$OA)
control_sf$OA <- as.character(control_sf$OA)

treated_sf$scheme <- as.character(treated_sf$scheme)
control_sf$scheme <- as.character(control_sf$scheme)

# popup safety
treated_sf$popup <- paste0(
  "<b>Treated OA:</b> ", treated_sf$OA
)

ctrl_lookup <- pairs_A |>
  group_by(control_OA) |>
  summarise(
    matched_to = paste(unique(treated_OA), collapse = ", "),
    .groups = "drop"
  )

control_sf <- control_sf |>
  left_join(ctrl_lookup, by = c("OA" = "control_OA"))

control_sf$popup <- paste0(
  "<b>Control OA:</b> ", control_sf$OA,
  "<br><b>Matched to:</b> ", control_sf$matched_to
)

# =============================================================================
# BUILD MAP
# =============================================================================

cat("[5/5] Building interactive map...\n")

m <- leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolylines(
    data = lad_sf,
    color = "#777777",
    weight = 0.7,
    opacity = 0.5,
    group = "LAD boundaries"
  )

# =============================================================================
# SCHEME LAYERS =====================================

schemes <- unique(na.omit(treated_sf$scheme))
valid_layers <- c()

for (sc in schemes) {
  
  cat("Processing scheme:", sc, "\n")
  
  # treated OAs in scheme
  tr_sf <- treated_sf |>
    filter(scheme == sc)
  
  tr_oas <- tr_sf$OA
  
  # controls linked via matching
  ctrl_oas <- pairs_A |>
    filter(treated_OA %in% tr_oas) |>
    pull(control_OA) |>
    unique()
  
  ctrl_sf <- control_sf |>
    filter(OA %in% ctrl_oas)
  
  if (nrow(tr_sf) == 0 && nrow(ctrl_sf) == 0) next
  
  layer_name <- paste0("Scheme: ", sc)
  valid_layers <- c(valid_layers, layer_name)
  
  # controls
  if (nrow(ctrl_sf) > 0) {
    m <- m |>
      addPolygons(
        data = ctrl_sf,
        fillColor = "#1F78B4",   
        fillOpacity = 0.75,     
        color = "#0B3C5D",      
        weight = 1.2,            
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
  
  # treated
  if (nrow(tr_sf) > 0) {
    m <- m |>
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

# =============================================================================
# FINAL CONTROLS
# =============================================================================

m <- m |>
  addLayersControl(
    overlayGroups = c("LAD boundaries", valid_layers),
    options = layersControlOptions(collapsed = FALSE)
  ) |>
  addLegend(
    position = "bottomright",
    colors = c("#D85A30", "#2E6FAB"),
    labels = c("Treated OA", "Matched controls"),
    title = "Matching groups",
    opacity = 0.9
  ) |>
  setView(lng = -2.0, lat = 53.5, zoom = 6)

# =============================================================================
# SAVE
# =============================================================================

output_file <- here("outputs", "interactive_matching_map.html")

dir.create(here("outputs"), showWarnings = FALSE)

saveWidget(
  widget = m,
  file = output_file,
  selfcontained = TRUE
)

cat("\nMAP SAVED:\n", output_file, "\n")