# =============================================================================
# EXPORT INTERACTIVE MATCHING MAP
# Recent pooled matching data, England only
# =============================================================================
#
# Inputs:
#   data/processed/OA_matched_treated_pooled.rds
#   data/processed/OA_matched_donors_pooled.rds
#   data/processed/OA_matching_pairs_pooled.rds
#   data/processed/OA_matching_census.rds
#   data/processed/shp_files/OA_subset.shp
#   data/processed/shp_files/LADs_filtered.shp
#
# Output:
#   outputs/interactive_matching_map_recent_england.html
#
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

cat("============================================================\n")
cat(" EXPORTING INTERACTIVE MATCHING MAP: RECENT POOLED, ENGLAND ONLY\n")
cat("============================================================\n\n")

# =============================================================================
# LOAD RECENT MATCHING DATA
# =============================================================================

cat("[1/5] Loading recent pooled matching data...\n")

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
    scheme = scheme,
    treated_OA = as.character(.data[[pair_treated_col]]),
    control_OA = as.character(.data[[pair_control_col]])
  )

treated_pooled <- treated_pooled %>%
  mutate(OA = as.character(OA))

controls_pooled <- controls_pooled %>%
  mutate(OA = as.character(OA))

full_data <- full_data %>%
  mutate(OA = as.character(OA))

if ("LAD24CD" %in% names(full_data)) {
  full_data <- full_data %>%
    filter(substr(LAD24CD, 1, 1) == "E")
}

treated_pooled <- treated_pooled %>%
  filter(OA %in% full_data$OA)

controls_pooled <- controls_pooled %>%
  filter(OA %in% full_data$OA)

pairs_pooled <- pairs_pooled %>%
  filter(treated_OA %in% treated_pooled$OA, control_OA %in% controls_pooled$OA)

cat("  Treated OAs:", n_distinct(treated_pooled$OA), "\n")
cat("  Control OAs:", n_distinct(controls_pooled$OA), "\n")
cat("  Matched pairs:", nrow(pairs_pooled), "\n\n")

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

if (is.na(oa_id)) {
  stop("No OA identifier found in OA_subset.shp.")
}

oa_sf <- oa_raw %>%
  rename(OA = all_of(oa_id)) %>%
  mutate(OA = as.character(OA)) %>%
  filter(OA %in% full_data$OA) %>%
  st_transform(4326) %>%
  st_simplify(dTolerance = 0.0001, preserveTopology = TRUE)

lad_raw <- st_read(
  here("data", "processed", "shp_files", "LADs_filtered.shp"),
  quiet = TRUE
)

lad_id <- intersect(
  c("LAD21CD", "LAD22CD", "LAD24CD", "lad_code"),
  names(lad_raw)
)[1]

lad_nm <- intersect(
  c("LAD21NM", "LAD22NM", "LAD24NM", "lad_name"),
  names(lad_raw)
)[1]

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

cat("[3/5] Building lookup tables...\n")

oa_meta <- full_data %>%
  select(OA, any_of(c("LAD24CD", "LAD24NM")), scheme) %>%
  distinct()

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
  oa_meta <- oa_meta %>%
    mutate(LAD24NM = NA_character_)
}

treated_scheme_lookup <- pairs_pooled %>%
  distinct(treated_OA, scheme) %>%
  group_by(treated_OA) %>%
  summarise(scheme = paste(sort(unique(scheme)), collapse = "; "), .groups = "drop")

control_scheme_lookup <- pairs_pooled %>%
  distinct(control_OA, scheme) %>%
  group_by(control_OA) %>%
  summarise(scheme = paste(sort(unique(scheme)), collapse = "; "), .groups = "drop")

# =============================================================================
# CREATE SPATIAL OBJECTS
# =============================================================================

cat("[4/5] Preparing spatial layers...\n")

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
  select(-any_of(c("scheme_meta", "scheme_matched")))

ctrl_lookup <- pairs_pooled %>%
  group_by(control_OA) %>%
  summarise(
    matched_to = paste(sort(unique(treated_OA)), collapse = ", "),
    matched_schemes = paste(sort(unique(scheme)), collapse = "; "),
    .groups = "drop"
  )

control_sf <- control_sf %>%
  left_join(ctrl_lookup, by = c("OA" = "control_OA"))

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
  control_lad_popup
)

# =============================================================================
# BUILD MAP
# =============================================================================

cat("[5/5] Building interactive map...\n")

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
  
  tr_sf <- treated_sf %>%
    filter(OA %in% tr_oas)
  
  ctrl_sf <- control_sf %>%
    filter(OA %in% ctrl_oas)
  
  if (nrow(tr_sf) == 0 && nrow(ctrl_sf) == 0) next
  
  layer_name <- paste0("Scheme: ", sc)
  valid_layers <- c(valid_layers, layer_name)
  
  if (nrow(ctrl_sf) > 0) {
    m <- m %>%
      addPolygons(
        data = ctrl_sf,
        fillColor = "#1F78B4",
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
    colors = c("#D85A30", "#1F78B4"),
    labels = c("Treated OA", "Matched control OA"),
    title = "Recent matching groups",
    opacity = 0.9
  ) %>%
  setView(lng = -2.0, lat = 53.2, zoom = 6)

# =============================================================================
# SAVE
# =============================================================================

output_file <- here("outputs", "interactive_matching_map_recent_england.html")
summary_file <- here("outputs", "interactive_matching_map_recent_england_summary.csv")

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

cat("\nMAP SAVED:\n", output_file, "\n")
cat("SUMMARY SAVED:\n", summary_file, "\n")
