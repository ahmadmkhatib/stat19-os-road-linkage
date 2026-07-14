# =============================================================================
# STATIC PUBLICATION MAP: ENGLISH CAZ SCHEMES
# =============================================================================
#
# Uses:
#   - pairs_pooled from the matching script
#   - LADs.shp for the correct England boundary
#   - CAZ_areas.shp for CAZ polygons
#
# The map:
#   - excludes roads and matched OA polygons
#   - displays Newcastle as Tyneside
#   - uses elbow connectors that do not cross labels
#   - provides sufficient margins for long labels
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(ggspatial)
  library(here)
})

map_crs <- 27700

# =============================================================================
# 1. Identify schemes included in the matched analysis
# =============================================================================

if (!exists("pairs_pooled")) {
  stop("pairs_pooled does not exist. Run the matching-data section first.")
}

if (!"scheme" %in% names(pairs_pooled)) {
  stop("The scheme column is missing from pairs_pooled.")
}

analysis_schemes <- pairs_pooled %>%
  filter(!is.na(scheme)) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  as.character()

cat("\nSchemes included in the static map:\n")
print(
  if_else(
    analysis_schemes == "Newcastle",
    "Tyneside",
    analysis_schemes
  )
)

# =============================================================================
# 2. Construct the correct England boundary
# =============================================================================

england_lads_raw <- st_read(
  here("data", "processed", "shp_files", "LADs.shp"),
  quiet = TRUE
)

lad_code_col <- intersect(
  c("LAD24CD", "LAD22CD", "LAD21CD", "lad_code"),
  names(england_lads_raw)
)[1]

if (is.na(lad_code_col)) {
  stop("Could not identify the LAD code column in LADs.shp.")
}

england_sf <- england_lads_raw %>%
  filter(substr(.data[[lad_code_col]], 1, 1) == "E") %>%
  st_transform(map_crs) %>%
  st_make_valid() %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) %>%
  st_make_valid()

if (nrow(england_sf) == 0) {
  stop("No English LAD geometries were retained.")
}

# =============================================================================
# 3. Load and prepare CAZ polygons
# =============================================================================

caz_raw <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_transform(map_crs) %>%
  st_make_valid()

scheme_col <- intersect(
  c("scheme", "Scheme", "SCHEME"),
  names(caz_raw)
)[1]

if (is.na(scheme_col)) {
  stop("Could not identify the scheme column in CAZ_areas.shp.")
}

caz_sf <- caz_raw %>%
  rename(scheme_internal = all_of(scheme_col)) %>%
  mutate(scheme_internal = as.character(scheme_internal)) %>%
  filter(scheme_internal %in% analysis_schemes) %>%
  group_by(scheme_internal) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) %>%
  mutate(
    display_name = if_else(
      scheme_internal == "Newcastle",
      "Tyneside",
      scheme_internal
    )
  ) %>%
  st_make_valid()

if (nrow(caz_sf) == 0) {
  stop("No CAZ geometries matched the schemes in pairs_pooled.")
}

missing_schemes <- setdiff(
  analysis_schemes,
  caz_sf$scheme_internal
)

if (length(missing_schemes) > 0) {
  warning(
    "No CAZ geometry found for: ",
    paste(missing_schemes, collapse = ", ")
  )
}

if ("Newcastle" %in% analysis_schemes &&
    !"Tyneside" %in% caz_sf$display_name) {
  stop("Newcastle geometry was not retained for display as Tyneside.")
}

# =============================================================================
# 4. Calculate CAZ points and label positions
# =============================================================================

caz_labels <- caz_sf %>%
  st_point_on_surface() %>%
  mutate(
    cx = st_coordinates(.)[, 1],
    cy = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  arrange(desc(cy))

england_bbox <- st_bbox(england_sf)

# Position scheme labels to the left of England.
label_x <- as.numeric(england_bbox["xmin"]) - 55000

# A small red point appears immediately to the right of each label.
label_dot_x <- label_x + 12000

# All horizontal connector lines end here before turning towards the CAZ.
elbow_x <- as.numeric(england_bbox["xmin"]) - 5000

caz_labels <- caz_labels %>%
  mutate(
    label_x = label_x,
    label_dot_x = label_dot_x,
    elbow_x = elbow_x,
    
    # Arrange labels from north to south.
    label_y = seq(
      from = as.numeric(england_bbox["ymax"]) - 20000,
      to = as.numeric(england_bbox["ymin"]) + 20000,
      length.out = n()
    )
  )

# Dynamic margins prevent Birmingham, Portsmouth or other long labels
# from being clipped.
left_label_space <- max(
  220000,
  max(nchar(caz_labels$display_name)) * 18000
)

# Extra space to the right for the north arrow and scale bar.
right_map_space <- 120000

# =============================================================================
# 5. Create the static map
# =============================================================================

england_caz_publication_map <- ggplot() +
  
  # England boundary
  geom_sf(
    data = england_sf,
    fill = "#F7F3EA",
    colour = "#555555",
    linewidth = 0.45
  ) +
  
  # CAZ polygons
  geom_sf(
    data = caz_sf,
    fill = "#D94841",
    colour = "#9B1D20",
    linewidth = 0.65,
    alpha = 0.95
  ) +
  
  # Horizontal part of each connector
  geom_segment(
    data = caz_labels,
    aes(
      x = label_dot_x,
      y = label_y,
      xend = elbow_x,
      yend = label_y
    ),
    colour = "#777777",
    linewidth = 0.35,
    lineend = "round"
  ) +
  
  # Angled part from the elbow to the CAZ
  geom_segment(
    data = caz_labels,
    aes(
      x = elbow_x,
      y = label_y,
      xend = cx,
      yend = cy
    ),
    colour = "#777777",
    linewidth = 0.35,
    lineend = "round"
  ) +
  
  # Small red point beside each scheme name
  geom_point(
    data = caz_labels,
    aes(
      x = label_dot_x,
      y = label_y
    ),
    colour = "#C83E3A",
    size = 1.8
  ) +
  
  # Right-aligned text ensures connector lines begin after the labels
  geom_text(
    data = caz_labels,
    aes(
      x = label_x,
      y = label_y,
      label = display_name
    ),
    hjust = 1,
    size = 4.1,
    colour = "#222222",
    family = "sans"
  ) +
  
  # Scale bar
  annotation_scale(
    location = "br",
    width_hint = 0.24,
    unit_category = "metric",
    pad_x = unit(0.55, "cm"),
    pad_y = unit(0.35, "cm"),
    text_cex = 0.8,
    line_width = 0.7
  ) +
  
  # North arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    
    # Smaller pad_x moves the arrow further right
    pad_x = unit(0.15, "cm"),
    pad_y = unit(1.25, "cm"),
    
    # Reduce arrow size
    height = unit(0.85, "cm"),
    width  = unit(0.85, "cm"),
    
    style = north_arrow_orienteering(
      fill = c("#333333", "white"),
      line_col = "#333333",
      text_size = 7
    )
  ) +
  
  # Plot extent, including label and map-furniture margins
  coord_sf(
    crs = map_crs,
    xlim = c(
      label_x - left_label_space,
      as.numeric(england_bbox["xmax"]) + right_map_space
    ),
    ylim = c(
      as.numeric(england_bbox["ymin"]) - 15000,
      as.numeric(england_bbox["ymax"]) + 15000
    ),
    expand = FALSE,
    clip = "on"
  ) +
  
  labs(
    title = "Clean Air Zones included in the pooled analysis",
    subtitle = "Seven English schemes",
    x = NULL,
    y = NULL
  ) +
  
  theme_void(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      colour = "#555555",
      margin = margin(b = 10)
    ),
    
    panel.background = element_rect(
      fill = "#EAF4FB",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  )

# Display the map in RStudio
print(england_caz_publication_map)

# =============================================================================
# 6. Save map
# =============================================================================

dir.create(
  here("outputs"),
  recursive = TRUE,
  showWarnings = FALSE
)

static_map_file <- here(
  "outputs",
  "England_CAZ_publication_map.png"
)

ggsave(
  filename = static_map_file,
  plot = england_caz_publication_map,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nStatic CAZ map saved:\n")
cat(" ", static_map_file, "\n")

