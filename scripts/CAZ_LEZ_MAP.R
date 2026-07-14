# =============================================================================
# England CAZ map: zones only
# =============================================================================
#
# The map deliberately excludes matched OA/road-sample layers. Matching data
# are read only to identify the seven schemes included in the pooled analysis.
#
# "Newcastle" remains the internal scheme identifier used for data joins, but
# it is displayed as "Tyneside" in the map and console summary.
#
# Inputs:
#   data/processed/OA_matched_full_pooled.rds
#   data/processed/shp_files/CAZ_areas.shp
#
# Output:
#   outputs/England_CAZ_zones_only_map.png
# =============================================================================

library(sf)
library(ggplot2)
library(here)
library(dplyr)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

target_crs <- 27700

# =============================================================================
# 1. Identify the schemes in the current matched analysis
# =============================================================================

matched_full <- readRDS(
  here("data", "processed", "OA_matched_full_pooled.rds")
)

if (!"scheme" %in% names(matched_full)) {
  stop("scheme column is missing from OA_matched_full_pooled.rds.")
}

matched_schemes <- matched_full %>%
  filter(!is.na(scheme)) %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("\nSchemes in the current matched analysis:\n")
print(if_else(matched_schemes == "Newcastle", "Tyneside", matched_schemes))

# =============================================================================
# 2. Load CAZ boundaries and construct the England base map
# =============================================================================

caz_raw <- st_read(
  here("data", "processed", "shp_files", "CAZ_areas.shp"),
  quiet = TRUE
) %>%
  st_make_valid() %>%
  st_transform(target_crs)

if (!"scheme" %in% names(caz_raw)) {
  stop("scheme column is missing from CAZ_areas.shp.")
}

caz_england_sf <- caz_raw %>%
  filter(
    scheme %in% matched_schemes,
    !str_detect(
      scheme,
      regex("Oxford|Glasgow|Aberdeen|Dundee|Edinburgh", ignore_case = TRUE)
    )
  ) %>%
  group_by(scheme) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(
    display_name = if_else(scheme == "Newcastle", "Tyneside", scheme)
  ) %>%
  st_make_valid()

if (nrow(caz_england_sf) == 0) {
  stop("No CAZ geometries matched the schemes in the current analysis.")
}

if ("Newcastle" %in% matched_schemes &&
    !any(caz_england_sf$display_name == "Tyneside")) {
  stop("The Newcastle CAZ geometry was not retained for display as Tyneside.")
}

uk <- ne_countries(
  country = "United Kingdom",
  scale = "medium",
  returnclass = "sf"
) %>%
  st_transform(target_crs)

# Clip the UK geometry to an England-focused plotting extent. The northern
# limit must extend beyond 550,000 metres to retain Newcastle/Tyneside.
england_box <- st_as_sfc(st_bbox(c(
  xmin = -200000,
  xmax = 700000,
  ymin = -100000,
  ymax = 660000
), crs = target_crs))

england <- st_intersection(uk, england_box)

# =============================================================================
# 3. Construct label positions
# =============================================================================

caz_centroids <- caz_england_sf %>%
  st_point_on_surface() %>%
  mutate(
    cx = st_coordinates(.)[, 1],
    cy = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  arrange(desc(cy))

england_bbox <- st_bbox(england)
caz_bbox <- st_bbox(caz_england_sf)
label_x_left <- as.numeric(england_bbox["xmin"]) - 90000

# Keep the first and last labels inside the plotting panel. Without this cap,
# the Tyneside label is positioned above the northern map boundary and clipped.
label_y_top <- min(
  as.numeric(caz_bbox["ymax"]) + 10000,
  as.numeric(england_bbox["ymax"]) - 25000
)
label_y_bottom <- max(
  as.numeric(caz_bbox["ymin"]) - 10000,
  as.numeric(england_bbox["ymin"]) + 25000
)

caz_labels <- caz_centroids %>%
  mutate(
    label_x = label_x_left,
    label_y = seq(
      from = label_y_top,
      to = label_y_bottom,
      length.out = n()
    ),
    line_end_x = label_x_left + 15000
  )

# =============================================================================
# 4. Draw the zones-only map
# =============================================================================

england_caz_map <- ggplot() +
  geom_sf(
    data = england,
    fill = "#f5f0e8",
    colour = "grey35",
    linewidth = 0.45
  ) +
  geom_sf(
    data = caz_england_sf,
    fill = "#E63946",
    colour = "#9B1D20",
    linewidth = 0.65,
    alpha = 0.88
  ) +
  geom_segment(
    data = caz_labels,
    aes(x = cx, y = cy, xend = line_end_x, yend = label_y),
    colour = "grey45",
    linewidth = 0.3
  ) +
  geom_point(
    data = caz_labels,
    aes(x = cx, y = cy),
    shape = 21,
    fill = "#E63946",
    colour = "#9B1D20",
    size = 2.2
  ) +
  # Use a borderless background label rather than plain text. Because the
  # connector segments are drawn first, this background masks the part of a
  # line that would otherwise run through a scheme name.
  geom_label(
    data = caz_labels,
    aes(x = label_x, y = label_y, label = display_name),
    hjust = 0,
    size = 4.2,
    colour = "grey10",
    fill = "#eaf4fb",
    linewidth = 0,
    label.padding = unit(0.08, "lines"),
    label.r = unit(0, "lines"),
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 380000,
    y = 300000,
    label = "England",
    size = 5.5,
    colour = "grey40",
    fontface = "italic",
    alpha = 0.85
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering(
      fill = c("grey30", "white"),
      line_col = "grey20",
      text_col = "grey20",
      text_size = 10
    ),
    height = unit(1.6, "cm"),
    width = unit(1.6, "cm")
  ) +
  coord_sf(
    crs = target_crs,
    xlim = c(
      label_x_left - 10000,
      as.numeric(england_bbox["xmax"]) + 50000
    ),
    ylim = c(
      as.numeric(england_bbox["ymin"]) - 10000,
      as.numeric(england_bbox["ymax"]) + 10000
    ),
    expand = FALSE
  ) +
  labs(
    title = "Clean Air Zones in England",
    subtitle = "Schemes included in the pooled analysis",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_line(colour = "grey92", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

print(england_caz_map)

ggsave(
  here("outputs", "England_CAZ_zones_only_map.png"),
  england_caz_map,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("\nSaved:\n")
cat("  outputs/England_CAZ_zones_only_map.png\n")
