# =============================================================================
# England-only CAZ map using the latest pooled matching outputs
# =============================================================================
#
# Inputs:
#   data/processed/OA_matched_full_pooled.rds
#   data/processed/shp_files/OA_subset.shp
#   data/processed/shp_files/CAZ_areas.shp
#
# Output:
#   outputs/England_CAZ_recent_matching_map.png
#
# Notes:
#   - Scotland is removed entirely.
#   - CAZ schemes and OA coverage are taken from the recent matched data.
#   - Treated and matched-control OAs are mapped separately.
# =============================================================================

library(sf)
library(ggplot2)
library(here)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

target_crs <- 27700

# =============================================================================
# 1. Load recent matching data
# =============================================================================

matched_full <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))

if (!"OA" %in% names(matched_full)) {
  stop("OA column is missing from OA_matched_full_pooled.rds.")
}

if (!"scheme" %in% names(matched_full)) {
  stop("scheme column is missing from OA_matched_full_pooled.rds.")
}

if (!"treat_indicator" %in% names(matched_full)) {
  if ("treated_OA" %in% names(matched_full)) {
    matched_full <- matched_full %>%
      mutate(treat_indicator = as.integer(treated_OA == 1))
  } else {
    stop("Need treat_indicator or treated_OA in OA_matched_full_pooled.rds.")
  }
}

matched_oa <- matched_full %>%
  mutate(
    group = if_else(treat_indicator == 1, "Treated OA", "Matched control OA")
  ) %>%
  distinct(scheme, OA, group)

matched_schemes <- matched_oa %>%
  distinct(scheme) %>%
  pull(scheme) %>%
  sort()

cat("\nMatched schemes in current data:\n")
print(matched_schemes)

# =============================================================================
# 2. Load spatial data and filter to England/current matched sample
# =============================================================================

oa_sf <- st_read(here("data", "processed", "shp_files", "OA_subset.shp"),
                 quiet = TRUE) %>%
  st_transform(target_crs)

if (!"OA" %in% names(oa_sf)) {
  stop("OA column is missing from OA_subset.shp.")
}

if ("LAD24CD" %in% names(oa_sf)) {
  oa_sf <- oa_sf %>%
    filter(substr(LAD24CD, 1, 1) == "E")
}

matched_oa_sf <- oa_sf %>%
  inner_join(matched_oa, by = "OA") %>%
  st_make_valid()

if (nrow(matched_oa_sf) == 0) {
  stop("No matched OAs joined to OA_subset.shp. Check OA identifiers.")
}

caz_raw <- st_read(here("data", "processed", "shp_files", "CAZ_areas.shp"),
                   quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(target_crs)

if (!"scheme" %in% names(caz_raw)) {
  stop("scheme column is missing from CAZ_areas.shp.")
}

caz_england_sf <- caz_raw %>%
  filter(
    scheme %in% matched_schemes,
    !str_detect(scheme, regex("Oxford|Glasgow|Aberdeen|Dundee|Edinburgh", ignore_case = TRUE))
  ) %>%
  rename(zone_name = scheme) %>%
  group_by(zone_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

if (nrow(caz_england_sf) == 0) {
  stop("No CAZ geometries matched to the current matched schemes.")
}

# England base map only. This avoids any Scotland geometry or label.
uk <- ne_countries(country = "United Kingdom", scale = "medium", returnclass = "sf") %>%
  st_transform(target_crs)

england_box <- st_as_sfc(st_bbox(c(
  xmin = -200000, xmax = 700000,
  ymin = -100000, ymax = 550000
), crs = target_crs))

england <- st_intersection(uk, england_box)

# =============================================================================
# 3. Label positions
# =============================================================================

caz_centroids <- caz_england_sf %>%
  st_centroid() %>%
  mutate(
    cx = st_coordinates(.)[, 1],
    cy = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  arrange(desc(cy))

england_bbox <- st_bbox(england)
matched_bbox <- st_bbox(matched_oa_sf)

label_x_left <- as.numeric(england_bbox["xmin"]) - 90000

caz_labels <- caz_centroids %>%
  mutate(
    label_x = label_x_left,
    label_y = seq(
      from = as.numeric(matched_bbox["ymax"]) + 20000,
      to = as.numeric(matched_bbox["ymin"]) - 20000,
      length.out = n()
    ),
    line_end_x = label_x_left + 15000
  )

# =============================================================================
# 4. Map
# =============================================================================

england_caz_map <- ggplot() +
  geom_sf(
    data = england,
    fill = "#f5f0e8",
    colour = "grey35",
    linewidth = 0.45
  ) +
  geom_sf(
    data = matched_oa_sf %>% filter(group == "Matched control OA"),
    fill = "#9ecae1",
    colour = NA,
    alpha = 0.45
  ) +
  geom_sf(
    data = matched_oa_sf %>% filter(group == "Treated OA"),
    fill = "#fdae6b",
    colour = NA,
    alpha = 0.75
  ) +
  geom_sf(
    data = caz_england_sf,
    fill = "#E63946",
    colour = "#9B1D20",
    linewidth = 0.55,
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
  geom_text(
    data = caz_labels,
    aes(x = label_x, y = label_y, label = zone_name),
    hjust = 0,
    size = 4.2,
    colour = "grey10",
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
    xlim = c(label_x_left - 10000, as.numeric(england_bbox["xmax"]) + 50000),
    ylim = c(as.numeric(england_bbox["ymin"]) - 10000,
             as.numeric(england_bbox["ymax"]) + 10000),
    expand = FALSE
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
  ) +
  labs(
    title = "Clean Air Zones and Recent Matched OA Sample",
    subtitle = "England only; treated and matched-control OAs from OA_matched_full_pooled.rds",
    x = NULL,
    y = NULL
  )

print(england_caz_map)

ggsave(
  here("outputs", "England_CAZ_recent_matching_map.png"),
  england_caz_map,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

matched_map_summary <- matched_oa %>%
  count(scheme, group, name = "n_oas") %>%
  arrange(scheme, group)

write_csv(
  matched_map_summary,
  here("outputs", "England_CAZ_recent_matching_map_summary.csv")
)

cat("\nSaved:\n")
cat("  outputs/England_CAZ_recent_matching_map.png\n")
cat("  outputs/England_CAZ_recent_matching_map_summary.csv\n")
