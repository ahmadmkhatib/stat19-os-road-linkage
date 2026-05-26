library(sf)
library(ggplot2)
library(here)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

target_crs <- 27700

# ── Load data ─────────────────────────────────────────────────────────────────
caz_raw <- st_read(here("data","processed","shp_files","CAZ_areas.shp")) %>%
  st_make_valid() %>%
  st_transform(target_crs)

oa_sub <- st_read(here("data","processed","shp_files","OA_subset.shp")) %>%
  st_transform(target_crs)

# ── UK base map split into England / Scotland ─────────────────────────────────
uk <- ne_countries(country = "United Kingdom", scale = "medium", returnclass = "sf") %>%
  st_transform(target_crs)

scotland_box <- st_as_sfc(st_bbox(c(
  xmin = -200000, xmax = 700000,
  ymin = 550000,  ymax = 1300000
), crs = target_crs))

england_box <- st_as_sfc(st_bbox(c(
  xmin = -200000, xmax = 700000,
  ymin = -100000, ymax = 550000
), crs = target_crs))

scotland <- st_intersection(uk, scotland_box)
england  <- st_intersection(uk, england_box)

# ── Clean up CAZ ──────────────────────────────────────────────────────────────
caz <- caz_raw %>%
  rename(zone_name = scheme) %>%
  filter(!grepl("Oxford", zone_name, ignore.case = TRUE)) %>%
  group_by(zone_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# ── Centroids ─────────────────────────────────────────────────────────────────
caz_centroids <- caz %>%
  st_centroid() %>%
  mutate(
    cx = st_coordinates(.)[,1],
    cy = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  arrange(desc(cy))

# ── Bounding boxes ────────────────────────────────────────────────────────────
map_bbox  <- st_bbox(st_union(england, scotland) %>% st_union())
data_bbox <- st_bbox(oa_sub)

map_xmin <- as.numeric(map_bbox["xmin"])
map_xmax <- as.numeric(map_bbox["xmax"])
map_ymin <- as.numeric(map_bbox["ymin"])
map_ymax <- as.numeric(map_bbox["ymax"])

# ── Split zones ───────────────────────────────────────────────────────────────
scotland_threshold <- 600000

caz_scotland <- caz_centroids %>% filter(cy >= scotland_threshold) %>% arrange(desc(cy))
caz_england  <- caz_centroids %>% filter(cy <  scotland_threshold) %>% arrange(desc(cy))

# ── LEFT labels ───────────────────────────────────────────────────────────────
label_x_left  <- map_xmin - 90000
line_gap_left <- 15000

n_eng <- nrow(caz_england)
label_y_left <- seq(
  from       = as.numeric(data_bbox["ymax"]) + 20000,
  to         = as.numeric(data_bbox["ymin"]) - 20000,
  length.out = n_eng
)

caz_labels_left <- caz_england %>%
  mutate(
    label_x    = label_x_left,
    label_y    = label_y_left,
    line_end_x = label_x_left + line_gap_left
  )

# ── RIGHT labels ──────────────────────────────────────────────────────────────
label_x_right  <- map_xmax + 40000
line_gap_right <- 15000

n_sco <- nrow(caz_scotland)
label_y_right <- seq(
  from       = map_ymax - 20000,
  to         = map_ymax - 20000 - (n_sco - 1) * 60000,
  length.out = n_sco
)

caz_labels_right <- caz_scotland %>%
  mutate(
    label_x    = label_x_right,
    label_y    = label_y_right,
    line_end_x = label_x_right - line_gap_right
  )

# ── Main map ──────────────────────────────────────────────────────────────────
main_map <- ggplot() +
  
  # England fill — no border colour so no internal E/S line
  geom_sf(data = england,
          fill = "#f5f0e8", colour = NA, linewidth = 0.5) +
  
  # Scotland fill — no border colour
  geom_sf(data = scotland,
          fill = "#e8f0e8", colour = NA, linewidth = 0.5) +
  
  # UK outline drawn on top — gives clean coastline without internal border
  geom_sf(data = uk,
          fill = NA, colour = "grey30", linewidth = 0.5) +
  
  geom_sf(data = oa_sub,
          fill = "#d4e6f1", colour = NA, alpha = 0.5) +
  
  geom_sf(data = caz,
          fill = "#E63946", colour = "#9B1D20", linewidth = 0.5, alpha = 0.9) +
  
  # LEFT leader lines
  geom_segment(
    data = caz_labels_left,
    aes(x = cx, y = cy, xend = line_end_x, yend = label_y),
    colour = "grey50", linewidth = 0.3
  ) +
  geom_point(
    data = caz_labels_left,
    aes(x = cx, y = cy),
    shape = 21, fill = "#E63946", colour = "#9B1D20", size = 2.2
  ) +
  geom_text(
    data = caz_labels_left,
    aes(x = label_x, y = label_y, label = zone_name),
    hjust = 0, size = 4.5, colour = "grey10", lineheight = 0.9
  ) +
  
  # RIGHT leader lines
  geom_segment(
    data = caz_labels_right,
    aes(x = cx, y = cy, xend = line_end_x, yend = label_y),
    colour = "grey50", linewidth = 0.3
  ) +
  geom_point(
    data = caz_labels_right,
    aes(x = cx, y = cy),
    shape = 21, fill = "#E63946", colour = "#9B1D20", size = 2.2
  ) +
  geom_text(
    data = caz_labels_right,
    aes(x = label_x, y = label_y, label = zone_name),
    hjust = 0, size = 4.5, colour = "grey10", lineheight = 0.9
  ) +
  
  # Country labels
  annotate("text", x = 380000, y = 300000,
           label = "England", size = 5.5, colour = "grey40",
           fontface = "italic", alpha = 0.8) +
  
  annotate("text", x = 230000, y = 720000,
           label = "Scotland", size = 5.5, colour = "grey40",
           fontface = "italic", alpha = 0.8) +
  
  # North arrow
  annotation_north_arrow(
    location    = "br",
    which_north = "true",
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering(
      fill     = c("grey30", "white"),
      line_col = "grey20",
      text_col = "grey20",
      text_size = 10
    ),
    height = unit(1.8, "cm"),
    width  = unit(1.8, "cm")
  ) +
  
  coord_sf(
    crs    = target_crs,
    xlim   = c(label_x_left - 10000,  label_x_right + 180000),
    ylim   = c(map_ymin - 10000,       map_ymax + 10000),
    expand = FALSE
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid       = element_line(colour = "grey92", linewidth = 0.3),
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, colour = "grey40"),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  ) +
  labs(
    title    = "Clean Air Zones (CAZ) and Low Emission Zones (LEZ)",
    subtitle = "England and Scotland",
    x = NULL, y = NULL
  )

main_map

ggsave(here("outputs","CAZ_LEZ_map.png"),
       main_map, width = 16, height = 13, dpi = 300, bg = "white")