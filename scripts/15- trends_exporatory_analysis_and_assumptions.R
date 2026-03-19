

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(zoo)
library(arrow)
library(units)
library(ggtext)   # for rich text labels
library(scales)   # for nicer axis formatting
library(patchwork)
# ── ────────────────────────────────────────────────────────────────
road_panel_model <- arrow::open_dataset(here("data", "processed", "road_panel_dataset")) %>%
  collect()


names(road_panel_model)
summary(road_panel_model$quarter_year)

# roads inside the CAZs  == the treatment @ road level 
road_caz_prop<- readRDS(here("data", "processed", "roads_caz_props.rds"))

table(road_caz_prop$scheme)
## remove oxford 
road_caz_prop<- road_caz_prop %>% filter(!scheme=="Oxford ZEZ Pilot")


names(road_panel_model)

road_panel_model <- road_panel_model %>%
  mutate(
    group = case_when(
      treated_group_50pct == 1  ~ "CAZ roads",
      control_group2 == 1 ~ "Non-CAZ city controls"
    ),
    scheme = if_else(is.na(scheme), "Outside_CAZ_areas", scheme)
  )

# ── ──────────────────────────────────────────────────────
# ── Quarterly aggregation ──────────────────────────────────────────────────────
quarterly_trends <- road_panel_model %>% filter(!scheme=="Oxford ZEZ Pilot") %>% 
  group_by(scheme, quarter_year, group) %>%
  summarise(
    total_injuries = sum(across(starts_with("total_inj_adj")), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # ensure quarter_year is a yearqtr object
    quarter_year = zoo::as.yearqtr(quarter_year),
    # convert to Date at start of quarter
    quarter_date = as.Date(quarter_year, frac = 0),
    total_injuries = as.numeric(total_injuries)
  )


# ── CAZ start dates ────────────────────────────────────────────────────────────
road_caz_prop <- road_caz_prop %>%
  mutate(
    caz_start_date = as.Date(zoo::as.yearqtr(caz_start_q, format = "%Y Q%q"))
  )

# ── Key period boundaries ──────────────────────────────────────────────────────
covid_start <- as.Date("2020-03-01")
covid_end   <- as.Date("2021-03-01")

# ── Join CAZ dates & assign 5-period labels ───────────────────────────────────
quarterly_trends <- quarterly_trends %>%
  left_join(road_caz_prop %>% select(scheme, caz_start_date), by = "scheme") %>%
  mutate(
    period = case_when(
      # Pre-COVID
      quarter_date < covid_start                                          ~ "1_pre_covid",
      # During COVID
      quarter_date >= covid_start & quarter_date < covid_end              ~ "2_covid",
      # Post-COVID, before CAZ (only meaningful where CAZ exists)
      quarter_date >= covid_end   & !is.na(caz_start_date) &
        quarter_date < caz_start_date                                     ~ "3_post_covid_pre_caz",
      # Post-CAZ
      !is.na(caz_start_date) & quarter_date >= caz_start_date            ~ "4_post_caz",
      # Post-COVID, no CAZ (Outside_CAZ_areas)
      quarter_date >= covid_end   & is.na(caz_start_date)                ~ "5_post_covid_no_caz",
      TRUE ~ NA_character_
    )
  )

# ── Shading rectangles ─────────────────────────────────────────────────────────
shading_periods <- data.frame(
  xmin  = c(covid_start),
  xmax  = c(covid_end),
  label = c("COVID-19\nPandemic"),
  fill  = c("#FF6B6B"),
  alpha = c(0.12)
)

# ── Colour palette ─────────────────────────────────────────────────────────────
group_colours <- c(
  "CAZ roads"              = "#1B4F8A",   # deep navy
  "Non-CAZ city controls"  = "#E07B39"    # warm amber
)

smooth_colours <- c(
  "CAZ roads"              = "#5B9BD5",
  "Non-CAZ city controls"  = "#F5A96A"
)

# ── Build the plot ─────────────────────────────────────────────────────────────
p <- ggplot(
  quarterly_trends %>% filter(!is.na(group)),
  aes(x = quarter_date, y = total_injuries, colour = group)
) +
  
  # ── 1. COVID shaded band ──────────────────────────────────────────────────
  geom_rect(
    data        = shading_periods,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
        fill = label, alpha = alpha),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("COVID-19\nPandemic" = "#FF6B6B")) +
  scale_alpha_identity() +
  
  # COVID annotation label
  annotate(
    "text",
    x     = covid_start + (covid_end - covid_start) / 2,
    y     = Inf,
    vjust = 1.4,
    label = "COVID-19",
    size  = 2.8,
    colour = "#CC3333",
    fontface = "bold",
    family = "sans"
  ) +
  
  # ── 2. LM trend lines per group × period segment ───────────────────────
  geom_smooth(
    aes(group = interaction(group, period)),
    method    = "lm",      # linear = straight line
    se        = FALSE,     # no ribbon
    colour    = "red",
    linewidth = 1.4,
    linetype  = "solid"
  ) +
  
  # ── 3. Observed injury lines ──────────────────────────────────────────────
  geom_line(
    aes(group = group),
    linewidth = 0.9,
    alpha     = 0.85
  ) +
  
  # ── 4. Points at each quarter ─────────────────────────────────────────────
  geom_point(
    aes(group = group),
    size  = 1.6,
    alpha = 0.7,
    shape = 19
  ) +
  
  # ── 5. CAZ vertical line (only for schemes with a CAZ date) ───────────────
  geom_vline(
    data      = road_caz_prop %>% filter(!is.na(caz_start_date)),
    aes(xintercept = as.numeric(caz_start_date)),
    colour    = "#2CA02C",
    linetype  = "solid",
    linewidth = 1.1
  ) +
  
  # CAZ label — sits at the top of the vertical green line
  geom_text(
    data        = road_caz_prop %>% filter(!is.na(caz_start_date)),
    inherit.aes = FALSE,
    aes(x = caz_start_date, y = Inf,
        label = paste0("CAZ\n", format(caz_start_date, "%b %Y"))),
    vjust    = 1.2,      # pulls text just below the top edge
    hjust    = -0.1,     # slight rightward nudge off the line
    size     = 2.6,
    colour   = "#2CA02C",
    fontface = "bold"
  ) +
  
  # ── 6. Period-break dotted lines at covid_start / covid_end ───────────────
  geom_vline(xintercept = as.numeric(covid_start),
             linetype = "dotted", colour = "#CC3333", linewidth = 0.6) +
  geom_vline(xintercept = as.numeric(covid_end),
             linetype = "dotted", colour = "#CC3333", linewidth = 0.6) +
  
  # ── 7. Scales & coordinates ───────────────────────────────────────────────
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand      = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.12))  # headroom for CAZ label
  ) +
  scale_colour_manual(values = group_colours) +
  
  # ── 8. Facets ─────────────────────────────────────────────────────────────
  facet_wrap(~ scheme, scales = "free_y", ncol = 2) +
  
  # ── 9. Annotations for period segments (strip at bottom) ──────────────────
  # Done via subtitle / caption instead to avoid per-facet clutter
  labs(
    title    = "Quarterly total Road Injury Trends by CAZ Area",
    subtitle = paste0(
      "<span style='color:#555'>Shading: </span>",
      "<span style='color:#CC3333;font-weight:bold'>COVID-19 pandemic</span>",
      "   |   ",
      "<span style='color:#2CA02C;font-weight:bold'>──  CAZ implementation date</span>",
      "   |   ",
      "Red lines = Linear trend per period"
    ),
    caption  = paste0(
      "Periods: ",
      "(1) Pre-COVID  |  ",
      "(2) During COVID  |  ",
      "(3) Post-COVID / Pre-CAZ  |  ",
      "(4) Post-CAZ  |  ",
      "(5) Post-COVID (no CAZ)"
    ),
    x      = NULL,
    y      = "Total Road Injuries (adjusted)",
    colour = NULL
  ) +
  
  # ── 10. Theme ─────────────────────────────────────────────────────────────
  theme_minimal(base_size = 11) +
  theme(
    # Title & subtitle
    plot.title         = element_text(face = "bold", size = 14, hjust = 0,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_markdown(size = 9, hjust = 0,
                                          margin = margin(b = 10)),
    plot.caption       = element_text(size = 7.5, colour = "#666666",
                                      hjust = 0, margin = margin(t = 8)),
    
    # Facet strip
    strip.text         = element_text(face = "bold", size = 10,
                                      colour = "#222222"),
    strip.background   = element_rect(fill = "#F0F4FA", colour = NA),
    
    # Legend
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 10),
    legend.key.width   = unit(1.8, "cm"),
    
    # Axes
    axis.text.x        = element_text(size = 7.5, colour = "#444444"),
    axis.text.y        = element_text(size = 8,   colour = "#444444"),
    axis.title.y       = element_text(size = 9,   colour = "#333333",
                                      margin = margin(r = 6)),
    
    # Grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.4),
    panel.border       = element_rect(colour = "#DDDDDD", fill = NA,
                                      linewidth = 0.5),
    
    # Spacing
    plot.margin        = margin(12, 16, 10, 12)
  )

print(p)

# ── Save ───────────────────────────────────────────────────────────────────────
ggsave(
  here("outputs", "quarterly_injury_trends_caz.png"),
  plot   = p,
  width  = 14,
  height = 9,
  dpi    = 300,
  bg     = "white"
)







# ── Quarterly aggregation — stratified by injury type ─────────────────────────
quarterly_trends <- road_panel_model %>%
  filter(!scheme == "Oxford ZEZ Pilot") %>%
  group_by(scheme, quarter_year, group) %>%
  summarise(
    KSI    = sum(across(starts_with("KSI_adj")),    na.rm = TRUE),
    Slight = sum(across(starts_with("Slight_adj")), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # pivot to long so injury_type becomes an aesthetic
  pivot_longer(
    cols      = c(KSI, Slight),
    names_to  = "injury_type",
    values_to = "injuries"
  ) %>%
  mutate(
    quarter_year = zoo::as.yearqtr(quarter_year),
    quarter_date = as.Date(quarter_year, frac = 0),
    injuries     = as.numeric(injuries)
  )

# ── Join CAZ dates & period labels (same as before) ───────────────────────────
quarterly_trends <- quarterly_trends %>%
  left_join(road_caz_prop %>% select(scheme, caz_start_date), by = "scheme") %>%
  mutate(
    period = case_when(
      quarter_date < covid_start                                         ~ "1_pre_covid",
      quarter_date >= covid_start & quarter_date < covid_end             ~ "2_covid",
      quarter_date >= covid_end   & !is.na(caz_start_date) &
        quarter_date < caz_start_date                                    ~ "3_post_covid_pre_caz",
      !is.na(caz_start_date) & quarter_date >= caz_start_date           ~ "4_post_caz",
      quarter_date >= covid_end   & is.na(caz_start_date)               ~ "5_post_covid_no_caz",
      TRUE ~ NA_character_
    )
  )

# ── Colour & linetype palettes ─────────────────────────────────────────────────
group_colours <- c(
  "CAZ roads"             = "#1B4F8A",
  "Non-CAZ city controls" = "#E07B39"
)

injury_linetypes <- c(
  "KSI"    = "solid",
  "Slight" = "dashed"
)

# ── Build the plot ─────────────────────────────────────────────────────────────
p1 <- ggplot(
  quarterly_trends %>% filter(!is.na(group)),
  aes(x = quarter_date, y = injuries,
      colour   = group,
      linetype = injury_type,
      shape    = injury_type)
) +
  
  # ── 1. COVID shaded band ────────────────────────────────────────────────
  geom_rect(
    data        = shading_periods,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
        fill = label, alpha = alpha),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("COVID-19\nPandemic" = "#FF6B6B")) +
  scale_alpha_identity() +
  
  annotate(
    "text",
    x     = covid_start + (covid_end - covid_start) / 2,
    y     = Inf, vjust = 1.4,
    label = "COVID-19", size = 2.8,
    colour = "#CC3333", fontface = "bold", family = "sans"
  ) +
  
  # ── 2. LM trend lines per group × period × injury_type ─────────────────
  geom_smooth(
    aes(group = interaction(group, period, injury_type)),
    method    = "lm",
    se        = FALSE,
    colour    = "red",
    linewidth = 1.1,
    linetype  = "solid"    # trend lines always solid red; injury type shown on observed lines
  ) +
  
  # ── 3. Observed lines ───────────────────────────────────────────────────
  geom_line(
    aes(group = interaction(group, injury_type)),
    linewidth = 0.9,
    alpha     = 0.85
  ) +
  
  # ── 4. Points ───────────────────────────────────────────────────────────
  geom_point(
    aes(group = interaction(group, injury_type)),
    size = 1.5, alpha = 0.7
  ) +
  
  # ── 5. CAZ vertical line ────────────────────────────────────────────────
  geom_vline(
    data      = road_caz_prop %>% filter(!is.na(caz_start_date)),
    aes(xintercept = as.numeric(caz_start_date)),
    colour    = "#2CA02C", linetype = "solid", linewidth = 1.1
  ) +
  
  geom_text(
    data        = road_caz_prop %>% filter(!is.na(caz_start_date)),
    inherit.aes = FALSE,
    aes(x = caz_start_date, y = Inf,
        label = paste0("CAZ\n", format(caz_start_date, "%b %Y"))),
    vjust = 1.2, hjust = -0.1, size = 2.6,
    colour = "#2CA02C", fontface = "bold"
  ) +
  
  # ── 6. COVID boundary lines ─────────────────────────────────────────────
  geom_vline(xintercept = as.numeric(covid_start),
             linetype = "dotted", colour = "#CC3333", linewidth = 0.6) +
  geom_vline(xintercept = as.numeric(covid_end),
             linetype = "dotted", colour = "#CC3333", linewidth = 0.6) +
  
  # ── 7. Scales ───────────────────────────────────────────────────────────
  scale_x_date(
    date_breaks = "6 months", date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  scale_colour_manual(values  = group_colours) +
  scale_linetype_manual(values = injury_linetypes) +
  scale_shape_manual(values   = c("KSI" = 17, "Slight" = 19)) +
  
  # ── 8. Facets ───────────────────────────────────────────────────────────
  facet_wrap(~ scheme, scales = "free_y", ncol = 2) +
  
  # ── 9. Labels ───────────────────────────────────────────────────────────
  labs(
    title    = "Quarterly Road Injury Trends by CAZ Area",
    subtitle = paste0(
      "<span style='color:#555'>Shading: </span>",
      "<span style='color:#CC3333;font-weight:bold'>COVID-19 pandemic</span>",
      "   |   ",
      "<span style='color:#2CA02C;font-weight:bold'>──  CAZ implementation date</span>",
      "   |   ",
      "Red lines = Linear trend per period"
    ),
    caption = paste0(
      "Periods: (1) Pre-COVID  |  (2) During COVID  |  ",
      "(3) Post-COVID / Pre-CAZ  |  (4) Post-CAZ  |  (5) Post-COVID (no CAZ)\n",
      "Injury type: solid line / triangle = KSI;  dashed line / circle = Slight"
    ),
    x        = NULL,
    y        = "Road Injuries (adjusted)",
    colour   = NULL,
    linetype = "Injury type",
    shape    = "Injury type"
  ) +
  
  # ── 10. Theme (unchanged from original) ─────────────────────────────────
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_markdown(size = 9, hjust = 0,
                                          margin = margin(b = 10)),
    plot.caption       = element_text(size = 7.5, colour = "#666666",
                                      hjust = 0, margin = margin(t = 8)),
    strip.text         = element_text(face = "bold", size = 10, colour = "#222222"),
    strip.background   = element_rect(fill = "#F0F4FA", colour = NA),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 10),
    legend.key.width   = unit(1.8, "cm"),
    axis.text.x        = element_text(size = 7.5, colour = "#444444"),
    axis.text.y        = element_text(size = 8,   colour = "#444444"),
    axis.title.y       = element_text(size = 9,   colour = "#333333",
                                      margin = margin(r = 6)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.4),
    panel.border       = element_rect(colour = "#DDDDDD", fill = NA,
                                      linewidth = 0.5),
    plot.margin        = margin(12, 16, 10, 12)
  )

print(p1)

ggsave(
  here("outputs", "quarterly_injury_trends_ksi_slight.png"),
  plot = p1, width = 14, height = 9, dpi = 300, bg = "white"
)

######################################
#########################################################################################################################
###########################    MAP 
library(sf)
library(ggplot2)
library(patchwork)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)

# ── 1. Load & reproject ──────────────────────────────────────────────────────
target_crs <- 27700

caz_raw <- st_read(here("data","processed","shp_files","CAZ_areas.shp")) %>%
  st_make_valid() %>%
  st_transform(target_crs)

oa_sub <- st_read(here("data","processed","shp_files","OA_subset.shp")) %>%
  st_transform(target_crs)

# England & Scotland only — drop Wales, NI, Ireland
# 
eng_scot <- ne_countries(
  country     = "United Kingdom",
  scale       = "medium",        # ← j
  returnclass = "sf"
) %>%
  st_transform(target_crs)



# ── 2. Dissolve CAZ polygons per scheme ──────────────────────────────────────
# names(caz_raw)   ← run to find your name column
name_col <- "scheme"   # ← CHANGE to your actual column

caz <- caz_raw %>%
  rename(zone_name = all_of(name_col)) %>%
  group_by(zone_name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# ── 3. Centroids ordered north → south ───────────────────────────────────────
caz_centroids <- caz %>%
  st_centroid() %>%
  mutate(
    cx = st_coordinates(.)[,1],
    cy = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  arrange(desc(cy))

n_zones <- nrow(caz_centroids)

# ── 4. Left-margin label positions ───────────────────────────────────────────
map_bbox  <- st_bbox(eng_scot)
data_bbox <- st_bbox(oa_sub)

label_x      <- as.numeric(map_bbox["xmin"]) - 80000
line_start_x <- label_x + 55000

label_y <- seq(
  from       = as.numeric(data_bbox["ymax"]) - 10000,
  to         = as.numeric(data_bbox["ymin"]) + 10000,
  length.out = n_zones
)

caz_labels <- caz_centroids %>%
  mutate(
    label_x  = label_x,
    label_y  = label_y,
    label_no = row_number()
  )

# ── 5. Pick 4 zones to zoom: northernmost, southernmost, 2 middle ────────────
n         <- nrow(caz_centroids)
mid1_idx  <- floor(n / 3)        # upper-middle
mid2_idx  <- floor(2 * n / 3)    # lower-middle

zoom_indices <- c(1, mid1_idx, mid2_idx, n)   # N, mid-N, mid-S, S
zoom_labels  <- c("top-left", "top-right", "bottom-left", "bottom-right")

make_zoom_bbox <- function(zone_sf, pad = 3000) {
  bb <- st_bbox(zone_sf)
  list(
    xmin = bb["xmin"] - pad, xmax = bb["xmax"] + pad,
    ymin = bb["ymin"] - pad, ymax = bb["ymax"] + pad,
    name = zone_sf$zone_name
  )
}

zoom_zones <- lapply(zoom_indices, function(i) {
  zone_name_i <- caz_centroids$zone_name[i]
  zone_sf     <- caz %>% filter(zone_name == zone_name_i)
  make_zoom_bbox(zone_sf, pad = 3000)
})

# ── 6. Helper: build one zoom inset ──────────────────────────────────────────
make_inset <- function(zb) {
  ggplot() +
    geom_sf(data = oa_sub,
            fill      = "grey90",
            colour    = "grey55",
            linewidth = 0.12) +
    geom_sf(data = caz,
            fill      = "#E63946",
            colour    = "#9B1D20",
            linewidth = 0.45) +
    coord_sf(
      crs    = target_crs,
      xlim   = c(zb$xmin, zb$xmax),
      ylim   = c(zb$ymin, zb$ymax),
      expand = FALSE
    ) +
    theme_void(base_size = 8) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.1),
      plot.title   = element_text(size = 7, hjust = 0.5, face = "bold",
                                  margin = margin(t = 3, b = 2))
    ) +
    labs(title = zb$name)
}

insets <- lapply(zoom_zones, make_inset)

# ── 7. Main map ───────────────────────────────────────────────────────────────
# Dashed boxes on main map for the 4 zoomed zones
rect_data <- do.call(rbind, lapply(zoom_zones, function(zb) {
  data.frame(xmin = zb$xmin, xmax = zb$xmax,
             ymin = zb$ymin, ymax = zb$ymax)
}))

main_map <- ggplot() +
  
  geom_sf(data = eng_scot,
          fill      = "grey97",
          colour    = "grey30",
          linewidth = 0.5) +
  
  geom_sf(data = oa_sub,
          fill   = "#d4e6f1",
          colour = NA,
          alpha  = 0.5) +
  
  geom_sf(data = caz,
          fill      = "#E63946",
          colour    = "#9B1D20",
          linewidth = 0.5,
          alpha     = 0.9) +
  
  # 4 dashed zoom boxes
  geom_rect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill      = NA,
    colour    = "black",
    linewidth = 0.7,
    linetype  = "dashed",
    inherit.aes = FALSE
  ) +
  
  # Leader lines to centroids
  geom_segment(
    data      = caz_labels,
    aes(x = line_start_x, y = label_y, xend = cx, yend = cy),
    colour    = "grey40",
    linewidth = 0.3
  ) +
  
  # Centroid dots
  geom_point(
    data  = caz_labels,
    aes(x = cx, y = cy),
    shape  = 21,
    fill   = "#E63946",
    colour = "#9B1D20",
    size   = 1.8
  ) +
  
  # Numbered scheme labels
  geom_text(
    data  = caz_labels,
    aes(x = label_x, y = label_y,
        label = paste0(label_no, ". ", zone_name)),
    hjust      = 0,
    size       = 2.5,
    colour     = "grey10",
    lineheight = 0.85
  ) +
  
  coord_sf(
    crs    = target_crs,
    xlim   = c(label_x - 5000,                      as.numeric(map_bbox["xmax"]) + 20000),
    ylim   = c(as.numeric(map_bbox["ymin"]) - 10000, as.numeric(map_bbox["ymax"]) + 10000),
    expand = FALSE
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_line(colour = "grey92", linewidth = 0.3),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 9, hjust = 0.5, colour = "grey40"),
    axis.text        = element_text(size = 8, colour = "grey50"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  ) +
  labs(
    title    = "Clean Air Zones (CAZ) and Low Emission Zones (LEZ)",
    subtitle = "England and Scotland",
    x = NULL, y = NULL
  )

# ── 8. Compose: 4 insets in corners around the main map ──────────────────────
# Positions: top-left, top-right, bottom-left, bottom-right
# Values are fractions of the main panel; > 1 or < 0 = outside

final_map <- main_map +
  
  # Top-left inset (northernmost zone)
  inset_element(insets[[1]],
                left = -0.42, bottom = 0.72,
                right = -0.02, top = 0.99) +
  
  # Top-right inset (upper-middle zone)
  inset_element(insets[[2]],
                left = 1.02, bottom = 0.72,
                right = 1.42, top = 0.99) +
  
  # Bottom-left inset (lower-middle zone)
  inset_element(insets[[3]],
                left = -0.42, bottom = 0.01,
                right = -0.02, top = 0.28) +
  
  # Bottom-right inset (southernmost zone)
  inset_element(insets[[4]],
                left = 1.02, bottom = 0.01,
                right = 1.42, top = 0.28)

final_map

 ggsave(here("outputs","CAZ_LEZ_map_4insets.png"),
        final_map, width = 16, height = 13, dpi = 300, bg = "white")

library(tidyverse)
library(here)
library(arrow)
library(sf)
library(zoo)
library(knitr)

# -----------------------------
# Load Relevant Data
# -----------------------------

# Road-level panel and classification
road_panel_model<-arrow::open_dataset(here("data","processed","road_panel_dataset")) %>% collect()
analysis_roads   <- readRDS(here("data", "processed", "analysis_roads.rds"))

# OA-level data
OA_analysis <- readRDS(here("data", "processed", "OA_level_from_polygons.rds"))

# Injuries data
injuries <- readRDS(here("data", "processed", "injuries_matched_final.rds"))


# -----------------------------
# Road-Level Counts
# -----------------------------
road_summary <- analysis_roads %>%
  summarise(
    n_roads_total = n(),
    n_treated_any = sum(treated_group_any),
    n_treated50 = sum(treated_group_50pct),
    n_control1 = sum(control_group1),
    n_control2 = sum(control_group2)
  )

# -----------------------------
#  Injury Summaries by Treatment/Control
inj_summary <- road_panel_model %>%
  group_by(treated_group_any,treated_group_50pct, control_group1, control_group2) %>%
  summarise(
    total_KSI = sum(across(starts_with("KSI_adj")), na.rm = TRUE),
    total_Slight = sum(across(starts_with("Slight_adj")), na.rm = TRUE),
    total_injuries = sum(across(starts_with("total_inj_adj")), na.rm = TRUE),
    n_roads = n_distinct(identifier),
    .groups = "drop"
  )
# -----------------------------
# Quarterly Injury Distribution
# -----------------------------
quarterly_summary <- road_panel_model %>%
  group_by(quarter_year) %>%
  summarise(
    total_KSI = sum(across(starts_with("KSI_adj")), na.rm = TRUE),
    total_Slight = sum(across(starts_with("Slight_adj")), na.rm = TRUE),
    total_injuries = sum(across(starts_with("total_inj_adj")), na.rm = TRUE)
  )

# -----------------------------
#  OA-Level Descriptive Statistics
# -----------------------------
OA_summary <- OA_analysis %>%
  summarise(
    n_OAs_total = n(),
    n_treated_any = sum(treated_OA_any),
    n_treated_50 =sum(treated_OA_50pct),
    n_buffer = sum(buffer_OA),
    n_control1 = sum(control_group1_OA),
    n_control2 = sum(control_group2_OA),
    median_dist_to_centre = median(dist_citycentre),
    mean_dist_to_centre = mean(dist_citycentre),
    min_dist_to_centre = min(dist_citycentre),
    max_dist_to_centre = max(dist_citycentre)
  )






road_summary %>% st_drop_geometry() %>%kable(
  caption = "Road-Level Counts by Treatment and Control Groups"
)


inj_summary %>% st_drop_geometry() %>%kable(
  caption = "Total Injuries by Road Group"
)

quarterly_summary %>% st_drop_geometry() %>% kable(
  caption = "Quarterly Injury Distribution"
)

OA_summary %>% st_drop_geometry() %>%kable(
  caption = "OA-Level Summary Statistics"
)





library(flextable)
library(officer)

doc <- read_docx()

doc <- body_add_flextable(doc,
                          flextable(st_drop_geometry(road_summary)))

doc <- body_add_flextable(doc,
                          flextable(st_drop_geometry(inj_summary)))

print(doc, target = "tables.docx")
















glimpse(road_panel_model)

road_panel_model %>% select(treated_any) %>% table(useNA = "always")

# check the dimensions
n_distinct(road_panel_model$identifier) * n_distinct(road_panel_model$quarter_year)


road_panel_model %>% 
  summarise(
    n_roads = n_distinct(identifier),
    n_quarters = n_distinct(quarter_year),
    total_rows = n()
  )



## balance - panel per road 

road_panel_model %>%
  count(identifier) %>%
  summarise(
    min_obs = min(n),
    max_obs = max(n)
  )


## ever treated roads 
road_panel_model %>%
  distinct(identifier, treated_group_any) %>%
  count(treated_group_any
  )





#ttt timmimg 
road_panel_model %>%
  filter(treated_group_any == 1) %>%
  distinct(identifier, caz_start_q) %>%
  count(caz_start_q)

#Visualise treatment rollout 
road_panel %>%
  group_by(quarter_year) %>%
  summarise(
    n_treated = sum(treated)
  ) %>%
  ggplot(aes(quarter_year, n_treated)) +
  geom_line() +
  labs(title = "Number of Treated Roads Over Time")


# distibution 
road_panel %>%
  summarise(
    mean_inj = mean(total_inj_adj),
    var_inj  = var(total_inj_adj),
    zero_prop = mean(total_inj_adj == 0)
  )




road_panel %>%
  ggplot(aes(total_inj_adj)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0, 10)) +
  labs(title = "Distribution of Quarterly Injuries per Road")

#road avarages 

road_panel %>%
  group_by(identifier) %>%
  summarise(avg_inj = mean(total_inj_adj)) %>%
  ggplot(aes(avg_inj)) +
  geom_histogram(bins = 40) +
  labs(title = "Average Injuries per Road")



road_panel %>%
  group_by(quarter_year, ever_treated) %>%
  summarise(mean_inj = mean(total_inj_adj), .groups = "drop") %>%
  ggplot(aes(quarter_year, mean_inj, colour = factor(ever_treated))) +
  geom_line(size = 1) +
  labs(
    title = "Pre-Treatment Trends: Treated vs Control Roads",
    colour = "Ever Treated"
  )


min_caz <- min(road_panel$caz_start_q, na.rm = TRUE)

pre_data <- road_panel %>%
  filter(quarter_year < min_caz)

pre_data %>%
  group_by(ever_treated) %>%
  summarise(
    mean_inj = mean(total_inj_adj),
    mean_length = mean(length, na.rm = TRUE),
    prop_A_road = mean(road_class == "A", na.rm = TRUE),
    n_roads = n_distinct(identifier)
  )

## hetrogenity byu calss 


road_panel %>%
  group_by(road_class, ever_treated) %>%
  summarise(mean_inj = mean(total_inj_adj), .groups = "drop") %>%
  ggplot(aes(road_class, mean_inj, fill = factor(ever_treated))) +
  geom_col(position = "dodge")



## serial corelation check 

top_road <- road_panel %>%
  group_by(identifier) %>%
  summarise(avg = mean(total_inj_adj)) %>%
  arrange(desc(avg)) %>%
  slice(1) %>%
  pull(identifier)

road_panel %>%
  filter(identifier == top_road) %>%
  ggplot(aes(quarter_year, total_inj_adj)) +
  geom_line()


#   event structure 


road_panel <- road_panel %>%
  mutate(
    event_time = as.numeric(quarter_year - caz_start_q)
  )


road_panel %>%
  filter(ever_treated == 1) %>%
  ggplot(aes(event_time)) +
  geom_histogram(bins = 30)


### zero inf 

road_panel %>%
  group_by(ever_treated) %>%
  summarise(
    zero_rate = mean(total_inj_adj == 0)
  )



## agregation 

road_panel %>%
  mutate(year = year(as.Date(quarter_year))) %>%
  group_by(identifier, year) %>%
  summarise(yearly_inj = sum(total_inj_adj), .groups = "drop") %>%
  summarise(
    mean_yearly = mean(yearly_inj),
    zero_prop_yearly = mean(yearly_inj == 0)
  )











