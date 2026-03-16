

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(zoo)
library(arrow)
library(units)
library(ggtext)   # for rich text labels
library(scales)   # for nicer axis formatting

# ── ────────────────────────────────────────────────────────────────
road_panel_model <- arrow::open_dataset(here("data", "processed", "road_panel_dataset")) %>%
  collect()



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
