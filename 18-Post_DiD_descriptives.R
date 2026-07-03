# =============================================================================
#  post-DiD descriptives:
# Indexed treated/control injury trends and Bradford quarterly counts
# =============================================================================
#
# Run after the matched road panel has been created, and after the DiD script
# if you want the same majority-quarter CAZ start adjustment used there.
#
# Outputs:
#   output/pooled/All_clean/descriptive_trends/
#     - indexed_2019_<scheme>.png
#     - indexed_2019_all_schemes_facets.png
#     - all_schemes_quarterly_injury_counts_indexed.csv
#     - bradford_quarterly_injury_counts.csv
#
# =============================================================================

library(tidyverse)
library(arrow)
library(here)
library(zoo)
library(lubridate)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename

outcome_var <- "total_inj_adj_All"

outdir <- here("output", "pooled", "All_clean", "descriptive_trends")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load matched road panel and align CAZ start quarter with the DiD script
# =============================================================================

road_panel <- arrow::read_parquet(
  here("data", "processed", "road_panel_matched_pooled.parquet")
) %>%
  mutate(
    quarter_year = as.yearqtr(quarter_year),
    caz_start_q = as.yearqtr(caz_start_q)
  )

if (!outcome_var %in% names(road_panel)) {
  stop("Outcome column not found in road_panel: ", outcome_var)
}

road_caz_props <- readRDS(here("data", "processed", "roads_caz_props.rds"))

scheme_start <- road_caz_props %>%
  distinct(scheme, startDt, caz_start_q) %>%
  filter(!is.na(startDt)) %>%
  mutate(
    start_date = dmy(startDt),
    raw_q = as.yearqtr(start_date),
    q_start = as.Date(raw_q),
    q_end = as.Date(raw_q + 0.25) - 1,
    q_mid = q_start + as.integer(difftime(q_end, q_start, units = "days")) / 2,
    caz_start_q_adj = if_else(start_date > q_mid, raw_q + 0.25, raw_q)
  ) %>%
  select(scheme, caz_start_q_adj)

road_panel <- road_panel %>%
  left_join(scheme_start, by = "scheme") %>%
  mutate(
    caz_start_q = coalesce(caz_start_q_adj, caz_start_q),
    group = factor(
      if_else(treat_group == 1, "Treated", "Control"),
      levels = c("Control", "Treated")
    )
  ) %>%
  select(-caz_start_q_adj)

scheme_timing <- road_panel %>%
  filter(treat_group == 1, !is.na(caz_start_q)) %>%
  distinct(scheme, caz_start_q) %>%
  arrange(caz_start_q)

min_qtr <- min(as.numeric(road_panel$quarter_year), na.rm = TRUE)

road_id_col <- case_when(
  "identifier" %in% names(road_panel) ~ "identifier",
  "panel_id" %in% names(road_panel) ~ "panel_id",
  TRUE ~ NA_character_
)

if (is.na(road_id_col)) {
  stop("Could not find a road-link identifier column: expected identifier or panel_id.")
}

# =============================================================================
# 2. Quarterly injury counts and 2019-indexed trends
# =============================================================================

quarterly_counts <- road_panel %>%
  left_join(scheme_timing %>% rename(scheme_start_q = caz_start_q), by = "scheme") %>%
  mutate(
    qtr_int = as.integer(round((as.numeric(quarter_year) - min_qtr) * 4)) + 1L,
    start_qtr_int = as.integer(round((as.numeric(scheme_start_q) - min_qtr) * 4)) + 1L,
    event_time = qtr_int - start_qtr_int
  ) %>%
  group_by(scheme, quarter_year, event_time, group, treat_group) %>%
  summarise(
    injury_count = sum(.data[[outcome_var]], na.rm = TRUE),
    n_road_links = n_distinct(.data[[road_id_col]]),
    n_oas = n_distinct(OA),
    .groups = "drop"
  )

baseline_2019 <- quarterly_counts %>%
  filter(
    quarter_year >= as.yearqtr("2019 Q1"),
    quarter_year <= as.yearqtr("2019 Q4")
  ) %>%
  group_by(scheme, group) %>%
  summarise(
    baseline_2019_mean = mean(injury_count, na.rm = TRUE),
    .groups = "drop"
  )

quarterly_indexed <- quarterly_counts %>%
  left_join(baseline_2019, by = c("scheme", "group")) %>%
  mutate(
    index_2019 = if_else(
      is.na(baseline_2019_mean) | baseline_2019_mean == 0,
      NA_real_,
      100 * injury_count / baseline_2019_mean
    ),
    quarter_date = as.Date(quarter_year),
    calendar_year = as.integer(format(quarter_date, "%Y")),
    calendar_quarter = paste0("Q", cycle(quarter_year))
  ) %>%
  arrange(scheme, group, quarter_year)

write_csv(
  quarterly_indexed,
  file.path(outdir, "all_schemes_quarterly_injury_counts_indexed.csv")
)

# =============================================================================
# 3. Plotting functions
# =============================================================================

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

plot_scheme_index <- function(sc) {
  plot_df <- quarterly_indexed %>%
    filter(scheme == sc)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  sc_start <- scheme_timing %>%
    filter(scheme == sc) %>%
    pull(caz_start_q) %>%
    first()
  
  dotted_dates <- plot_df %>%
    filter(event_time %in% c(-6, -5)) %>%
    distinct(quarter_date) %>%
    pull(quarter_date)
  
  ggplot(plot_df, aes(x = quarter_date, y = index_2019, colour = group)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_vline(
      xintercept = as.numeric(dotted_dates),
      linetype = "dotted",
      colour = "black",
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = as.numeric(as.Date(sc_start)),
      linetype = "dashed",
      colour = "black",
      linewidth = 0.5
    ) +
    scale_colour_manual(
      values = c("Control" = "#F8766D", "Treated" = "#00BFC4"),
      drop = FALSE
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(
      title = paste0(sc, ": treated vs. control, indexed to 2019 average"),
      subtitle = "Dotted = event times -6/-5; dashed = actual CAZ start",
      x = NULL,
      y = "Index (2019 = 100)",
      colour = "group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

schemes_all <- sort(unique(quarterly_indexed$scheme))

scheme_plots <- set_names(map(schemes_all, plot_scheme_index), schemes_all)

walk2(scheme_plots, names(scheme_plots), function(p, sc) {
  if (is.null(p)) return(NULL)
  ggsave(
    filename = file.path(outdir, paste0("indexed_2019_", safe_filename(sc), ".png")),
    plot = p,
    width = 9,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
})

facet_lines <- scheme_timing %>%
  transmute(
    scheme,
    line_type = "CAZ start",
    quarter_date = as.Date(caz_start_q)
  ) %>%
  bind_rows(
    quarterly_indexed %>%
      filter(event_time %in% c(-6, -5)) %>%
      distinct(scheme, quarter_date) %>%
      mutate(line_type = "Event time -6/-5")
  )

facet_plot <- ggplot(
  quarterly_indexed,
  aes(x = quarter_date, y = index_2019, colour = group)
) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  geom_vline(
    data = facet_lines,
    aes(xintercept = as.numeric(quarter_date), linetype = line_type),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.35
  ) +
  facet_wrap(~ scheme, scales = "free_y") +
  scale_colour_manual(
    values = c("Control" = "#F8766D", "Treated" = "#00BFC4"),
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = c("Event time -6/-5" = "dotted", "CAZ start" = "dashed")
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Treated vs. matched-control injury trends by CAZ",
    subtitle = "Quarterly injury counts indexed to each group's 2019 average",
    x = NULL,
    y = "Index (2019 = 100)",
    colour = "group",
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(outdir, "indexed_2019_all_schemes_facets.png"),
  plot = facet_plot,
  width = 14,
  height = 9,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 4. Bradford CAZ vs. control: injury counts by quarter/year
# =============================================================================

bradford_quarterly_counts <- quarterly_indexed %>%
  filter(scheme == "Bradford") %>%
  select(
    scheme,
    calendar_year,
    calendar_quarter,
    quarter_year,
    event_time,
    group,
    injury_count,
    n_road_links,
    n_oas,
    baseline_2019_mean,
    index_2019
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(
      injury_count,
      n_road_links,
      n_oas,
      baseline_2019_mean,
      index_2019
    ),
    names_glue = "{group}_{.value}"
  ) %>%
  mutate(
    treated_control_ratio =
      if_else(Control_injury_count == 0, NA_real_,
              Treated_injury_count / Control_injury_count)
  ) %>%
  arrange(quarter_year)

write_csv(
  bradford_quarterly_counts,
  file.path(outdir, "bradford_quarterly_injury_counts.csv")
)

cat("\nSaved indexed CAZ trend plots and Bradford table to:\n")
cat(outdir, "\n")
cat("\nBradford quarterly injury-count table:\n")
print(bradford_quarterly_counts, n = Inf)
