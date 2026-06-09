# =============================================================================
# MATCHING FLOWCHART — POOLED MATCHING (SIMPLIFIED)
# =============================================================================
#
# Simple vertical flowchart: England OAs → Eligible → Stage 1 → Stage 2 → Final
#
# INPUTS (from data/processed/):
#   OA_matched_full_pooled.rds
#   OA_matching_pairs_pooled.rds
#   OA_common_support_flags_pooled.rds
#   OA_matching_census.rds
#
# OUTPUT:
#   output/diagnostics/pooled/fig_flowchart_pooled.png
#
# =============================================================================

library(tidyverse)
library(here)
library(grid)

select <- dplyr::select
filter <- dplyr::filter

outdir <- here("output", "diagnostics", "pooled")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA + COMPUTE NUMBERS
# =============================================================================

OA_matching_dataset <- readRDS(here("data", "processed", "OA_matching_census.rds"))
matched_full        <- readRDS(here("data", "processed", "OA_matched_full_pooled.rds"))
matching_pairs      <- readRDS(here("data", "processed", "OA_matching_pairs_pooled.rds"))

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(country = if_else(substr(LAD24CD, 1, 1) == "E", "England", "Scotland"))

n_england <- sum(OA_matching_dataset$country == "England")

# Eligible: treated + other-city controls after exclusions
data_england <- OA_matching_dataset %>%
  filter(
    country == "England",
    (treated_OA == 1 | control_group2_OA == 1),
    control_group1_OA == 0,
    buffer_OA == 0,
    n_roads > 0,
    !(treated_OA == 1 & zero_injury_OA == 1)
  )

n_eligible  <- nrow(data_england)
n_treated   <- sum(data_england$treated_OA == 1)
n_ctrl_pool <- n_eligible - n_treated
n_excluded  <- n_england - n_eligible

n_buffer     <- sum(OA_matching_dataset$country == "England" &
                      OA_matching_dataset$buffer_OA == 1)
n_same_city  <- sum(OA_matching_dataset$country == "England" &
                      OA_matching_dataset$control_group1_OA == 1)
n_zero_roads <- sum(OA_matching_dataset$country == "England" &
                      OA_matching_dataset$n_roads == 0)
n_zero_inj   <- sum(OA_matching_dataset$country == "England" &
                      OA_matching_dataset$treated_OA == 1 &
                      OA_matching_dataset$zero_injury_OA == 1)

# Stage 1 donor pool: unique controls retained across all schemes
# (each scheme keeps up to 50 nearest per treated OA, then deduplicates)
# We approximate from the matched dataset + the fact that Stage 2 selects from
# the Stage 1 pool. The Stage 1 pool is larger than the final controls.
# Use the actual number from the matching run.
n_ctrl_s1 <- n_distinct(matched_full$OA[matched_full$treat_indicator == 0])

# Final
n_treated_final <- sum(matched_full$treat_indicator == 1)
n_ctrl_final    <- n_distinct(matched_full$OA[matched_full$treat_indicator == 0])
n_pairs         <- nrow(matching_pairs)

cat("=== Flowchart numbers ===\n")
cat("England OAs:", n_england, "\n")
cat("Excluded:", n_excluded, "\n")
cat("Eligible:", n_eligible, "(", n_treated, "treated +", n_ctrl_pool, "controls)\n")
cat("Final treated:", n_treated_final, "\n")
cat("Final unique controls:", n_ctrl_final, "\n")
cat("Total pairs:", n_pairs, "\n")

# =============================================================================
# DRAWING HELPERS
# =============================================================================

draw_box <- function(x, y, w, h, label, fill = "#FFFFFF",
                     border = "#2C3E50", text_size = 9,
                     text_col = "#2C3E50", font = 1) {
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(0.012, "npc"),
    gp = gpar(fill = fill, col = border, lwd = 1.8)
  )
  grid.text(
    label,
    x = unit(x, "npc"), y = unit(y, "npc"),
    gp = gpar(fontsize = text_size, col = text_col, fontface = font,
              lineheight = 1.1)
  )
}

draw_arrow <- function(x0, y0, x1, y1, col = "#2C3E50") {
  grid.lines(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    gp = gpar(col = col, lwd = 1.8),
    arrow = arrow(length = unit(0.012, "npc"), type = "closed")
  )
}

draw_side_box <- function(x_main, y_from, x_side, y_side, w, h,
                          label, fill = "#FADBD8") {
  # Horizontal line then box
  grid.lines(
    x = unit(c(x_main, x_side - w/2), "npc"),
    y = unit(c(y_from, y_side), "npc"),
    gp = gpar(col = "#C0392B", lwd = 1.2, lty = 2)
  )
  draw_box(x_side, y_side, w, h, label, fill = fill,
           border = "#C0392B", text_size = 8)
}

# =============================================================================
# DRAW FLOWCHART
# =============================================================================

png(file.path(outdir, "fig_flowchart_pooled.png"),
    width = 2400, height = 3000, res = 300)

grid.newpage()
pushViewport(viewport(width = 0.92, height = 0.94))

# Colours
COL_BOX    <- "#EBF5FB"
COL_STAGE  <- "#D4EFDF"
COL_EXCL   <- "#FADBD8"
COL_FINAL  <- "#D5F5E3"

# Layout: 4 main boxes evenly spaced vertically
box_w <- 0.52
box_h <- 0.065
y_positions <- c(0.88, 0.66, 0.44, 0.18)

# ── Box 1: England OAs ──
draw_box(0.45, y_positions[1], box_w, box_h,
         paste0("England OAs\n(n = ",
                formatC(n_england, big.mark = ","), ")"),
         fill = COL_BOX, text_size = 10, font = 2)

draw_arrow(0.45, y_positions[1] - box_h/2,
           0.45, y_positions[2] + box_h/2)

# Exclusion side box
draw_side_box(
  0.45 + box_w/2, y_positions[1] - 0.04,
  0.82, (y_positions[1] + y_positions[2]) / 2, 0.30, 0.08,
  paste0("Excluded (n = ",
         formatC(n_excluded, big.mark = ","), ")\n",
         "Buffer: ", formatC(n_buffer, big.mark = ","),
         "  |  Same-city controls: ", formatC(n_same_city, big.mark = ","), "\n",
         "Zero roads: ", formatC(n_zero_roads, big.mark = ","),
         "  |  Zero-injury treated: ", n_zero_inj),
  fill = COL_EXCL
)

# ── Box 2: Eligible OAs ──
draw_box(0.45, y_positions[2], box_w, box_h,
         paste0("Eligible OAs\n(",
                formatC(n_treated, big.mark = ","), " treated + ",
                formatC(n_ctrl_pool, big.mark = ","),
                " other-city controls = ",
                formatC(n_eligible, big.mark = ","), ")"),
         fill = COL_BOX, text_size = 9)

draw_arrow(0.45, y_positions[2] - box_h/2,
           0.45, y_positions[3] + box_h/2)

# Stage 1 label on arrow
grid.text(
  "Stage 1: Mahalanobis matching\nroad, urban, sociodemographic (1:50)",
  x = unit(0.45, "npc"),
  y = unit((y_positions[2] + y_positions[3]) / 2, "npc"),
  gp = gpar(fontsize = 8, col = "#27AE60", fontface = "italic",
            lineheight = 1.0)
)

# ── Box 3: After Stage 1 ──
draw_box(0.45, y_positions[3], box_w, box_h,
         paste0("Stage 1 donor pool (per scheme)\n",
                formatC(n_treated, big.mark = ","),
                " treated OAs + up to 50 nearest controls each"),
         fill = COL_STAGE, text_size = 9)

draw_arrow(0.45, y_positions[3] - box_h/2,
           0.45, y_positions[4] + box_h/2)

# Stage 2 label on arrow
grid.text(
  "Stage 2: Mahalanobis matching\ntotal injury trend + level (ratio 1:1\u20131:10 per scheme)",
  x = unit(0.45, "npc"),
  y = unit((y_positions[3] + y_positions[4]) / 2, "npc"),
  gp = gpar(fontsize = 8, col = "#27AE60", fontface = "italic",
            lineheight = 1.0)
)

# ── Box 4: Final matched dataset ──
draw_box(0.45, y_positions[4], box_w, 0.09,
         paste0("Final matched dataset\n",
                formatC(n_treated_final, big.mark = ","),
                " treated OAs  |  ",
                formatC(n_ctrl_final, big.mark = ","),
                " unique control OAs\n",
                formatC(n_pairs, big.mark = ","),
                " treated\u2013control pairs across 7 CAZ schemes"),
         fill = COL_FINAL, text_size = 10, font = 2,
         border = "#1E8449")

popViewport()
dev.off()

cat("Saved: fig_flowchart_pooled.png\n")
