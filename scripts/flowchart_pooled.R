# =============================================================================
# MATCHING FLOWCHART — POOLED MATCHING (SIMPLIFIED)
# =============================================================================
#
# Simple vertical flowchart:
# England OAs → Eligible → Stage 1 → Stage 2 → Final
#
# INPUTS (from data/processed/):
#   OA_matched_full_pooled.rds
#   OA_matching_pairs_pooled.rds
#   OA_common_support_flags_pooled.rds
#   OA_matching_census.rds
#   OA_stage1_summary_pooled.rds
#   OA_stage1_controls_pooled.rds
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

# New Stage 1 outputs from matching script
stage1_summary_path  <- here("data", "processed", "OA_stage1_summary_pooled.rds")
stage1_controls_path <- here("data", "processed", "OA_stage1_controls_pooled.rds")

stage1_summary_available  <- file.exists(stage1_summary_path)
stage1_controls_available <- file.exists(stage1_controls_path)

if (stage1_summary_available) {
  stage1_summary_combined <- readRDS(stage1_summary_path)
} else {
  warning("OA_stage1_summary_pooled.rds not found. Stage 1 treated count will use eligible treated count.")
  stage1_summary_combined <- NULL
}

if (stage1_controls_available) {
  stage1_controls_combined <- readRDS(stage1_controls_path)
} else {
  warning("OA_stage1_controls_pooled.rds not found. Stage 1 control count will fall back to final controls, which understates the Stage 1 donor pool.")
  stage1_controls_combined <- NULL
}

OA_matching_dataset <- OA_matching_dataset %>%
  mutate(
    country = case_when(
      substr(LAD24CD, 1, 1) == "E" ~ "England",
      substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
      TRUE                          ~ "Other"
    )
  )

n_england <- sum(OA_matching_dataset$country == "England", na.rm = TRUE)

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
n_treated   <- sum(data_england$treated_OA == 1, na.rm = TRUE)
n_ctrl_pool <- n_eligible - n_treated
n_excluded  <- n_england - n_eligible

n_buffer <- sum(
  OA_matching_dataset$country == "England" &
    OA_matching_dataset$buffer_OA == 1,
  na.rm = TRUE
)

n_same_city <- sum(
  OA_matching_dataset$country == "England" &
    OA_matching_dataset$control_group1_OA == 1,
  na.rm = TRUE
)

n_zero_roads <- sum(
  OA_matching_dataset$country == "England" &
    OA_matching_dataset$n_roads == 0,
  na.rm = TRUE
)

n_zero_inj <- sum(
  OA_matching_dataset$country == "England" &
    OA_matching_dataset$treated_OA == 1 &
    OA_matching_dataset$zero_injury_OA == 1,
  na.rm = TRUE
)

# =============================================================================
# STAGE 1 NUMBERS — UNIQUE OA COUNTS ONLY
# =============================================================================

stage1_controls_path <- here("data", "processed", "OA_stage1_controls_pooled.rds")

if (!file.exists(stage1_controls_path)) {
  stop(
    "OA_stage1_controls_pooled.rds not found. ",
    "Run the updated matching script first so Stage 1 retained controls are saved."
  )
}

stage1_controls_combined <- readRDS(stage1_controls_path)

# All eligible treated OAs pass into Stage 1.
# Stage 1 is used to reduce the control donor pool.
n_treated_s1_retained <- n_treated
n_treated_s1_excluded <- 0

# Unique eligible controls before Stage 1
eligible_control_oas <- data_england %>%
  filter(treated_OA == 0) %>%
  distinct(OA)

n_ctrl_pool <- nrow(eligible_control_oas)

# Unique control OAs retained by Stage 1 in at least one scheme
stage1_control_oas <- stage1_controls_combined %>%
  distinct(OA)

n_ctrl_s1_retained <- nrow(stage1_control_oas)

# Controls excluded at Stage 1
n_ctrl_s1_excluded <- n_ctrl_pool - n_ctrl_s1_retained

# Optional: identify excluded control OAs
stage1_excluded_control_oas <- eligible_control_oas %>%
  anti_join(stage1_control_oas, by = "OA")

# Total unique OAs after Stage 1
n_stage1_retained_unique <- n_treated_s1_retained + n_ctrl_s1_retained
n_stage1_excluded_unique <- n_treated_s1_excluded + n_ctrl_s1_excluded

cat("\n=== Stage 1 numbers — unique OA counts ===\n")
cat("Eligible treated OAs:", formatC(n_treated, big.mark = ","), "\n")
cat("Eligible control OAs:", formatC(n_ctrl_pool, big.mark = ","), "\n")
cat("Stage 1 retained treated OAs:", formatC(n_treated_s1_retained, big.mark = ","), "\n")
cat("Stage 1 excluded treated OAs:", formatC(n_treated_s1_excluded, big.mark = ","), "\n")
cat("Stage 1 retained unique control OAs:", formatC(n_ctrl_s1_retained, big.mark = ","), "\n")
cat("Stage 1 excluded control OAs:", formatC(n_ctrl_s1_excluded, big.mark = ","), "\n")
cat("Stage 1 retained unique OAs total:", formatC(n_stage1_retained_unique, big.mark = ","), "\n")
cat("Stage 1 excluded unique OAs total:", formatC(n_stage1_excluded_unique, big.mark = ","), "\n")

# =============================================================================
# FINAL MATCHED DATASET NUMBERS
# =============================================================================

n_treated_final <- matched_full %>%
  filter(treat_indicator == 1) %>%
  summarise(n = n_distinct(OA)) %>%
  pull(n)

n_ctrl_final <- matched_full %>%
  filter(treat_indicator == 0) %>%
  summarise(n = n_distinct(OA)) %>%
  pull(n)

n_pairs <- nrow(matching_pairs)

cat("=== Flowchart numbers ===\n")
cat("England OAs:", n_england, "\n")
cat("Excluded before matching:", n_excluded, "\n")
cat("Eligible:", n_eligible, "(", n_treated, "treated +", n_ctrl_pool, "controls)\n")
cat("Stage 1 retained:", n_stage1_retained, "(", n_treated_s1, "treated +", n_ctrl_s1, "controls)\n")
cat("Stage 1 excluded:", n_stage1_excluded, "(", n_treated_s1_excluded, "treated +", n_ctrl_s1_excluded, "controls)\n")
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
  grid.lines(
    x = unit(c(x_main, x_side - w / 2), "npc"),
    y = unit(c(y_from, y_side), "npc"),
    gp = gpar(col = "#C0392B", lwd = 1.2, lty = 2)
  )
  draw_box(
    x_side, y_side, w, h, label,
    fill = fill,
    border = "#C0392B",
    text_size = 8
  )
}

# =============================================================================
# DRAW FLOWCHART
# =============================================================================

png(file.path(outdir, "fig_flowchart_pooled.png"),
    width = 2400, height = 3200, res = 300)

grid.newpage()
pushViewport(viewport(width = 0.92, height = 0.94))

# Colours
COL_BOX    <- "#EBF5FB"
COL_STAGE  <- "#D4EFDF"
COL_EXCL   <- "#FADBD8"
COL_FINAL  <- "#D5F5E3"

# Layout
box_w <- 0.52
box_h <- 0.070
y_positions <- c(0.88, 0.66, 0.43, 0.18)

# ── Box 1: England OAs ──
draw_box(
  0.45, y_positions[1], box_w, box_h,
  paste0(
    "England OAs\n",
    "(n = ", formatC(n_england, big.mark = ","), ")"
  ),
  fill = COL_BOX,
  text_size = 10,
  font = 2
)

draw_arrow(
  0.45, y_positions[1] - box_h / 2,
  0.45, y_positions[2] + box_h / 2
)

# Exclusion side box before eligibility
draw_side_box(
  0.45 + box_w / 2, y_positions[1] - 0.04,
  0.82, (y_positions[1] + y_positions[2]) / 2,
  0.30, 0.085,
  paste0(
    "Excluded before matching\n",
    "(n = ", formatC(n_excluded, big.mark = ","), ")\n",
    "Buffer: ", formatC(n_buffer, big.mark = ","), "\n",
    "Same-city controls: ", formatC(n_same_city, big.mark = ","), "\n",
    "Zero roads: ", formatC(n_zero_roads, big.mark = ","), "\n",
    "Zero-injury treated: ", formatC(n_zero_inj, big.mark = ",")
  ),
  fill = COL_EXCL
)

# ── Box 2: Eligible OAs ──
draw_box(
  0.45, y_positions[2], box_w, box_h,
  paste0(
    "Eligible OAs for matching\n",
    formatC(n_treated, big.mark = ","), " treated + ",
    formatC(n_ctrl_pool, big.mark = ","), " other-city controls\n",
    "Total eligible = ", formatC(n_eligible, big.mark = ",")
  ),
  fill = COL_BOX,
  text_size = 9
)

draw_arrow(
  0.45, y_positions[2] - box_h / 2,
  0.45, y_positions[3] + box_h / 2
)

# Stage 1 label on arrow
grid.text(
  "Stage 1: Mahalanobis matching\nroad, urban, business and sociodemographic variables",
  x = unit(0.45, "npc"),
  y = unit((y_positions[2] + y_positions[3]) / 2 + 0.025, "npc"),
  gp = gpar(
    fontsize = 8,
    col = "#27AE60",
    fontface = "italic",
    lineheight = 1.0
  )
)

# Stage 1 excluded side box
draw_side_box(
  0.45 + box_w / 2,
  (y_positions[2] + y_positions[3]) / 2 - 0.025,
  0.82,
  (y_positions[2] + y_positions[3]) / 2 - 0.055,
  0.30,
  0.075,
  paste0(
    "Excluded at Stage 1\n",
    "(n = ", formatC(n_stage1_excluded, big.mark = ","), ")\n",
    "Treated excluded: ", formatC(n_treated_s1_excluded, big.mark = ","), "\n",
    "Controls not selected: ", formatC(n_ctrl_s1_excluded, big.mark = ",")
  ),
  fill = COL_EXCL
)

# ── Box 3: After Stage 1 ──
draw_box(
  0.45, y_positions[3], box_w, 0.085,
  paste0(
    "After Stage 1 matching\n",
    formatC(n_treated_s1, big.mark = ","), " treated OAs retained\n",
    formatC(n_ctrl_s1, big.mark = ","), " unique control OAs retained\n",
    "Stage 1 retained total = ", formatC(n_stage1_retained, big.mark = ",")
  ),
  fill = COL_STAGE,
  text_size = 9,
  border = "#239B56"
)

draw_arrow(
  0.45, y_positions[3] - 0.085 / 2,
  0.45, y_positions[4] + 0.09 / 2
)

# Stage 2 label on arrow
grid.text(
  "Stage 2: Mahalanobis matching\ntotal injury trend + level; ratio selected per scheme",
  x = unit(0.45, "npc"),
  y = unit((y_positions[3] + y_positions[4]) / 2, "npc"),
  gp = gpar(
    fontsize = 8,
    col = "#27AE60",
    fontface = "italic",
    lineheight = 1.0
  )
)

# ── Box 4: Final matched dataset ──
draw_box(
  0.45, y_positions[4], box_w, 0.095,
  paste0(
    "Final matched dataset\n",
    formatC(n_treated_final, big.mark = ","), " treated OAs  |  ",
    formatC(n_ctrl_final, big.mark = ","), " unique control OAs\n",
    formatC(n_pairs, big.mark = ","), " treated–control pairs across 7 CAZ schemes"
  ),
  fill = COL_FINAL,
  text_size = 10,
  font = 2,
  border = "#1E8449"
)

popViewport()
dev.off()

cat("Saved: fig_flowchart_pooled.png\n")