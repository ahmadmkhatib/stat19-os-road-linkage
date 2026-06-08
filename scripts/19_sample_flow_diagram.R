library(tidyverse)
library(here)

outdir <- here("output", "diagnostics")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# ENGLAND-ONLY FLOW COUNTS
# =============================================================================

counts <- list(
  # Box 2: Large Urban Areas (England)
  matching_eng = 57805,
  treated_eng = 702, buffer_eng = 1597,
  ctrl1_eng = 8421, ctrl2_eng = 47085,

  # Exclusion cascade
  excl_total = 10974,
  excl_buffer = 1597,
  excl_samecity = 8421,
  excl_zero_roads = 879,
  excl_zero_inj_treated = 77,
  post_buffer = 56208,
  post_samecity = 47787,
  post_roads = 46908,

  # Eligible study population (England)
  eligible = 46831,
  elig_treated = 613,
  elig_ctrl = 46218,

  # Stage 1: Structural restriction (per-scheme, ratio 1:50)
  s1_total = 2887, s1_treated = 613, s1_ctrl = 2274,
  s1_ctrl_excl = 43944, s1_ctrl_excl_pct = 95.1,

  # Stage 2 / Final analytical sample
  s2_ctrl_excl = 1440, s2_ctrl_excl_pct = 63.3,
  final_total = 1447, final_treated = 613, final_ctrl = 834
)
fN <- function(x) format(x, big.mark = ",", scientific = FALSE)

# =============================================================================
# TEXT
# =============================================================================

txt <- list(
  b2_head = "Large Urban Areas (population \u2265100,000)\nin England",
  b2_body = paste0("N = ", fN(counts$matching_eng), " OAs\n",
                   "Treated: ", fN(counts$treated_eng),
                   "   Buffer: ", fN(counts$buffer_eng), "\n",
                   "Same-city controls: ", fN(counts$ctrl1_eng),
                   "   Other-city: ", fN(counts$ctrl2_eng)),

  b3_head = "Eligible Study Population",
  b3_body = paste0("N = ", fN(counts$eligible), " OAs\n",
                   "Treated: ", fN(counts$elig_treated),
                   "   Controls: ", fN(counts$elig_ctrl)),

  b4_head = "Stage 1: Structural Control Pool Restriction",
  b4_body = paste0("N = ", fN(counts$s1_total), " OAs\n",
                   "Treated: ", fN(counts$s1_treated),
                   "   Controls: ", fN(counts$s1_ctrl)),

  b4b_head = "Stage 2: Pretreatment Trends and Levels Matching",
  b4b_body = paste0("N = ", fN(counts$final_total), " OAs\n",
                    "Treated: ", fN(counts$final_treated),
                    "   Controls: ", fN(counts$final_ctrl)),
  b4b_foot = "Analytical Sample",

  eA_head = paste0("Excluded:  ", fN(counts$excl_total), " OAs"),
  eA_body = paste0(
    "Buffer zone: ", fN(counts$excl_buffer),
    "  \u2192  ", fN(counts$post_buffer), "\n",
    "Same-city controls (CAZ): ", fN(counts$excl_samecity),
    "  \u2192  ", fN(counts$post_samecity), "\n",
    "Zero roads: ", fN(counts$excl_zero_roads),
    "  \u2192  ", fN(counts$post_roads), "\n",
    "Zero-injury treated: ", fN(counts$excl_zero_inj_treated),
    "  \u2192  ", fN(counts$eligible)),

  eB_head = paste0("OAs excluded: ", fN(counts$s1_ctrl_excl)),
  eB_body = paste0(counts$s1_ctrl_excl_pct, "% of eligible controls\n",
                   "No Stage 1 structural match found"),

  eC_head = paste0("OAs excluded: ", fN(counts$s2_ctrl_excl)),
  eC_body = paste0(counts$s2_ctrl_excl_pct, "% of Stage 1 pool")
)

# =============================================================================
# AUTO-SIZING — boxes scale to text content
# =============================================================================
count_lines <- function(s) length(strsplit(s, "\n")[[1]])

LH  <- 4.0
PAD <- 3.5
GAP <- 2.5

box_h <- function(...) {
  items <- list(...)
  lines <- sum(sapply(items, count_lines))
  n     <- length(items)
  lines * LH + (n - 1) * 1.5 + 2 * PAD
}

h2   <- box_h(txt$b2_head,  txt$b2_body)
h3   <- box_h(txt$b3_head,  txt$b3_body)
h4   <- box_h(txt$b4_head,  txt$b4_body)
h4b  <- box_h(txt$b4b_head, txt$b4b_body, txt$b4b_foot)
eA_h <- box_h(txt$eA_head,  txt$eA_body)
eB_h <- box_h(txt$eB_head,  txt$eB_body)
eC_h <- box_h(txt$eC_head,  txt$eC_body)

# =============================================================================
# LAYOUT
# =============================================================================
CX  <- 34;  BX1 <- 1;  BX2 <- 67
ECX <- 84;  EX1 <- 70; EX2 <- 98

gap_A <- eA_h + 2 * GAP
gap_B <- eB_h + 2 * GAP

total_h <- h2 + gap_A + h3 + GAP + h4 + gap_B + h4b

y       <- total_h + 3
b2_top  <- y;    b2_bot  <- b2_top - h2;   b2_cy  <- (b2_top + b2_bot) / 2
b3_top  <- b2_bot - gap_A; b3_bot <- b3_top - h3;  b3_cy  <- (b3_top + b3_bot) / 2
b4_top  <- b3_bot - GAP;   b4_bot <- b4_top - h4;  b4_cy  <- (b4_top + b4_bot) / 2
b4b_top <- b4_bot - gap_B; b4b_bot <- b4b_top - h4b; b4b_cy <- (b4b_top + b4b_bot) / 2

eA_top <- b2_bot  - GAP;  eA_bot <- eA_top - eA_h;  eA_cy <- (eA_top + eA_bot) / 2
eB_top <- b4_bot  - GAP;  eB_bot <- eB_top - eB_h;  eB_cy <- (eB_top + eB_bot) / 2
eC_top <- b4b_bot - GAP;  eC_bot <- eC_top - eC_h;  eC_cy <- (eC_top + eC_bot) / 2

canvas_top <- b2_top + 4
canvas_bot <- eC_bot - 3

# =============================================================================
# COLOURS
# =============================================================================
fill_data  <- "#F0F4FA";  fill_elig  <- "#E6EFF9";  fill_stage <- "#F0EAF8"
fill_final <- "#E4F2E4";  fill_excl  <- "#FDF0ED"
brd_data   <- "#4A7FB5";  brd_elig   <- "#2E6FAB";  brd_stage  <- "#6B3FA0"
brd_final  <- "#2A7A2A";  brd_excl   <- "#C0392B"
col_arrow  <- "#444444";  col_title  <- "#1A2E5A"

fH <- 9;  fB <- 8;  fI <- 7;  fS <- 8.5;  fE <- 8;  fe <- 7

arr <- arrow(length = unit(0.3, "cm"), type = "closed")

# =============================================================================
# PLOT
# =============================================================================

p <- ggplot() +

  ## -- BOX 2 -----------------------------------------------------------------
  annotate("rect", xmin = BX1, xmax = BX2, ymin = b2_bot, ymax = b2_top,
           fill = fill_data, colour = brd_data, linewidth = 0.7) +
  annotate("text", x = CX, y = b2_top - PAD - LH * count_lines(txt$b2_head) / 2,
           label = txt$b2_head, size = fH, fontface = "bold", colour = col_title,
           lineheight = 1.15) +
  annotate("text", x = CX, y = b2_cy - LH * (count_lines(txt$b2_head) - 1) / 2,
           label = txt$b2_body, size = fB, colour = "#333333", lineheight = 1.15) +

  annotate("segment", x = CX,  xend = EX1, y = b2_bot, yend = b2_bot,
           colour = brd_excl, linewidth = 0.5, arrow = arr) +
  annotate("segment", x = CX,  xend = CX,  y = b2_bot, yend = b3_top,
           colour = col_arrow, linewidth = 0.5, arrow = arr) +

  ## -- EXCLUSION A -----------------------------------------------------------
  annotate("rect", xmin = EX1, xmax = EX2, ymin = eA_bot, ymax = eA_top,
           fill = fill_excl, colour = brd_excl, linewidth = 0.7) +
  annotate("text", x = ECX, y = eA_top - PAD - LH * 0.5,
           label = txt$eA_head, size = fE, fontface = "bold", colour = brd_excl) +
  annotate("text", x = ECX,
           y = eA_top - PAD - LH * 1.5 - LH * count_lines(txt$eA_body) / 2,
           label = txt$eA_body, size = fe, colour = "#333333", lineheight = 1.15) +

  ## -- BOX 3 -----------------------------------------------------------------
  annotate("rect", xmin = BX1, xmax = BX2, ymin = b3_bot, ymax = b3_top,
           fill = fill_elig, colour = brd_elig, linewidth = 0.7) +
  annotate("text", x = CX, y = b3_top - PAD - LH * 0.5,
           label = txt$b3_head, size = fH, fontface = "bold", colour = col_title) +
  annotate("text", x = CX, y = b3_cy,
           label = txt$b3_body, size = fB, colour = "#333333", lineheight = 1.15) +

  annotate("segment", x = CX, xend = CX, y = b3_bot, yend = b4_top,
           colour = col_arrow, linewidth = 0.5, arrow = arr) +

  ## -- BOX 4 -----------------------------------------------------------------
  annotate("rect", xmin = BX1, xmax = BX2, ymin = b4_bot, ymax = b4_top,
           fill = fill_stage, colour = brd_stage, linewidth = 0.7) +
  annotate("text", x = CX, y = b4_top - PAD - LH * 0.5,
           label = txt$b4_head, size = fH, fontface = "bold", colour = "#3A1A6E") +
  annotate("text", x = CX, y = b4_cy,
           label = txt$b4_body, size = fB, colour = "#333333", lineheight = 1.15) +

  annotate("segment", x = CX,  xend = EX1, y = b4_bot, yend = b4_bot,
           colour = brd_excl, linewidth = 0.5, arrow = arr) +
  annotate("segment", x = CX,  xend = CX,  y = b4_bot, yend = b4b_top,
           colour = col_arrow, linewidth = 0.5, arrow = arr) +

  ## -- EXCLUSION B -----------------------------------------------------------
  annotate("rect", xmin = EX1, xmax = EX2, ymin = eB_bot, ymax = eB_top,
           fill = fill_excl, colour = brd_excl, linewidth = 0.7) +
  annotate("text", x = ECX, y = eB_top - PAD - LH * 0.5,
           label = txt$eB_head, size = fE, fontface = "bold", colour = brd_excl) +
  annotate("text", x = ECX,
           y = eB_top - PAD - LH * 1.5 - LH * count_lines(txt$eB_body) / 2,
           label = txt$eB_body, size = fe, colour = "#333333", lineheight = 1.15) +

  ## -- BOX 4b: Stage 2 / Analytical Sample -----------------------------------
  annotate("rect", xmin = BX1, xmax = BX2, ymin = b4b_bot, ymax = b4b_top,
           fill = fill_final, colour = brd_final, linewidth = 0.9) +
  annotate("text", x = CX, y = b4b_top - PAD - LH * 0.5,
           label = txt$b4b_head, size = fH, fontface = "bold", colour = "#3A1A6E",
           lineheight = 1.15) +
  annotate("text", x = CX, y = b4b_cy,
           label = txt$b4b_body, size = fB, fontface = "bold",
           colour = "#1A4D1A", lineheight = 1.3) +
  annotate("text", x = CX,
           y = b4b_bot + PAD + LH * count_lines(txt$b4b_foot) / 2,
           label = txt$b4b_foot, size = fH + 1, fontface = "bold",
           colour = "#1A4D1A") +

  annotate("segment", x = CX, xend = EX1, y = b4b_bot, yend = b4b_bot,
           colour = brd_excl, linewidth = 0.5, arrow = arr) +

  ## -- EXCLUSION C -----------------------------------------------------------
  annotate("rect", xmin = EX1, xmax = EX2, ymin = eC_bot, ymax = eC_top,
           fill = fill_excl, colour = brd_excl, linewidth = 0.7) +
  annotate("text", x = ECX, y = eC_top - PAD - LH * 0.5,
           label = txt$eC_head, size = fE, fontface = "bold", colour = brd_excl) +
  annotate("text", x = ECX,
           y = eC_top - PAD - LH * 1.5 - LH * count_lines(txt$eC_body) / 2,
           label = txt$eC_body, size = fe, colour = "#333333", lineheight = 1.15) +

  ## -- TITLES ----------------------------------------------------------------
  labs(
    title    = "Flow Diagram for OA-level Sample Construction",
    subtitle = "Output Areas (OAs) \u00b7 England only",
    caption  = "CAZ = Clean Air Zone  \u00b7  7 English schemes"
  ) +

  coord_cartesian(xlim = c(0, 100),
                  ylim = c(canvas_bot, canvas_top), expand = FALSE) +
  theme_void(base_size = 18) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(18, 18, 18, 18),
    plot.title      = element_text(size = 26, face = "bold", colour = col_title,
                                   hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle   = element_text(size = 18, colour = "#555555",
                                   hjust = 0.5, margin = margin(b = 12)),
    plot.caption    = element_text(size = 15, colour = "#888888",
                                   hjust = 0.5, margin = margin(t = 12))
  )

ggsave(file.path(outdir, "fig15_sample_flow_diagram_england.png"),
       plot = p, width = 22, height = 32, dpi = 300, bg = "white")
cat("Saved:", file.path(outdir, "fig15_sample_flow_diagram_england.png"), "\n")
