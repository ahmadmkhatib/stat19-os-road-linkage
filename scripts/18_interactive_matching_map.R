# =============================================================================
# SCRIPT 18 — INTERACTIVE MAP: TREATED OAs & MATCHED CONTROLS
# =============================================================================
#
# PURPOSE:
#   Shiny + Leaflet interactive map. Select any treated OA from the dropdown
#   (or click it on the map) and all matched controls assigned to it
#   light up in blue. The selected treated OA is highlighted in red.
#
# HOW TO RUN:
#   source this file in RStudio, or:
#   shiny::runApp("scripts/18_interactive_matching_map.R")
#
# INPUTS (data/processed/):
#   OA_matched_full_mixed.rds      — full matched dataset (treated + controls + covariates)
#   OA_matched_treated_mixed.rds   — treated OA IDs, weights, stratum
#   OA_matched_donors_mixed.rds    — control OA IDs, weights
#   OA_matching_census.rds         — full pool (LAD codes/names, country, etc.)
#   OA_matching_pairs_mixed.rds    — exact treated→control pairs (saved by script 16)
#
# SOURCE: 16_Matching_England_othercityControlsScotland_mix.R
#   England: other-city controls only
#   Scotland: other-city + same-city controls
#
# INPUTS (data/spatial/):
#   OA_boundaries_GB.gpkg
#   LAD_boundaries_GB.gpkg
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(leaflet)
  library(sf)
  library(tidyverse)
  library(here)
  library(MASS)        # ginv fallback for singular covariance
})

select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate

cat("================================================================\n")
cat("  Script 18 — Interactive Matching Map\n")
cat("================================================================\n\n")

# =============================================================================
# SECTION 1 — LOAD PROCESSED DATA
# =============================================================================

cat("[1/5] Loading processed data...\n")

matched_A  <- readRDS(here("data", "processed", "OA_matched_full_mixed.rds"))
treated_A  <- readRDS(here("data", "processed", "OA_matched_treated_mixed.rds"))
controls_A <- readRDS(here("data", "processed", "OA_matched_donors_mixed.rds"))
full_data  <- readRDS(here("data", "processed", "OA_matching_census.rds"))

# Helper: add country from LAD24CD prefix if not already present
add_country <- function(df) {
  if ("country" %in% names(df)) return(df)
  if ("LAD24CD" %in% names(df)) {
    return(df |> mutate(
      country = case_when(
        substr(LAD24CD, 1, 1) == "E" ~ "England",
        substr(LAD24CD, 1, 1) == "S" ~ "Scotland",
        TRUE                         ~ "Unknown"
      )
    ))
  }
  mutate(df, country = "Unknown")
}

matched_A <- add_country(matched_A)
full_data  <- add_country(full_data)

# Recode "EnglandWales" → "England" (data artefact — no Wales in this study)
recode_country <- function(df) {
  if ("country" %in% names(df))
    df$country <- ifelse(df$country == "EnglandWales", "England", df$country)
  df
}
matched_A <- recode_country(matched_A)
full_data  <- recode_country(full_data)

cat("  Treated OAs :", nrow(treated_A), "\n")
cat("  Control OAs :", nrow(controls_A), "\n\n")

# =============================================================================
# SECTION 2 — MATCHING PAIRS (load cached or reconstruct)
# =============================================================================

cat("[2/5] Matching pairs...\n")

pairs_path <- here("data", "processed", "OA_matching_pairs_mixed.rds")

if (file.exists(pairs_path)) {

  pairs_A <- readRDS(pairs_path)
  cat("  Loaded exact pairs from script 16:", nrow(pairs_A),
      "treated-control assignments\n")

} else {

  # Fallback: reconstruct approximate pairs if the file is missing.
  # This should not normally occur — re-run script 16 to regenerate.
  cat("  WARNING: OA_matching_pairs_mixed.rds not found.\n")
  cat("  Re-run 16_Matching_England_othercityControlsScotland_mix.R to generate it.\n")
  cat("  Falling back to Mahalanobis reconstruction (approximate)...\n\n")

  s2_trends <- c(
    "trend_car_KSI_pkm",   "trend_car_slight_pkm",
    "trend_cyc_KSI_pkm",   "trend_cyc_slight_pkm",
    "trend_ped_KSI_pkm",   "trend_ped_slight_pkm",
    "trend_other_KSI_pkm", "trend_other_slight_pkm",
    "trend_total_pkm"
  )
  s2_levels_log <- paste0("log1p_", c(
    "mean_car_KSI_pkm",   "mean_car_slight_pkm",
    "mean_cyc_KSI_pkm",   "mean_cyc_slight_pkm",
    "mean_ped_KSI_pkm",   "mean_ped_slight_pkm",
    "mean_other_KSI_pkm", "mean_other_slight_pkm",
    "mean_total_pkm"
  ))
  s2_vars <- intersect(c(s2_trends, s2_levels_log), names(matched_A))

  if (length(s2_vars) < 3)
    stop("Too few Stage 2 variables found in OA_matched_full_mixed.rds (",
         length(s2_vars), "). Re-run script 16.")

  cat("  Stage 2 vars available:", length(s2_vars), "/",
      length(c(s2_trends, s2_levels_log)), "\n")

  K <- max(1L, min(round(sum(controls_A$weights) / nrow(treated_A)), 10L))
  cat("  Estimated ratio K =", K, "\n\n")

  treated_rows <- matched_A |>
    filter(treat_indicator == 1) |>
    select(OA, country, all_of(s2_vars)) |>
    drop_na(all_of(s2_vars))

  control_rows <- matched_A |>
    filter(treat_indicator == 0) |>
    select(OA, country, all_of(s2_vars)) |>
    drop_na(all_of(s2_vars))

  cat("  Treated (complete):", nrow(treated_rows),
      "| Controls (complete):", nrow(control_rows), "\n")

  S     <- cov(as.matrix(treated_rows[, s2_vars]), use = "pairwise.complete.obs")
  S     <- S + diag(1e-8 * max(diag(S), na.rm = TRUE), ncol(S))
  S_inv <- tryCatch(solve(S), error = function(e) {
    message("  Covariance singular — using pseudo-inverse (MASS::ginv)")
    MASS::ginv(S)
  })

  t_mat <- as.matrix(treated_rows[, s2_vars])
  c_mat <- as.matrix(control_rows[, s2_vars])
  n_t   <- nrow(treated_rows)

  cat("  Computing distances for", n_t, "treated OAs...\n")
  step <- max(1L, floor(n_t / 10L))

  pairs_A <- purrr::map_df(seq_len(n_t), function(i) {
    if (i %% step == 0)
      cat(sprintf("    %d / %d (%.0f%%)\n", i, n_t, 100 * i / n_t))
    t_id      <- treated_rows$OA[i]
    t_country <- treated_rows$country[i]
    t_vec     <- t_mat[i, ]
    c_idx     <- which(control_rows$country == t_country)
    if (length(c_idx) == 0L) return(tibble())
    D     <- sweep(c_mat[c_idx, , drop = FALSE], 2, t_vec, "-")
    dist  <- rowSums((D %*% S_inv) * D)
    k_here <- min(K, length(c_idx))
    top   <- order(dist)[seq_len(k_here)]
    tibble(treated_OA = t_id,
           control_OA = control_rows$OA[c_idx[top]],
           mdist      = round(dist[top], 4))
  })

  cat("  Pairs reconstructed (approximate):", nrow(pairs_A), "\n")
  saveRDS(pairs_A, pairs_path)
  cat("  Cached to:", pairs_path, "\n\n")
}

stopifnot(
  "pairs_A needs treated_OA" = "treated_OA" %in% names(pairs_A),
  "pairs_A needs control_OA" = "control_OA" %in% names(pairs_A)
)

# =============================================================================
# SECTION 3 — LOAD & PREPARE SPATIAL DATA
# =============================================================================

cat("[3/5] Loading spatial data...\n")

oa_path  <- here("data", "processed", "shp_files", "OA_subset.shp")
lad_path <- here("data", "processed", "shp_files", "LADs_filtered.shp")


# --- OA boundaries ---
oa_raw <- st_read(oa_path, quiet = TRUE)
oa_id  <- intersect(c("OA21CD", "OA11CD", "geo_code", "OA"), names(oa_raw))[1]

oa_raw <- oa_raw |> rename(OA = all_of(oa_id))

analysis_oas <- union(treated_A$OA, controls_A$OA)

oa_sf <- oa_raw |>
  filter(OA %in% analysis_oas) |>
  st_transform(4326) |>
  st_simplify(preserveTopology = TRUE, dTolerance = 0.0001) |>
  st_make_valid()

# --- LAD boundaries ---
lad_raw <- st_read(lad_path, quiet = TRUE)
lad_id  <- intersect(c("LAD21CD","LAD22CD","LAD24CD","lad_code"), names(lad_raw))[1]
lad_nm  <- intersect(c("LAD21NM","LAD22NM","LAD24NM","lad_name"),  names(lad_raw))[1]


lad_sf <- lad_raw |>
  rename(LAD_CODE = all_of(lad_id),
         LAD_NAME = all_of(lad_nm)) |>
  st_transform(4326) |>
  st_simplify(preserveTopology = TRUE, dTolerance = 0.0005) |>
  st_make_valid()

cat("  OA polygons loaded  :", nrow(oa_sf), "/", length(analysis_oas), "\n")
cat("  LAD polygons loaded :", nrow(lad_sf), "\n\n")

# =============================================================================
# SECTION 4 — LOOKUP TABLES & ANNOTATED LAYERS
# =============================================================================

cat("[4/5] Building lookup tables and map layers...\n")

# OA → LAD name / country lookup
oa_meta <- full_data |>
  select(OA, LAD24CD, any_of("LAD24NM"), country) |>
  distinct()

# If LAD name missing from full_data, pull from LAD boundaries
if (!"LAD24NM" %in% names(oa_meta)) {
  lad_names_df <- lad_sf |>
    st_drop_geometry() |>
    select(LAD_CODE, LAD_NAME) |>
    distinct()
  oa_meta <- oa_meta |>
    left_join(lad_names_df, by = c("LAD24CD" = "LAD_CODE")) |>
    rename(LAD24NM = LAD_NAME)
}

oa_meta <- oa_meta |>
  mutate(lad_label = coalesce(LAD24NM, LAD24CD, "Unknown LAD"))

# Controls per treated OA (for info panel)
ctrl_count <- pairs_A |>
  count(treated_OA, name = "n_controls")

# Dropdown choices: labelled by OA + city + country
treated_info <- treated_A |>
  left_join(oa_meta, by = "OA") |>
  left_join(ctrl_count, by = c("OA" = "treated_OA")) |>
  mutate(
    n_controls     = replace_na(n_controls, 0L),
    dropdown_label = paste0(OA, "  [", lad_label, ", ", coalesce(country, "?"), "]")
  ) |>
  arrange(lad_label, OA)

# Annotated treated OA geometry layer
oa_treated <- oa_sf |>
  filter(OA %in% treated_A$OA) |>
  left_join(oa_meta, by = "OA") |>
  left_join(
    treated_A |> select(OA, baseline_injury_stratum),
    by = "OA"
  ) |>
  mutate(
    popup = paste0(
      "<b>Treated OA:</b> ", OA, "<br>",
      "<b>LAD:</b> ",  coalesce(lad_label, "—"), "<br>",
      "<b>Country:</b> ", coalesce(country, "—"), "<br>",
      "<b>Injury stratum:</b> ",
      coalesce(as.character(baseline_injury_stratum), "—")
    )
  )

# Annotated control OA geometry layer
oa_control <- oa_sf |>
  filter(OA %in% controls_A$OA) |>
  left_join(oa_meta, by = "OA") |>
  left_join(controls_A |> select(OA, weights), by = "OA") |>
  mutate(
    popup = paste0(
      "<b>Control OA:</b> ", OA, "<br>",
      "<b>LAD:</b> ", coalesce(lad_label, "—"), "<br>",
      "<b>Country:</b> ", coalesce(country, "—"), "<br>",
      "<b>Match weight:</b> ", round(coalesce(weights, 1), 2)
    )
  )

cat("  Treated OAs with geometry:", nrow(oa_treated), "\n")
cat("  Control OAs with geometry:", nrow(oa_control), "\n\n")

# =============================================================================
# SECTION 5 — COLOUR CONSTANTS
# =============================================================================

COL_T_SEL     <- "#D85A30"   # Selected treated OA            (rust red)
COL_T_DIM     <- "#F2C4B0"   # All other treated OAs          (pale orange)
COL_C_HI      <- "#2E6FAB"   # Matched controls for selection  (blue)
COL_C_DIM     <- "#AECCE8"   # All other controls              (pale blue)
COL_BORDER_T  <- "#8B2010"   # Border: selected treated
COL_BORDER_C  <- "#1A4F8A"   # Border: highlighted controls
COL_LAD       <- "#888888"   # LAD boundary lines


cat("[5/5] Launching Shiny app...\n\n")

# =============================================================================
# SCHEME LOOKUP
# =============================================================================

scheme_lookup <- full_data[full_data$treated_OA == 1,
                            c("OA", "scheme", "LAD24NM", "LAD24CD")] |>
  as_tibble() |>
  dplyr::filter(!is.na(scheme), scheme != "") |>
  dplyr::distinct()

scheme_oa_counts <- scheme_lookup |>
  dplyr::count(scheme, name = "n_treated") |>
  dplyr::arrange(scheme)

scheme_choices <- setNames(
  scheme_oa_counts$scheme,
  paste0(scheme_oa_counts$scheme, "  (", scheme_oa_counts$n_treated, " treated OAs)")
)

cat("  Schemes available:", length(scheme_choices), "\n\n")

# =============================================================================
# SHINY UI
# =============================================================================

ui <- page_sidebar(

  title = "Interactive Matching Map — Treated OAs & Matched Controls",

  theme = bs_theme(bootswatch = "flatly", primary = "#1A2E5A"),

  sidebar = sidebar(
    width = 340,
    open  = "open",

    # ---- Scheme selector ----
    card(
      card_header(class = "bg-primary text-white",
                  icon("city"), " View by Scheme"),
      card_body(
        p(class = "text-muted small mb-2",
          "Select a CAZ/LEZ scheme to highlight all treated OAs and
           their combined matched controls."),
        selectInput(
          inputId  = "selected_scheme",
          label    = NULL,
          choices  = c("--- none ---" = "", scheme_choices),
          selected = "",
          width    = "100%"
        ),
        uiOutput("scheme_info_panel"),
        actionButton("btn_zoom_scheme", "Zoom to scheme",
                     class = "btn-secondary btn-sm w-100 mt-1")
      )
    ),

    # ---- Individual OA selector (restored single-select) ----
    card(
      card_header(class = "bg-primary text-white",
                  icon("map-pin"), " Select Individual OA"),
      card_body(
        p(class = "text-muted small mb-2",
          "Search by OA code or city name. Click a red polygon on the map
           to jump to it."),
        selectizeInput(
          inputId  = "selected_oa",
          label    = "Treated OA:",
          choices  = setNames(treated_info$OA, treated_info$dropdown_label),
          selected = treated_info$OA[1],
          width    = "100%",
          options  = list(placeholder = "Search OA code or city...",
                          maxOptions  = 150)
        ),
        actionButton("btn_zoom", "Zoom to selected OA",
                     class = "btn-primary btn-sm w-100 mt-1"),
        hr(style = "margin: 8px 0;"),
        p(class = "small text-muted mb-1", style = "font-weight:600;",
          icon("circle-dot"), " Zoom to a matched control:"),
        uiOutput("ctrl_zoom_buttons")
      )
    ),

    # ---- OA info panel ----
    card(
      card_header("Selected OA details"),
      card_body(uiOutput("oa_info_panel"))
    ),

    # ---- Display toggles ----
    card(
      card_header("Display options"),
      card_body(
        checkboxInput("show_all_treated",
                      "Show all other treated OAs (dimmed)", value = TRUE),
        checkboxInput("show_all_controls",
                      "Show all other controls (dimmed)",   value = FALSE)
      )
    ),

    # ---- Legend ----
    card(
      card_header("Legend"),
      card_body(
        HTML(paste0(
          "<div style='font-size:0.82rem; line-height:2.0;'>",
          "<span style='display:inline-block;width:14px;height:14px;background:#D85A30;",
          "border:2px solid #8B2010;border-radius:2px;margin-right:6px;vertical-align:middle;'></span>",
          "Selected / scheme treated OAs<br>",
          "<span style='display:inline-block;width:14px;height:14px;background:#2E6FAB;",
          "border:2px solid #1A4F8A;border-radius:2px;margin-right:6px;vertical-align:middle;'></span>",
          "Matched controls<br>",
          "<span style='display:inline-block;width:14px;height:14px;background:#F2C4B0;",
          "border:1px solid #C07050;border-radius:2px;margin-right:6px;vertical-align:middle;'></span>",
          "Other treated OAs (dimmed)<br>",
          "<span style='display:inline-block;width:14px;height:14px;background:#AECCE8;",
          "border:1px solid #7AA8CC;border-radius:2px;margin-right:6px;vertical-align:middle;'></span>",
          "Other controls (dimmed)</div>"
        ))
      )
    ),

    p(class = "text-muted mt-2", style = "font-size:0.75rem;",
      "Mixed analysis | Two-stage Mahalanobis matching",
      br(), "England: other-city controls | Scotland: other-city + same-city",
      br(), "Pairs: ", if (file.exists(pairs_path)) "exact (from script 16)" else "reconstructed (approx)",
      br(), paste0("Treated n=", nrow(treated_A), " | Controls n=", nrow(controls_A)))
  ),

  # ---- MAIN PANEL: MAP ----
  card(
    full_screen = TRUE,
    card_body(padding = 0,
              leafletOutput("map", width = "100%", height = "calc(100vh - 80px)"))
  )
)

# =============================================================================
# SHINY SERVER
# =============================================================================

server <- function(input, output, session) {

  # ---- Individual OA reactives --------------------------------------------

  sel_oa <- reactive({
    req(input$selected_oa)
    input$selected_oa
  })

  sel_controls <- reactive({
    pairs_A |>
      dplyr::filter(treated_OA == sel_oa()) |>
      dplyr::pull(control_OA)
  })

  sel_meta <- reactive({
    treated_info |> dplyr::filter(OA == sel_oa())
  })

  # ---- Filter OA dropdown to selected scheme --------------------------------
  observe({
    s <- input$selected_scheme

    if (is.null(s) || s == "") {
      # No scheme selected — restore full list
      choices <- setNames(treated_info$OA, treated_info$dropdown_label)
      updateSelectizeInput(session, "selected_oa",
                           choices  = choices,
                           selected = treated_info$OA[1])
    } else {
      # Filter to OAs in this scheme
      scheme_oas  <- scheme_lookup$OA[scheme_lookup$scheme == s]
      scheme_info <- treated_info |> dplyr::filter(OA %in% scheme_oas)
      choices     <- setNames(scheme_info$OA, scheme_info$dropdown_label)
      updateSelectizeInput(session, "selected_oa",
                           choices  = choices,
                           selected = scheme_info$OA[1])
    }
  })

  # ---- Scheme reactives ---------------------------------------------------

  sel_scheme <- reactive({
    s <- input$selected_scheme
    if (is.null(s) || s == "") return(NULL)
    s
  })

  scheme_treated_oas <- reactive({
    s <- sel_scheme()
    if (is.null(s)) return(character(0))
    scheme_lookup$OA[scheme_lookup$scheme == s]
  })

  scheme_control_oas <- reactive({
    t_oas <- scheme_treated_oas()
    if (length(t_oas) == 0) return(character(0))
    pairs_A |>
      dplyr::filter(treated_OA %in% t_oas) |>
      dplyr::pull(control_OA) |>
      unique()
  })

  # ---- Scheme info panel --------------------------------------------------
  output$scheme_info_panel <- renderUI({
    s <- sel_scheme()
    if (is.null(s)) return(NULL)
    t_oas <- scheme_treated_oas()
    c_oas <- scheme_control_oas()
    tags$table(
      class = "table table-sm table-borderless mb-1 mt-1",
      style = "font-size:0.82rem;",
      tags$tbody(
        tags$tr(tags$th("Scheme:"),         tags$td(tags$b(s))),
        tags$tr(tags$th("Treated OAs:"),    tags$td(tags$b(length(t_oas)))),
        tags$tr(tags$th("Matched controls:"), tags$td(tags$b(length(c_oas))))
      )
    )
  })

  # ---- Individual OA info panel -------------------------------------------
  output$oa_info_panel <- renderUI({
    info  <- sel_meta()
    ctrls <- sel_controls()
    if (nrow(info) == 0)
      return(p(class = "text-muted small", "No OA selected."))
    tags$table(
      class = "table table-sm table-borderless mb-0",
      style = "font-size:0.82rem;",
      tags$tbody(
        tags$tr(tags$th("OA Code:"),
                tags$td(tags$code(info$OA))),
        tags$tr(tags$th("City / LAD:"),
                tags$td(coalesce(info$lad_label, "—"))),
        tags$tr(tags$th("Country:"),
                tags$td(coalesce(info$country, "—"))),
        tags$tr(tags$th("Inj. stratum:"),
                tags$td(coalesce(as.character(info$baseline_injury_stratum), "—"))),
        tags$tr(tags$th("Controls:"),
                tags$td(tags$b(length(ctrls))))
      )
    )
  })

  # ---- Per-control zoom buttons -------------------------------------------
  output$ctrl_zoom_buttons <- renderUI({
    ctrl_oas <- sel_controls()
    if (length(ctrl_oas) == 0)
      return(p(class = "text-muted small", "No controls assigned."))
    lapply(seq_along(ctrl_oas), function(i) {
      oa  <- ctrl_oas[i]
      lad <- oa_meta$lad_label[oa_meta$OA == oa][1]
      lbl <- if (length(lad) > 0 && !is.na(lad))
        paste0("Control ", i, "  - ", lad) else paste0("Control ", i, "  - ", oa)
      actionButton(
        inputId = paste0("zoom_ctrl_slot_", i),
        label   = lbl,
        class   = "btn-outline-secondary btn-sm w-100 mb-1",
        onclick = sprintf(
          "Shiny.setInputValue('zoom_ctrl_clicked', '%s', {priority: 'event'})", oa)
      )
    })
  })

  # ---- Base map (rendered once) -------------------------------------------
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$CartoDB.Positron,
                       options = tileOptions(opacity = 0.9)) |>
      addPolylines(data = lad_sf, color = COL_LAD, weight = 0.7,
                   opacity = 0.5, group = "LAD boundaries") |>
      addLayersControl(overlayGroups = "LAD boundaries",
                       options = layersControlOptions(collapsed = TRUE)) |>
      setView(lng = -2.0, lat = 53.5, zoom = 6)
  })

  # ---- Map: individual OA layers ------------------------------------------
  observe({
    oa_sel <- sel_oa()
    c_oas  <- sel_controls()
    proxy  <- leafletProxy("map")
    proxy |>
      clearGroup("bg_treated") |>
      clearGroup("bg_controls") |>
      clearGroup("hi_controls") |>
      clearGroup("sel_treated")

    if (isTRUE(input$show_all_treated)) {
      bg_t <- oa_treated |> dplyr::filter(OA != oa_sel)
      if (nrow(bg_t) > 0)
        proxy |> addPolygons(
          data = bg_t, fillColor = COL_T_DIM, fillOpacity = 0.35,
          color = "#C07050", weight = 0.4, popup = ~popup, label = ~OA,
          layerId = ~paste0("t_", OA), group = "bg_treated",
          highlightOptions = highlightOptions(fillOpacity = 0.60, bringToFront = FALSE))
    }

    if (isTRUE(input$show_all_controls)) {
      bg_c <- oa_control |> dplyr::filter(!OA %in% c_oas)
      if (nrow(bg_c) > 0)
        proxy |> addPolygons(
          data = bg_c, fillColor = COL_C_DIM, fillOpacity = 0.18,
          color = "#7AA8CC", weight = 0.3, popup = ~popup, label = ~OA,
          group = "bg_controls")
    }

    hi_c <- oa_control |> dplyr::filter(OA %in% c_oas)
    if (nrow(hi_c) > 0) {
      hi_c <- hi_c |>
        dplyr::mutate(popup = paste0(popup, "<br><b style='color:", COL_C_HI,
                                     ";'>&#9679; Matched to: ", oa_sel, "</b>"))
      proxy |> addPolygons(
        data = hi_c, fillColor = COL_C_HI, fillOpacity = 0.75,
        color = COL_BORDER_C, weight = 1.8, popup = ~popup, label = ~OA,
        layerId = ~paste0("c_", OA), group = "hi_controls",
        highlightOptions = highlightOptions(fillOpacity = 0.95, weight = 2.5, bringToFront = TRUE))
    }

    sel_t <- oa_treated |> dplyr::filter(OA == oa_sel)
    if (nrow(sel_t) > 0) {
      sel_t <- sel_t |>
        dplyr::mutate(popup = paste0(popup, "<br><b style='color:", COL_T_SEL,
                                     ";'>&#9679; SELECTED | ", length(c_oas), " controls</b>"))
      proxy |> addPolygons(
        data = sel_t, fillColor = COL_T_SEL, fillOpacity = 0.85,
        color = COL_BORDER_T, weight = 2.5, popup = ~popup, label = ~OA,
        layerId = ~paste0("t_", OA), group = "sel_treated",
        highlightOptions = highlightOptions(fillOpacity = 1.0, weight = 3.5, bringToFront = TRUE))
    }
  })

  # ---- Map: scheme layers -------------------------------------------------
  observe({
    s     <- sel_scheme()
    proxy <- leafletProxy("map")
    proxy |> clearGroup("scheme_treated") |> clearGroup("scheme_controls")
    if (is.null(s) || s == "") return()

    t_oas <- scheme_treated_oas()
    c_oas <- scheme_control_oas()

    sc_c <- oa_control |> dplyr::filter(OA %in% c_oas)
    if (nrow(sc_c) > 0) {
      sc_c <- sc_c |>
        dplyr::mutate(popup = paste0(popup, "<br><b style='color:", COL_C_HI,
                                     ";'>&#9679; Scheme: ", s, "</b>"))
      proxy |> addPolygons(
        data = sc_c, fillColor = COL_C_HI, fillOpacity = 0.70,
        color = COL_BORDER_C, weight = 1.5, popup = ~popup, label = ~OA,
        layerId = ~paste0("c_", OA), group = "scheme_controls",
        highlightOptions = highlightOptions(fillOpacity = 0.95, weight = 2.5, bringToFront = TRUE))
    }

    sc_t <- oa_treated |> dplyr::filter(OA %in% t_oas)
    if (nrow(sc_t) > 0) {
      n_ctrl_per <- pairs_A |>
        dplyr::filter(treated_OA %in% t_oas) |>
        dplyr::count(treated_OA, name = "n_ctrl")
      sc_t <- sc_t |>
        dplyr::left_join(n_ctrl_per, by = c("OA" = "treated_OA")) |>
        dplyr::mutate(
          n_ctrl = replace_na(n_ctrl, 0L),
          popup  = paste0(popup, "<br><b style='color:", COL_T_SEL,
                          ";'>&#9679; Scheme: ", s, " | ", n_ctrl, " controls</b>"))
      proxy |> addPolygons(
        data = sc_t, fillColor = COL_T_SEL, fillOpacity = 0.85,
        color = COL_BORDER_T, weight = 2.0, popup = ~popup, label = ~OA,
        group = "scheme_treated",
        highlightOptions = highlightOptions(fillOpacity = 1.0, weight = 3.0, bringToFront = TRUE))
    }
  })

  # ---- Click handler: treated polygon or control polygon -------------------
  observeEvent(input$map_shape_click, {
    click_id <- input$map_shape_click$id
    if (is.null(click_id)) return()

    # --- Clicked a TREATED OA: update dropdown ---
    if (startsWith(click_id, "t_")) {
      oa_code <- sub("^t_", "", click_id)
      if (oa_code %in% treated_A$OA)
        updateSelectizeInput(session, "selected_oa", selected = oa_code)
      return()
    }

    # --- Clicked a CONTROL OA: look up which treated OA(s) it was matched to ---
    if (startsWith(click_id, "c_")) {
      ctrl_code <- sub("^c_", "", click_id)

      parents <- pairs_A |>
        dplyr::filter(control_OA == ctrl_code) |>
        dplyr::pull(treated_OA) |>
        unique()

      if (length(parents) == 0) {
        showNotification(
          paste0("Control OA ", ctrl_code, " — no parent found in pairs table."),
          type = "warning", duration = 5
        )
        return()
      }

      # Build LAD labels for each parent
      parent_labels <- vapply(parents, function(p) {
        lad <- treated_info$lad_label[treated_info$OA == p][1]
        if (length(lad) > 0 && !is.na(lad))
          paste0(p, " [", lad, "]")
        else p
      }, character(1))

      if (length(parents) == 1) {
        # Single parent: select it in the dropdown and notify
        msg <- paste0("Control ", ctrl_code,
                      " \u2192 matched to: ", parent_labels)
        showNotification(msg, type = "message", duration = 6)
        updateSelectizeInput(session, "selected_oa", selected = parents)

      } else {
        # Multiple parents: notify with full list (replace=TRUE means shared controls)
        msg <- paste0(
          "Control ", ctrl_code, " is shared across ",
          length(parents), " treated OAs:\n",
          paste(parent_labels, collapse = "\n")
        )
        showNotification(
          tags$div(
            tags$b(paste0("Control ", ctrl_code, " — matched to ", length(parents), " treated OAs:")),
            tags$ul(lapply(parent_labels, tags$li))
          ),
          type = "message", duration = 10
        )
        # Select the first parent in the dropdown
        updateSelectizeInput(session, "selected_oa", selected = parents[1])
      }
    }
  })

  # ---- Zoom: individual control OA ----------------------------------------
  observeEvent(input$zoom_ctrl_clicked, {
    geom <- oa_sf |> dplyr::filter(OA == input$zoom_ctrl_clicked)
    if (nrow(geom) == 0) return()
    bb <- st_bbox(geom)
    leafletProxy("map") |>
      fitBounds(as.numeric(bb["xmin"]) - 0.01, as.numeric(bb["ymin"]) - 0.01,
                as.numeric(bb["xmax"]) + 0.01, as.numeric(bb["ymax"]) + 0.01)
  }, ignoreInit = TRUE)

  # ---- Zoom: individual treated OA ----------------------------------------
  observeEvent(input$btn_zoom, {
    geom <- oa_sf |> dplyr::filter(OA == sel_oa())
    if (nrow(geom) == 0) return()
    bb <- st_bbox(geom)
    leafletProxy("map") |>
      fitBounds(as.numeric(bb["xmin"]) - 0.01, as.numeric(bb["ymin"]) - 0.01,
                as.numeric(bb["xmax"]) + 0.01, as.numeric(bb["ymax"]) + 0.01)
  })

  # ---- Zoom: scheme (all treated + all controls) --------------------------
  observeEvent(input$btn_zoom_scheme, {
    t_oas <- scheme_treated_oas()
    c_oas <- scheme_control_oas()
    if (length(t_oas) == 0) return()
    zoom_sf <- oa_sf |> dplyr::filter(OA %in% c(t_oas, c_oas))
    if (nrow(zoom_sf) == 0) return()
    bb      <- st_bbox(zoom_sf)
    pad_lng <- max(0.01, (as.numeric(bb["xmax"]) - as.numeric(bb["xmin"])) * 0.04)
    pad_lat <- max(0.01, (as.numeric(bb["ymax"]) - as.numeric(bb["ymin"])) * 0.04)
    leafletProxy("map") |>
      fitBounds(
        lng1 = as.numeric(bb["xmin"]) - pad_lng,
        lat1 = as.numeric(bb["ymin"]) - pad_lat,
        lng2 = as.numeric(bb["xmax"]) + pad_lng,
        lat2 = as.numeric(bb["ymax"]) + pad_lat
      )
  })

}

# =============================================================================
# LAUNCH
# =============================================================================

shinyApp(ui, server)
