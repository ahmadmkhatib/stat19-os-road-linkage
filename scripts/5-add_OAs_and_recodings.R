# ============================================================
# Script 04: Recode matched RTIs and assign Output Areas (OA)
# ============================================================

library(sf)
library(tmap)
library(mapview)
library(tidyverse)
library(lubridate)
library(here)

# ------------------------------------------------------------
# Load matched RTI data
# -------------------------------
injuries_matched <- read_rds(here("data", "processed", "injuries_matched.rds"))

# Ensure RTIs are rows (each row = one RTI)
# Convert to sf using injury coordinates (British National Grid)

injuries_matched_sf <- injuries_matched %>%
  st_as_sf(
    coords = c("injury_x", "injury_y"),
    crs = 27700,
    remove = FALSE
  )

# Check CRS
st_crs(injuries_matched_sf)

# ------------------------------------------------------------
# Load Output Area boundaries (England & Wales and scotland 2011) ### scotland does not have 2021 
# ------------------------------------------------------------

oa_2011 <- st_read("../stat19-os-road-linkage-data/infuse_oa_lyr_2011_clipped.shp") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  rename(OA_CODE = geo_code) 

# ------------------------------------------------------------
# Spatial join – assign 1 OA per RTI
# ------------------------------------------------------------

# st_within ensures one OA per point
injuries_with_oa <- injuries_matched_sf %>%
  st_join(
    oa_2011,
    join = st_within,
    left = TRUE
  )

any(duplicated(injuries_with_oa$injury_id))
# ---------------------
# some quality Assurance
# ---------------------

injuries_with_oa %>%
  summarise(
    total = n(),
    missing_oa = sum(is.na(OA_CODE)),
    pct_missing = 100 * missing_oa / total
  )


# fix missing using nearest feature
missing_idx <- which(is.na(injuries_with_oa$OA_CODE))

injuries_with_oa$OA_CODE[missing_idx] <-
  oa_2011$OA_CODE[
    st_nearest_feature(injuries_with_oa[missing_idx, ], oa_2011)
  ]


# -------------------------------
#  time variables
# --------------------------------

injuries_with_oa <- injuries_with_oa %>%
  mutate(
    date = as.Date(date),
    month_year = floor_date(date, "month"),
    quarter_year = floor_date(date, "quarter"),
    year = year(date)
  )


# Recode casualty_type1
# -----------------------------
injuries_with_oa <- injuries_with_oa %>%
  st_drop_geometry() %>%
  rename(identifier = matched_roadID) %>%
  mutate(
    casualty_type1 = if_else(
      casualty_type1 == "Car or van driver or occupant",
      "Car/Van",
      casualty_type1
    )
  )



write_rds(
  injuries_with_oa,
  here("data", "processed", "injuries_matched_OA.rds")
)



# ----------------------------- create a shp file 
dir.create("data/processed/shp_files", recursive = TRUE)

names(injuries_with_oa) <- substr(names(injuries_with_oa), 1, 10)



st_write(injuries_with_oa,
         here("data","processed","shp_files","injuries_final.gpkg"),
         layer = "injuries_final",
         delete_layer = TRUE)


inj <- st_read("data/processed/shp_files/injuries_final.gpkg")




