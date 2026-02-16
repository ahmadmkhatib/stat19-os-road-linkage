# ============================================================
# Script 04: Recode matched RTIs and assign Output Areas (OA)
# ============================================================

library(sf)
library(tidyverse)
library(lubridate)
library(here)

# ------------------------------------------------------------
# 1. Load matched RTI data
# ------------------------------------------------------------

matched <- read_rds(here("data", "processed", "RTIs_final.rds"))

# Ensure RTIs are rows (each row = one RTI)
# Convert to sf using injury coordinates (British National Grid)

matched_sf <- matched %>%
  st_as_sf(
    coords = c("injury_x", "injury_y"),
    crs = 27700,
    remove = FALSE
  )

# Check CRS
st_crs(matched_sf)

# ------------------------------------------------------------
# 2. Load Output Area boundaries (England & Wales and scotland 2011) ### scotland does not have 2021 
# ------------------------------------------------------------

oa_2011 <- st_read("data/raw/infuse_oa_lyr_2011_clipped.shp") %>%
  st_transform(27700) %>%
  st_make_valid() %>%
  rename(OA_CODE = geo_code) 

# ------------------------------------------------------------
# Spatial join – assign 1 OA per RTI
# ------------------------------------------------------------

# st_within ensures one OA per point
RTI_DATA <- matched_sf %>%
  st_join(
    oa_2011,
    join = st_within,
    left = TRUE
  )


# ------------------------------------------------------------
# Quality Assurance
# ------------------------------------------------------------

RTI_DATA %>%
  summarise(
    total = n(),
    missing_oa = sum(is.na(OA_CODE)),
    pct_missing = 100 * missing_oa / total
  )


# fix missing using nearest feature
missing_idx <- which(is.na(RTI_DATA$OA_CODE))

RTI_DATA$OA_CODE[missing_idx] <-
  oa_2011$OA_CODE[
    st_nearest_feature(RTI_DATA[missing_idx, ], oa_2011)
  ]


# ------------------------------------------------------------
# additional time variables
# ------------------------------------------------------------

RTI_DATA <- RTI_DATA %>%
  mutate(
    date = as.Date(date),
    month_year = floor_date(date, "month"),
    quarter_year = floor_date(date, "quarter"),
    year = year(date)
  )

# ----------------------------------------------------------

# ------------------------------------------------------------
# --------------------------------------

write_rds(
  RTI_DATA,
  here("data", "processed", "RTIs_final_with_OA.rds")
)

#
