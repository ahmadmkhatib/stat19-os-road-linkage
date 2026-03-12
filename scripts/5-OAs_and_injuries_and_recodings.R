# ============================================================
#Recode matched RTIs and assign Output Areas (OA)
# ============================================================

library(sf)
library(tmap)
library(mapview)
library(tidyverse)
library(lubridate)
library(here)



injuries_matched <- read_rds(here("data", "processed", "injuries_matched.rds"))


# Convert to sf using injury coordinates (British National Grid)
injuries_matched_sf <- injuries_matched %>%
  st_as_sf(
    coords = c("injury_x", "injury_y"),
    crs = 27700,
    remove = FALSE
  )

# Load Output Area boundaries (England and scotland 2011) ### scotland does not have 2021 
# ------------------------------------------------------------

oa_scot <- st_read("../stat19-os-road-linkage-data/OutputArea2022_EoR.shp") %>%
  st_transform(27700) %>%
  st_make_valid() 

oa_eng <- st_read("../stat19-os-road-linkage-data/OA_2021_EW_BFE_V9.shp") %>%
  st_transform(27700) %>%
  st_make_valid()


glimpse(oa_scot)
glimpse(oa_eng)

oa_scot <-oa_scot %>%   select(OA = code, geometry)
oa_eng<-oa_eng %>% select(OA = OA21CD, geometry)

oa<-bind_rows(oa_eng, oa_scot)

# ------------------------------------------------------------
#   join –    
# ------------------------------------------------------------
injuries_with_oa <- injuries_matched_sf %>%
  st_join(
    oa,
    join = st_within,
    left = TRUE
  )

any(duplicated(injuries_with_oa$injury_id))
# ---------------------
# checks

injuries_with_oa %>%
  summarise(
    total = n(),
    missing_oa = sum(is.na(OA)),
    pct_missing = 100 * missing_oa / total
  )


# fix missing using nearest feature
missing_idx <- which(is.na(injuries_with_oa$OA))

injuries_with_oa$OA[missing_idx] <-
  oa$OA[
    st_nearest_feature(injuries_with_oa[missing_idx, ], oa)
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


### save OA combined England and Scot
write_rds(
  oa,
  here("data", "processed", "OAs_comb.rds")
)

st_write(oa,
         here("data","processed","shp_files","OAs_comb.shp"),
         layer = "injuries_final",
         delete_layer = TRUE)
