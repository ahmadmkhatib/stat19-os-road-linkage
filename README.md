# STAT19 – OS Open Roads Spatial Linkage

This repository provides a spatial linkage pipeline for matching UK STATS19 road traffic collision data to OS Open Roads network links.

The framework is designed to support transport research and policy evaluations.

---

# Overview

This project spatially links individual STATS19 road traffic injury (RTI) records to the most plausible road link in the Ordnance Survey (OS) Open Roads network.

The pipeline includes:

This project:

1. Defines a reproducible set of large UK cities (≥100,000 population)
2. Create LADs subset
3. Constructs a control LAD pool from the remaining large cities
4. Prepares OS Open Roads within selected LADs
5. Prepares STATS19 injury data
6. Performs hierarchy-aware spatial road linkage
7. Assigns Output Areas (OA)
8. Validation checks

The objective is to create a clean, reproducible GB injury–road level dataset
---

# Repository Structure

```
├── README.md
├── RS_project.Rproj
├── scripts/
│   ├── 00_create_large_city_LAD_subset.R
|   ├── 01_prepare_roads.R
│   ├── 02_prepare_stats19.R
│   └── 03_match_injuries_to_roads_by_type.R
│   ├── 04_Recoding_the_matched_RTI_Data_and_add_OAs.R
│   └── 05_validation_checks_final_data.R

├── data/
│   ├── raw/
│   └── processed/
```

---

# Data Requirements

Due to licensing restrictions, raw datasets are **not included** in this repository.

Users must independently obtain:

- **STATS19 collision, vehicle, and casualty data**  
  Source: UK Department for Transport  

- **OS Open Roads network shapefile**  
  Source: Ordnance Survey  

- **Local Authority District boundaries shapefile**  
 Source: geoportal.statistics.gov.uk

---

# Data Setup

Create the following folder structure:

data/
  raw/
    OS highways all.shp
    LAD_DEC_24_UK_BGC.shp
    dft-road-casualty-statistics-collision-1979-latest.csv
    dft-road-casualty-statistics-vehicle-1979-latest.csv
    dft-road-casualty-statistics-casualty-1979-latest.csv
  processed/
```

All processed outputs will be written to:

```
data/processed/
```

---

# Coordinate Reference System

All spatial analysis is performed in:

British National Grid (EPSG:27700)

This ensures:

- Distance calculations are in metres  
- Accurate nearest-road matching  
- Consistent spatial operations  

---

# Workflow

The analysis consists of  scripts that must be run in sequence.

---

---

# Script 00 – Create Large City LAD Subset

Creates a dataset of UK cities with population ≥100,000 and assigns LAD codes.

### Inputs

- LAD boundaries (ONS)
- Built-Up Areas (England & Wales)
- Scotland settlement population file
- City latitude/longitude file

### Output
data/processed/big_cities_with_LADs.rds

Contains:

- City name  
- LAD code (LAD24CD)  
- LAD name (LAD24NM)  
- Population  
- Country  

This file is used downstream to define control LADs.

---

# Script 01 – Prepare OS Roads

Purpose: Prepare OS Open Roads data within selected LADs.

### Steps

- Load LAD boundaries
- Define treatment LADs (England CAZ + Scotland LEZ)
- Construct Control LADs dynamically
- Subset LAD geometries
- Load OS Open Roads
- Spatially filter roads to selected LADs
- Recode road hierarchy
- Save processed roads

### Output



# 01_prepare_roads.R

## Purpose

Prepare and standardise the OS road network for spatial linkage.

## Steps

- Load OS roads shapefile  
- Clean attribute fields  
- Standardise road classification  
- Create simplified road class variable:
  - Motorway  
  - A  
  - B  
  - minor (C + Unclassified)  
- Transform to EPSG:27700  
- Save the cleaned road network  

## Output

```
data/processed/roads_final.rds
```

---

# 02_prepare_stats19.R

## Purpose

Prepare STATS19 injury data for spatial matching.
Each row in the final dataset represents one casualty (one injury).
The dataset is constructed starting from the casualty table.

Joins are performed as:

Casualty → Collision (by collision_index)

Casualty → Vehicle (by collision_index + vehicle_reference)

Steps

Import collision, vehicle, and casualty datasets

Parse collision dates and filter to 2015 onwards

Select required variables only (memory efficient)

Join:

Collision data to casualties

Vehicle data to casualties using collision_index + vehicle_reference

Create a unique injury identifier:

injury_id = collision_index + "_" + casualty_reference

Recode:

Road class hierarchy

Casualty severity (KSI / Slight)

Vehicle type

Propulsion type

Convert to spatial points using:

location_easting_osgr

location_northing_osgr

Set CRS to EPSG:27700 (British National Grid)

Attach Local Authority District (LAD)

Filter to the selected LAD subset only
## Output

```
data/processed/RTIs_final.rds
```

This file contains one spatial record per injury.

---

# 03_match_injuries_to_roads.R

## Purpose

Link each injury to the most plausible road link using hierarchy-aware spatial matching.

---

## Matching Rules

### 1. Primary Match (Same Class ≤ 50 metres)

- A-road injuries → A roads  
- B-road injuries → B roads  
- Motorway injuries → Motorways  
- C and Unclassified injuries → OS minor roads  

If a same-class road is found within 50 metres, the injury is matched.

---

### 2. Fallback Match (≤ 100 metres, Any Class)

If no same-class road is found within 50 metres,  
the injury is matched to the nearest road of any class within 100 metres.

---

### 3. Drop Rule

Injuries more than 100 metres from any road link are excluded.

---

## Distance Variables Recorded

- `dist_same` – distance to nearest same-class road  
- `dist_any` – distance to nearest road of any class  
- `match_type`:
  - `same_class`
  - `fallback_any`

Each injury is matched to one road link only.

---

## Output

```
data/processed/injuries_matched.rds
```

This dataset contains:

- Injury attributes  
- Matched road identifier  
- Match distance  
- Match type  

---

# Script 04: Recode matched RTIs and assign Output Areas (OA)

This script processes the previously matched RTI dataset (`injuries_matched.rds`) to:

1. Convert each RTI into a spatial point (`sf`) using British National Grid coordinates.
2. Assign **one Output Area (OA)** per RTI using 2011 OA boundaries covering England, Wales, and Scotland.
3. Recode and create additional temporal variables (`month_year`, `quarter_year`, `year`) for time-based analyses.
4. Perform quality assurance checks to ensure all RTIs are assigned an OA, with fallback nearest-neighbour matching for missing points.

---
#  Script 05: validation checks on the injuries_with_oa dataset
1. Duplicats 
2. Spatial geometries
3. CRS 
4. OA coverage  
5. Distance-to-road matches
6. Valid temporal range 


## Workflow

1. **Load matched RTI data**  
   - Each row represents one RTI.
   - Converts injury coordinates into `sf` points (EPSG:27700).

2. **Load Output Area boundaries**  
   - Uses `infuse_oa_lyr_2011_clipped.shp`, which covers England, Wales, and Scotland.
   - CRS harmonisation and geometry validation applied.

3. **Spatial linkage**  
   - Primary: `st_within` ensures each RTI is assigned to the OA polygon it falls inside.
   - Secondary: `st_nearest_feature` assigns the nearest OA for points outside any polygon.

4. **Quality Assurance (QA)**  
   - Checks total RTIs, number of missing OA assignments, and percentage missing.
   - After nearest-neighbour fallback, the dataset should have minimal or no missing OA codes.

5. **Add temporal variables**  
   - `month_year` → first day of the month of the collision.
   - `quarter_year` → first day of the quarter of the collision.
   - `year` → collision year.

6. **Output**  
   - Final dataset: `injuries_with_OA.rds`  
   - Each row = one RTI, including spatial coordinates, matched road link, OA code, and recoded variables.

---

## Notes

- OA boundaries are based on **2011 UK census geography**, as 2021 coverage does not include Scotland.  
- The dataset is ready for further spatial analysis, including linkage to local authority or small-area indicators.
- All operations are performed in **EPSG:27700 (British National Grid)** to maintain consistency with the road network and injury locations.

---

---

# Software Requirements

R (>= 4.3 recommended)

Required packages:

- sf  
- tidyverse  
- data.table  
- stats19  
- stringr  
- lubridate  
- here  

Install using:

```r
install.packages(c("sf","tidyverse","data.table","stats19","stringr","lubridate","here"))
```

---

# Reproducing the Full Pipeline

Run scripts in order:

```r
source("scripts/01_prepare_roads.R")
source("scripts/02_prepare_stats19.R")
source("scripts/03_match_injuries_to_roads.R")
```

Processed outputs will appear in:

```
data/processed/
```

---

# Intended Use

This repository provides a generalisable UK-wide linkage framework.

It is intended for:

- Academic research  
- Road safety evaluations


---

# Author

Ahmad Alkhatib  
Public Health and Policy Evaluation Unit  
Imperial College London  

---

# License

This project is licensed under the MIT License.

Raw datasets are subject to their respective licences:

- Department for Transport STATS19 licence  
- Ordnance Survey Open Government Licence  
