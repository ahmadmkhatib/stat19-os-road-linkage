# STAT19 – OS Open Roads Spatial Linkage & Area Assignment (UK)

![Reproducible Research](https://img.shields.io/badge/Reproducible-Yes-brightgreen)
![R Version](https://img.shields.io/badge/R-%3E%3D4.3-blue)
![Spatial CRS](https://img.shields.io/badge/CRS-EPSG%3A27700-lightgrey)



This repository provides a spatial linkage framework connecting the UK STATS19 road traffic injuries and collision data to:
UK Road network geometry from the OS Open Roads network links.

Local Authority Districts (LADs)

Output Areas (OA, 2011 boundaries)
 

The framework is designed to create an injury-level dataset matched to road segments and small-area geography to support transport research and policy evaluations.

---

# Overview

This project spatially links individual STATS19 road traffic injury (RTI) records to the most plausible road link in the Ordnance Survey (OS) Open Roads network.

The pipeline:

Selects large UK cities (≥100,000 population)

Identifies intervention and control areas

Filters OS Open Roads to selected LADs

Constructs injury-level STATS19 dataset (2015+)

Data Quality Assurance 

Spatially matches injuries to plausible road segments

Assigns Output Areas

Runs validation checks

All spatial operations use British National Grid (EPSG:27700).

The objective is to create a clean, reproducible GB injury–road level dataset
---

# Repository Structure

```
data/
 ├── raw/                # Raw input datasets (not versioned)
 ├── processed/          # Derived RDS and spatial outputs
 │    ├── injuries_final.rds
 │    ├── injuries_matched.rds
 │    ├── injuries_with_oa.rds
 │    ├── roads_filtered.rds
 │    ├── LADs_sub.rds
 │    └── validation_summary.rds
scripts/
 ├── 00_create_large_city_LAD_subset.R
 ├── 01_prepare_os_roads.R
 ├── 02_prepare_stats19_injuries.R
 ├── 03_match_injuries_to_roads_by_type.R
 ├── 04_assign_output_areas.R
 └── 05_validation_checks.R

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

- Load the OS roads shapefile  
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
data/processed/injuries_final.rds
```

This file contains one spatial record per injury.

---

## 03_match_injuries_to_roads_by_type.R

## Purpose

Link each injury to the nearest most plausible road link using hierarchy-aware spatial matching.

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
## Matching uses:

st_nearest_feature()

st_distance()

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


---
#  Script 05: validation checks on the injuries_with_oa dataset
1. Purpose:
Formal QA and validation of linkage pipeline.

Checks Performed:

Duplicate injury_id

Empty geometries

Invalid geometries

CRS verification

OA assignment completeness

Road distance plausibility (>100m flagged)

Date range consistency

Future-dated injuries

Missingness summary

Output:

validation_summary.rds

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

sf
tidyverse
lubridate
readxl
stringr
here
tmap


Install using:

```r
install.packages(c("sf","tidyverse","tmap","stats19","stringr","lubridate","here"))
```

---

# Reproducing the Full Pipeline
Place raw datasets into:
```r
data/raw/
```
Run scripts in order:

```r
00_create_large_city_LAD_subset.R
01_prepare_os_roads.R
02_prepare_stats19_injuries.R
03_match_injuries_to_roads_by_type.R
04_assign_output_areas.R
05_validation_checks.R

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


# Data Dictionary

Injury-Level Analytical Dataset (injuries_OA.rds)


| Variable                | Description                                                     | Source                     |
| ----------------------- | --------------------------------------------------------------- | -------------------------- |
| injury_id               | Unique injury identifier (collision_index + casualty_reference) | Derived                    |
| collision_index         | Unique collision identifier                                     | STATS19                    |
| casualty_reference      | Casualty ID within collision                                    | STATS19                    |
| date                    | Date of collision                                               | STATS19                    |
| year                    | Calendar year                                                   | Derived                    |
| month                   | Calendar month (1–12)                                           | Derived                    |
| month_year              | Year–month identifier                                           | Derived                    |
| quarter                 | Calendar quarter (1–4)                                          | Derived                    |
| quarter_year            | Year–quarter identifier                                         | Derived                    |
| severity                | Binary indicator: KSI (Killed/Seriously Injured) vs Slight      | Derived (STATS19)          |
| casualty_severity       | Original STATS19 severity code                                  | STATS19                    |
| casualty_type           | Pedestrian / Cyclist / Vehicle occupant (harmonised)            | STATS19                    |
| age                     | Casualty age                                                    | STATS19                    |
| sex                     | Casualty sex                                                    | STATS19                    |
| vehicle_type            | Simplified vehicle class                                        | STATS19 (harmonised)       |
| propulsion_type         | Fuel/propulsion category                                        | STATS19                    |
| junction                | Binary indicator of junction involvement                        | STATS19                    |
| road_class_police       | Road class recorded in STATS19                                  | STATS19                    |
| speed_limit             | Posted speed limit (mph)                                        | STATS19                    |
| easting                 | British National Grid easting                                   | STATS19                    |
| northing                | British National Grid northing                                  | STATS19                    |
| geometry                | Injury point geometry (EPSG:27700)                              | Derived                    |
| road_link_id            | OS Open Roads unique segment identifier                         | OS Open Roads              |
| road_class_matched      | Road class of matched segment                                   | OS Open Roads              |
| road_name               | Road name (if available)                                        | OS Open Roads              |
| trunk_road              | Indicator for trunk road (if available)                         | OS Open Roads              |
| road_distance_m         | Distance from injury to matched road segment (metres)           | Derived (spatial matching) |
| match_type              | Matching type (same-class ≤50m / fallback ≤100m)                | Derived                    |
| LAD24CD                 | Local Authority District code                                   | Spatial join               |
| LAD24NM                 | Local Authority District name                                   | Spatial join               |
| OA11CD                  | Output Area code (2011 boundaries)                              | Spatial join               |
| OA_assignment_method    | Within-OA polygon or nearest-OA centroid                        | Derived                    |
| intervention_area       | Indicator for CAZ/LEZ intervention LAD (if applicable)          | Derived                    |
| intervention_start_year | Policy implementation year (if merged)                          | Derived                    |
| post_policy             | Post-intervention indicator                                     | Derived                    |


```mermaid
flowchart TD

A[STATS19 Raw Data] --> B[Filter 2015+ Collisions]
B --> C[Construct Injury-Level Dataset]
C --> D[Spatial Join to LADs]

E[OS Open Roads] --> F[Filter to Selected LADs]
F --> G[Recode Road Classes]

D --> H[Class-Based Nearest Matching]
G --> H

H --> I[Matched Injury-Road Dataset]
I --> J[Assign Output Areas]

J --> K[Final Analytical Dataset]
K --> L[Validation & QA Report]



