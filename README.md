# STAT19 – OS Open Roads Spatial Linkage

This repository provides a spatial linkage pipeline for matching UK STATS19 road traffic collision data to OS Open Roads network links.

The framework is designed to support road safety research, transport policy evaluation, and spatial injury analysis.

---

# Overview

This project spatially links individual STATS19 road traffic injury (RTI) records to the most plausible road link in the Ordnance Survey (OS) Open Roads network.

The pipeline includes:

- Cleaning and preparing STATS19 collision-level data  
- Cleaning and preparing the OS Open Roads network data  
- Coordinate reference system harmonisation  
- Spatial linkage using:
  - Hierarchy-aware nearest-road matching  
  - Fixed distance thresholds (50m / 100m)  
- Validation and sensitivity checks (e.g., road class agreement)

The objective is to create a clean, reproducible injury–road level dataset suitable for:

- Road safety analysis  
- Road hierarchy comparisons  
- Clean Air Zone (CAZ) and control city evaluations  
- Spatial modelling and mapping  

---

# Repository Structure

```
├── README.md
├── RS_project.Rproj
├── scripts/
│   ├── 01_prepare_roads.R
│   ├── 02_prepare_stats19.R
│   └── 03_match_injuries_to_roads.R
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

- **OS Open Roads network data**  
  Source: Ordnance Survey  

- **Local Authority District boundaries**  
  LAD_DEC_24_UK_BGC.shp  

---

# Data Setup

Create the following folder structure:

```
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

The analysis consists of three scripts that must be run in sequence.

---

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
- Save cleaned road network  

## Output

```
data/processed/roads_final.rds
```

---

# 02_prepare_stats19.R

## Purpose

Prepare STATS19 injury data for spatial matching.

## Steps

1. Import collision, vehicle, and casualty datasets  
2. Join datasets using `collision_index`  
3. Parse dates and filter to 2015 onwards  
4. Recode road hierarchy variables  
5. Convert injury locations to spatial points  
6. Transform to EPSG:27700  
7. Attach Local Authority District (LAD)  
8. Exclude London boroughs (codes starting with E09000)  
9. Create:
   - `injury_id`
   - `injury_class` (Motorway, A, B, minor)

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
data/processed/matchedRTI.rds
```

This dataset contains:

- Injury attributes  
- Matched road identifier  
- Match distance  
- Match type  

---

# Quality Assurance

The pipeline includes checks for:

- Date parsing validity  
- Road class distribution  
- Matching type distribution  
- Proportion of unmatched injuries  
- Distance summaries  

Typical results:

- Majority matched via same-class rule  
- Small proportion matched via fallback  
- Very small percentage (>100m) excluded  

This reflects expected positional imprecision in STATS19 reporting.

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
- Road safety evaluation  
- Transport policy impact analysis  
- Spatial epidemiology studies  

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
