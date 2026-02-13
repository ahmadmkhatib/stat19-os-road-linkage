# STAT19 – OS Open Roads Spatial Linkage Framework

This repository provides a reproducible spatial linkage pipeline for matching UK STATS19 road traffic collision data to OS Open Roads network links.

The framework is designed to support road safety research, transport policy evaluation, and spatial injury analysis.

---

## Overview

The pipeline includes:

- Cleaning and preparing STATS19 collision-level data
- Cleaning and preparing the OS Open Roads network data
- Coordinate reference system harmonisation
- Spatial linkage using:
  - Fixed-distance buffer matching
  - Nearest-road matching
- Validation and sensitivity checks (e.g., road class agreement)

---

## Data Requirements

Due to licensing restrictions, raw datasets are **not included** in this repository.

Users must obtain:

- **STATS19 collision data** from the UK Department for Transport
- **OS Open Roads** from Ordnance Survey

Both datasets must be downloaded independently and stored locally.

---

## Repository Structure

---

## Software Requirements

This project uses R (>= 4.3) and requires packages such as:

- sf
- tidyverse
- data.table
- units
- here

---

## Intended Use

This repository provides a generalisable UK-wide linkage framework.

It is intended for:

- Academic research
- Road safety evaluation
- Policy impact analysis
- Spatial epidemiology studies

---

## Author

Ahmad Alkhatib  
Imperial College London

---

## License

This project is licensed under the MIT License.

