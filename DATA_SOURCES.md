# Data Sources

This repository does not redistribute raw data.

Users must download the following datasets directly from official sources.

---

## 1. STATS19 Road Safety Data

Source: Department for Transport  
Licence: Open Government Licence v3.0  

Download page:  
https://www.gov.uk/government/statistical-data-sets/road-safety-open-data

Required tables:
- collisions
- casualties
- vehicles

Place extracted CSV files in:
data/raw/

---

## 2. OS Open Roads

Source: Ordnance Survey  
Licence: Open Government Licence  

Download:
https://www.ordnancesurvey.co.uk/products/os-open-roads

Place GeoPackage or Shapefile in:
data/raw/

---

## 3. Output Area Boundaries

Source: Office for National Statistics  
Licence: Open Government Licence  

Download:
https://geoportal.statistics.gov.uk/

Use 2011 Output Area boundaries (EPSG:27700).

Place file in:
data/raw/
