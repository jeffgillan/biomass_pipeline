# itcSegment v1.0 API Reference for Claude Code

> **Purpose**: This document provides the exact function signatures, parameter codes,
> and usage patterns for the itcSegment R package. This package is the ONLY actively
> maintained CRAN package that provides built-in allometric DBH and AGB estimation
> from ALS-derived metrics using the Jucker et al. (2017) global equations.
>
> **CRITICAL NOTE ON PARAMETER NAMING**: The `dbh()` and `agb()` functions use `CA`
> as the parameter name, but this refers to **Crown Diameter in meters**, NOT Crown Area.
> This is a known source of confusion. Always pass crown diameter, not crown area.

---

## 1. `dbh()` — Diameter at Breast Height Prediction

Predicts DBH from tree height and crown diameter using Jucker et al. (2017) allometric equations.

### Signature

```r
dbh(H = NULL, CA = NULL, biome = 0)
```

### Parameters

| Parameter | Type    | Description                               | Units    |
|-----------|---------|-------------------------------------------|----------|
| `H`       | Numeric | Tree height                               | **meters** |
| `CA`      | Numeric | Crown diameter (**NOT crown area**)       | **meters** |
| `biome`   | Integer | Biome code (0–24), see table below        | —        |

### Returns

Numeric vector of DBH values in **centimeters**.

### Biome Codes

| Code | Biome Description                                          |
|------|------------------------------------------------------------|
| 0    | **Global** (general-purpose, use when unsure)              |
| 1    | Afrotropic – Tropical forests – Angiosperm                 |
| 2    | Afrotropic – Woodlands and savannas – Angiosperm           |
| 3    | Australasia – Temperate mixed forests – Angiosperm         |
| 4    | Australasia – Temperate mixed forests – Gymnosperm         |
| 5    | Australasia – Woodlands and savannas – Angiosperm          |
| 6    | Indo-Malaya – Tropical forests – Angiosperm                |
| 7    | **Nearctic – Boreal forests – Angiosperm**                 |
| 8    | **Nearctic – Boreal forests – Gymnosperm**                 |
| 9    | **Nearctic – Temperate coniferous forests – Angiosperm**   |
| 10   | **Nearctic – Temperate coniferous forests – Gymnosperm**   |
| 11   | **Nearctic – Temperate mixed forests – Angiosperm**        |
| 12   | **Nearctic – Temperate mixed forests – Gymnosperm**        |
| 13   | **Nearctic – Woodlands and savannas – Angiosperm**         |
| 14   | **Nearctic – Woodlands and savannas – Gymnosperm**         |
| 15   | Neotropic – Tropical forests – Angiosperm                  |
| 16   | Palearctic – Boreal forests – Angiosperm                   |
| 17   | Palearctic – Boreal forests – Gymnosperm                   |
| 18   | Palearctic – Temperate coniferous forests – Angiosperm     |
| 19   | Palearctic – Temperate coniferous forests – Gymnosperm     |
| 20   | Palearctic – Temperate mixed forests – Angiosperm          |
| 21   | Palearctic – Temperate mixed forests – Gymnosperm          |
| 22   | Palearctic – Woodlands and savannas – Angiosperm           |
| 23   | Palearctic – Woodlands and savannas – Gymnosperm           |
| 24   | Palearctic – Woodlands and savannas – Gymnosperm           |

**Biome selection guidance for US forests:**
- Codes 7–14 (bolded) cover North American (Nearctic) biomes
- For the US Southwest (e.g., Arizona): code 13 or 14 (Woodlands and savannas)
- For Pacific Northwest conifers: code 10 (Temperate coniferous – Gymnosperm)
- For eastern deciduous forests: code 11 (Temperate mixed – Angiosperm)
- For southeastern pine plantations: code 10 or 12
- When forest type is mixed or unknown: code 0 (Global)

### Internal Model

The function implements: `DBH = exp(a + b × ln(H) + c × ln(CD))`

Where a, b, c are biome-specific coefficients from the Jucker et al. (2017) lookup table.
For biome 0 (Global): a = 0.557, b = 0.809, c = 0.056.

### Example

```r
library(itcSegment)

# Single tree: 25m height, 8m crown diameter
dbh(H = 25, CA = 8, biome = 0)
# Returns DBH in cm

# Vectorized: works on vectors from lidR crown_metrics output
metrics$DBH_cm <- dbh(
  H = metrics$Z,                           # Tree height from lidR (meters)
  CA = 2 * sqrt(metrics$convhull_area / pi), # Crown diameter from convex hull area
  biome = 10                                # Nearctic temperate coniferous gymnosperm
)
```

---

## 2. `agb()` — Aboveground Biomass Prediction

Predicts AGB from tree height and crown diameter using Jucker et al. (2017) equations.

### Signature

```r
agb(H = NULL, CA = NULL, species = 1)
```

### Parameters

| Parameter | Type    | Description                             | Units    |
|-----------|---------|---------------------------------------- |----------|
| `H`       | Numeric | Tree height                             | **meters** |
| `CA`      | Numeric | Crown diameter (**NOT crown area**)     | **meters** |
| `species` | Integer | Species group code                      | —        |

### Species Codes

| Code | Species Group  |
|------|----------------|
| 1    | **Gymnosperm** (conifers: pine, spruce, fir, etc.) |
| 2    | **Angiosperm** (broadleaf: oak, maple, etc.)       |

### Returns

Numeric vector of AGB values in **kilograms**.

### Example

```r
# Predict AGB for conifer trees
metrics$AGB_kg <- agb(
  H = metrics$Z,
  CA = 2 * sqrt(metrics$convhull_area / pi),
  species = 1  # gymnosperm
)

# Convert to Mg (metric tons)
metrics$AGB_Mg <- metrics$AGB_kg / 1000
```

---

## 3. `itcLiDAR()` — Individual Tree Crown Segmentation

Performs ITC segmentation directly on point cloud data. Note: this is the original
Dalponte & Coomes (2016) algorithm. lidR's `dalponte2016()` is a C++ reimplementation
that runs MUCH faster. Use `itcLiDAR()` only if you specifically need itcSegment's
integrated workflow.

### Signature

```r
itcLiDAR(
  X = NULL,      # X coordinates
  Y = NULL,      # Y coordinates
  Z = NULL,      # Z coordinates (height-normalized)
  epsg = NULL,   # EPSG code for CRS
  resolution = 0.5,  # CHM resolution
  MinSearchFilSize = 3,   # Min search filter size
  MaxSearchFilSize = 7,   # Max search filter size
  TRESHSeed = 0.45,       # Seed threshold
  TRESHCrown = 0.55,      # Crown threshold
  minDist = 10,            # Min distance between seeds
  maxDist = 40,            # Max crown diameter
  HesightMin = 2           # Min tree height
)
```

### Returns

`SpatialPolygonsDataFrame` with columns:
- `Height_m`: Tree height in meters
- `CA_m2`: Crown area in square meters
- `X`, `Y`: Tree top coordinates

### Usage Note

For most workflows, it is better to use `lidR::segment_trees(dalponte2016(...))` for
segmentation and then pass the results to `itcSegment::dbh()` and `itcSegment::agb()`
for allometric estimation. This gives you the speed of lidR's C++ engine with
itcSegment's allometric equations.

---

## 4. `itcLiDARallo()` — ITC Segmentation with Allometric Output

Combined function that performs segmentation AND allometric estimation in one step.

### Signature

```r
itcLiDARallo(
  X = NULL,
  Y = NULL,
  Z = NULL,
  epsg = NULL,
  resolution = 0.5,
  MinSearchFilSize = 3,
  MaxSearchFilSize = 7,
  TRESHSeed = 0.45,
  TRESHCrown = 0.55,
  minDist = 10,
  maxDist = 40,
  HeightMin = 2,
  cw = 1,        # Crown width allometric parameter
  ch = 1,        # Crown height allometric parameter
  biession = 1   # Biome/species parameter
)
```

---

## 5. Common Integration Pattern: lidR + itcSegment

The recommended workflow uses lidR for processing and segmentation, then itcSegment
for allometric estimation:

```r
library(lidR)
library(itcSegment)

# --- lidR handles point cloud processing ---
las <- readLAS("forest.laz") |>
  classify_ground(csf()) |>
  normalize_height(tin())

# Generate CHM
chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.2))
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)

# Detect trees and segment crowns
ttops <- locate_trees(chm, lmf(ws = function(x) x * 0.07 + 3, hmin = 2))
las <- segment_trees(las, dalponte2016(chm, ttops))

# Extract per-tree metrics
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")

# --- itcSegment handles allometric estimation ---
# CRITICAL: Convert crown AREA to crown DIAMETER before passing to dbh()/agb()
metrics$crown_diameter <- 2 * sqrt(metrics$convhull_area / pi)

# Predict DBH (cm) and AGB (kg)
metrics$DBH_cm <- dbh(H = metrics$Z, CA = metrics$crown_diameter, biome = 0)
metrics$AGB_kg <- agb(H = metrics$Z, CA = metrics$crown_diameter, species = 1)
```

---

## Reference

Jucker, T., Caspersen, J., Chave, J., et al. (2017). Allometric equations for
integrating remote sensing imagery into forest monitoring programs. *Global Change
Biology*, 23(1), 177–190. https://doi.org/10.1111/gcb.13388
