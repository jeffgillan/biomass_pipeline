# ForestTools v1.0.x API Reference for Claude Code

> **Purpose**: ForestTools provides alternative tree detection (`vwf`) and crown
> segmentation (`mcws`) algorithms to lidR. It operates on CHM rasters and is fully
> compatible with `terra`/`sf` objects. Use ForestTools when you want to compare
> results against lidR's built-in algorithms, or when marker-controlled watershed
> segmentation outperforms region-growing approaches for your forest type.

---

## 1. `vwf()` — Variable Window Filter (Tree Detection)

Detects treetops using a local maximum filter with variable window size based on
the Popescu & Wynne method. This is an alternative to `lidR::locate_trees(lmf())`.

### Signature

```r
vwf(
  CHM,                  # SpatRaster canopy height model
  winFun,               # Function: height -> window RADIUS (not diameter)
  minHeight = NULL,     # Minimum tree height to detect
  warnings = TRUE,      # Show warnings
  minWinNeib = "queen", # Neighborhood type for smallest windows: "queen" or "rook"
  IDfield = "treeID",   # Name of ID field in output
  resolution_round = 5  # Decimal rounding for resolution
)
```

### Parameters

| Parameter    | Type       | Description                                          |
|-------------|------------|------------------------------------------------------|
| `CHM`        | SpatRaster | Canopy height model (terra raster). Also accepts legacy RasterLayer |
| `winFun`     | Function   | Maps pixel height to search window **radius** (not diameter!) |
| `minHeight`  | Numeric    | Trees shorter than this are excluded. Set >= 2m      |
| `minWinNeib` | Character  | `"queen"` (8-connected) or `"rook"` (4-connected) for 3x3 windows |
| `IDfield`    | Character  | Name of the tree ID column in output. Default `"treeID"` |

### CRITICAL: winFun defines RADIUS, not diameter

In `vwf()`, the window function returns the search **radius** in map units.
In `lidR::lmf()`, `ws` defines the full window **size** (diameter).

```r
# ForestTools vwf: function returns RADIUS
winFun <- function(x) { x * 0.06 + 0.5 }  # 30m tree -> 2.3m radius -> 4.6m window

# lidR lmf: ws defines full window SIZE (diameter)
ws_fun <- function(x) { x * 0.12 + 1.0 }  # 30m tree -> 4.6m window
```

### Returns

`sf` POINT object with columns:
- `treeID`: Integer tree identifier
- `height`: Tree height value from CHM
- geometry: Point location

### Example

```r
library(ForestTools)
library(terra)

chm <- rast("chm.tif")

# Detect treetops
ttops <- vwf(chm, winFun = function(x) { x * 0.06 + 0.5 }, minHeight = 2)
```

---

## 2. `mcws()` — Marker-Controlled Watershed Segmentation

Segments a CHM into individual crown outlines using watershed segmentation guided
by treetop locations. This is an alternative to `lidR::segment_trees(dalponte2016())`.

### Signature

```r
mcws(
  treetops,          # sf POINT object from vwf() or locate_trees()
  CHM,               # SpatRaster canopy height model
  minHeight = 0,     # Minimum CHM height for crown membership
  format = "raster", # Output format: "raster" or "polygons"
  OSGeoPath = NULL,  # DEPRECATED, no longer needed
  IDfield = "treeID" # Name of tree ID field in treetops
)
```

### Parameters

| Parameter  | Type       | Description                                         |
|-----------|------------|-----------------------------------------------------|
| `treetops` | sf POINT   | Treetop locations from `vwf()` or `locate_trees()`  |
| `CHM`      | SpatRaster | Canopy height model                                 |
| `minHeight`| Numeric    | CHM pixels below this height are excluded from crowns. Set LOWER than `vwf(minHeight)` |
| `format`   | Character  | `"raster"` (SpatRaster with crown IDs) or `"polygons"` (sf POLYGON). Raster is faster; polygons remove orphan segments |
| `IDfield`  | Character  | Name of tree ID field. Must match the field in `treetops` |

### Returns

- If `format = "raster"`: `SpatRaster` where cell values are tree IDs
- If `format = "polygons"`: `sf` POLYGON object with tree ID column

### Example

```r
# Detect treetops
ttops <- vwf(chm, winFun = function(x) { x * 0.06 + 0.5 }, minHeight = 2)

# Segment crowns as raster (faster)
crown_raster <- mcws(ttops, chm, minHeight = 1, format = "raster")

# Segment crowns as polygons (for area calculations)
crown_polys <- mcws(ttops, chm, minHeight = 1, format = "polygons")
```

---

## 3. `glcm()` — Gray-Level Co-occurrence Matrix Texture Metrics

Computes textural metrics per crown segment. Useful for species classification.

### Signature

```r
glcm(
  image,           # SpatRaster (e.g., orthophoto band, CHM, or intensity raster)
  segs,            # SpatRaster of crown segments from mcws()
  n_grey = 32      # Number of grey levels for quantization
)
```

### Returns

`data.frame` with texture metrics per segment (mean, variance, homogeneity,
contrast, dissimilarity, entropy, second_moment, correlation).

---

## 4. `sp_summarise()` — Spatial Summarization

Computes summary statistics for trees within areas of interest or on a grid.
Useful for generating rasterized stand-level metrics.

### Signature

```r
sp_summarise(
  trees,              # sf POINT (treetops) or sf POLYGON (crowns)
  areas = NULL,       # Optional sf POLYGON areas of interest
  variables = NULL,   # Column names to summarize (e.g., "height", "crownArea")
  grid = NULL,        # Numeric grid size for rasterized output
  statFuns = list(min = min, max = max, mean = mean, median = median, sd = sd)
)
```

### Parameters

| Parameter  | Type       | Description                                         |
|-----------|------------|-----------------------------------------------------|
| `trees`    | sf         | Tree locations (POINT) or crown outlines (POLYGON)  |
| `areas`    | sf POLYGON | Areas of interest to summarize within               |
| `variables`| Character  | Column names to compute statistics for              |
| `grid`     | Numeric    | Grid cell size for spatial grid output              |
| `statFuns` | List       | Named list of summary functions                     |

### Returns

- If `areas` provided: `sf` with statistics per area
- If `grid` provided: `SpatRaster` stack with TreeCount + stats per grid cell

### Example (rasterized tree count and height stats)

```r
# 20m grid of tree density and max height
grid_stats <- sp_summarise(ttops, grid = 20, variables = "height")
plot(grid_stats$TreeCount)
plot(grid_stats$heightMax)
```

---

## 5. Interoperability with lidR

ForestTools functions accept both `terra` and legacy `raster`/`sp` objects.
When working with lidR v4.x output:

```r
library(lidR)
library(ForestTools)

# lidR generates the CHM
chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.2))
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)

# ForestTools detects trees and segments crowns
ttops <- vwf(chm, winFun = function(x) { x * 0.06 + 0.5 }, minHeight = 2)
crowns <- mcws(ttops, chm, minHeight = 1, format = "polygons")

# Back to lidR/terra for rasterization
crown_areas <- sf::st_area(crowns)
```

**Key compatibility notes:**
- `vwf()` output is `sf` POINT — directly usable as input to `lidR::segment_trees(dalponte2016(chm, ttops))`
- `mcws()` with `format = "raster"` returns `SpatRaster` — compatible with `terra` operations
- The `height` column from `vwf()` is named differently than `lidR::locate_trees()` which uses `Z`. Rename if cross-using:
  ```r
  # If passing ForestTools ttops to lidR functions:
  names(ttops)[names(ttops) == "height"] <- "Z"
  ```
