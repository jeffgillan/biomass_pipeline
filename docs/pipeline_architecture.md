# Pipeline Architecture: ALS → Individual Trees → DBH → Biomass Raster

> **Purpose**: This document defines the design decisions, workflow structure, and
> rationale for the biomass estimation pipeline. Claude Code should follow this
> architecture when generating R scripts. Do not improvise alternative approaches
> without explicit user direction.

---

## Project Goal

Process airborne LiDAR (ALS) point cloud data through an individual tree detection
pipeline to produce a rasterized biomass map suitable for web display via MapLibre GL JS.

**Constraints:**
- No field plot data available (prototype phase)
- DBH and biomass estimated entirely via Jucker et al. (2017) global allometric equations
- Output format: Cloud Optimized GeoTIFF (COG) for web serving

---

## Three-Tier Package Architecture

### Tier 1: Point Cloud Processing & Tree Segmentation → `lidR`

lidR handles all point cloud I/O, preprocessing, CHM generation, tree detection,
segmentation, and per-tree metric extraction.

### Tier 2: DBH Estimation → `itcSegment`

itcSegment's `dbh()` function applies the Jucker global allometry with biome-specific
parameters to convert tree height + crown diameter into DBH estimates.

### Tier 3: Biomass Estimation → `itcSegment` (for prototype)

itcSegment's `agb()` function provides direct AGB estimation from height and crown
diameter. For this prototype without field data, this is simpler and introduces
fewer error propagation steps than computing DBH first, then applying separate
biomass equations.

---

## Modular Script Structure

The pipeline should be implemented as separate R scripts, each handling one
processing phase. This allows isolated debugging and reprocessing of individual steps.

```
R/
├── 01_read_and_classify.R    # Read LAS, classify ground, noise removal
├── 02_normalize_and_chm.R    # Height normalization, CHM generation, smoothing
├── 03_detect_and_segment.R   # Tree detection, crown segmentation, metrics
├── 04_estimate_dbh_biomass.R # Allometric DBH and AGB estimation
├── 05_rasterize_output.R     # Rasterize tree-level AGB to grid, write COG
└── run_pipeline.R            # Orchestration script that sources all steps
```

Each script should:
1. Read its input from a defined location (file or R object)
2. Perform its processing step
3. Write intermediate output for the next step
4. Include logging/messaging for progress tracking

---

## Step-by-Step Processing Logic

### Step 1: Read and Classify Ground

```
Input:  Raw .las/.laz file(s)
Output: LAS object with ground points classified (class 2)
```

- Use `readLAS()` with appropriate select/filter parameters
- Run `classify_noise()` with `sor()` or `ivf()` to flag outliers
- Run `classify_ground()` with `csf()` algorithm
- CSF parameters: `sloop_smooth = TRUE` for variable terrain,
  `class_threshold = 0.5`, `cloth_resolution = 0.5`

### Step 2: Normalize Heights and Generate CHM

```
Input:  Ground-classified LAS
Output: Height-normalized LAS + smoothed CHM raster
```

- Run `normalize_height()` with `tin()` interpolation
- Filter out any remaining negative heights: `filter_poi(las, Z >= 0)`
- Generate CHM: `rasterize_canopy()` with `pitfree()` algorithm
  - Resolution: 0.5m
  - `subcircle = 0.2` (critical for ALS data)
  - `thresholds = c(0, 2, 5, 10, 15, 20, 25, 30)`
- Smooth CHM: 3×3 median filter via `terra::focal()`
- Save CHM as intermediate GeoTIFF

### Step 3: Detect Trees and Segment Crowns

```
Input:  Height-normalized LAS + smoothed CHM
Output: Segmented LAS + sf dataframe of tree metrics with crown geometries
```

- Detect treetops: `locate_trees()` with `lmf()`
  - Window size function: `ws = function(x) { x * 0.07 + 3 }`
  - Minimum height: `hmin = 2` (adjust based on forest type)
- Segment trees: `segment_trees()` with `dalponte2016(chm, ttops)`
  - Default parameters are generally good: `th_tree = 2, th_seed = 0.45, th_cr = 0.55`
- Extract metrics: `crown_metrics()` with `.stdtreemetrics` and `geom = "convex"`
- Compute crown diameter: `CD = 2 * sqrt(convhull_area / pi)`

**Algorithm selection rationale:**
- `dalponte2016()` chosen over `silva2016()` because it produces more natural
  crown boundaries and handles heterogeneous canopies better
- `dalponte2016()` chosen over `li2012()` because the point-cloud method is
  orders of magnitude slower and does not improve results at typical ALS densities
- `lmf()` chosen over ForestTools `vwf()` for simplicity; both produce comparable
  results. ForestTools can be used as an alternative for comparison.

### Step 4: Estimate DBH and Biomass

```
Input:  Tree metrics sf dataframe with height and crown diameter
Output: Same dataframe with DBH (cm) and AGB (kg) columns added
```

- Predict DBH: `itcSegment::dbh(H, CA, biome)`
  - H = tree height from `metrics$Z` (meters)
  - CA = crown diameter (meters) — **NOT crown area**
  - biome = select based on study area (see itcSegment reference)
- Predict AGB: `itcSegment::agb(H, CA, species)`
  - species = 1 (gymnosperm) or 2 (angiosperm)
  - Returns AGB in kilograms

**Handling unknown species:**
Without species information, two options:
1. Use `species = 1` for conifer-dominated areas, `species = 2` for broadleaf
2. Compute both and average (rough approximation)
3. Use `dbh()` + a generic allometric equation instead of `agb()` directly

### Step 5: Rasterize Biomass Output

```
Input:  Tree metrics with AGB values + crown polygon geometries
Output: Cloud Optimized GeoTIFF of biomass density (Mg/ha)
```

- Choose output resolution: 10–20m pixels recommended
- Rasterize using `terra::rasterize()`:
  - Rasterize crown polygons, summing AGB per pixel
  - Convert to density: multiply by `(10000 / pixel_area)` to get Mg/ha
- Optional: apply focal smoothing to fill gaps from undetected understory
- Write as COG: `terra::writeRaster()` with `filetype = "COG"`

**Unit conversion chain:**
```
itcSegment::agb() returns kg per tree
Sum per pixel → kg per pixel
÷ 1000 → Mg per pixel
× (10000 / pixel_area_m2) → Mg/ha
```

For a 20m pixel: `AGB_Mg_ha = (sum_kg / 1000) * (10000 / 400)`

---

## Output Specifications for Web Display

The final biomass raster should be:

- **Format**: Cloud Optimized GeoTIFF (COG)
- **CRS**: EPSG:4326 (WGS 84) for web mapping, or original projected CRS if using
  server-side reprojection
- **Resolution**: 10–20m
- **Data type**: Float32
- **Compression**: DEFLATE
- **NoData value**: Explicitly set (e.g., -9999 or NaN)
- **Color ramp metadata**: Optional, but helps for quick visualization

```r
terra::writeRaster(
  biomass_raster,
  "biomass.tif",
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=DEFLATE", "TILED=YES", "OVERVIEW_RESAMPLING=AVERAGE"),
  filetype = "COG",
  NAflag = -9999
)
```

---

## LAScatalog Processing (for large datasets)

If the study area spans multiple LAS tiles, wrap the pipeline in LAScatalog
processing rather than loading everything into memory:

```r
ctg <- readLAScatalog("path/to/tiles/")
opt_chunk_buffer(ctg) <- 30  # 30m buffer prevents edge artifacts in segmentation
opt_output_files(ctg) <- "output/{*}_processed"

# Each lidR function works transparently on the catalog
```

The 30m buffer is critical: tree crowns that span tile boundaries will be
incorrectly segmented without sufficient buffer overlap.

---

## Quality Checks

After each step, include basic validation:

1. **After ground classification**: Check that ~20-40% of points are classified as ground (varies by canopy density)
2. **After normalization**: Verify no negative Z values; max Z should be reasonable (< 60m for most forests)
3. **After CHM**: Visual check that CHM looks reasonable; no large data voids
4. **After tree detection**: Count of trees should be plausible for the area (e.g., 100-500 trees/ha for typical forests)
5. **After DBH estimation**: Distribution should be roughly log-normal; typical range 10-80cm for mature forests
6. **After biomass**: Total AGB should be in plausible range for the biome (e.g., 50-400 Mg/ha for temperate forests)
