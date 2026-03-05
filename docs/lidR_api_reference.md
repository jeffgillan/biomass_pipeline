# lidR v4.2.x API Reference for Claude Code

> **Purpose**: This document provides accurate, current function signatures and usage patterns
> for the lidR R package (v4.x). Use this as the authoritative reference when writing R code
> that processes ALS point cloud data. Do NOT rely on web searches — many online examples
> use deprecated v2/v3 syntax.

## Critical: Deprecated Functions (DO NOT USE)

These functions existed in lidR v2/v3 and appear frequently in online tutorials.
**They will produce errors or warnings in v4.x.**

| Deprecated (v2/v3)          | Current (v4.x)                          |
|-----------------------------|------------------------------------------|
| `lastrees()`                | `segment_trees()`                        |
| `tree_metrics()`            | `crown_metrics()`                        |
| `delineate_crowns()`        | `crown_metrics(..., geom = "convex")`    |
| `find_trees()`              | `locate_trees()`                         |
| `lasground()`               | `classify_ground()`                      |
| `lasnormalize()`            | `normalize_height()`                     |
| `grid_canopy()`             | `rasterize_canopy()`                     |
| `grid_terrain()`            | `rasterize_terrain()`                    |
| `grid_metrics()`            | `pixel_metrics()`                        |
| `lasfilter*()`              | `filter_poi()`                           |
| `lasclip*()`                | `clip_*()`                               |

Also note: lidR v4.x returns `terra::SpatRaster` objects by default (not `raster::RasterLayer`).
The spatial framework is `sf`/`terra`, not `sp`/`raster`.

---

## 1. Reading Point Cloud Data

### `readLAS()`

```r
readLAS(
  files,
  select = "*",      # Which attributes to load: "xyzicrna" etc.
  filter = ""         # LASlib filter string applied at read time
)
```

**Parameters:**
- `files`: Character. Path to .las or .laz file(s)
- `select`: Character. Attributes to load. Default `"*"` loads all. Use `"xyz"` for coordinates only, `"xyzi"` for coords + intensity, `"xyzc"` for coords + classification, `"*"` for everything
- `filter`: Character. LASlib filter applied during reading (memory efficient). Examples:
  - `"-keep_first"` — keep only first returns
  - `"-drop_z_below 0"` — remove negative heights
  - `"-keep_class 2"` — keep only ground points
  - `"-keep_random_fraction 0.5"` — subsample 50%

**Returns:** `LAS` object

**Example:**
```r
las <- readLAS("forest.laz", select = "xyzicr", filter = "-drop_z_below 0")
```

### `readLAScatalog()`

```r
readLAScatalog(
  folder,
  select = "*",
  filter = "",
  chunk_size = 0,     # 0 = use original file tiling
  chunk_buffer = 30   # Buffer in CRS units (meters)
)
```

**Returns:** `LAScatalog` object for processing tiled datasets

---

## 2. Ground Classification

### `classify_ground()`

```r
classify_ground(las, algorithm, last_returns = TRUE)
```

**Parameters:**
- `las`: LAS or LAScatalog object
- `algorithm`: Ground classification algorithm (see below)
- `last_returns`: Logical. Use only last returns for classification. Default `TRUE`

**Algorithm options:**

```r
# Cloth Simulation Filter (recommended for most terrain types)
csf(
  sloop_smooth = FALSE,   # Handle steep slopes
  class_threshold = 0.5,  # Distance threshold (meters)
  cloth_resolution = 0.5, # Cloth grid resolution (meters)
  rigidness = 1,          # 1=flat, 2=gentle slope, 3=steep
  iterations = 500,       # Max iterations
  time_step = 0.65        # Time step
)

# Progressive Morphological Filter
pmf(
  ws = seq(3, 21, 3),     # Window sizes (increasing sequence)
  th = seq(0.1, 2, length.out = length(ws))  # Height thresholds
)
```

**Returns:** LAS object with Classification attribute updated (class 2 = ground)

**Example:**
```r
las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 0.5))
```

---

## 3. Height Normalization

### `normalize_height()`

```r
normalize_height(las, algorithm, use_class = c(2L, 9L), dtm = NULL)
```

**Parameters:**
- `las`: LAS or LAScatalog object (must have classified ground points)
- `algorithm`: Interpolation algorithm for DTM
- `use_class`: Integer vector. Ground point classes to use. Default `c(2L, 9L)`
- `dtm`: Optional pre-computed SpatRaster DTM. If provided, `algorithm` is ignored

**Algorithm options:**

```r
# Triangulated Irregular Network (recommended)
tin()

# Inverse Distance Weighting
knnidw(k = 10, p = 2)

# Kriging
kriging(k = 10)
```

**Returns:** LAS object with Z values replaced by height above ground

**Example:**
```r
las <- normalize_height(las, tin())
# Verify: all Z values should now be >= 0 (approximately)
```

---

## 4. Canopy Height Model Generation

### `rasterize_canopy()`

```r
rasterize_canopy(las, res = 1, algorithm = p2r())
```

**Parameters:**
- `las`: LAS or LAScatalog object (should be height-normalized)
- `res`: Numeric. Output resolution in CRS units (meters). Typical: 0.5 for tree detection
- `algorithm`: CHM algorithm (see below)

**Algorithm options:**

```r
# Point-to-raster (simple, fast)
p2r(
  subcircle = 0,       # Expand points to circles of this radius
  na.fill = NULL        # Optional fill algorithm for NA pixels
)

# Pit-free algorithm (RECOMMENDED for tree detection)
pitfree(
  thresholds = c(0, 2, 5, 10, 15, 20, 25, 30),  # Height thresholds
  max_edge = c(0, 1.5),  # Max Delaunay edge lengths [below lowest thresh, above]
  subcircle = 0.2         # IMPORTANT: Set > 0 for ALS data to fill gaps
)

# Delaunay triangulation
dsmtin(max_edge = 0)
```

**Returns:** `SpatRaster` (terra)

**IMPORTANT**: The CHM typically needs smoothing before tree detection:

```r
chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.2))
# Smooth with 3x3 median filter to reduce noise
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)
```

---

## 5. Individual Tree Detection

### `locate_trees()`

```r
locate_trees(las, algorithm, uniqueness = "incremental")
```

**Parameters:**
- `las`: LAS, LAScatalog, or SpatRaster (CHM)
- `algorithm`: Tree detection algorithm
- `uniqueness`: Character. How to assign tree IDs. `"incremental"`, `"gpstime"`, or `"bitmerge"`

**Algorithm options:**

```r
# Local Maximum Filter (recommended)
lmf(
  ws,             # Window size: numeric OR function of height
  hmin = 2,       # Minimum tree height to detect
  shape = "circular"  # "circular" or "square"
)
```

**CRITICAL**: `ws` can be a function for variable window size (recommended for heterogeneous forests):

```r
# Window size as a function of pixel/point height
# This means a 30m tree uses a ~5m search window, a 10m tree uses ~3.7m
ws_fun <- function(x) { x * 0.06 + 3 }

ttops <- locate_trees(chm, lmf(ws = ws_fun, hmin = 2))
```

If `ws` is a fixed number, it defines a constant window radius in map units.

**Returns:** `sf` POINT object with columns:
- `treeID`: Integer tree identifier
- `Z`: Tree height (from CHM or point cloud)
- geometry: Point location of treetop

**Example:**
```r
ttops <- locate_trees(chm, lmf(ws = function(x) { x * 0.07 + 3 }, hmin = 2))
```

---

## 6. Tree Crown Segmentation

### `segment_trees()`

```r
segment_trees(las, algorithm, attribute = "treeID", uniqueness = "incremental")
```

**Parameters:**
- `las`: LAS or LAScatalog object
- `algorithm`: Segmentation algorithm (see below)
- `attribute`: Character. Name of the new attribute. Default `"treeID"`

**Algorithm options:**

```r
# Dalponte 2016 - CHM-based region growing (RECOMMENDED for most cases)
# Requires both CHM and treetops as input
dalponte2016(
  chm,            # SpatRaster CHM
  treetops,       # sf POINT from locate_trees()
  th_tree = 2,    # Min height for a pixel to be a tree
  th_seed = 0.45, # Seed threshold (fraction of treetop height)
  th_cr = 0.55,   # Crown threshold (fraction of treetop height)
  max_cr = 10     # Max crown radius in pixels
)

# Silva 2016 - Voronoi-based (fast, good for plantations)
silva2016(
  chm,            # SpatRaster CHM
  treetops,       # sf POINT from locate_trees()
  max_cr_factor = 0.6,  # Max crown factor
  exclusion = 0.3       # Exclusion zone as fraction of height
)

# Li 2012 - Point cloud-based (NO CHM needed, but slow for large datasets)
li2012(
  dt1 = 1.5,      # Threshold 1
  dt2 = 2,         # Threshold 2
  R = 2,           # Search radius
  Zu = 15,         # Height threshold
  hmin = 2,        # Min tree height
  speed_up = 10    # Decimation factor for speed
)

# Watershed segmentation
watershed(
  chm,            # SpatRaster CHM
  th_tree = 2,    # Min height
  tol = 1,        # Tolerance
  ext = 1          # Extension
)
```

**Returns:** LAS object with `treeID` attribute added to each point

**Example (recommended workflow):**
```r
las <- segment_trees(las, dalponte2016(chm, ttops))
```

---

## 7. Crown Metrics Extraction

### `crown_metrics()`

This is the primary function for extracting per-tree attributes. It replaces the
deprecated `tree_metrics()` and `delineate_crowns()`.

```r
crown_metrics(
  las,
  func = ~list(z_max = max(Z)),  # Metric function (formula or function)
  geom = "point",       # Output geometry: "point", "convex", "concave", "bbox"
  concaveman = c(3, 0), # Parameters for concave hull (if geom="concave")
  attribute = "treeID"  # Name of the tree ID attribute
)
```

**Parameters:**
- `las`: LAS object with `treeID` attribute (from `segment_trees()`)
- `func`: Formula or function defining metrics to compute. Uses tidy evaluation with `~`
- `geom`: Character. Output geometry type:
  - `"point"` — centroid of each tree's points
  - `"convex"` — convex hull polygon of each crown
  - `"concave"` — concave hull polygon (requires concaveman package)
  - `"bbox"` — bounding box of each crown
- `attribute`: Character. Name of the tree segmentation attribute

**Built-in metric sets:**
```r
# Basic tree metrics (height, point count, convex hull area)
.stdtreemetrics
# Returns: treeID, Z (max height), npoints, convhull_area

# Height metrics (percentiles, moments, etc.)
.stdmetrics_z
# Returns: zmax, zmean, zsd, zskew, zkurt, zentropy, zq5-zq95, etc.

# Intensity metrics
.stdmetrics_i

# 3D shape metrics (eigenvalue-based)
.stdshapemetrics
```

**Custom metrics example:**
```r
custom_func <- ~list(
  height = max(Z),
  mean_z = mean(Z),
  sd_z = sd(Z),
  n_pts = length(Z),
  mean_intensity = mean(Intensity),
  crown_area = st_area(st_convex_hull(st_combine(st_as_sf(data.frame(X=X, Y=Y), coords=c("X","Y")))))
)

metrics <- crown_metrics(las, func = custom_func, geom = "convex")
```

**Simpler and more common usage:**
```r
# Get standard tree metrics with convex hull polygons
metrics <- crown_metrics(las, func = .stdtreemetrics, geom = "convex")
# Result is sf object with columns: treeID, Z, npoints, convhull_area, geometry

# Derive crown diameter from convex hull area
metrics$crown_diameter <- 2 * sqrt(metrics$convhull_area / pi)
```

**Returns:** `sf` object (point or polygon depending on `geom`)

---

## 8. Pixel-Level Metrics (Area-Based Approach)

### `pixel_metrics()`

```r
pixel_metrics(
  las,
  func,            # Metric function
  res = 20,        # Pixel size in CRS units
  start = c(0, 0), # Origin alignment
  filter = NULL,    # Point filter formula
  by_echo = NULL    # Compute by echo type
)
```

**Example (for wall-to-wall structural metrics):**
```r
metrics_raster <- pixel_metrics(las, .stdmetrics_z, res = 20)
# Returns a multi-band SpatRaster with one band per metric
```

---

## 9. Terrain Rasterization

### `rasterize_terrain()`

```r
rasterize_terrain(
  las,
  res = 1,
  algorithm = tin(),
  use_class = c(2L, 9L)
)
```

**Returns:** `SpatRaster` DTM

---

## 10. LAScatalog Processing Engine

For large datasets split across multiple tiles:

```r
ctg <- readLAScatalog("path/to/tiles/")

# Configure processing options
opt_chunk_size(ctg) <- 500       # Chunk size in meters (0 = use file boundaries)
opt_chunk_buffer(ctg) <- 30      # Buffer to avoid edge effects
opt_output_files(ctg) <- "output/{ORIGINALFILENAME}_{ID}"  # Output template
opt_laz_compression(ctg) <- TRUE # Compress output

# Enable parallel processing
library(future)
plan(multisession, workers = 4)

# All lidR functions work transparently on LAScatalog
nlas <- normalize_height(ctg, tin())
```

**Important catalog options:**
```r
opt_select(ctg) <- "xyzicr"     # Attributes to read
opt_filter(ctg) <- "-drop_z_below 0"  # Read filter
opt_chunk_size(ctg) <- 0        # 0 = process by original tile boundaries
opt_chunk_buffer(ctg) <- 30     # Buffer size (critical for tree segmentation)
opt_progress(ctg) <- TRUE       # Show progress bar
```

---

## 11. Utility Functions

```r
# Check LAS file for issues
las_check(las)

# Filter points
filter_poi(las, Classification != 7)  # Remove noise points
filter_poi(las, Z >= 0 & Z <= 50)     # Height range filter

# Clip to area of interest
clip_rectangle(las, xmin, ymin, xmax, ymax)
clip_circle(las, xcenter, ycenter, radius)
clip_roi(las, sf_polygon)

# Decimate points
decimate_points(las, random_per_voxel(res = 1, n = 1))
decimate_points(las, highest(res = 1))

# Classify noise
classify_noise(las, sor(k = 10, m = 3))  # Statistical Outlier Removal
classify_noise(las, ivf(res = 5, n = 6)) # Isolated Voxel Filter
```

---

## 12. Writing Output

```r
# Write LAS/LAZ
writeLAS(las, "output.laz")

# Write raster (CHM, DTM, metrics) — use terra
terra::writeRaster(chm, "chm.tif", overwrite = TRUE)

# Write as Cloud Optimized GeoTIFF (COG)
terra::writeRaster(
  biomass_raster,
  "biomass_cog.tif",
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "TILED=YES"),
  filetype = "COG"
)
```

---

## Key Dependencies

```r
# Core dependencies (installed automatically with lidR)
library(lidR)    # v4.2.3
library(terra)   # Raster operations
library(sf)      # Vector/spatial operations
library(data.table)  # Internal data handling

# The following are NOT needed (legacy):
# library(raster)  # Replaced by terra
# library(sp)      # Replaced by sf
# library(rgdal)   # Removed from CRAN
# library(rgeos)   # Removed from CRAN
```
