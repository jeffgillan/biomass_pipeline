# Known Pitfalls: ALS Biomass Pipeline in R

> **Purpose**: This document catalogs known issues, gotchas, and failure modes
> encountered when building LiDAR-based forest biomass pipelines in R. Claude Code
> should consult this before writing or debugging pipeline code.

---

## 1. lidR API Pitfalls

### 1.1 CHM MUST be smoothed before tree detection

**Problem**: Raw CHMs from `rasterize_canopy()` contain noise spikes and micro-pits
that cause `locate_trees()` to detect hundreds of false treetops.

**Solution**: Always apply at least one pass of median filtering:
```r
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)
```
For noisier data, two passes may be needed. Do NOT use mean filtering — it blurs
treetops and shifts their detected locations.

### 1.2 `pitfree()` requires `subcircle > 0` for ALS data

**Problem**: Without the `subcircle` parameter, the pitfree algorithm produces CHMs
with extensive NA gaps between discrete LiDAR returns, especially in sparser canopy areas.

**Solution**: Set `subcircle = 0.2` (or 0.15–0.25 depending on point density):
```r
chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.2))
```

### 1.3 `dalponte2016()` requires BOTH the CHM and treetops

**Problem**: Unlike some algorithms, `dalponte2016()` is a two-step process that
needs the CHM raster AND the detected treetop points.

**Correct:**
```r
las <- segment_trees(las, dalponte2016(chm, ttops))
```

**Wrong:**
```r
las <- segment_trees(las, dalponte2016(chm))  # Missing treetops!
```

### 1.4 `crown_metrics()` replaced `tree_metrics()` AND `delineate_crowns()`

**Problem**: Old tutorials reference `tree_metrics()` and `delineate_crowns()` which
no longer exist in lidR v4.x.

**Correct v4.x usage:**
```r
# Get metrics as points
metrics <- crown_metrics(las, .stdtreemetrics)

# Get metrics with crown polygon geometries
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")
```

### 1.5 CRS must be projected, not geographic

**Problem**: All spatial operations in lidR assume projected coordinates (meters).
Geographic CRS (lat/lon in degrees) will produce nonsensical results for distances,
areas, and window sizes.

**Check:** `sf::st_crs(las)` should show a projected CRS (UTM, State Plane, etc.)
If data is in EPSG:4326, reproject first using external tools or `sf::st_transform()`.

### 1.6 The `Z` column in crown_metrics output is MAX height, not mean

When using `.stdtreemetrics`, the `Z` column is the **maximum** height of points
in that tree's crown — i.e., the treetop height. This is what you want for
allometric equations. Do not confuse it with mean crown height.

### 1.7 `segment_trees()` modifies the LAS object in place

The function adds a `treeID` attribute to each point. Points not assigned to any
tree will have `treeID = NA`. Filter these out when needed:
```r
# Only points assigned to trees
tree_points <- filter_poi(las, !is.na(treeID))
```

### 1.8 LAScatalog buffer size matters for segmentation

**Problem**: If the chunk buffer is too small, trees near tile edges get truncated
or duplicated.

**Solution**: Set buffer to at least the maximum expected crown diameter:
```r
opt_chunk_buffer(ctg) <- 30  # 30m is safe for most forests
```

### 1.9 `terra` vs `raster` confusion

lidR v4.x defaults to `terra::SpatRaster` output. Old code using `raster::` functions
will fail. Common translations:

```r
# OLD (will fail)
raster::focal(chm, ...)
raster::writeRaster(...)
sp::spplot(...)

# NEW (correct)
terra::focal(chm, ...)
terra::writeRaster(...)
plot(metrics["Z"])  # sf objects plot directly
```

---

## 2. itcSegment Pitfalls

### 2.1 The `CA` parameter in `dbh()` is Crown DIAMETER, not Crown Area

**This is the single most common error when using itcSegment.**

The parameter is named `CA` but it expects **crown diameter in meters**, not crown area.
The source code documentation says "Crown diameter in meters."

**Correct:**
```r
# Derive crown diameter from convex hull area
CD <- 2 * sqrt(convhull_area / pi)
dbh(H = height, CA = CD, biome = 0)
```

**Wrong:**
```r
dbh(H = height, CA = convhull_area, biome = 0)  # WRONG! This passes area, not diameter
```

### 2.2 `dbh()` returns centimeters, `agb()` returns kilograms

Be careful with unit chains:
- `dbh()` → output in **cm**
- `agb()` → output in **kg**
- Most biomass maps report in **Mg/ha** (metric tons per hectare)

Conversion: `Mg/ha = (sum_of_kg_per_pixel / 1000) * (10000 / pixel_area_m2)`

### 2.3 Biome code 0 (Global) may not be optimal

While biome = 0 works everywhere, using a biome-specific code improves accuracy
significantly. For US forests, codes 7–14 are the Nearctic biomes. Choose based on:
- Forest type (boreal, temperate coniferous, temperate mixed, woodland)
- Species type (angiosperm vs gymnosperm)

### 2.4 `agb()` species parameter: 1 = gymnosperm, 2 = angiosperm

This is counterintuitive (most people expect 1 = broadleaf). Double-check:
- Pine, spruce, fir, cedar, hemlock → `species = 1`
- Oak, maple, birch, poplar, eucalyptus → `species = 2`

### 2.5 Both `dbh()` and `agb()` accept vectors

They are vectorized — pass entire columns, not individual values:
```r
metrics$DBH <- dbh(H = metrics$Z, CA = metrics$CD, biome = 10)  # Correct
# No need for apply() or loops
```

---

## 3. ForestTools Pitfalls

### 3.1 `vwf()` winFun returns RADIUS, not diameter

Unlike `lidR::lmf()` where `ws` is the full window size, `vwf()`'s `winFun`
returns the search **radius**. A function returning `x * 0.06 + 0.5` for a 30m
tree gives radius = 2.3m (window diameter = 4.6m).

### 3.2 `vwf()` output column is `height`, not `Z`

`lidR::locate_trees()` returns a column named `Z`.
`ForestTools::vwf()` returns a column named `height`.

If passing `vwf()` output into `lidR::segment_trees(dalponte2016(chm, ttops))`,
lidR may expect a `Z` column. Rename if needed:
```r
names(ttops)[names(ttops) == "height"] <- "Z"
```

### 3.3 `mcws()` minHeight should be LOWER than vwf() minHeight

`mcws(minHeight)` defines the minimum CHM value for a pixel to be assigned to a
crown. This should be lower than the minimum treetop detection height, or you'll
get crowns with no interior pixels.

```r
ttops <- vwf(chm, winFun, minHeight = 2)
crowns <- mcws(ttops, chm, minHeight = 1)  # Lower than 2!
```

---

## 4. General R/Spatial Pitfalls

### 4.1 Memory management for large point clouds

ALS datasets can be multi-GB. Key strategies:
```r
# Read only needed attributes
las <- readLAS("file.laz", select = "xyzicr")

# Use LAScatalog for tiled processing
ctg <- readLAScatalog("tiles/")

# Decimate if needed for testing
las_test <- decimate_points(las, random_per_voxel(res = 1, n = 1))
```

### 4.2 `sf::st_area()` returns units objects

When computing crown areas from sf geometries, the result has units attached:
```r
area <- sf::st_area(crown_polygon)
# Returns: 45.2 [m^2]  (with units)

# Strip units for numeric operations:
area_numeric <- as.numeric(sf::st_area(crown_polygon))
```

### 4.3 COG creation requires GDAL ≥ 3.1

Cloud Optimized GeoTIFF support requires a recent GDAL version. Check with:
```r
terra::gdal()  # Shows GDAL version
```
If `filetype = "COG"` fails, fall back to regular GeoTIFF with tiling:
```r
terra::writeRaster(r, "out.tif", gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
```

### 4.4 Rasterizing polygons: use `terra::rasterize()`, not `fasterize()`

`fasterize` depends on the deprecated `raster`/`sp` stack. Use `terra::rasterize()`:
```r
template <- terra::rast(ext(las), res = 20, crs = terra::crs(las))
biomass_raster <- terra::rasterize(
  terra::vect(crown_polygons),  # Convert sf to SpatVector
  template,
  field = "AGB_kg",
  fun = "sum"
)
```

### 4.5 `terra::vect()` vs `sf` objects

`terra::rasterize()` expects a `SpatVector`, not an `sf` object. Convert:
```r
spat_vec <- terra::vect(sf_object)
```

### 4.6 NaN/NA handling in allometric estimation

Trees with very small crowns (< 1m² area) or very short heights (< 2m) can
produce unrealistic DBH/AGB values. Filter before allometric estimation:
```r
metrics <- metrics[metrics$Z >= 2 & metrics$convhull_area >= 1, ]
```

### 4.7 Docker/containerization: use rocker/geospatial base image

The `rocker/geospatial` image includes pre-compiled GDAL, PROJ, GEOS, and the
R spatial packages (sf, terra, stars). This avoids the most common failure mode
in containerized R geospatial workflows — system library compilation errors.

```dockerfile
FROM rocker/geospatial:4.4
RUN R -e "install.packages(c('lidR', 'itcSegment', 'ForestTools'))"
```

### 4.8 `renv` for reproducibility

Use `renv` to lock package versions:
```r
renv::init()
renv::install(c("lidR", "itcSegment", "ForestTools", "terra", "sf"))
renv::snapshot()
```
Commit `renv.lock` to the repo so Claude Code can restore the exact environment.
