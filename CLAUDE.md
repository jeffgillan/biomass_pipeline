# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ALS (Airborne LiDAR) to biomass estimation pipeline in R. Processes point cloud data through individual tree detection to produce rasterized biomass maps as Cloud Optimized GeoTIFFs for web display via MapLibre GL JS. No field plot data — DBH and biomass are estimated entirely via Jucker et al. (2017) global allometric equations.

## Three-Tier Package Architecture

1. **lidR** — Point cloud I/O, ground classification, height normalization, CHM generation, tree detection, segmentation, crown metrics
2. **itcSegment** — Allometric DBH (`dbh()`) and AGB (`agb()`) estimation from height + crown diameter
3. **itcSegment** (prototype) / **BIOMASS** or **allodb** (production) — Biomass conversion

## Pipeline Structure

Scripts in `R/` are modular and sequential:
```
01_read_and_classify.R    → Read LAS, classify ground (CSF), noise removal
02_normalize_and_chm.R    → Height normalization (TIN), CHM (pitfree), 3x3 median smooth
03_detect_and_segment.R   → locate_trees(lmf), segment_trees(dalponte2016), crown_metrics
04_estimate_dbh_biomass.R → itcSegment::dbh() and agb()
05_rasterize_output.R     → Rasterize AGB to grid, write COG
run_pipeline.R            → Orchestration script sourcing all steps
```

## Critical Gotchas

**itcSegment `CA` parameter is Crown DIAMETER, not Crown Area.** This is the most common error. Always convert:
```r
CD <- 2 * sqrt(convhull_area / pi)
dbh(H = height, CA = CD, biome = 0)
```

**Units:** `dbh()` returns cm, `agb()` returns kg. Biomass maps use Mg/ha: `(sum_kg / 1000) * (10000 / pixel_area_m2)`

**`agb()` species codes:** 1 = gymnosperm (conifers), 2 = angiosperm (broadleaf). This is counterintuitive.

**CHM must be smoothed** before tree detection — use `terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)`. Never use mean filtering.

**`pitfree()` requires `subcircle = 0.2`** for ALS data or the CHM will have extensive NA gaps.

**`dalponte2016()` requires both CHM and treetops** — `segment_trees(las, dalponte2016(chm, ttops))`.

**lidR v4.x deprecations:** Use `segment_trees()` not `lastrees()`, `crown_metrics()` not `tree_metrics()`/`delineate_crowns()`, `rasterize_canopy()` not `grid_canopy()`. Output is `terra::SpatRaster` / `sf`, not `raster`/`sp`.

**CRS must be projected** (UTM, State Plane), not geographic (EPSG:4326). Reproject before processing.

**ForestTools `vwf()` winFun returns RADIUS**, while lidR `lmf()` ws is full window SIZE. Also `vwf()` output column is `height`, not `Z` — rename if passing to lidR functions.

## Key Algorithm Parameters

- CHM resolution: 0.5m with `pitfree(subcircle = 0.2)`
- Tree detection: `lmf(ws = function(x) { x * 0.07 + 3 }, hmin = 2)`
- Segmentation: `dalponte2016(chm, ttops)` with defaults `th_tree = 2, th_seed = 0.45, th_cr = 0.55`
- LAScatalog buffer: 30m minimum for tile-edge segmentation
- Output raster: 10-20m resolution, Float32, DEFLATE compression, COG format

## Biome Codes for US Forests (itcSegment)

Codes 7-14 are Nearctic. For the sample data (New Mexico): code 13 or 14 (Woodlands and savannas). Use code 0 (Global) when unsure.

## Data Files

- `USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.laz` — Sample ALS point cloud (South Central NM)
- `peloncillo_east_chm.tif` — Pre-computed CHM raster
- `USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.xml` — LAS metadata

## Reference Documentation

- `docs/pipeline_architecture.md` — Authoritative pipeline design and step-by-step logic
- `docs/lidR_api_reference.md` — lidR v4.x function signatures (use instead of web searches)
- `docs/itcSegment_api_reference.md` — itcSegment API with biome/species code tables
- `docs/ForestTools_api_reference.md` — ForestTools vwf/mcws API
- `docs/known_pitfalls.md` — Comprehensive catalog of failure modes and fixes

Consult these docs before writing or debugging pipeline code. Do not improvise alternative approaches without explicit user direction.
