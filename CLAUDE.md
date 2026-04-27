# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ALS (Airborne LiDAR) to biomass estimation pipeline in R. Processes point cloud data through individual tree detection to produce rasterized biomass maps as Cloud Optimized GeoTIFFs for web display via MapLibre GL JS. No field plot data — DBH and biomass are estimated entirely via Jucker et al. (2017) global allometric equations.

## Running the Pipeline

```bash
# Native (Conda)
conda env create -f environment.yml
conda activate biomass_pipeline
R -e "install.packages(c('itcSegment', 'BIOMASS', 'RCSF'), repos='https://cloud.r-project.org')"
Rscript R/run_pipeline.R /path/to/input.laz

# Docker
docker build -t biomass-pipeline .
docker run -v $(pwd):/pipeline/data -v $(pwd)/output:/pipeline/output biomass-pipeline
```

The pipeline creates `output/intermediate/` (per-step scratch files) and `output/final/` (deliverables).

## Pipeline Structure

Scripts in `R/` are modular and sequential:
```
01_read_and_classify.R    → Read LAS, classify ground (CSF), noise removal (SOR)
02_normalize_and_chm.R    → Height normalization (TIN), CHM (pitfree 0.5m), 3x3 median smooth
03_detect_and_segment.R   → locate_trees(lmf), segment_trees(dalponte2016), crown_metrics
04_estimate_dbh_biomass.R → itcSegment point estimates + 500-iteration Monte Carlo (Chave 2014)
05_rasterize_output.R     → Rasterize AGB to 10m grid (Mg/ha), write COG
run_pipeline.R            → Orchestration: CLI arg parsing, package checks, step timing
```

## Three-Tier Package Architecture

1. **lidR** — Point cloud I/O, ground classification, height normalization, CHM generation, tree detection, segmentation, crown metrics
2. **itcSegment** — Allometric DBH (`dbh()`) and AGB (`agb()`) point estimates from height + crown diameter
3. **BIOMASS** / **allodb** — Production biomass conversion (Monte Carlo cross-check via `BIOMASS::AGBmonteCarlo()`)

## Monte Carlo Uncertainty (Step 04)

Step 04 does more than point estimates. It runs 500 iterations perturbing five independent error sources, then applies the **Chave 2014 pantropical allometric equation** per iteration:

```r
AGB_kg = 0.0673 * (WD * DBH_cm^2 * H_m)^0.976
```

Error sources perturbed per iteration:
- Allometric coefficients: `a=0.557±0.05`, `b=0.809±0.03`, `c=0.056±0.02` (Jucker RSE=0.40 on ln scale)
- ALS height: `N(0, 1.0m)` RMSE
- Crown diameter: lognormal with CV=15%
- Wood density: `N(0.50, 0.08)` g/cm³

Per-tree outputs: `AGB_mc_mean`, `AGB_mc_sd`, `AGB_mc_q05/q95`, `AGB_mc_cv`.

## Output Files

```
output/final/
  biomass_agb_mgha.tif           — AGB density (Mg/ha, Float32, COG, EPSG:4326)
  biomass_uncertainty_mgha.tif   — 1 SD uncertainty (Mg/ha, COG)
  biomass_uncertainty_cv_pct.tif — Coefficient of Variation (%, COG)
  tree_crowns_biomass.gpkg       — Individual tree crown polygons with attributes
```

COG specs: 10m resolution (100 m² pixels), Float32, DEFLATE compression, -9999 NoData, EPSG:4326.

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

**CRS must be projected** (UTM, State Plane), not geographic (EPSG:4326). Reproject before processing. Input sample data is EPSG:6341 (compound CRS) — step 01 auto-reprojects to EPSG:32612 (UTM 12N). Final rasters are reprojected to EPSG:4326 in step 05 for web display.

**ForestTools `vwf()` winFun returns RADIUS**, while lidR `lmf()` ws is full window SIZE. Also `vwf()` output column is `height`, not `Z` — rename if passing to lidR functions.

**`terra::rasterize()` expects SpatVector, not sf** — convert with `terra::vect()` before rasterizing.

## Key Algorithm Parameters

- CHM resolution: 0.5m with `pitfree(subcircle = 0.2)`
- Tree detection: `lmf(ws = function(x) { x * 0.07 + 3 }, hmin = 2)`
- Segmentation: `dalponte2016(chm, ttops)` with defaults `th_tree = 2, th_seed = 0.45, th_cr = 0.55`
- Tree filtering: Z ≥ 2m AND convhull_area ≥ 1 m² (removes spurious detections)
- LAScatalog buffer: 30m minimum for tile-edge segmentation
- Output raster: 10m resolution (100 m² pixels), Float32, DEFLATE compression, COG format

## Biome Codes for US Forests (itcSegment)

Codes 7-14 are Nearctic. For the sample data (New Mexico): code 13 or 14 (Woodlands and savannas). Use code 0 (Global) when unsure.

## Validation

There is no formal test suite. Each step prints inline quality checks:
- Step 01: ground point % (expects 10–60%)
- Step 02: max height sanity (warns if >30m for semi-arid woodland)
- Step 03: tree density (expects 5–500 trees/ha), crown diameter and height ranges
- Step 04: DBH/AGB range checks across all methods
- Step 05: pixel counts, AGB/CV statistics, output file sizes

## Data Files

- `USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.laz` — Sample ALS point cloud (South Central NM, ~10.8M points, 226 ha, EPSG:6341)
- `peloncillo_east_chm.tif` — Pre-computed CHM raster
- `USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.xml` — LAS metadata

## Reference Documentation

- `docs/pipeline_architecture.md` — Authoritative pipeline design and step-by-step logic
- `docs/lidR_api_reference.md` — lidR v4.x function signatures (use instead of web searches)
- `docs/itcSegment_api_reference.md` — itcSegment API with biome/species code tables
- `docs/ForestTools_api_reference.md` — ForestTools vwf/mcws API
- `docs/known_pitfalls.md` — Comprehensive catalog of failure modes and fixes

Consult these docs before writing or debugging pipeline code. Do not improvise alternative approaches without explicit user direction.
