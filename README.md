# Aerial LiDAR Biomass Pipeline

Airborne LiDAR (ALS) to aboveground biomass (AGB) estimation pipeline for the Peloncillo Mountains, New Mexico. Processes a USGS 3DEP point cloud through individual tree detection and allometric estimation to produce rasterized biomass maps as Cloud Optimized GeoTIFFs. DBH and biomass are estimated via Jucker et al. (2017) global allometric equations with Monte Carlo uncertainty propagation.

<br>

## Pipeline Overview

```
.laz point cloud
    |
    v
01_read_and_classify.R    -- Read, noise removal (SOR), ground classification (CSF)
    |
    v
02_normalize_and_chm.R    -- TIN height normalization, pitfree CHM (0.5m), 3x3 median smooth
    |
    v
03_detect_and_segment.R   -- Tree detection (lmf), crown segmentation (dalponte2016), metrics
    |
    v
04_estimate_dbh_biomass.R -- itcSegment DBH/AGB, Jucker MC error propagation, BIOMASS cross-check
    |
    v
05_rasterize_output.R     -- 20m raster grid, Mg/ha conversion, COG output
    |
    v
biomass_agb_mgha.tif + biomass_uncertainty_mgha.tif (COG)
```
<br>
<br>

## Quick Start -- Local (Pathway 1)

```bash
# Create conda environment
conda env create -f environment.yml
conda activate biomass_pipeline

# Install CRAN-only packages
R -e "install.packages(c('itcSegment', 'BIOMASS', 'RCSF'), repos='https://cloud.r-project.org')"

# Place .laz file in project root, then run
Rscript R/run_pipeline.R
```

## Quick Start -- Docker (Pathway 2)

```bash
docker build -t biomass-pipeline .
docker run \
  -v $(pwd):/pipeline/data \
  -v $(pwd)/output:/pipeline/output \
  biomass-pipeline
```
<br>
<br>

## Input Data

- **File:** `USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.laz`
- Download from [Cyverse](https://data.cyverse.org/dav-anon/iplant/home/jgillan/living_carbon_demo/USGS_LPC_NM_SouthCentral_2018_D19_12RXV885860.laz)
- **Source:** USGS 3DEP (South Central NM, 2018, D19 collection)
- **CRS:** NAD83(2011) / UTM 12N (EPSG:6341 compound)
- **Total Points:** ~10.8 million
- **Point Density:** 5-10 pts/SqMeter
- **Area:** 226 hectares; 559 acres; 2.2 sqKM; 0.87 sqMiles
- **Size:** 69 MiB
- Download other lidar data from [USGS 3DEP](https://apps.nationalmap.gov/lidar-explorer/).

<br>
<br>

## Outputs

| File | Description | Units |
|------|-------------|-------|
| `output/final/biomass_agb_mgha.tif` | Aboveground biomass density | Mg/ha (megagrams per hectare) |
| `output/final/biomass_uncertainty_mgha.tif` | AGB uncertainty (1 SD from MC) | Mg/ha |
| `output/intermediate/tree_crowns.geojson` | Polygon of every detected tree | kg per tree |

Both final TIFs are Cloud Optimized GeoTIFFs (Float32, DEFLATE compression, -9999 NoData).

<br>
<br>

## Above Ground Biomass Calculation

This pipeline uses the [Jucker et al. 2017](https://onlinelibrary.wiley.com/doi/10.1111/gcb.13388) allometric model which is a log-linear regression: it predicts the natural log of DBH from the natural log of height and crown diameter. Like any regression, the model was fitted to a dataset — in this case, about 108,000 trees from around the world — and the predictions carry residual error from that fit.

## Error Estimation

Two independent uncertainty methods:

1. **Jucker Monte Carlo (500 iterations):** Propagates uncertainty in allometric coefficients, height measurement (1m RMSE), crown delineation (15% CV), allometric residual (RSE=0.40 on ln scale), and wood density (0.50 +/- 0.08 g/cm3) through Jucker 2017 + Chave 2014 equations. Per-tree outputs: mean, SD, 5th/95th percentiles, CV.

2. **BIOMASS::AGBmonteCarlo():** Independent cross-check using the BIOMASS R package's built-in error propagation with Chave 2004 DBH error model.

Pixel-level uncertainty: sum-of-variances across trees within each 20m pixel, converted to SD in Mg/ha.


---                                                                                                              
  Path 1: itcSegment (point estimate only)                                                                         
                                                                                                                   
  - itcSegment::dbh(H, CA, biome=14) → itcSegment::agb(D, species=1)                                               
  - Produces a single AGB_kg value per tree — no error bounds                                                      
  - Serves as a comparison baseline                                                                                
                                                                                                                   
  Path 2: Jucker 2017 Monte Carlo (primary uncertainty)

  500 MC iterations per tree, propagating 5 independent error sources each iteration:

  ┌───────────────────────────────────┬───────────────────────────────┬───────────────────────┐
  │           Error source            │         How perturbed         │       Magnitude       │
  ├───────────────────────────────────┼───────────────────────────────┼───────────────────────┤
  │ Allometric coefficients (a, b, c) │ N(coeff, SE)                  │ SE = 0.05, 0.03, 0.02 │
  ├───────────────────────────────────┼───────────────────────────────┼───────────────────────┤
  │ ALS height measurement            │ H + N(0, 1.0)                 │ 1.0 m RMSE            │
  ├───────────────────────────────────┼───────────────────────────────┼───────────────────────┤
  │ Crown diameter delineation        │ Lognormal, CV=15%             │ multiplicative        │
  ├───────────────────────────────────┼───────────────────────────────┼───────────────────────┤
  │ Allometric residual               │ + N(0, 0.40) on ln(DBH) scale │ RSE = 0.40            │
  ├───────────────────────────────────┼───────────────────────────────┼───────────────────────┤
  │ Wood density                      │ N(0.50, 0.08), floor 0.2      │ species variation     │
  └───────────────────────────────────┴───────────────────────────────┴───────────────────────┘

  Each iteration: perturbed H,CD → Jucker ln(DBH) equation → Chave 2014 AGB formula → one AGB sample. From 500
  samples: mean, SD, 5th/95th percentiles, CV per tree.

  Path 3: BIOMASS::AGBmonteCarlo (cross-check)

  - Uses itcSegment DBH as input (not Jucker DBH)
  - Propagates height error (1m), wood density error (0.08), and DBH error (Chave 2004 model)
  - 500 iterations → mean + SD per tree
  - Not used in final rasters — purely for validation against Path 2

  Raster-level propagation (Step 05)

  - Only Path 2 (Jucker MC) feeds the uncertainty raster
  - Per-tree variance: AGB_mc_var = AGB_mc_sd²
  - Per-pixel: sum of variances across all trees in pixel (assumes independence)
  - Pixel SD: √(sum_var), then converted to Mg/ha
  - Output: biomass_uncertainty_mgha.tif = 1 SD band

  Key assumption: tree-level errors are independent. Systematic biases (CHM method, allometric model choice, biome
  code) aren't captured.



## License

MIT
