# From ALS point clouds to tree-level DBH: the R package landscape

**The `lidR` package dominates individual tree segmentation from airborne LiDAR, but DBH estimation requires a separate allometric step that no single package handles end-to-end.** The critical bottleneck in the ALS-to-biomass pipeline is not segmentation — where mature, well-validated algorithms exist — but rather the indirect estimation of stem diameter from canopy-level metrics. The `itcSegment` package is the only actively maintained CRAN package that provides built-in allometric DBH prediction from ALS-derived height and crown diameter, using the Jucker et al. (2017) global equations. A practical workflow combining `lidR` for point cloud processing and tree segmentation, `itcSegment::dbh()` for diameter estimation, and `BIOMASS` or `allodb` for biomass conversion currently represents the most robust R-based pipeline. Several historically important packages (`rLiDAR`, `TreeLS`) are now archived from CRAN, while newer tools like `lasR` and `cloud2trees` signal where the ecosystem is heading.

## lidR is the gravitational center of ALS processing in R

The `lidR` package (v4.2.3, January 2026) remains the most comprehensive and actively maintained R package for airborne LiDAR analysis. Maintained by Jean-Romain Roussel — now through r-lidar.com, a commercial consulting entity that keeps the software free and open-source — the package has accumulated **over 1,000 citations** and forms the foundation of virtually every R-based ALS workflow. Since v4.x, it has fully transitioned to the modern `sf`/`terra` spatial framework.

For individual tree detection, `lidR` provides `locate_trees()` with a **local maximum filter** (`lmf()`) that operates on either raster CHMs or raw point clouds. The variable window size capability is essential for heterogeneous forests — users supply a function relating tree height to expected crown width (e.g., `ws = function(x) { x * 0.07 + 3 }`), allowing the search window to adapt to canopy structure. For segmentation, four algorithms ship with the core package. **`dalponte2016()`** implements seed-based region growing on the CHM and is the most widely validated choice for conifer and mixed forests. **`silva2016()`** uses Voronoi tessellation from treetop seeds — fast and effective in uniform plantations but producing geometric rather than natural crown boundaries. **`li2012()`** is the only point-cloud-based algorithm, bypassing CHM generation entirely, though it suffers from computational complexity worse than O(n²) and can struggle with large datasets. **`watershed()`** applies morphological watershed segmentation without requiring pre-detected treetops but tends to over-segment in heterogeneous canopies.

The `crown_metrics()` function (which replaced the deprecated `tree_metrics()` and `delineate_crowns()`) extracts per-tree attributes with full geometric flexibility — outputting point, convex hull, concave hull, or bounding box geometries. Built-in metric sets include `.stdtreemetrics` (height, point count, convex hull area), `stdmetrics_z()` (height percentiles, L-moments, entropy), `stdmetrics_i()` (intensity statistics), and `stdshapemetrics()` (eigenvalue-based 3D shape features). Custom metric functions can compute anything derivable from the per-tree point cloud. Critically, **lidR does not estimate DBH** — this is a data limitation, not a software one. It outputs the tree height and crown area values that feed downstream allometric models.

For production-scale work, the `LAScatalog` processing engine handles tiled datasets with configurable chunk sizes, spatial indexing via `.lax` files, and parallelism through the `future` package. The `lidRplugins` extension (GitHub only, last updated February 2023) adds experimental algorithms including `ptree()` and `hamraz2016()` for deciduous forest segmentation, though it may require maintenance due to the `rgeos` deprecation. A newer package, **`lasR`**, from the same development team, targets terabyte-scale processing with **10–100× speed improvements** over lidR through single-pass C++ pipelining — worth monitoring for operational deployments.

## The supporting cast: ForestTools, itcSegment, and the archived packages

**ForestTools** (v1.0.3, February 2025) is the strongest complement to lidR among actively maintained packages. Its `vwf()` function implements the Popescu & Wynne variable window filter for treetop detection, and `mcws()` provides marker-controlled watershed segmentation — a different algorithm family than lidR's options. It also computes GLCM texture metrics per crown segment, which can feed species classification models. ForestTools operates on CHM rasters and is fully compatible with `terra`/`sf` objects produced by lidR. It does not estimate DBH.

**`itcSegment`** (v1.0, August 2023) occupies a unique niche because it is the **only CRAN package providing built-in allometric DBH and AGB estimation from ALS-derived metrics**. Its `dbh()` function implements the Jucker et al. (2017) global allometric equations with **25 biome-specific parameter sets**, predicting stem diameter from tree height and crown diameter. The companion `agb()` function estimates aboveground biomass with separate coefficients for gymnosperms and angiosperms. The package also provides its own ITC delineation via `itcLiDAR()` (the original Dalponte & Coomes 2016 algorithm), though lidR's C++ reimplementation runs hundreds to thousands of times faster. The allometric functions, however, remain uniquely valuable even when segmentation is performed in lidR. Author Michele Dalponte is a co-author on the Jucker et al. paper, lending credibility to the implementation.

Several historically important packages are now **archived or deprecated**. `rLiDAR` was removed from CRAN in October 2023 due to its dependency on the obsolete `rgeos` package; its core Silva et al. (2016) algorithm lives on natively in lidR. **`TreeLS`** was archived from CRAN in February 2021 but remains installable from GitHub. It provides direct stem-level DBH estimation via Hough Transform detection and RANSAC/IRLS circle fitting at 1.3 m height — powerful capabilities, but designed for TLS data requiring dense stem point representation. A 2020 study demonstrated TreeLS with UAV-LiDAR at **>100 pts/m²** achieving ~11% relative RMSE for DBH, but standard ALS densities of 1–10 pts/m² are entirely insufficient for circle fitting. The `treetop` package (v0.0.5) wraps lidR's `silva2016` in a Shiny GUI and adds no new algorithms.

Two TLS-focused packages merit brief mention for users who may integrate terrestrial and airborne data. **FORTLS** (v1.5.3, actively maintained) provides RANSAC-based DBH estimation and DBSCAN tree detection from TLS/MLS scans. **ITSMe** (GitHub only) offers direct DBH measurement from segmented tree point clouds with explicit support for UAV-LiDAR. Neither is suitable for standard ALS data.

A newer entry, **`crownsegmentr`**, implements the Adaptive Mean Shift 3D (AMS3D) algorithm from Ferraz et al. (2016) — a genuinely different approach to 3D point-cloud-based crown delineation that uses adaptive kernel sizes scaled by crown-to-height ratios. It integrates directly with lidR's preprocessing functions and is worth evaluating for complex multi-layered canopies where CHM-based methods struggle.

## The DBH estimation problem and how to solve it

ALS sensors view forests from above and generate point clouds dominated by upper canopy returns. At typical acquisition densities of **2–8 pts/m²**, virtually no stem-level detail exists at 1.3 m height. DBH estimation from ALS is therefore fundamentally indirect — a modeling problem, not a measurement problem.

The dominant approach uses **allometric equations relating tree height and crown diameter to stem diameter**. The Jucker et al. (2017) universal model, built from a global database of 108,753 trees, takes the form `ln(DBH) = a + b × ln(H) + c × ln(CD)` and has been validated across forest types worldwide. A 2025 study in southeastern Queensland confirmed this model outperformed alternatives for ALS-derived predictions. Height alone predicts DBH with R² of roughly **0.5–0.7** for generalized models; adding crown diameter pushes this to **0.7–0.9** depending on forest type and species homogeneity. Species-specific or region-specific allometries substantially improve accuracy — using local allometries can shift AGB estimates by **11–13%** compared to pantropical equations.

Machine learning approaches represent an alternative when field calibration data is available. Random Forest, XGBoost, and partial least squares regression trained on lidR-derived crown metrics (height percentiles, crown area, intensity statistics, voxel metrics) have achieved R² values of **0.82–0.90** for DBH prediction. SHAP analysis from a 2025 UAV-LiDAR study confirmed that height metrics and voxel-based features are the most influential predictors. However, these models require matched field-measured DBH for training and are inherently site-specific.

A critical workflow decision concerns whether to estimate DBH as an intermediate step or predict biomass directly from height and crown diameter. Each intermediate variable introduces additional error propagation — Popescu (2007) reported **47% RMSE** at the individual tree level for AGB via DBH as an intermediate, while direct H+CD→AGB approaches achieved roughly 33% lower error. For users whose downstream application specifically requires DBH (e.g., for stand structure analysis or diameter distributions), the intermediate step is unavoidable. For pure biomass estimation, direct allometric models from `itcSegment::agb()` or custom equations may be preferable.

## Package interoperability and recommended workflow architecture

The most practical R pipeline chains three tiers of packages, each handling a distinct phase:

**Tier 1 — Point cloud processing and tree segmentation:** `lidR` handles everything from reading LAS/LAZ files through ground classification (`classify_ground()` with CSF or PMF), height normalization, CHM generation (the `pitfree()` algorithm with subcircle adjustment plus median smoothing is recommended), tree detection, segmentation, and per-tree metric extraction. ForestTools provides alternative detection (`vwf()`) and segmentation (`mcws()`) algorithms that may outperform lidR's built-in options in certain forest types — testing both is advisable.

**Tier 2 — DBH estimation:** `itcSegment::dbh()` applies the Jucker global allometry with biome-specific parameters. For US forests, `rFIA` can query the FIA database for regional height-diameter relationships by species group, enabling local allometric calibration. Custom statistical models trained on field plots offer the highest accuracy but require matched ground truth.

**Tier 3 — Biomass conversion:** `BIOMASS` (Réjou-Méchain et al. 2017) provides the Chave et al. (2014) pantropical equation with taxonomy-based wood density lookup and **full Bayesian error propagation** via `AGBmonteCarlo()`. For temperate and boreal forests, `allodb` (rOpenSci, GitHub) weights 570 parsed allometric equations by taxonomic and climatic similarity to the study site. Both require DBH as input.

In code, the core pipeline is:

```r
library(lidR); library(itcSegment); library(BIOMASS)

# Process point cloud
las <- readLAS("file.laz") |>
  classify_ground(csf()) |>
  normalize_height(tin())

# Generate CHM and detect trees
chm <- rasterize_canopy(las, res = 0.5, pitfree(subcircle = 0.2))
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)
ttops <- locate_trees(chm, lmf(ws = function(x) x * 0.07 + 3))

# Segment and extract metrics
las <- segment_trees(las, dalponte2016(chm, ttops))
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")
metrics$CD <- 2 * sqrt(metrics$convhull_area / pi)

# Estimate DBH (Jucker global model, biome = 0)
metrics$DBH <- itcSegment::dbh(H = metrics$Z, CA = metrics$CD, biome = 0)

# Estimate biomass (tropical example)
WD <- getWoodDensity(genus = species_data$genus, species = species_data$species)
metrics$AGB <- computeAGB(D = metrics$DBH, WD = WD$meanWD, H = metrics$Z)
```

## Accuracy expectations and practical limitations to plan for

Individual tree detection rates from ALS vary dramatically with forest structure. **Overstory detection exceeds 90%** in most conditions, but **understory detection drops to roughly 60%** because canopy occlusion blocks returns from lower strata. The Kaartinen et al. (2012) benchmark found detection rates ranging from 25–102% across algorithms and forest types (values above 100% indicate commission errors — single crowns split into multiple segments). Dense, multi-layered broadleaf forests present the greatest challenge; open conifer plantations are easiest. Point density matters: CHM-based methods function at ≥2 pts/m² but improve substantially above **8 pts/m²**, while understory segmentation requires ~170 pts/m² — achievable only with UAV-LiDAR, not standard manned ALS.

Species classification from ALS structural data alone achieves **76–87% overall accuracy** for broad groups (conifer vs. broadleaf), but species-level identification generally requires spectral data fusion. This matters because allometric equation selection depends on species — misidentification propagates directly to DBH and biomass errors.

A **hybrid approach** combining individual tree detection for dominant canopy trees with area-based methods (`lidR::pixel_metrics()`) for total stand-level estimates is increasingly recommended in the literature. This captures biomass from undetected understory trees while retaining individual-level precision where detection succeeds. For the user's allometric biomass modeling goal, this hybrid framing may prove more robust than relying solely on individual tree detection in structurally complex forests.

## Conclusion

The R ecosystem for ALS-to-biomass analysis is mature but fragmented across purpose-built packages. `lidR` is non-negotiable as the processing backbone — no other package approaches its breadth, performance, or maintenance cadence. The key insight for pipeline design is that **DBH estimation is a modeling decision, not a software feature**: users must choose between the Jucker global allometry (via `itcSegment::dbh()`), locally calibrated statistical models, or machine learning approaches trained on field data. The Jucker equations provide a defensible starting point, but local calibration against field plots will always yield better accuracy. For users processing large datasets, the transition path from `lidR` to `lasR` (same development team, 10–100× faster) is worth planning for. The experimental `cloud2trees` package on GitHub attempts to wrap the entire pipeline into a single function call, but its dependence on the archived `TreeLS` package makes it fragile for production use. The most resilient approach remains assembling the pipeline from well-maintained components: `lidR` → `itcSegment` → `BIOMASS`/`allodb`, with ForestTools as an alternative segmentation engine and `rFIA` for allometric calibration data.