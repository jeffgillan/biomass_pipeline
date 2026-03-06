##############################################################################
## 03_detect_and_segment.R
## Tree detection (lmf), crown segmentation (dalponte2016), crown metrics
##############################################################################

library(lidR)
library(terra)
library(sf)

cat("=== Step 3: Detect Trees and Segment Crowns ===\n")

# --- Paths ---
input_las  <- file.path("output", "intermediate", "las_normalized.laz")
input_chm  <- file.path("output", "intermediate", "chm_smooth.tif")
output_rds <- file.path("output", "intermediate", "tree_metrics.rds")

if (!file.exists(input_las)) stop("Input not found: ", input_las)
if (!file.exists(input_chm)) stop("Input not found: ", input_chm)

# --- Read ---
cat("Reading normalized LAS and smoothed CHM...\n")
las <- readLAS(input_las)
chm <- terra::rast(input_chm)

# Ensure CHM CRS matches LAS CRS (compound CRS can be lost on disk round-trip)
terra::crs(chm) <- sf::st_crs(las)$wkt

# --- Detect treetops ---
cat("Detecting treetops (lmf, variable window, hmin=2)...\n")
ttops <- locate_trees(chm, lmf(
  ws   = function(x) { x * 0.07 + 3 },
  hmin = 2
))
cat(sprintf("  Treetops detected: %d\n", nrow(ttops)))

# --- Segment crowns ---
cat("Segmenting trees (dalponte2016)...\n")
las <- segment_trees(las, dalponte2016(chm, ttops))

# --- Extract crown metrics with convex hull geometries ---
cat("Extracting crown metrics...\n")
metrics <- crown_metrics(las, .stdtreemetrics, geom = "convex")
cat(sprintf("  Trees with metrics: %d\n", nrow(metrics)))

# --- Compute crown diameter from convex hull area ---
metrics$crown_diameter <- 2 * sqrt(as.numeric(metrics$convhull_area) / pi)

# --- Filter: keep trees with Z >= 2m AND convhull_area >= 1 m^2 ---
n_before <- nrow(metrics)
metrics <- metrics[metrics$Z >= 2 & as.numeric(metrics$convhull_area) >= 1, ]
cat(sprintf("  Trees after filtering (Z>=2, area>=1): %d (removed %d)\n",
            nrow(metrics), n_before - nrow(metrics)))

# --- Quality checks ---
extent_area_ha <- as.numeric(st_area(st_as_sfc(st_bbox(metrics)))) / 10000
tree_density <- nrow(metrics) / extent_area_ha
cat(sprintf("  Approximate tree density: %.0f trees/ha\n", tree_density))
cat(sprintf("  Crown diameter range: %.1f - %.1f m\n",
            min(metrics$crown_diameter), max(metrics$crown_diameter)))
cat(sprintf("  Height range: %.1f - %.1f m\n", min(metrics$Z), max(metrics$Z)))

if (tree_density < 5 || tree_density > 500) {
  warning(sprintf("Tree density %.0f trees/ha outside typical range", tree_density))
}

# --- Save ---
cat(sprintf("Writing tree metrics to %s...\n", output_rds))
saveRDS(metrics, output_rds)

cat("Step 3 complete.\n\n")
