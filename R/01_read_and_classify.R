##############################################################################
## 01_read_and_classify.R
## Read LAS, reproject to UTM, classify noise and ground points
##############################################################################

library(lidR)
library(sf)

cat("=== Step 1: Read and Classify ===\n")

# --- Paths ---
input_laz  <- input_las_file  # set by run_pipeline.R from CLI argument
output_laz <- file.path("output", "intermediate", "las_classified.laz")

if (!file.exists(input_laz)) stop("Input LAZ not found: ", input_laz)

# --- Read ---
cat("Reading LAZ file...\n")
las <- readLAS(input_laz, select = "xyzicr")
cat(sprintf("  Points read: %d\n", npoints(las)))
cat(sprintf("  CRS: %s\n", st_crs(las)$input))

# --- Verify CRS is projected (UTM) ---
current_crs <- st_crs(las)
if (is.na(current_crs)) {
  stop("LAS file has no CRS defined")
}
# Check if geographic (needs reprojection) vs projected (already good)
is_geographic <- isTRUE(current_crs$IsGeographic)
if (is_geographic) {
  cat("Reprojecting to EPSG:32612 (UTM 12N)...\n")
  las <- st_transform(las, st_crs(32612))
  cat(sprintf("  New CRS: %s\n", st_crs(las)$input))
} else {
  cat("  CRS is already projected (UTM). No reprojection needed.\n")
}

# --- Classify noise (Statistical Outlier Removal) ---
cat("Classifying noise (SOR k=10, m=3)...\n")
las <- classify_noise(las, sor(k = 10, m = 3))
n_noise <- sum(las$Classification == 18)
cat(sprintf("  Noise points flagged: %d (%.1f%%)\n", n_noise, 100 * n_noise / npoints(las)))

# Remove noise points
las <- filter_poi(las, Classification != 18)
cat(sprintf("  Points after noise removal: %d\n", npoints(las)))

# --- Classify ground (CSF) ---
cat("Classifying ground (CSF)...\n")
las <- classify_ground(las, csf(
  sloop_smooth    = TRUE,
  class_threshold = 0.5,
  cloth_resolution = 0.5
))

n_ground <- sum(las$Classification == 2)
pct_ground <- 100 * n_ground / npoints(las)
cat(sprintf("  Ground points: %d (%.1f%%)\n", n_ground, pct_ground))

# Quality check
if (pct_ground < 10 || pct_ground > 60) {
  warning(sprintf("Ground point percentage (%.1f%%) outside expected range 10-60%%", pct_ground))
}

# --- Write output ---
cat(sprintf("Writing classified LAS to %s...\n", output_laz))
writeLAS(las, output_laz)

cat("Step 1 complete.\n\n")
