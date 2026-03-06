##############################################################################
## 02_normalize_and_chm.R
## Height normalization (TIN), CHM generation (pitfree), 3x3 median smooth
##############################################################################

library(lidR)
library(terra)

cat("=== Step 2: Normalize Heights and Generate CHM ===\n")

# --- Paths ---
input_laz     <- file.path("output", "intermediate", "las_classified.laz")
output_las    <- file.path("output", "intermediate", "las_normalized.laz")
output_chm    <- file.path("output", "intermediate", "chm_smooth.tif")

if (!file.exists(input_laz)) stop("Input not found: ", input_laz)

# --- Read ---
cat("Reading classified LAS...\n")
las <- readLAS(input_laz)

# --- Normalize height using TIN interpolation ---
cat("Normalizing heights (TIN)...\n")
las <- normalize_height(las, tin())

# Remove residual negative Z values
n_neg <- sum(las$Z < 0)
if (n_neg > 0) {
  cat(sprintf("  Removing %d points with Z < 0\n", n_neg))
  las <- filter_poi(las, Z >= 0)
}

max_z <- max(las$Z)
cat(sprintf("  Max normalized height: %.1f m\n", max_z))

if (max_z > 30) {
  warning(sprintf("Max height %.1f m is high for semi-arid woodland; verify data", max_z))
}

# --- Write normalized LAS ---
cat(sprintf("Writing normalized LAS to %s...\n", output_las))
writeLAS(las, output_las)

# --- Generate CHM (pitfree algorithm) ---
cat("Generating CHM (pitfree, 0.5m resolution, subcircle=0.2)...\n")
chm <- rasterize_canopy(las, res = 0.5, pitfree(
  subcircle  = 0.2,
  thresholds = c(0, 2, 5, 10, 15, 20, 25, 30)
))

# --- 3x3 median smooth ---
cat("Applying 3x3 median filter...\n")
chm <- terra::focal(chm, matrix(1, 3, 3), fun = median, na.rm = TRUE)

# --- Write CHM ---
cat(sprintf("Writing smoothed CHM to %s...\n", output_chm))
terra::writeRaster(chm, output_chm, overwrite = TRUE)

cat(sprintf("  CHM dimensions: %d x %d\n", nrow(chm), ncol(chm)))
cat(sprintf("  CHM value range: %.1f - %.1f m\n", minmax(chm)[1], minmax(chm)[2]))

cat("Step 2 complete.\n\n")
