##############################################################################
## 05_rasterize_output.R
## Rasterize tree-level AGB to 10m grid, convert to Mg/ha, write COG
##############################################################################

library(terra)
library(sf)

cat("=== Step 5: Rasterize Biomass Output ===\n")

# --- Paths ---
input_rds       <- file.path("output", "intermediate", "tree_biomass.rds")
output_agb      <- file.path("output", "final", "biomass_agb_mgha.tif")
output_uncert   <- file.path("output", "final", "biomass_uncertainty_mgha.tif")

if (!file.exists(input_rds)) stop("Input not found: ", input_rds)

# --- Read tree biomass ---
metrics <- readRDS(input_rds)
cat(sprintf("  Trees loaded: %d\n", nrow(metrics)))

# --- Create 10m raster template from tree extent ---
res_m <- 10
pixel_area_m2 <- res_m^2  # 100 m^2

# Need per-tree variance for uncertainty propagation
# Variance of each tree's AGB estimate (from MC)
metrics$AGB_mc_var <- metrics$AGB_mc_sd^2

# Convert sf polygons to centroids, then to SpatVector for terra::rasterize
# Using centroids ensures every tree lands in exactly one pixel (no missed trees
# from polygon-center rule when crowns are smaller than pixel size)
tree_centroids <- sf::st_centroid(metrics)
tree_vect <- terra::vect(tree_centroids)

# Build template raster
template <- terra::rast(tree_vect, resolution = res_m)

# --- Rasterize AGB (sum per pixel) ---
cat("Rasterizing AGB (sum per pixel)...\n")
agb_sum <- terra::rasterize(tree_vect, template, field = "AGB_kg", fun = "sum")

# Convert to Mg/ha: (sum_kg / 1000) * (10000 / pixel_area_m2)
agb_mgha <- (agb_sum / 1000) * (10000 / pixel_area_m2)

cat(sprintf("  AGB range: %.1f - %.1f Mg/ha\n",
            minmax(agb_mgha)[1], minmax(agb_mgha)[2]))

# --- Rasterize uncertainty (sum-of-variances, then sqrt, same unit conversion) ---
cat("Rasterizing uncertainty...\n")
var_sum <- terra::rasterize(tree_vect, template, field = "AGB_mc_var", fun = "sum")

# SD in kg, then convert to Mg/ha
sd_kg <- sqrt(var_sum)
uncert_mgha <- (sd_kg / 1000) * (10000 / pixel_area_m2)

cat(sprintf("  Uncertainty range: %.1f - %.1f Mg/ha\n",
            minmax(uncert_mgha)[1], minmax(uncert_mgha)[2]))

# --- Reproject to EPSG:4326 for web map display ---
if (Sys.getenv("PROJ_DATA") == "" && dir.exists("/opt/conda/share/proj")) {
  Sys.setenv(PROJ_DATA = "/opt/conda/share/proj")
}
cat("Reprojecting rasters to EPSG:4326 (WGS84)...\n")
agb_mgha    <- terra::project(agb_mgha, "EPSG:4326", method = "bilinear")
uncert_mgha <- terra::project(uncert_mgha, "EPSG:4326", method = "bilinear")

# --- Write as COG ---
# Check GDAL version for COG support
gdal_version <- terra::gdal(lib = "gdal")
gdal_major <- as.numeric(strsplit(gdal_version, "\\.")[[1]][1])
use_cog <- gdal_major >= 3

write_cog <- function(raster, path) {
  if (use_cog) {
    cat(sprintf("Writing COG: %s\n", path))
    terra::writeRaster(raster, path,
      overwrite = TRUE,
      datatype  = "FLT4S",
      gdal      = c("COMPRESS=DEFLATE", "OVERVIEW_RESAMPLING=AVERAGE"),
      filetype  = "COG",
      NAflag    = -9999
    )
  } else {
    cat(sprintf("GDAL < 3.1 — writing tiled GeoTIFF: %s\n", path))
    terra::writeRaster(raster, path,
      overwrite = TRUE,
      datatype  = "FLT4S",
      gdal      = c("COMPRESS=DEFLATE", "TILED=YES"),
      NAflag    = -9999
    )
  }
}

write_cog(agb_mgha, output_agb)
write_cog(uncert_mgha, output_uncert)

# --- Write tree crown polygons as GeoPackage for QGIS ---
output_gpkg <- file.path("output", "final", "tree_crowns_biomass.gpkg")
cat(sprintf("Writing crown polygons to %s...\n", output_gpkg))
metrics_4326 <- sf::st_transform(metrics, 4326)
sf::st_write(metrics_4326, output_gpkg, delete_dsn = TRUE, quiet = TRUE)

# --- Quality check ---
agb_vals <- terra::values(agb_mgha, na.rm = TRUE)
if (length(agb_vals) > 0) {
  cat(sprintf("\n--- Output Summary ---\n"))
  cat(sprintf("  Pixels with data: %d\n", length(agb_vals)))
  cat(sprintf("  AGB median: %.1f Mg/ha\n", median(agb_vals)))
  cat(sprintf("  AGB mean:   %.1f Mg/ha\n", mean(agb_vals)))
  cat(sprintf("  AGB range:  %.1f - %.1f Mg/ha\n", min(agb_vals), max(agb_vals)))
}

cat("Step 5 complete.\n\n")
