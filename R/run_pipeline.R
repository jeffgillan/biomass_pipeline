##############################################################################
## run_pipeline.R
## Orchestration script: runs the full ALS-to-biomass pipeline
##############################################################################

cat("============================================================\n")
cat("  ALS Biomass Pipeline\n")
cat("============================================================\n\n")

# Set PROJ_DATA if not already set (conda environments)
if (Sys.getenv("PROJ_DATA") == "" && dir.exists("/opt/conda/share/proj")) {
  Sys.setenv(PROJ_DATA = "/opt/conda/share/proj")
}

t_start <- Sys.time()

# --- Create output directories ---
dir.create("output/intermediate", recursive = TRUE, showWarnings = FALSE)
dir.create("output/final", recursive = TRUE, showWarnings = FALSE)

# --- Check required packages ---
cat("Checking required packages...\n")
required_pkgs <- c("lidR", "terra", "sf", "itcSegment")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c('", paste(missing, collapse = "','"), "'))")
}

optional_pkgs <- c("BIOMASS")
missing_opt <- optional_pkgs[!sapply(optional_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_opt) > 0) {
  cat(sprintf("  Optional packages not available: %s (will skip cross-check)\n",
              paste(missing_opt, collapse = ", ")))
}

gdal_ver <- terra::gdal(lib = "gdal")
cat(sprintf("  GDAL version: %s\n", gdal_ver))
cat("  All required packages found.\n\n")

# --- Run pipeline steps ---
run_step <- function(script, step_name) {
  cat(sprintf("--- %s ---\n", step_name))
  t0 <- Sys.time()
  source(script, local = FALSE)
  elapsed <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("  [%s completed in %.1f seconds]\n\n", step_name, as.numeric(elapsed)))
}

run_step("R/01_read_and_classify.R",    "Step 1: Read & Classify")
run_step("R/02_normalize_and_chm.R",    "Step 2: Normalize & CHM")
run_step("R/03_detect_and_segment.R",   "Step 3: Detect & Segment")
run_step("R/04_estimate_dbh_biomass.R", "Step 4: DBH & Biomass")
run_step("R/05_rasterize_output.R",     "Step 5: Rasterize Output")

# --- Summary ---
total_elapsed <- difftime(Sys.time(), t_start, units = "secs")
cat("============================================================\n")
cat(sprintf("  Pipeline complete in %.1f seconds\n", as.numeric(total_elapsed)))
cat("============================================================\n")

# List outputs
cat("\nOutputs:\n")
final_files <- list.files("output/final", full.names = TRUE)
for (f in final_files) {
  cat(sprintf("  %s (%.1f KB)\n", f, file.size(f) / 1024))
}
