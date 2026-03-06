##############################################################################
## 04_estimate_dbh_biomass.R
## Point estimates (itcSegment), Monte Carlo error propagation (Jucker 2017),
## BIOMASS::AGBmonteCarlo() cross-check
##############################################################################

library(itcSegment)
library(sf)

cat("=== Step 4: Estimate DBH and Biomass with Uncertainty ===\n")

# --- Paths ---
input_rds  <- file.path("output", "intermediate", "tree_metrics.rds")
output_rds <- file.path("output", "intermediate", "tree_biomass.rds")

if (!file.exists(input_rds)) stop("Input not found: ", input_rds)

# --- Read tree metrics ---
metrics <- readRDS(input_rds)
n_trees <- nrow(metrics)
cat(sprintf("  Trees loaded: %d\n", n_trees))

# ========================================================================
# 4A. Point estimates via itcSegment
# ========================================================================
cat("Computing point estimates (itcSegment)...\n")

# CRITICAL: CA = crown DIAMETER, not crown area
metrics$DBH_cm <- itcSegment::dbh(
  H     = metrics$Z,
  CA    = metrics$crown_diameter,
  biome = 14  # Nearctic Woodlands/savannas - Gymnosperm
)

metrics$AGB_kg <- itcSegment::agb(
  H       = metrics$Z,
  CA      = metrics$crown_diameter,
  species = 1  # gymnosperm
)

cat(sprintf("  DBH range: %.1f - %.1f cm\n", min(metrics$DBH_cm), max(metrics$DBH_cm)))
cat(sprintf("  AGB range: %.1f - %.1f kg\n", min(metrics$AGB_kg), max(metrics$AGB_kg)))

# ========================================================================
# 4B. Monte Carlo error propagation on Jucker allometric model
# ========================================================================
cat("Running Monte Carlo error propagation (500 iterations)...\n")

set.seed(42)
n_mc <- 500

# Jucker et al. (2017) global allometric coefficients for ln(DBH) ~ ln(H) + ln(CD)
# From Table S4: Global model
jucker_a   <- 0.557
jucker_b   <- 0.809
jucker_c   <- 0.056
jucker_rse <- 0.40  # residual standard error on ln scale

# Coefficient uncertainties (approximate SE from Jucker 2017)
jucker_a_se <- 0.05
jucker_b_se <- 0.03
jucker_c_se <- 0.02

# Wood density for gymnosperm semi-arid
wd_mean <- 0.50
wd_sd   <- 0.08

# Store MC results: n_trees x n_mc matrix
agb_mc_matrix <- matrix(NA_real_, nrow = n_trees, ncol = n_mc)

for (i in seq_len(n_mc)) {
  # Sample allometric coefficients
  a_i <- rnorm(1, jucker_a, jucker_a_se)
  b_i <- rnorm(1, jucker_b, jucker_b_se)
  c_i <- rnorm(1, jucker_c, jucker_c_se)

  # Perturb height (additive, SD = 1.0 m ALS RMSE)
  H_i <- metrics$Z + rnorm(n_trees, 0, 1.0)
  H_i <- pmax(H_i, 0.5)  # floor to avoid log of zero/negative

  # Perturb crown diameter (lognormal, CV = 15%)
  cd_cv <- 0.15
  cd_sigma <- sqrt(log(1 + cd_cv^2))
  cd_mu    <- log(metrics$crown_diameter) - 0.5 * cd_sigma^2
  CD_i <- exp(rnorm(n_trees, cd_mu, cd_sigma))

  # Jucker model: ln(DBH) = a + b*ln(H) + c*ln(CD) + epsilon
  ln_dbh_i <- a_i + b_i * log(H_i) + c_i * log(CD_i) + rnorm(n_trees, 0, jucker_rse)
  DBH_cm_i <- exp(ln_dbh_i)

  # Sample wood density
  WD_i <- rnorm(n_trees, wd_mean, wd_sd)
  WD_i <- pmax(WD_i, 0.2)  # floor

  # Chave 2014 pantropical: AGB_kg = 0.0673 * (WD * D_cm^2 * H_m)^0.976
  # D is in cm, H in m, AGB in kg
  agb_mc_matrix[, i] <- 0.0673 * (WD_i * DBH_cm_i^2 * H_i)^0.976
}

# Summarize MC results per tree
metrics$AGB_mc_mean <- rowMeans(agb_mc_matrix)
metrics$AGB_mc_sd   <- apply(agb_mc_matrix, 1, sd)
metrics$AGB_mc_q05  <- apply(agb_mc_matrix, 1, quantile, probs = 0.05)
metrics$AGB_mc_q95  <- apply(agb_mc_matrix, 1, quantile, probs = 0.95)
metrics$AGB_mc_cv   <- metrics$AGB_mc_sd / metrics$AGB_mc_mean

cat(sprintf("  MC AGB mean range: %.1f - %.1f kg\n",
            min(metrics$AGB_mc_mean), max(metrics$AGB_mc_mean)))
cat(sprintf("  MC AGB CV range: %.2f - %.2f\n",
            min(metrics$AGB_mc_cv), max(metrics$AGB_mc_cv)))

# ========================================================================
# 4C. BIOMASS::AGBmonteCarlo() as independent check
# ========================================================================
cat("Running BIOMASS::AGBmonteCarlo() cross-check...\n")

has_biomass <- requireNamespace("BIOMASS", quietly = TRUE)

if (has_biomass) {
  MC_result <- BIOMASS::AGBmonteCarlo(
    D       = metrics$DBH_cm,
    WD      = rep(wd_mean, n_trees),
    errWD   = rep(wd_sd, n_trees),
    H       = metrics$Z,
    errH    = rep(1.0, n_trees),
    Dpropag = "chave2004"
  )
  metrics$AGB_BIOMASS_mean <- rowMeans(MC_result$AGB_simu)
  metrics$AGB_BIOMASS_sd   <- apply(MC_result$AGB_simu, 1, sd)

  cat(sprintf("  BIOMASS AGB mean range: %.1f - %.1f kg\n",
              min(metrics$AGB_BIOMASS_mean), max(metrics$AGB_BIOMASS_mean)))
} else {
  warning("BIOMASS package not available; skipping cross-check")
  metrics$AGB_BIOMASS_mean <- NA_real_
  metrics$AGB_BIOMASS_sd   <- NA_real_
}

# ========================================================================
# Quality checks
# ========================================================================
cat("\n--- Quality Summary ---\n")
cat(sprintf("  Trees: %d\n", n_trees))
cat(sprintf("  DBH: median=%.1f cm, range=%.1f-%.1f cm\n",
            median(metrics$DBH_cm), min(metrics$DBH_cm), max(metrics$DBH_cm)))
cat(sprintf("  itcSegment AGB: median=%.1f kg, range=%.1f-%.1f kg\n",
            median(metrics$AGB_kg), min(metrics$AGB_kg), max(metrics$AGB_kg)))
cat(sprintf("  Jucker MC AGB: median=%.1f kg, range=%.1f-%.1f kg\n",
            median(metrics$AGB_mc_mean), min(metrics$AGB_mc_mean), max(metrics$AGB_mc_mean)))
if (has_biomass) {
  cat(sprintf("  BIOMASS AGB: median=%.1f kg, range=%.1f-%.1f kg\n",
              median(metrics$AGB_BIOMASS_mean), min(metrics$AGB_BIOMASS_mean), max(metrics$AGB_BIOMASS_mean)))
}

# --- Save ---
cat(sprintf("\nWriting tree biomass to %s...\n", output_rds))
saveRDS(metrics, output_rds)

cat("Step 4 complete.\n\n")
