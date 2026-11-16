# =============================================================================
# A00_appendix_setup.R  --  shared settings for every appendix script
# -----------------------------------------------------------------------------
# Sourced at the top of A02-A10. Run every appendix script from the REPOSITORY
# ROOT (the folder holding the .Rproj), e.g.
#     Rscript appendix_salinity_temperature_wateryear/A03_tps_salinity_wateryear_fits.R
#
# Appendix scope (manuscript appendix only; the main-text figures are salinity-free):
#   * salinity    -- TPS fits 14-16, EDSM fits 22, 24, 25, 26
#   * temperature -- EDSM fit 23 (TPS temperature tracks the seasonal spline,
#                    r = 0.75-0.95 with water_year_day in every year, so it is not fitted)
#   * water year  -- TPS fit 17 (= fit 16 + factor(water_year)),
#                    EDSM fit 27 (= fit 24 + factor(water_year)); diagnostics only, NOT
#                    members of either reported model set and given no stacking weight.
#   * supporting  -- CDEC outflow / salinity / turbidity time series (Fig. S1), salinity
#                    contour maps with catch (Fig. S2), coefficient and ELPD tables.
# =============================================================================

if (!length(list.files(".", pattern = "\\.Rproj$")))
  stop("Run the appendix scripts from the repository root (the folder with the .Rproj file).")

source("R/source_rmd_chunks.R")

# Main-text notebooks whose functions the appendix re-uses (single source of truth).
NB_TPS_WQ  <- "01_tps_catch_models/tps_catch_wq_models.Rmd"
NB_EDSM    <- "02_edsm_catch_models/edsm_catch_kfold_models.Rmd"
NB_MCF     <- "05_figures/model_comparison_forest.Rmd"
NB_FIGS    <- "05_figures/main_text_ensemble_figures.Rmd"

# Fit folders (shared with the main text: appendix fits sit next to the reported ones,
# numbered after them, so the main-text code never picks them up by accident).
FIT_TPS    <- "fits_tps_wq_Vfix"
FIT_EDSM   <- "fits_edsm_Vfix"

# Data
EDSM_SALFIX <- "data/tow_12_6_edsm_wq_salfix.csv"   # corrected EDSM salinity (A02)
WEB         <- "data/web_downloads"                 # 00_data/00_download_web_data.R
DERIVED     <- "data/derived_appendix"              # intermediate tables built here

# Appendix figures (all prefixed sal_)
OUT_APPX    <- "output_figs_appendix"

for (d in c(DERIVED, OUT_APPX)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# PDF device that embeds fonts on macOS (quartz) and falls back to cairo / pdf elsewhere.
open_pdf_device <- function(path, w, h) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = path, type = "pdf", width = w, height = h, bg = "white")
  } else if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename = path, width = w, height = h, bg = "white")
  } else {
    grDevices::pdf(file = path, width = w, height = h, bg = "white")
  }
}
