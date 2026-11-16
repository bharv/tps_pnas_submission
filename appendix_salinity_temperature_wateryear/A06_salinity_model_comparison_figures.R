# =============================================================================
# A06_salinity_model_comparison_figures.R
# -----------------------------------------------------------------------------
# Appendix versions of the two catch model-comparison figures (ELPD forest +
# stacking weight + covariate table) with the salinity / temperature models
# included: EDSM 27-model set (fits 0-26) and TPS 17-model set (fits 0-16).
# Water-year fits 17 (TPS) and 27 (EDSM) are diagnostics and are not shown.
#
# Every function (weights, covariate parsing, row assembly, figure builder, saver)
# is sourced from 05_figures/model_comparison_forest.Rmd; only the model lists and
# the covariate vocabulary are extended here.
#
# Inputs : fits_edsm_Vfix/kfold5_station_criteria_bundle_27.rds and
#          fit_edsm0-26_catch_kfold5_station_comparisons.csv (A05);
#          fits_tps_wq_Vfix/fit_tps_wq0..16.rds (their $criteria$loo) and
#          fit_tps_wq0-16_catch_looic_comparisons.csv (A03)
# Outputs: output_figs_appendix/sal_model_comparison_edsm_kfold_lines.{pdf,tiff,png}
#          output_figs_appendix/sal_model_comparison_tps_wq_looic_lines.{pdf,tiff,png}
# Cost   : ~1 min (reads the 17 TPS fits one at a time for their LOO objects).
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

# definitions only -- the data-loading chunks (read-comparisons-data, assemble-data)
# and the main-text figure chunks are NOT run
source_rmd_chunks(NB_MCF, c("setup", "weight-helpers", "covariate-defs", "model-cores",
                            "read-comparisons", "assemble", "figure-builder"))
OUT_DIR <- OUT_APPX                                   # redirect save_figure() output

# ---- extend the covariate vocabulary ------------------------------------------------
COV_LABELS <- c(COV_LABELS, "salinity_z" = "Salinity", "temperature_z" = "Temperature")
COV_ORDER  <- unique(c(COV_ORDER, "Salinity", "Temperature"))

# ---- EDSM: 27 models ------------------------------------------------------------------
cores_edsm_sal <- c(cores_edsm,
  fit_edsm_22 = "out30_z * (distance_m_z + distance_along_channel_z) + salinity_z",
  fit_edsm_23 = "out30_z * (distance_m_z + distance_along_channel_z) + temperature_z",
  fit_edsm_24 = "out30_z + distance_along_channel_z + secchi_z + salinity_z",
  fit_edsm_25 = "out30_z * distance_along_channel_z + secchi_z + salinity_z",
  fit_edsm_26 = "out30_z * (distance_m_z + distance_along_channel_z) + secchi_z + salinity_z")

edsm_crit_sal <- readRDS(file.path(FIT_EDSM, "kfold5_station_criteria_bundle_27.rds"))[names(cores_edsm_sal)]
w_edsm_sal    <- stacking_weights_from_criteria(edsm_crit_sal)
comp_edsm_sal <- read_comparison(file.path(FIT_EDSM, "fit_edsm0-26_catch_kfold5_station_comparisons.csv"))

# ---- TPS water quality: 17 models -----------------------------------------------------
cores_tps_wq_sal <- c(cores_tps_wq,
  fit_tps_wq14 = "out30_z * distance_m_z + distance_along_channel_z + salinity_z",
  fit_tps_wq15 = "out30_z * distance_m_z + distance_along_channel_z + Edge_Type + salinity_z",
  fit_tps_wq16 = "out30_z * distance_m_z + distance_along_channel_z + Edge_Type + secchi_z + salinity_z")

wq_crit_sal <- setNames(lapply(0:16, function(i)
  readRDS(file.path(FIT_TPS, sprintf("fit_tps_wq%d.rds", i)))$criteria$loo),
  paste0("fit_tps_wq", 0:16))
w_wq_sal    <- stacking_weights_from_criteria(wq_crit_sal)
comp_wq_sal <- read_comparison(file.path(FIT_TPS, "fit_tps_wq0-16_catch_looic_comparisons.csv"))

# ---- assemble & draw ------------------------------------------------------------------
g_edsm_sal <- assemble_group(comp_edsm_sal, w_edsm_sal, cores_edsm_sal)
g_wq_sal   <- assemble_group(comp_wq_sal,   w_wq_sal,   cores_tps_wq_sal)

fig_edsm_sal <- make_comparison_figure(
  g_edsm_sal, title = "EDSM catch model comparison (salinity models included)",
  criterion_lab = "grouped K-fold CV on (1 | StationCode)",
  cell_aspect = 1.40, v_border = 0.06, style = "lines")
save_figure(fig_edsm_sal, "sal_model_comparison_edsm_kfold_lines.png",
            n_rows = nrow(g_edsm_sal$df), n_cov = length(g_edsm_sal$cov_cols))

fig_wq_sal <- make_comparison_figure(
  g_wq_sal, title = "TPS catch model comparison (salinity models included)",
  criterion_lab = "PSIS-LOO with moment matching",
  cell_aspect = 1.51, v_border = 0.06, style = "lines")
save_figure(fig_wq_sal, "sal_model_comparison_tps_wq_looic_lines.png",
            n_rows = nrow(g_wq_sal$df), n_cov = length(g_wq_sal$cov_cols))
