# =============================================================================
# A05_edsm_kfold_comparison_27.R
# -----------------------------------------------------------------------------
# Grouped 5-fold ELPD comparison and stacking / pseudo-BMA+ weights over the
# 27-model EDSM set (reported fits 0-21 + appendix fits 22-26). Fit 27 (the
# water-year diagnostic) is deliberately NOT included.
#
# Reads only small criteria objects: the reported-set bundle written by
# 02_edsm_catch_models/edsm_catch_kfold_models.Rmd (`harvest-bundle` chunk) plus
# the slim sidecars written by A04. No fit is loaded.
#
# Outputs (fits_edsm_Vfix/):
#   kfold5_station_criteria_bundle_27.rds         criteria for fits 0-26 (used by A06, A07, A10)
#   fit_edsm0-26_kfold5_station_weights.csv       ELPD, SE, p_kfold, stacking, pBMA+ per model
#   fit_edsm0-26_catch_kfold5_station_comparisons.csv   loo_compare table
# Cost: seconds.
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")
suppressMessages(library(loo))

b <- readRDS(file.path(FIT_EDSM, "kfold5_station_criteria_bundle.rds"))   # fits 0-21
for (n in 22:26)
  b[[sprintf("fit_edsm_%d", n)]] <- readRDS(file.path(FIT_EDSM, sprintf("kfold_fit_edsm_%d.rds", n)))
b <- b[paste0("fit_edsm_", 0:26)]
stopifnot(!any(vapply(b, is.null, logical(1))))

# ---- ELPD ranking ---------------------------------------------------------------------
cmp <- loo_compare(b)
print(round(cmp[, c("elpd_diff", "se_diff", "elpd_kfold", "se_elpd_kfold", "p_kfold")], 2))
comp <- data.frame(fit = rownames(cmp), elpd_diff = cmp[, "elpd_diff"], se_diff = cmp[, "se_diff"],
                   elpd_kfold = cmp[, "elpd_kfold"], se_elpd_kfold = cmp[, "se_elpd_kfold"])
write.csv(comp, file.path(FIT_EDSM, "fit_edsm0-26_catch_kfold5_station_comparisons.csv"),
          row.names = FALSE)

# ---- weights from the pointwise K-fold ELPD -------------------------------------------
e  <- sapply(b, function(x) x$pointwise[, "elpd_kfold"])
sw <- loo::stacking_weights(e)
pw <- loo::pseudobma_weights(e, BB = TRUE)

sal <- paste0("fit_edsm_", c(22, 24, 25, 26))
sec <- paste0("fit_edsm_", c(10, 11, 13, 16:21, 24:26))
message(sprintf("salinity-bearing total stacking weight: %.4f (pBMA+ %.4f)",
                sum(sw[colnames(e) %in% sal]), sum(pw[colnames(e) %in% sal])))
message(sprintf("Secchi-bearing total stacking weight:   %.4f (pBMA+ %.4f)",
                sum(sw[colnames(e) %in% sec]), sum(pw[colnames(e) %in% sec])))

saveRDS(b, file.path(FIT_EDSM, "kfold5_station_criteria_bundle_27.rds"))
w <- data.frame(model      = colnames(e),
                elpd_kfold = sapply(b, function(x) x$estimates["elpd_kfold", "Estimate"]),
                se         = sapply(b, function(x) x$estimates["elpd_kfold", "SE"]),
                p_kfold    = sapply(b, function(x) x$estimates["p_kfold", "Estimate"]),
                stacking   = as.numeric(sw),
                pbma_bb    = as.numeric(pw))
write.csv(w[order(-w$elpd_kfold), ], file.path(FIT_EDSM, "fit_edsm0-26_kfold5_station_weights.csv"),
          row.names = FALSE)
