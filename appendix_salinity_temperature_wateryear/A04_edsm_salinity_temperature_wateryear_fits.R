# =============================================================================
# A04_edsm_salinity_temperature_wateryear_fits.R
# -----------------------------------------------------------------------------
# EDSM appendix fits, on the same 4,241 tows, row order, priors, sampler settings,
# (1 | StationCode) intercept and the SAME cached 5-fold station partition as the
# reported fits 0-21:
#   fit 22 = fit-9 backbone + salinity        (screen, the way Secchi/DO/habitat were screened)
#   fit 23 = fit-9 backbone + temperature     (screen; temperature did not pass)
#   fit 24 = out30 + dist along + Secchi + salinity             } salinity crossed back onto
#   fit 25 = out30 * dist along + Secchi + salinity             } the Secchi-bearing backbones,
#   fit 26 = out30 * (dist off + dist along) + Secchi + salinity} as fits 17-21 did for Secchi
#   fit 27 = fit 24 + factor(water_year)      (water-year DIAGNOSTIC only, never stacked)
# fit-9 backbone = out30_z * (distance_m_z + distance_along_channel_z)
#
# All six use data/tow_12_6_edsm_wq_salfix.csv (A02): corrected salinity_z and a new
# temperature_z; every other column is identical to the reported frame (asserted).
#
# The runner, priors (b ~ normal(0, 2.5), shape ~ gamma(2, 0.1), Intercept ~
# normal(-5.5, 2), sds ~ normal(20, 5)), K-fold helper and fold vector are sourced
# from 02_edsm_catch_models/edsm_catch_kfold_models.Rmd.
#
# Outputs: fits_edsm_Vfix/fit_edsm_22..27.rds and slim K-fold sidecars
#          kfold_fit_edsm_22..27.rds (so A05 never reloads the ~100 MB fits);
#          summaries appended to fits_edsm_Vfix/edsm_fits_Vfix_output.txt
# Cost   : ~1-2 h per model (fit + 5 K-fold refits). Existing fits are skipped.
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

source_rmd_chunks(NB_EDSM, c("libraries", "load-data", "priors", "kfold-helper",
                             "station-folds", "fit-defs-runner"))
suppressMessages({library(readr); library(loo)})

# ---- data + hard alignment guards ------------------------------------------------
edsm_salfix <- read_csv(EDSM_SALFIX, show_col_types = FALSE)
stopifnot(nrow(edsm_salfix) == nrow(tow_12_6_edsm_wq_nona),
          identical(edsm_salfix$StationCode, tow_12_6_edsm_wq_nona$StationCode),
          identical(edsm_salfix$n_l_75,      tow_12_6_edsm_wq_nona$n_l_75),
          identical(edsm_salfix$out30_z,     tow_12_6_edsm_wq_nona$out30_z),
          identical(edsm_salfix$distance_m_z, tow_12_6_edsm_wq_nona$distance_m_z),
          identical(edsm_salfix$distance_along_channel_z,
                    tow_12_6_edsm_wq_nona$distance_along_channel_z),
          length(station_folds) == nrow(edsm_salfix),
          all(tapply(station_folds, edsm_salfix$StationCode,
                     function(x) length(unique(x))) == 1))

# ---- 1. salinity / temperature fits 22-26 -----------------------------------------
fit_defs_appendix <- tibble::tribble(
  ~num, ~fixef,
  22L,  "out30_z * (distance_m_z + distance_along_channel_z) + salinity_z",
  23L,  "out30_z * (distance_m_z + distance_along_channel_z) + temperature_z",
  24L,  "out30_z + distance_along_channel_z + secchi_z + salinity_z",
  25L,  "out30_z * distance_along_channel_z + secchi_z + salinity_z",
  26L,  "out30_z * (distance_m_z + distance_along_channel_z) + secchi_z + salinity_z"
)

setup_kfold_plan(chains = 4)
for (i in seq_len(nrow(fit_defs_appendix))) {
  n <- fit_defs_appendix$num[i]
  f <- file.path(FIT_EDSM, sprintf("fit_edsm_%d.rds", n))
  if (file.exists(f)) { message("SKIP fit ", n, " - already exists"); next }
  fit_station_model(num = n, fixef = fit_defs_appendix$fixef[i],
                    folds = station_folds, data = edsm_salfix)
  fit <- readRDS(f)
  saveRDS(fit$criteria[["kfold5_StationCode"]],
          file.path(FIT_EDSM, sprintf("kfold_fit_edsm_%d.rds", n)))
  rm(fit); gc(verbose = FALSE)
}

# ---- 2. water-year diagnostic fit 27 = fit 24 + factor(water_year) -----------------
# Reference level is the first water year (2017).
f27 <- file.path(FIT_EDSM, "fit_edsm_27.rds")
if (!file.exists(f27)) {
  edsm <- edsm_salfix %>% mutate(wy_f = factor(water_year))
  form27 <- n_l_75 ~ out30_z + distance_along_channel_z + secchi_z + salinity_z + wy_f +
    s(water_year_day) + offset(log(Volume)) + (1 | StationCode)

  fit <- brm(formula = form27, family = negbinomial(link = "log"), data = edsm,
             prior = catch_priors, iter = 3000, warmup = 2000, chains = 4,
             control = list(max_treedepth = 12, adapt_delta = 0.99),
             save_pars = save_pars(all = FALSE), seed = 123)

  set.seed(123)
  kf <- brms::kfold(fit, folds = station_folds, compare = FALSE, save_fits = FALSE,
                    chains = 4, cores = 4, iter = 3000, warmup = 2000,
                    control = list(max_treedepth = 12, adapt_delta = 0.99))
  fit$criteria[["kfold5_StationCode"]] <- kf

  saveRDS(fit, f27)
  saveRDS(kf, file.path(FIT_EDSM, "kfold_fit_edsm_27.rds"))
  fileConn <- file(file.path(FIT_EDSM, "edsm_fits_Vfix_output.txt"), open = "a")
  cat("\n\n\nSummary of fit_edsm_27 (APPENDIX DIAGNOSTIC: fit 24 + factor(water_year)):\n\n",
      file = fileConn)
  cat(capture.output(summary(fit)), file = fileConn, sep = "\n")
  close(fileConn)
  rm(fit); gc()
}

# ---- 3. fit 27 vs fit 24 (paired ELPD) ------------------------------------------------
k24 <- readRDS(file.path(FIT_EDSM, "kfold_fit_edsm_24.rds"))
k27 <- readRDS(file.path(FIT_EDSM, "kfold_fit_edsm_27.rds"))
d   <- k27$pointwise[, "elpd_kfold"] - k24$pointwise[, "elpd_kfold"]
message(sprintf("fit 27 - fit 24: dELPD = %.2f (SE %.2f)", sum(d), sqrt(length(d)) * sd(d)))
