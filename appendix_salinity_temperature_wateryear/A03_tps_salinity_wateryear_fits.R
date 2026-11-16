# =============================================================================
# A03_tps_salinity_wateryear_fits.R
# -----------------------------------------------------------------------------
# TPS appendix fits, on the same 818-tow frame, priors, sampler settings and
# moment-matched PSIS-LOO as the reported fits 0-13:
#   fit 14 = base + salinity                     (one-at-a-time screen, like fits 3-6)
#   fit 15 = base + edge type + salinity
#   fit 16 = base + edge type + Secchi + salinity
#   fit 17 = fit 16 + factor(water_year)         (water-year DIAGNOSTIC only)
# base = out30_z * distance_m_z + distance_along_channel_z
# every model also carries s(water_year_day) + offset(log(Volume)).
#
# The model runner, priors, data loading and LOO helpers are NOT copied here: they
# are sourced from the main-text notebook (01_tps_catch_models/tps_catch_wq_models.Rmd)
# so the appendix fits are guaranteed to be built exactly like the reported ones.
#
# Then:
#   * loo_compare over fits 0-16  -> fit_tps_wq0-16_catch_looic_comparisons.csv
#   * stacking & pseudo-BMA+ weights computed on the SAME loo objects over 0-13
#     (reported set) and over 0-16 (salinity added), so the difference isolates
#     what salinity does to the ensemble
#                                  -> fit_tps_wq0-16_stacking_weights_salinity_compare.csv
#   Fit 17 is reported against fit 16 only (paired ELPD); it is never stacked.
#
# Inputs : data/tow_12_6_tps_wq_nona.csv, fits_tps_wq_Vfix/fit_tps_wq0..13.rds
# Outputs: fits_tps_wq_Vfix/fit_tps_wq14..17.rds + slim LOO sidecars loo_tps_wq14..17.rds,
#          summaries appended to fits_tps_wq_Vfix/tps_wq_fits_Vfix_output.txt
# Cost   : ~5-10 min per fit. Existing fits are never overwritten.
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

# Runner, priors (b ~ normal(0, 5), shape ~ gamma(2, 0.1)), data and LOO helpers
source_rmd_chunks(NB_TPS_WQ, c("libraries", "config", "load-data", "priors",
                               "loo-helpers", "fit-defs-runner", "load-loo-fn"))
stopifnot(nrow(tow_12_6_tps_wq_nona) == 818)

# ---- 1. salinity fits 14-16 ---------------------------------------------------------
fit_defs_appendix <- tibble::tribble(
  ~num, ~round,  ~fixef,                                                    ~note,
  14L,  "appx",  paste(BASE_FIXEF, "+ salinity_z"),                         "base + salinity",
  15L,  "appx",  paste(BASE_FIXEF, "+ Edge_Type + salinity_z"),             "base + edge + salinity",
  16L,  "appx",  paste(BASE_FIXEF, "+ Edge_Type + secchi_z + salinity_z"),  "base + edge + Secchi + salinity"
)
fit_defs_all <- dplyr::bind_rows(fit_defs, fit_defs_appendix)

# Skips any fit already on disk (overwrite = FALSE is the runner's default).
run_tps_wq_fits(nums = 14:16, defs = fit_defs_all, save_sidecar = TRUE)

# ---- 2. water-year diagnostic fit 17 = fit 16 + factor(water_year) ------------------
# Reference level is the first water year (2018). Not added to any stacked set.
fit17_file <- file.path(FIT_DIR, "fit_tps_wq17.rds")
if (!file.exists(fit17_file)) {
  tps <- tow_12_6_tps_wq_nona %>% mutate(wy_f = factor(water_year))
  form17 <- n_l_75 ~ out30_z * distance_m_z + distance_along_channel_z + Edge_Type +
    secchi_z + salinity_z + wy_f + s(water_year_day) + offset(log(Volume))

  fit <- brm(formula = form17, family = negbinomial(link = "log"), data = tps,
             prior = catch_priors, iter = 3000, warmup = 2000, chains = 4,
             control = list(max_treedepth = 12, adapt_delta = 0.99),
             save_pars = save_pars(all = TRUE), seed = 123)

  # same LOO treatment as fits 0-16: PSIS-LOO, then moment matching
  l <- loo(fit, save_psis = TRUE)
  message("fit 17: Pareto k > 0.7 before moment matching: ", sum(l$diagnostics$pareto_k > 0.7))
  fit$criteria$loo <- loo_moment_match(fit, loo = l, cores = 1)

  saveRDS(fit, fit17_file)
  saveRDS(slim_loo(fit$criteria$loo), file.path(FIT_DIR, "loo_tps_wq17.rds"))
  fileConn <- file(OUT_TXT, open = "a")
  cat("\n\n\nSummary of fit_tps_wq17 (APPENDIX DIAGNOSTIC: fit 16 + factor(water_year)):\n\n",
      file = fileConn)
  cat(capture.output(summary(fit)), file = fileConn, sep = "\n")
  close(fileConn)
  rm(fit); gc()
}

# ---- 3. ELPD ranking 0-16 and stacking with vs without salinity ---------------------
loo_all <- load_loo_tps_wq(0:16, prefer_sidecar = TRUE)
stopifnot(length(loo_all) == 17)

cmp <- loo_compare(loo_all)
print(cmp)
write.csv(as.data.frame(cmp), file.path(FIT_DIR, "fit_tps_wq0-16_catch_looic_comparisons.csv"))

set.seed(1)
base_set <- paste0("fit_tps_wq", 0:13)
w_base <- loo_model_weights(loo_all[base_set], method = "stacking")
w_sal  <- loo_model_weights(loo_all,           method = "stacking")
p_base <- loo_model_weights(loo_all[base_set], method = "pseudobma")
p_sal  <- loo_model_weights(loo_all,           method = "pseudobma")

pad <- function(w) c(as.numeric(w), rep(NA, length(loo_all) - length(w)))
stack_cmp <- tibble::tibble(
  model          = names(loo_all),
  num            = as.integer(sub("^fit_tps_wq", "", names(loo_all))),
  stack_no_sal   = pad(w_base),          # fits 14-16 absent from the reported set
  stack_with_sal = as.numeric(w_sal),
  pbma_no_sal    = pad(p_base),
  pbma_with_sal  = as.numeric(p_sal)) %>%
  mutate(stack_change = round(stack_with_sal - stack_no_sal, 4)) %>%
  left_join(fit_defs_all %>% select(num, fixef, note), by = "num") %>%
  arrange(desc(stack_with_sal))
print(stack_cmp, n = Inf, width = Inf)
write.csv(stack_cmp, file.path(FIT_DIR, "fit_tps_wq0-16_stacking_weights_salinity_compare.csv"),
          row.names = FALSE)

# ---- 4. fit 17 vs fit 16: what a year factor does (paired ELPD) ----------------------
l16 <- load_loo_tps_wq(16, prefer_sidecar = TRUE)[[1]]
l17 <- readRDS(file.path(FIT_DIR, "loo_tps_wq17.rds"))
d   <- l17$pointwise[, "elpd_loo"] - l16$pointwise[, "elpd_loo"]
message(sprintf("fit 17 - fit 16: dELPD = %.2f (SE %.2f)", sum(d), sqrt(length(d)) * sd(d)))
