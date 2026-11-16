# =============================================================================
# A07_salinity_ensemble_figures.R
# -----------------------------------------------------------------------------
# Appendix versions of the two stacked-ensemble figures with the salinity /
# temperature models included (TPS fits 0-16, EDSM fits 0-26; water-year fits
# 17 / 27 excluded):
#   * combined TPS + EDSM ensemble conditional effects, original units, with a
#     salinity row added
#   * the "butterfly" coefficient forest (alt4 layout) with a salinity row
#
# Two stages, both re-using the function and input chunks of
# 05_figures/main_text_ensemble_figures.Rmd:
#   1. CACHE (slow, ~20-30 min, once): stream every fit, keep its b_* draws +
#      predictive criterion and accumulate the stacked conditional effects
#      -> cache/sal_fast_figure_cache_kfold.rds
#   2. RENDER (~1 min): draw both figures from the cache.
# Set REBUILD_CACHE <- TRUE to force stage 1 after refitting.
#
# EDSM salinity is back-transformed to psu with the mean/SD of the CORRECTED
# salinity (data/tow_12_6_edsm_wq_salfix.csv), never the original column.
#
# Outputs: output_figs_appendix/sal_combined_TPS_EDSM_conditional_effects_origscale_kfold.*
#          output_figs_appendix/sal_combined_TPS_EDSM_butterfly_forest_alt4_kfold.*
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

source_rmd_chunks(NB_FIGS, c("setup", "raw-data", "fit-registry", "helpers", "compute-funcs",
                             "ccep-cached", "combined-inputs", "facet-cond-inputs",
                             "origscale-inputs", "butterfly-alt4-fn", "alt4-inputs"))

REBUILD_CACHE <- FALSE
SAL_CACHE     <- "cache/sal_fast_figure_cache_kfold.rds"

# ---- expanded model sets ---------------------------------------------------------------
edsm_nums            <- 0:26
fit_names_edsm_catch <- paste0("fit_edsm_catch", edsm_nums)
fit_paths_edsm       <- setNames(file.path(FIT_EDSM, paste0("fit_edsm_", edsm_nums, ".rds")),
                                 fit_names_edsm_catch)
tps_nums             <- 0:16
fit_names_tps_catch  <- paste0("fit_tps_catch", tps_nums)
fit_paths_tps        <- setNames(file.path(FIT_TPS, paste0("fit_tps_wq", tps_nums, ".rds")),
                                 fit_names_tps_catch)
fit_paths_all <- c(fit_paths_tps, fit_paths_edsm)

# ---- add a salinity row to the conditional-effects layout ------------------------------
edsm_salfix <- read.csv(EDSM_SALFIX, stringsAsFactors = FALSE)

facet_spec <- dplyr::bind_rows(facet_spec, tibble::tribble(
  ~row, ~survey, ~covariate,
  7L,   "tps",   "b_salinity_z",
  7L,   "edsm",  "b_salinity_z"))
combined_cond_plot_inputs$facet_spec           <- facet_spec
combined_cond_plot_inputs$tps_selected_params  <- dplyr::pull(dplyr::filter(facet_spec, survey == "tps"),  covariate)
combined_cond_plot_inputs$edsm_selected_params <- dplyr::pull(dplyr::filter(facet_spec, survey == "edsm"), covariate)
combined_cond_plot_inputs$param_labels <- c(combined_cond_plot_inputs$param_labels,
                                            "b_salinity_z" = "Salinity (z-score)")
combined_cond_plot_inputs$main_title    <- "TPS vs EDSM: ensemble conditional effects (salinity models included)"
combined_cond_plot_inputs$output_height <- 13

tps_z_stats_origscale$salinity_z  <- list(mean = mean(tow_12_6_tps_wq_nona$Salinity, na.rm = TRUE),
                                          sd   = sd(tow_12_6_tps_wq_nona$Salinity,   na.rm = TRUE), digits = 1)
edsm_z_stats_origscale$salinity_z <- list(mean = mean(edsm_salfix$Salinity, na.rm = TRUE),
                                          sd   = sd(edsm_salfix$Salinity,   na.rm = TRUE), digits = 1)

combined_cond_plot_inputs_origscale <- modifyList(combined_cond_plot_inputs, list(
  param_labels = c(combined_cond_plot_inputs_origscale$param_labels, "b_salinity_z" = "Salinity (PSU)"),
  output_file  = file.path(OUT_APPX, "sal_combined_TPS_EDSM_conditional_effects_origscale_kfold")))

# ---- stage 1: cache -------------------------------------------------------------------
if (REBUILD_CACHE || !file.exists(SAL_CACHE)) {
  stopifnot(all(file.exists(fit_paths_all)))
  edsm_bundle   <- readRDS(file.path(FIT_EDSM, "kfold5_station_criteria_bundle_27.rds"))
  edsm_criteria <- setNames(edsm_bundle[paste0("fit_edsm_", edsm_nums)], fit_names_edsm_catch)

  catch_be <- extract_beta_and_criterion(c(fit_names_tps_catch, fit_names_edsm_catch), fit_paths_all)
  beta_draws_cache   <- catch_be$draws
  tps_criteria_cache <- catch_be$criteria[fit_names_tps_catch]

  ce_args <- combined_cond_plot_inputs
  cond_catch <- list(
    tps_result  = compute_survey_ce(fit_names_tps_catch, fit_paths_all, ce_args$tps_selected_params,
                    tps_criteria_cache, resolution = ce_args$resolution,
                    re_formula = ce_args$re_formula, probs = ce_args$probs,
                    epred_fn = epred_marginal_lognormal),
    edsm_result = compute_survey_ce(fit_names_edsm_catch, fit_paths_all, ce_args$edsm_selected_params,
                    edsm_criteria, resolution = ce_args$resolution,
                    re_formula = ce_args$re_formula, probs = ce_args$probs,
                    epred_fn = epred_marginal_lognormal))

  saveRDS(list(beta_draws_cache = beta_draws_cache, tps_criteria_cache = tps_criteria_cache,
               edsm_criteria = edsm_criteria, cond_catch = cond_catch,
               inputs = combined_cond_plot_inputs, inputs_os = combined_cond_plot_inputs_origscale,
               tps_z = tps_z_stats_origscale, edsm_z = edsm_z_stats_origscale,
               meta = list(built_at = Sys.time(), edsm_nums = edsm_nums, tps_nums = tps_nums)),
          SAL_CACHE)
}

# ---- stage 2: render ------------------------------------------------------------------
fc <- readRDS(SAL_CACHE)
cond_catch             <- fc$cond_catch
beta_draws_cache       <- fc$beta_draws_cache
tps_criteria_cache     <- fc$tps_criteria_cache
edsm_criteria          <- fc$edsm_criteria
tps_z_stats_origscale  <- fc$tps_z
edsm_z_stats_origscale <- fc$edsm_z
ci <- combined_cond_plot_inputs_origscale    # output path above, layout from the cache build

# (a) combined conditional effects, original units
ens <- combined_cond_effects_plot_cached(
  tps_result = cond_catch$tps_result, edsm_result = cond_catch$edsm_result,
  wide_y_params = ci$wide_y_params, narrow_y_params = ci$narrow_y_params,
  param_labels = ci$param_labels, y_axis_label = ci$y_axis_label,
  tps_title = ci$tps_title, edsm_title = ci$edsm_title, main_title = ci$main_title,
  show_facet_subtitles = ci$show_facet_subtitles, legend_position = ci$legend_position,
  tps_z_stats = tps_z_stats_origscale, edsm_z_stats = edsm_z_stats_origscale,
  facet_spec = ci$facet_spec)
of <- ci$output_file
ggsave(paste0(of, ".pdf"),  ens, device = "pdf",  width = 6.5, height = 13)
ggsave(paste0(of, ".tiff"), ens, device = "tiff", dpi = 600, compression = "lzw", width = 6.5, height = 13)
ggsave(paste0(of, "_preview.png"), ens, device = "png", dpi = 150, width = 6.5, height = 13)

# (b) butterfly coefficient forest, alt4 layout, salinity row inserted after Secchi
combined_param_labels <- c(combined_param_labels, "b_salinity_z" = "Salinity")
combined_preferred_param_order <- append(combined_preferred_param_order, "b_salinity_z",
  after = which(combined_preferred_param_order == "b_secchi_z"))

lab4 <- combined_param_labels
for (w in c("Distance", "Dissolved", "Emergent", "Terminal"))
  lab4 <- sub(paste0(w, " "), paste0(w, "\n"), lab4, fixed = TRUE)
lab4 <- center_md_lines(lab4, family = "Helvetica", size = 8.5)
lab4 <- gsub("\n", "  \n", lab4, fixed = TRUE)     # markdown hard breaks for marquee
names(lab4) <- names(combined_param_labels)

bi <- modifyList(butterfly_forest_alt4_inputs, list(
  param_labels = lab4, preferred_param_order = combined_preferred_param_order,
  output_file  = file.path(OUT_APPX, "sal_combined_TPS_EDSM_butterfly_forest_alt4_kfold")))

bf <- combined_cluster_butterfly_plot_alt4(
  tps_fit_names = fit_names_tps_catch, edsm_fit_names = fit_names_edsm_catch,
  fit_paths = fit_paths_all, tps_criteria = tps_criteria_cache,
  edsm_criteria = edsm_criteria, draws_cache = beta_draws_cache,
  patterns_to_exclude = bi$patterns_to_exclude,
  preferred_param_order = bi$preferred_param_order, param_labels = bi$param_labels,
  survey_colors = bi$survey_colors, declutter_threshold = bi$declutter_threshold,
  min_weight = bi$min_weight, x_axis_label = bi$x_axis_label, x_limits = bi$x_limits,
  slab_scale = bi$slab_scale, band_half = bi$band_half, slab_n = bi$slab_n,
  alpha_range = bi$alpha_range, median_tick_len = bi$median_tick_len,
  outer_pad_frac_top = bi$outer_pad_frac_top, outer_pad_frac_bottom = bi$outer_pad_frac_bottom,
  legend_on_param = bi$legend_on_param, legend_x = bi$legend_x,
  legend_y_nudge = bi$legend_y_nudge, base_family = bi$base_family,
  annot_size = bi$annot_size, axis_text_y_size = bi$axis_text_y_size,
  axis_text_x_size = bi$axis_text_x_size, axis_title_size = bi$axis_title_size,
  legend_text_size = bi$legend_text_size, annot_x_gap = bi$annot_x_gap,
  tps_annot_angle = bi$tps_annot_angle, edsm_annot_angle = bi$edsm_annot_angle,
  tps_annot_hjust = bi$tps_annot_hjust, edsm_annot_hjust = bi$edsm_annot_hjust,
  tps_annot_x = bi$tps_annot_x, edsm_annot_x = bi$edsm_annot_x,
  tps_annot_y = bi$tps_annot_y, edsm_annot_y = bi$edsm_annot_y, band_gap = bi$band_gap)

ggsave(paste0(bi$output_file, ".pdf"),  bf, device = "pdf", width = bi$output_width, height = bi$output_height)
ggsave(paste0(bi$output_file, ".tiff"), bf, device = "tiff", dpi = bi$output_dpi, compression = "lzw",
       width = bi$output_width, height = bi$output_height)
ggsave(paste0(bi$output_file, "_preview.png"), bf, device = "png", dpi = 200,
       width = bi$output_width, height = bi$output_height)
