# Appendix: salinity, temperature and water year

Everything here supports the manuscript **appendix** only. The main-text figures, model sets and
stacked ensembles contain no salinity, temperature or water-year term.

**Why these are in an appendix.** Salinity cannot be identified separately from outflow in these
data. Outflow carries parr seaward, and it also freshens the water they are carried into. So low
salinity at positive tows is what both hypotheses predict. Temperature tracks the seasonal spline
(r = 0.75–0.95 with day of water year in every year) and was a null in EDSM. A water-year factor
would absorb the between-year outflow contrast that the models are built to estimate. The
water-year fits (TPS 17, EDSM 27) are therefore diagnostics: they are reported, never stacked.

Run the scripts from the **repository root**, in order. Each one skips work whose output already
exists (fits are never overwritten).

| script | what it does | cost |
|---|---|---|
| `A00_appendix_setup.R` | shared paths/helpers; sourced by every script below | – |
| `A01_download_appendix_web_data.R` | CDEC + USGS NWIS conductance, CDEC outflow & turbidity, TIGERweb water | 10–20 min (network) |
| `A02_build_edsm_salfix_data.R` | EDSM frame with **corrected** salinity (µS/cm → PSS-78) + `temperature_z` | seconds |
| `A03_tps_salinity_wateryear_fits.R` | TPS fits 14–16 (salinity) + 17 (water year); LOO ranking 0–16; stacking with/without salinity | ~30 min |
| `A04_edsm_salinity_temperature_wateryear_fits.R` | EDSM fits 22, 24–26 (salinity), 23 (temperature), 27 (water year), same folds as 0–21 | ~10 h |
| `A05_edsm_kfold_comparison_27.R` | grouped K-fold ranking + stacking / pseudo-BMA+ weights over EDSM 0–26 | seconds |
| `A06_salinity_model_comparison_figures.R` | `sal_model_comparison_*` figures | ~1 min |
| `A07_salinity_ensemble_figures.R` | `sal_combined_TPS_EDSM_*` conditional-effects and butterfly figures (cache-backed) | ~25 min first run, then ~1 min |
| `A08_cdec_outflow_salinity_turbidity.R` | Fig. S1 (outflow / salinity / turbidity by year) + memory and recovery statistics | <1 min |
| `A09_salinity_contour_maps.R` | Fig. S2 (salinity field at the freshest date of each year, with catch) + envelope statistics | a few min |
| `A10_appendix_tables_and_stats.R` | paired ELPD contrasts, outflow by year (m³ s⁻¹), stacked-coefficient tables | a few min |

The model runners, priors, fold vector, figure functions and plotting inputs are **not copied**
into these scripts. They are sourced, by chunk label, from the main-text notebooks (via
`R/source_rmd_chunks.R`), so the appendix fits and figures are built exactly like the reported
ones.

Figures go to `output_figs_appendix/` (prefix `sal_`). Intermediate tables go to
`data/derived_appendix/`.
