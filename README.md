# tps_pnas_submission

Analysis code for the Tidal Parr Survey (TPS) and Enhanced Delta Smelt Monitoring (EDSM)
juvenile Chinook salmon analysis submitted to *PNAS*. It covers data preparation, every reported
Bayesian model set, the model diagnostics, and the code that draws the final figures. Also covered is the code for the
salinity / temperature / water-year analysis appendix.

**This repository holds code only.** Survey data, fitted models (several GB of `.rds`), caches and
rendered figures are not included. All data needed for analysis and figure building is available from public web services and can be
re-downloaded with `00_data/00_download_web_data.R`.

---

## Layout and run order

Run everything from the **repository root**: open `tps_pnas_submission.Rproj` in RStudio, or
`setwd()` to this folder. Each notebook finds the root itself, so knitting from a sub-folder
also works.

```
00_data/
  00_download_web_data.R                 web data: CDEC, USGS NWIS, USGS/DWR DEM, TIGERweb
  01_build_model_datasets.Rmd            the four modelling frames (TPS x2, EDSM, fork length)
  02_tidal_stage_depth_estimation.Rmd    NOAA tides + VDatum + DEM -> per-station depth check
01_tps_catch_models/
  tps_catch_wq_models.Rmd                TPS catch, fits 0-13, PSIS-LOO (moment matched)   [main]
  tps_catch_landscape_models.Rmd         TPS catch, landscape fits 1-10, PSIS-LOO       [support]
02_edsm_catch_models/
  edsm_catch_kfold_models.Rmd            EDSM catch, fits 0-21, grouped 5-fold CV        [main]
03_tps_forklength_models/
  tps_forklength_kfold_models.Rmd        TPS fork length, fits 1-19, leave-one-station-out [main]
  fit9_family_comparison.Rmd             response-family check of the best FL model     [support]
  fit9_bimodality_diagnosis.Rmd          why the FL residual QQ is S-shaped (bimodality)  [support]
  gear_check_mamou_vs_otter.Rmd          Mamou vs otter trawl: catch and fork length     [support]
04_diagnostics/
  covariate_correlations_by_water_year.Rmd   correlation matrices, by survey x water year
  dharma_residual_diagnostics.Rmd            DHARMa residual tests (knit once per model set)
  pareto_k_moment_match_audit.Rmd            PSIS-LOO Pareto-k audit
05_figures/
  main_text_ensemble_figures.Rmd         stacked-ensemble effects + coefficient "butterfly" forests
  model_comparison_forest.Rmd            ELPD forest + stacking weights + covariate table
  habitat_edge_maps.Rmd                  edge-type and water-body maps
  pub_table_tows_all.Rmd                 effort / catch summary table
appendix_salinity_temperature_wateryear/ salinity, temperature & water-year appendix (own README)
appendix_seasonal_smooth_sensitivity/    s(water-year day) basis-dimension sensitivity appendix
  appendix_spline_sensitivity.R          refits TPS catch at k = 10, 6, 4, 3 (+ standardized) -> cache + figure
  appendix_spline_sensitivity.Rmd        appendix text, tables and Fig. S1 from the cache
R/
  source_rmd_chunks.R                    lets the appendix re-use main-text code by chunk label
```

Steps 01–03 fit the models. The fitting chunks are `eval = FALSE`, and the TPS runner refuses to
overwrite an existing fit, so knitting a model notebook only prints the model set. Run the fitting
chunks interactively when you want to refit. Steps 04–05 read the saved fits.

## Model sets reported in the main text

All catch models are negative binomial (log link) with `offset(log(Volume))` and a
`s(water_year_day)` seasonal spline. The response is the number of Chinook < 75 mm fork length per
tow. Continuous covariates are z-scored within each modelling frame.

| set | folder | n | models | random effect | predictive criterion |
|---|---|---|---|---|---|
| TPS catch (water quality, edge type) | `fits_tps_wq_Vfix/` | 818 tows | 0–13 | none | PSIS-LOO, moment matched |
| TPS catch (landscape) | `fits_tps_landscape/` | 838 tows | 1–10 | none | PSIS-LOO, moment matched |
| EDSM catch | `fits_edsm_Vfix/` | 4,241 tows | 0–21 | `(1 | StationCode)` | grouped 5-fold CV by station |
| TPS fork length (gaussian) | `fits_tps_fl_kfold/` | 134 fish | 1–19 | `(1 | StationCode)` | leave-one-station-out CV (K = 69) |

Ensemble effects are stacking-weighted mixtures across each set (`loo::stacking_weights`).
Appendix-only fits share these folders and are numbered after the reported ones: TPS 14–17 and
EDSM 22–27.

## Data

**Supplied by the survey programs (not downloaded by code)** — place in `data/`:

| file | content |
|---|---|
| `All_Surveys_Combined_v7.3.csv` | tow-level TPS + EDSM catch, effort, water quality, distances, habitat/edge class, 30-day mean Delta outflow per tow |
| `TPSCatchWithCSLength.csv` | TPS fish-level fork lengths |
| `Tidal_Parr_Complete_tides.csv`, `station_code_timestamp.csv` | TPS sample date-times |
| `Tide_Stations.csv` | NOAA CO-OPS gauges considered for water-surface elevation |

**Pulled from public web services** (code in this repository):

| source | used for | where |
|---|---|---|
| CDEC (CA DWR) — Delta Outflow Index, conductance, turbidity | appendix Fig. S1, salinity maps | `00_data/00_download_web_data.R` |
| USGS NWIS — specific conductance (00095) at 11455820, 11455780, 11455508 | appendix salinity maps | `00_data/00_download_web_data.R` |
| USGS/DWR 2016 SF Bay-Delta 10 m DEM (ScienceBase 58599681e4b01224f329b484) | bed elevation at TPS stations | `00_data/00_download_web_data.R` |
| US Census TIGERweb hydrography | water basemap for every map | `00_data/00_download_web_data.R` (also `05_figures/habitat_edge_maps.Rmd`) |
| NOAA CO-OPS tides & currents API; NOAA VDatum API | water-surface elevation, MLLW→NAVD88 | `00_data/02_tidal_stage_depth_estimation.Rmd` |

`00_data/01_build_model_datasets.Rmd` rebuilds the four modelling frames into `data/rebuilt/`. It
then checks every modelling column against the files the fits were built from. Last run: all four
frames identical (838 / 818 / 4,241 / 134 rows).

## Software

R 4.5.3. Main packages and the versions used: brms 2.23.0, rstan 2.32.7, loo 2.9.0, DHARMa 0.5.0,
posterior 1.7.0, tidyverse 2.0.0, ggplot2 4.0.3, patchwork 1.3.2, ggdist 3.3.3, marquee 1.2.1,
sf 1.1.0, terra 1.9.11, oce 1.8.3, EflowStats 5.2.0, future 1.70.0, kableExtra 1.4.0, ragg 1.5.2,
zoo 1.8.15, jsonlite 2.0.0, PerformanceAnalytics 2.1.0, lubridate 1.9.5.

Fitting is slow. A TPS catch model takes ~5 min. An EDSM model with grouped K-fold takes ~1–2 h.
A fork-length model with leave-one-station-out CV takes hours. The figure notebooks read a cache
built once from the fits, so a normal render takes about a minute.
