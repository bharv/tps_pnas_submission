# =============================================================================
# 00_download_web_data.R
# -----------------------------------------------------------------------------
# Re-pulls EVERY dataset this analysis obtained from a public web service, so the
# repository can stay code-only. Nothing here is needed to re-render figures from
# existing fits; it is needed to rebuild the derived inputs from scratch.
#
# Run from the repository root (the folder holding the .Rproj):
#     source("00_data/00_download_web_data.R")
#     download_all()                 # everything below, skipping files already on disk
#   or call one section at a time, e.g. download_cdec_timeseries().
#
# What comes from where
# ---------------------
#  (1) CDEC (California Data Exchange Center, CA Dept. of Water Resources)
#        Delta Outflow Index (DTO), specific conductance and turbidity series.
#        Used ONLY by the salinity / temperature / water-year APPENDIX
#        (appendix figure S1 and the salinity contour maps).
#  (2) USGS NWIS instantaneous values, parameter 00095 (specific conductance)
#        at Carquinez Br. (11455820), Benicia Br. (11455780), Van Sickle Is.
#        (11455508). Appendix salinity contour maps only.
#  (3) USGS / DWR 2016 San Francisco Bay-Delta 10 m bathymetric/topographic DEM
#        (ScienceBase item 58599681e4b01224f329b484, NAVD88 metres, UTM 10N).
#        Used to estimate bed elevation at every TPS station ->
#        data/tps_bed_navd88.csv (read by 00_data/02_tidal_stage_depth_estimation.Rmd).
#  (4) US Census TIGERweb areal hydrography -> data/spatial/delta_hydro_tiger.rds
#        Vector water basemap for every map (habitat maps, salinity maps).
#  (5) NOAA CO-OPS tides & currents API + NOAA VDatum API
#        Water-surface elevation at the sampled date-times and MLLW->NAVD88
#        offsets. These calls are made (and cached in memory per station/day)
#        INSIDE 00_data/02_tidal_stage_depth_estimation.Rmd, which is the single
#        place that code lives (API endpoints in its `setup` chunk, request
#        functions in `helpers`, the pull itself in `compute`). Nothing to do here.
#
# NOT downloaded by code (supplied by the survey programs, see README):
#   data/All_Surveys_Combined_v7.3.csv  -- tow-level TPS + EDSM catch/effort/
#        water quality, with OUT_Moving_Mean_30d (30-d trailing mean of the
#        DAYFLOW/DTO Delta outflow) already attached per tow;
#   data/TPSCatchWithCSLength.csv, data/EDSM_KDTR.csv -- fish-level lengths;
#   data/Tidal_Parr_Complete_tides.csv, data/station_code_timestamp.csv,
#   data/Tide_Stations.csv -- sample date-times and the NOAA gauge list.
#
# Service notes (as of 2026-11): CDEC's CSVDataServlet and USGS's legacy
# waterservices.usgs.gov "iv" endpoint both worked when this was written. The
# USGS endpoint returns HTTP 503 for long windows, hence the month-by-month
# requests with retries. If USGS retires the legacy service, the same data are
# available from the USGS Water Data API (https://api.waterdata.usgs.gov).
# =============================================================================

suppressPackageStartupMessages({
  library(utils)
})

# All paths are relative to the repository root.
ROOT <- getwd()   # run from the repository root (the folder holding the .Rproj)
WEB  <- file.path(ROOT, "data", "web_downloads")

# Small helper: download `url` to `dest` unless a non-empty copy already exists.
# Retries a few times because both CDEC and USGS occasionally drop requests.
fetch <- function(url, dest, tries = 3, pause = 4, timeout = 300) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(timeout, getOption("timeout"))); on.exit(options(old))
  for (i in seq_len(tries)) {
    ok <- tryCatch({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok && file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
    Sys.sleep(pause)
  }
  warning("download failed after ", tries, " tries: ", url)
  invisible(NULL)
}

# CDEC CSV servlet URL for one station / sensor / duration over a date window.
# dur_code: D = daily, H = hourly, E = event (15-min).
cdec_url <- function(station, sensor, dur, start, end) {
  sprintf(paste0("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?",
                 "Stations=%s&SensorNums=%s&dur_code=%s&Start=%s&End=%s"),
          station, sensor, dur, start, end)
}

# -----------------------------------------------------------------------------
# (1a) CDEC series for appendix figure S1 (outflow / salinity / turbidity by year)
#      Sensor 23 = Delta Outflow Index (cfs); 100 = EC (uS/cm); 324 = specific
#      conductance @25 C; 27 = turbidity (NTU).
# -----------------------------------------------------------------------------
download_cdec_timeseries <- function(start = "2016-10-01", end = "2021-07-15") {
  specs <- c("DTO:23:D",                                   # Delta outflow
             "MAL:100:D", "BDL:324:D", "PCT:100:D",        # conductance -> salinity
             "MRZ:100:D", "CLL:100:D", "GZL:100:E", "GZL:100:H",
             "MAL:27:H",  "GZL:27:E",  "GZL:27:H",         # turbidity
             "MRZ:27:E",  "MRZ:27:H")
  out_dir <- file.path(WEB, "cdec_timeseries")
  for (s in specs) {
    p <- strsplit(s, ":")[[1]]
    fetch(cdec_url(p[1], p[2], p[3], start, end),
          file.path(out_dir, sprintf("%s_%s_%s.csv", p[1], p[2], p[3])))
  }
  message("CDEC time series -> ", out_dir)
}

# -----------------------------------------------------------------------------
# (1b) CDEC conductance stations for the appendix salinity contour maps.
#      The Suisun Marsh interior stations (BDL, MSL, GYS, CYG, TEA, IBS, FLT, GZS,
#      GZM) are downloaded for completeness but EXCLUDED from the interpolation in
#      appendix A10 (they sit behind the Suisun Marsh Salinity Control Gates).
# -----------------------------------------------------------------------------
download_cdec_salinity_field <- function(start = "2016-10-15", end = "2021-07-05") {
  specs <- c("MRZ:100:D", "PCT:100:D", "MAL:100:D", "CLL:100:D", "BDL:324:D",
             "GZL:100:E", "GZM:100:E", "HUN:100:E", "GZS:100:E", "RYC:100:E",
             "HON:100:E", "MSL:100:E", "ANC:100:E", "SSI:100:E", "RIV:100:E",
             "GYS:100:E", "CYG:100:E", "FLT:100:E", "TEA:100:E", "IBS:100:E")
  out_dir <- file.path(WEB, "salinity_field", "cdec")
  for (s in specs) {
    p <- strsplit(s, ":")[[1]]
    fetch(cdec_url(p[1], p[2], p[3], start, end),
          file.path(out_dir, sprintf("cdec_%s_%s_%s.csv", p[1], p[2], p[3])))
  }
  message("CDEC salinity-field series -> ", out_dir)
}

# -----------------------------------------------------------------------------
# (2) USGS NWIS instantaneous specific conductance (00095), Nov-Jun of WY2017-2021,
#     one request per calendar month (the service 503s on long windows).
# -----------------------------------------------------------------------------
download_usgs_conductance <- function(sites = c("11455820", "11455780", "11455508"),
                                      water_years = 2017:2021) {
  out_dir <- file.path(WEB, "salinity_field", "usgs_monthly")
  # Nov-Dec of the previous calendar year + Jan-Jun of the water year
  months <- unlist(lapply(water_years, function(wy)
    c(sprintf("%d-%02d", wy - 1L, 11:12), sprintf("%d-%02d", wy, 1:6))))
  for (site in sites) for (m in months) {
    first <- as.Date(paste0(m, "-01"))
    last  <- seq(first, by = "month", length.out = 2)[2] - 1
    url <- sprintf(paste0("https://waterservices.usgs.gov/nwis/iv/?format=rdb&sites=%s",
                          "&parameterCd=00095&startDT=%s&endDT=%s"), site, first, last)
    fetch(url, file.path(out_dir, sprintf("%s_%s.rdb", site, m)), timeout = 180)
    Sys.sleep(1)                       # be polite to the service
  }
  message("USGS conductance -> ", out_dir)
}

# -----------------------------------------------------------------------------
# (3) 2016 SF Bay-Delta 10 m DEM -> bed elevation (NAVD88 m) at each TPS station.
#     The DEM is ~1.4 GB, so it goes to a temporary folder and is deleted after
#     the point extraction unless keep_dem = TRUE. Only the small per-station
#     table is kept: data/tps_bed_navd88.csv (StationCode, bed_navd88_m).
#     bed_navd88_m = the EXACT DEM cell under the station coordinate (no buffer;
#     the GPS coordinates are precise, so the cell is used even on a bank).
# -----------------------------------------------------------------------------
download_dem_and_extract_beds <- function(stations_csv = file.path(ROOT, "data", "tps_forklength_data.csv"),
                                          out_csv      = file.path(ROOT, "data", "tps_bed_navd88.csv"),
                                          keep_dem     = FALSE) {
  if (file.exists(out_csv)) { message("exists: ", out_csv); return(invisible(out_csv)) }
  stopifnot(requireNamespace("terra", quietly = TRUE),
            requireNamespace("jsonlite", quietly = TRUE))

  # Resolve the current download URLs from the ScienceBase item record rather
  # than hard-coding them (ScienceBase file URLs are content hashes).
  item  <- jsonlite::fromJSON(paste0("https://www.sciencebase.gov/catalog/item/",
                                     "58599681e4b01224f329b484?format=json&fields=files"))
  files <- item$files
  asc_url <- files$url[files$name == "sfbaydeltadem10m2016.asc"]
  prj_url <- files$url[files$name == "sfbaydeltadem10m2016.prj"]
  stopifnot(length(asc_url) == 1)

  dem_dir <- if (keep_dem) file.path(WEB, "dem") else file.path(tempdir(), "dem")
  asc <- fetch(asc_url, file.path(dem_dir, "sfbaydeltadem10m2016.asc"), timeout = 3600)
  fetch(prj_url, file.path(dem_dir, "sfbaydeltadem10m2016.prj"))

  dem <- terra::rast(asc)
  terra::crs(dem) <- "EPSG:26910"                 # NAD83 / UTM zone 10N, metres

  st <- utils::read.csv(stations_csv)
  st <- unique(st[, c("StationCode", "Station_Longitude", "Station_Latitude")])
  pts <- terra::vect(st, geom = c("Station_Longitude", "Station_Latitude"), crs = "EPSG:4326")
  pts <- terra::project(pts, "EPSG:26910")
  bed <- terra::extract(dem, pts)[, 2]            # exact cell, no buffer

  utils::write.csv(data.frame(StationCode = st$StationCode, bed_navd88_m = round(bed, 2)),
                   out_csv, row.names = FALSE)
  if (!keep_dem) unlink(dem_dir, recursive = TRUE)
  message("bed elevations -> ", out_csv)
}

# -----------------------------------------------------------------------------
# (4) TIGERweb areal hydrography for the Suisun Bay / west Delta study area.
#     Identical query to 05_figures/habitat_edge_maps.Rmd (which also fetches it
#     on first run). Geometry only; attributes are not used.
# -----------------------------------------------------------------------------
download_tiger_hydro <- function(out_rds = file.path(ROOT, "data", "spatial", "delta_hydro_tiger.rds")) {
  if (file.exists(out_rds)) { message("exists: ", out_rds); return(invisible(out_rds)) }
  stopifnot(requireNamespace("sf", quietly = TRUE))
  bbox <- c(xmin = -122.60, ymin = 37.92, xmax = -121.68, ymax = 38.32)
  url <- paste0(
    "https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/Hydro/",
    "MapServer/1/query?",
    "geometry=", paste(bbox[c("xmin", "ymin", "xmax", "ymax")], collapse = ","),
    "&geometryType=esriGeometryEnvelope&inSR=4326&outSR=4326",
    "&spatialRel=esriSpatialRelIntersects&where=1%3D1&outFields=NAME",
    "&returnGeometry=true&f=geojson")
  tmp <- tempfile(fileext = ".geojson")
  fetch(url, tmp)
  h <- sf::st_make_valid(sf::st_read(tmp, quiet = TRUE))
  h <- suppressWarnings(sf::st_crop(h, bbox))
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sf::st_geometry(h), out_rds)
  message("TIGER hydrography -> ", out_rds)
}

# -----------------------------------------------------------------------------
# Everything, in dependency order. Each step skips files already present.
# -----------------------------------------------------------------------------
download_all <- function() {
  download_tiger_hydro()
  download_dem_and_extract_beds()
  download_cdec_timeseries()
  download_cdec_salinity_field()
  download_usgs_conductance()
  message("\nNOAA tide / VDatum data are pulled (and cached) by ",
          "00_data/02_tidal_stage_depth_estimation.Rmd when it is run.")
  invisible(TRUE)
}
