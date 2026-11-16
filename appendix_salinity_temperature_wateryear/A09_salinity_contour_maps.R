# =============================================================================
# A09_salinity_contour_maps.R                       (appendix Figure S2 + stats)
# -----------------------------------------------------------------------------
# Maps the interpolated surface-salinity field across the study area at the
# freshest moment of each water year (the date of minimum salinity at the
# westernmost station, Carquinez Br.), with every TPS and EDSM tow of that water
# year overlaid (catch sized by fish per tow, zero-catch tows as small crosses) and
# the three "catch-enclosing" contours (0.1, 1 and 6 psu, labelled with the share
# of all fish caught at or below each).
#
# Interpolation: inverse-distance weighting with IN-WATER distance (shortest path
# through a 100 m water raster built from TIGER hydrography, land as a barrier),
# weights 1 / (d + 50 m)^2.5, from 14 open-estuary conductance stations (CDEC +
# USGS), PSS-78 at 25 C. Suisun Marsh interior stations are excluded (behind the
# Salinity Control Gates). West of Carquinez Br. the field is extrapolated and is
# shaded as such.
#
# Parts
#   1. daily salinity from the CDEC conductance downloads  -> cdec_daily_salinity.rds
#   2. station table, mapped dates, water mask, in-water distance surfaces
#                                                          -> sal_field_inputs.rds
#   3. Figure S2, pooled symbols (clean + full-contour "topo" variant)
#   4. Figure S2, TPS vs EDSM symbols (clean + topo)
#   5. envelope statistics (area / westward reach of each contour, % tows & fish inside)
#   6. % of fish at or below each threshold, by survey
# (intermediates in data/derived_appendix/, figures in output_figs_appendix/sal_*)
#
# Inputs: data/web_downloads/salinity_field/{cdec,usgs_monthly}/ (A01),
#         data/spatial/delta_hydro_tiger.rds, data/tow_12_6_tps_wq_nona.csv,
#         data/tow_12_6_edsm_wq_salfix.csv (A02)
# Cost  : a few minutes (the 14 in-water distance surfaces dominate)
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

# =============================================================================
# Part 1 - daily salinity at the CDEC stations
# =============================================================================
# Build a daily salinity table (PSS-78) from the downloaded CDEC conductance series.
suppressMessages({library(tidyverse); library(oce)})
CD <- file.path(WEB, "salinity_field", "cdec")

files <- list.files(CD, pattern = "^cdec_.*\\.csv$", full.names = TRUE)
cat("cdec files:", length(files), "\n")

read_one <- function(f) {
  x <- suppressWarnings(read_csv(f, col_types = cols(.default = col_character())))
  if (!"VALUE" %in% names(x) || nrow(x) == 0) return(NULL)
  tibble(stn  = x$STATION_ID,
         date = as.Date(substr(x$`DATE TIME`, 1, 8), format = "%Y%m%d"),
         ec   = suppressWarnings(as.numeric(x$VALUE))) %>%
    filter(!is.na(date), !is.na(ec), ec > 0, ec < 1e6)
}

raw <- map_dfr(files, read_one)
cat("rows:", nrow(raw), " stations:", paste(sort(unique(raw$stn)), collapse = " "), "\n")

# daily MEDIAN, not mean: several CDEC EC series carry isolated sensor spikes
# (e.g. Rio Vista, a freshwater station, has single readings implying ~50 psu).
daily <- raw %>% group_by(stn, date) %>%
  summarise(ec = median(ec, na.rm = TRUE), n_obs = n(), .groups = "drop") %>%
  mutate(salinity = oce::swSCTp(conductivity = ec/1000,
                                temperature  = rep(25, length(ec)),
                                pressure     = rep(0,  length(ec)),
                                conductivityUnit = "mS/cm", eos = "unesco")) %>%
  # drop anything above seawater: physically impossible anywhere in this estuary
  filter(salinity <= 35)

wy_of <- function(d) ifelse(as.integer(format(d, "%m")) >= 10,
                            as.integer(format(d, "%Y")) + 1L, as.integer(format(d, "%Y")))
inwin <- function(d) { m <- as.integer(format(d, "%m")); m >= 11 | m <= 6 }

daily <- daily %>% mutate(wy = wy_of(date)) %>% filter(wy %in% 2017:2021, inwin(date))

cat("\n--- days of data per station per water year ---\n")
print(daily %>% count(stn, wy) %>% pivot_wider(names_from = wy, values_from = n, values_fill = 0) %>%
        as.data.frame())

cat("\n--- salinity range by station ---\n")
print(daily %>% group_by(stn) %>%
        summarise(n = n(), min = round(min(salinity),2), med = round(median(salinity),2),
                  max = round(max(salinity),2)) %>% as.data.frame())

saveRDS(daily, file.path(DERIVED, "cdec_daily_salinity.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 2 - salinity field inputs (stations, dates, water mask, in-water distances)
# =============================================================================
# Build the per-water-year salinity field (in-water IDW) used by the contour map.
# Saves a tidy raster data frame + station table + chosen dates; the plotting
# script reads only this, so the slow part runs once.
suppressMessages({library(tidyverse); library(sf); library(terra); library(oce)})

## ============================ 1. station salinity ==========================
cdec <- readRDS(file.path(DERIVED, "cdec_daily_salinity.rds"))      # stn, date, salinity, wy

read_usgs_chunks <- function(site) {
  fs <- list.files(file.path(WEB, "salinity_field", "usgs_monthly"), pattern = paste0("^", site, "_.*\\.rdb$"),
                   full.names = TRUE)
  if (!length(fs)) return(NULL)
  map_dfr(fs, function(f) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[!startsWith(ln, "#")]
    if (length(ln) < 3) return(NULL)
    hdr <- strsplit(ln[1], "\t")[[1]]
    # prefer the UPPER (near-surface) conductance series
    cols <- grep("_00095$", hdr, value = TRUE)
    if (!length(cols)) return(NULL)
    pick <- cols[1]
    df <- read_tsv(paste(ln[-2], collapse = "\n"), col_types = cols(.default = col_character()),
                   progress = FALSE)
    tibble(date = as.Date(substr(df$datetime, 1, 10)),
           ec   = suppressWarnings(as.numeric(df[[pick]])))
  }) %>% filter(!is.na(date), !is.na(ec), ec > 0, ec < 1e6)
}

USGS_UPPER_TS <- c("11455820" = "16238", "11455780" = "16223", "11455508" = "164314")
usgs <- map_dfr(names(USGS_UPPER_TS), function(site) {
  x <- read_usgs_chunks(site)
  if (is.null(x) || !nrow(x)) { message("no USGS data for ", site); return(NULL) }
  x %>% group_by(date) %>% summarise(ec = median(ec), n_obs = n(), .groups = "drop") %>%
    mutate(stn = site)
})

if (nrow(usgs)) {
  usgs <- usgs %>%
    mutate(salinity = oce::swSCTp(conductivity = ec/1000,
                                  temperature = rep(25, length(ec)),
                                  pressure = rep(0, length(ec)),
                                  conductivityUnit = "mS/cm", eos = "unesco")) %>%
    filter(salinity <= 35) %>%
    mutate(wy = ifelse(as.integer(format(date, "%m")) >= 10,
                       as.integer(format(date, "%Y")) + 1L, as.integer(format(date, "%Y")))) %>%
    filter(wy %in% 2017:2021,
           { m <- as.integer(format(date, "%m")); m >= 11 | m <= 6 })
}

sal_daily <- bind_rows(cdec %>% select(stn, date, wy, salinity, n_obs),
                       usgs %>% select(stn, date, wy, salinity, n_obs))

## ============================ 2. station coordinates =======================
STN <- tribble(
  ~stn,       ~lon,         ~lat,       ~label,                      ~src,
  "11455820", -122.22552,   38.06131,   "Carquinez Br (Crockett)",   "USGS",
  "MRZ",      -122.140488,  38.027637,  "Martinez",                  "CDEC",
  "11455780", -122.12663,   38.04492,   "Benicia Br",                "USGS",
  "GZL",      -122.038117,  38.124250,  "Grizzly Bay",               "CDEC",
  "HUN",      -122.05270,   38.15600,   "Hunter Cut",                "CDEC",
  "PCT",      -121.999971,  38.056788,  "Port Chicago",              "CDEC",
  "RYC",      -121.995878,  38.083971,  "Ryer cutoff",               "CDEC",
  "HON",      -121.939200,  38.072400,  "Honker Bay",                "CDEC",
  "MAL",      -121.920090,  38.042805,  "Mallard Is.",               "CDEC",
  "11455508", -121.88755,   38.04953,   "Van Sickle Is.",            "USGS",
  "CLL",      -121.850123,  38.073981,  "Collinsville",              "CDEC",
  "ANC",      -121.802965,  38.017817,  "Antioch",                   "CDEC",
  "SSI",      -121.761736,  38.074097,  "Sherman Is.",               "CDEC",
  "RIV",      -121.684543,  38.159203,  "Rio Vista",                 "CDEC"
)
# Suisun Marsh interior stations (BDL, MSL, GYS, CYG, TEA, IBS, FLT) are excluded:
# they sit behind the Suisun Marsh Salinity Control Gates, are not resolved by the
# TIGER hydrography, and do not represent the open-estuary gradient.

sal_daily <- sal_daily %>% filter(stn %in% STN$stn)
cat("stations with data:", paste(sort(unique(sal_daily$stn)), collapse = " "), "\n")
cat("\ndays per station per wy:\n")
print(sal_daily %>% count(stn, wy) %>% pivot_wider(names_from = wy, values_from = n, values_fill = 0) %>%
        left_join(STN %>% select(stn, label, lon), "stn") %>% arrange(lon) %>% as.data.frame())

## ============================ 3. the five dates ============================
WEST <- "11455820"                      # westernmost station with data in these years
if (!WEST %in% sal_daily$stn) stop("westernmost station has no data; cannot pick dates")

dates <- sal_daily %>% filter(stn == WEST) %>% group_by(wy) %>%
  slice_min(salinity, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(wy, date, west_sal = round(salinity, 2))
cat("\n--- minimum-salinity date at the westernmost station (", WEST, ") ---\n")
print(as.data.frame(dates))

## values at every station on those dates (allow +/-2 d gap-fill)
snap_val <- function(st, d) {
  x <- sal_daily %>% filter(stn == st, abs(as.numeric(date - d)) <= 2)
  if (!nrow(x)) return(NA_real_)
  x$salinity[which.min(abs(as.numeric(x$date - d)))]
}
field <- expand_grid(dates, STN %>% select(stn, lon, lat, label)) %>%
  mutate(salinity = map2_dbl(stn, date, snap_val)) %>% filter(!is.na(salinity))
cat("\n--- station salinity on each mapped date ---\n")
print(field %>% select(wy, date, stn, label, lon, salinity) %>%
        mutate(salinity = round(salinity, 2)) %>% arrange(wy, lon) %>% as.data.frame())

## ============================ 4. water mask ================================
UTM <- 32610; RES <- 100
hydro <- st_sf(geometry = readRDS("data/spatial/delta_hydro_tiger.rds"))
st_crs(hydro) <- 4326
h <- st_transform(hydro, UTM) |> st_make_valid()
INTERP <- st_bbox(c(xmin = -122.56, ymin = 37.95, xmax = -121.68, ymax = 38.31),
                  crs = st_crs(4326)) |> st_as_sfc() |> st_transform(UTM)
h <- suppressWarnings(st_crop(h, INTERP))

r0    <- rast(vect(INTERP), resolution = RES, crs = paste0("EPSG:", UTM))
water <- rasterize(vect(h), r0, field = 1, background = NA)
pch   <- patches(water, directions = 8, zeroAsNA = TRUE)
fr    <- freq(pch); fr <- fr[order(-fr$count), ]
main  <- ifel(pch == fr$value[1], 1, NA)
cat("\nmain water body cells:", global(main, "notNA")[1,1], "\n")

## ============================ 5. in-water distance per station =============
pts  <- st_as_sf(STN, coords = c("lon","lat"), crs = 4326) |> st_transform(UTM)
wxy  <- crds(as.points(main))
sxy  <- st_coordinates(pts)
snapxy <- t(sapply(seq_len(nrow(sxy)), function(i) {
  d <- (wxy[,1]-sxy[i,1])^2 + (wxy[,2]-sxy[i,2])^2
  wxy[which.min(d), ]
}))
STN$snap_x <- snapxy[,1]; STN$snap_y <- snapxy[,2]

cat("computing in-water distance surfaces...\n")
dist_list <- list()
for (i in seq_len(nrow(STN))) {
  # gridDist traverses non-NA, non-target cells and treats NA as a barrier:
  # water = 0 (traversable), land = NA (barrier), the station cell = 1 (target).
  seed <- ifel(is.na(main), NA, 0)
  cell <- cellFromXY(main, cbind(STN$snap_x[i], STN$snap_y[i]))
  seed[cell] <- 1
  dist_list[[STN$stn[i]]] <- gridDist(seed, target = 1) |> mask(main)
  cat("  ", STN$stn[i], "")
}
cat("\n")
dstack <- rast(dist_list)
names(dstack) <- STN$stn
saveRDS(list(dstack = wrap(dstack), main = wrap(main), STN = STN,
             dates = dates, field = field, sal_daily = sal_daily, UTM = UTM, RES = RES),
        file.path(DERIVED, "sal_field_inputs.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 3 - Figure S2, pooled symbols
# =============================================================================
# Figure: salinity field at the freshest moment of each water year, with catch overlaid.
#   version "clean" : heat map + the three catch-enclosing contours only
#   version "topo"  : heat map + full contour set, catch-enclosing ones in bold
suppressMessages({library(tidyverse); library(sf); library(terra); library(patchwork); library(scales)})

inp    <- readRDS(file.path(DERIVED, "sal_field_inputs.rds"))
dstack <- rast(inp$dstack); main <- rast(inp$main)
STN <- inp$STN; dates <- inp$dates; field <- inp$field

## ----------------------------- tows ---------------------------------------
rd <- function(f, lab) read_csv(f, show_col_types = FALSE) %>%
  transmute(survey = lab, wy = water_year,
            lon = Station_Longitude, lat = Station_Latitude,
            salinity = Salinity, catch = n_l_75)
tows <- bind_rows(rd("data/tow_12_6_tps_wq_nona.csv",   "TPS"),
                  rd("data/tow_12_6_edsm_wq_salfix.csv","EDSM"))
pos  <- tows %>% filter(catch > 0); zero <- tows %>% filter(catch == 0)
fish_sal  <- rep(pos$salinity, pos$catch)
pct_below <- function(x) 100 * mean(fish_sal <= x)

## catch-enclosing salinity thresholds (pooled across years, fish-weighted)
THR  <- c(0.1, 1, 6)
TLAB <- sprintf("%s psu – %.0f%% of catch", c("0.1", "1", "6"), sapply(THR, pct_below))
TCOL <- c("#F9A825", "#E65100", "#7B1FA2")
names(TCOL) <- as.character(THR); names(TLAB) <- as.character(THR)
cat("catch-enclosing contours:\n"); for (i in seq_along(THR)) cat("  ", TLAB[i], "\n")

## ----------------------------- fields --------------------------------------
IDW_POWER <- 2.5
make_field <- function(w) {
  f   <- field %>% filter(wy == w)
  d   <- dstack[[f$stn]]
  wts <- 1 / (d + 50)^IDW_POWER
  num <- wts
  for (i in seq_len(nlyr(wts))) num[[i]] <- wts[[i]] * f$salinity[i]
  (sum(num) / sum(wts)) |> mask(main) |> project("EPSG:4326", method = "bilinear")
}
flds <- setNames(lapply(dates$wy, make_field), dates$wy)

## ----------------------------- constants -----------------------------------
hydro <- st_sf(geometry = readRDS("data/spatial/delta_hydro_tiger.rds")); st_crs(hydro) <- 4326
XLIM <- c(-122.510, -121.750)        # TPS sampling extent + margin
YLIM <- c(  38.000,   38.245)
WEST_LIM <- -122.22552               # westernmost salinity station (Carquinez Br)
LAND <- "#e6e6e6"
TOPO_LV <- c(0.06, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.5, 2.5, 4, 8, 12, 16, 20)
BASE_PT <- 7

contour_sf <- function(r, levels) {
  mm <- minmax(r)
  levels <- levels[levels > mm[1] & levels < mm[2]]
  if (!length(levels)) return(NULL)
  cl <- suppressWarnings(as.contour(r, levels = levels))
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  st_as_sf(cl) %>% mutate(level = as.numeric(as.character(level)))
}

thm <- theme_void(base_size = BASE_PT, base_family = "Helvetica") +
  theme(plot.title = element_text(size = BASE_PT, face = "bold", hjust = 0,
                                  margin = margin(b = 1)),
        plot.margin = margin(1, 2, 1, 2),
        legend.position = "right",
        legend.title = element_text(size = BASE_PT),
        legend.text  = element_text(size = BASE_PT - 1),
        legend.key.height = unit(8, "pt"), legend.key.width = unit(8, "pt"),
        legend.margin = margin(0, 0, 0, 0))

panel <- function(w, version) {
  r   <- flds[[as.character(w)]]
  rdf <- as.data.frame(r, xy = TRUE, na.rm = TRUE) |> setNames(c("lon","lat","salinity"))
  d   <- dates %>% filter(wy == w)
  pz  <- zero %>% filter(wy == w); pp <- pos %>% filter(wy == w)
  thr <- contour_sf(r, THR)

  g <- ggplot() +
    annotate("rect", xmin = XLIM[1], xmax = XLIM[2], ymin = YLIM[1], ymax = YLIM[2],
             fill = LAND, colour = NA) +
    geom_sf(data = hydro, fill = "white", colour = NA) +
    geom_raster(data = rdf, aes(lon, lat, fill = salinity))

  if (version == "topo") {
    tp <- contour_sf(r, setdiff(TOPO_LV, THR))
    if (!is.null(tp)) g <- g +
      geom_sf(data = tp, colour = "white",  linewidth = 0.38, alpha = 0.60) +
      geom_sf(data = tp, colour = "grey10", linewidth = 0.15, alpha = 0.80)
  }
  if (!is.null(thr)) g <- g +
    geom_sf(data = thr, colour = "white", linewidth = 1.15, alpha = 0.80, show.legend = FALSE) +
    geom_sf(data = thr, aes(colour = factor(level)), linewidth = 0.62, show.legend = FALSE)

  g +
    # dummy layer so every panel contributes an identical contour legend
    geom_line(data = tibble(lon = NA_real_, lat = NA_real_, lvl = as.character(THR)),
              aes(lon, lat, colour = lvl), key_glyph = "path") +
    # west of the westernmost salinity station the field is extrapolated
    annotate("rect", xmin = XLIM[1], xmax = WEST_LIM, ymin = YLIM[1], ymax = YLIM[2],
             fill = "white", alpha = 0.42, colour = NA) +
    annotate("segment", x = WEST_LIM, xend = WEST_LIM, y = YLIM[1], yend = YLIM[2],
             linetype = "22", linewidth = 0.22, colour = "grey20") +
    geom_point(data = pz, aes(lon, lat), shape = 4, size = 0.34, stroke = 0.34,
               colour = "white", alpha = 0.55) +
    geom_point(data = pz, aes(lon, lat), shape = 4, size = 0.30, stroke = 0.15,
               colour = "grey15", alpha = 0.60) +
    geom_point(data = pp, aes(lon, lat, size = catch), shape = 21, fill = "#D81B60",
               colour = "black", stroke = 0.20, alpha = 0.95) +
    geom_point(data = STN, aes(lon, lat), shape = 24, size = 0.7,
               fill = "white", colour = "black", stroke = 0.20) +
    annotate("segment", x = XLIM[1] + 0.025, xend = XLIM[1] + 0.025 + 10/87.6,
             y = YLIM[1] + 0.018, yend = YLIM[1] + 0.018, linewidth = 0.5, colour = "grey10") +
    annotate("text", x = XLIM[1] + 0.025 + 5/87.6, y = YLIM[1] + 0.033,
             label = "10 km", size = (BASE_PT - 1.5)/.pt, colour = "grey10") +
    scale_fill_gradientn(
      colours = c("#f7fcf5","#c7e9c0","#74c476","#31a354","#1f6ea8","#2b3f72","#141f3d"),
      trans = "log10", limits = c(0.05, 20),
      breaks = c(0.05, 0.2, 1, 5, 20), labels = c("0.05","0.2","1","5","20"),
      oob = scales::squish, name = "Salinity (psu)") +
    scale_colour_manual(values = TCOL, labels = TLAB, limits = as.character(THR),
                        drop = FALSE, name = NULL) +
    scale_size_area(max_size = 2.6, breaks = c(1, 3, 10, 18),
                    limits = c(1, max(pos$catch)), name = "Fish per tow") +
    guides(size   = guide_legend(order = 1),
           colour = guide_legend(order = 2, override.aes = list(linewidth = 0.8)),
           fill   = guide_colourbar(order = 3, barheight = unit(40, "pt"),
                                    barwidth = unit(7, "pt"))) +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = sprintf("WY%d  —  %s   (max %.1f psu in view)",
                         w, format(d$date, "%d %b %Y"), max(rdf$salinity))) +
    thm
}


build <- function(version, stem) {
  ps  <- lapply(dates$wy, panel, version = version)
  fig <- (ps[[1]] | ps[[2]]) / (ps[[3]] | ps[[4]]) / (ps[[5]] | patchwork::guide_area()) +
    patchwork::plot_layout(guides = "collect")
  W <- 17.8/2.54; H <- 13.5/2.54
  base <- file.path(OUT_APPX, stem)
  open_pdf_device(paste0(base, ".pdf"), W, H); print(fig); dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = W, height = H, units = "in",
                 res = 600, compression = "lzw", background = "white"); print(fig); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = W, height = H, units = "in",
                res = 200, background = "white"); print(fig); dev.off()
  cat("wrote", stem, "\n")
}

build("clean", "sal_salinity_contour_catch_maps")
build("topo",  "sal_salinity_contour_catch_maps_topo")
cat("DONE\n")

# =============================================================================
# Part 4 - Figure S2, TPS and EDSM distinguished
# =============================================================================
# Figure: salinity field at the freshest moment of each water year, with catch
# overlaid, TPS and EDSM distinguished by symbol.
#   version "clean" : heat map + the three catch-enclosing contours only
#   version "topo"  : heat map + full contour set, catch-enclosing ones in bold
suppressMessages({library(tidyverse); library(sf); library(terra); library(patchwork); library(scales)})

inp    <- readRDS(file.path(DERIVED, "sal_field_inputs.rds"))
dstack <- rast(inp$dstack); main <- rast(inp$main)
STN <- inp$STN; dates <- inp$dates; field <- inp$field

## ----------------------------- tows ---------------------------------------
rd <- function(f, lab) read_csv(f, show_col_types = FALSE) %>%
  transmute(survey = lab, wy = water_year,
            lon = Station_Longitude, lat = Station_Latitude,
            salinity = Salinity, catch = n_l_75)
tows <- bind_rows(rd("data/tow_12_6_tps_wq_nona.csv",   "TPS"),
                  rd("data/tow_12_6_edsm_wq_salfix.csv","EDSM"))

## symbol groups: survey x (catch / no catch)
GRP <- c("TPS, catch", "EDSM, catch", "TPS, no catch", "EDSM, no catch")
GSH <- c("TPS, catch" = 23, "EDSM, catch" = 21, "TPS, no catch" = 3, "EDSM, no catch" = 4)
tows <- tows %>%
  mutate(grp = factor(paste0(survey, if_else(catch > 0, ", catch", ", no catch")),
                      levels = GRP))
pos  <- tows %>% filter(catch > 0)
zero <- tows %>% filter(catch == 0)

## catch-enclosing salinity thresholds (pooled across years AND surveys,
## fish-weighted).  Printed per survey so the reader can see the same
## thresholds hold in each survey separately.
THR <- c(0.1, 1, 6)
pct_of <- function(d, x) 100 * mean(rep(d$salinity, d$catch) <= x)
cat("=== % of fish at or below each contour ===\n")
print(map_dfr(c("TPS","EDSM","pooled"), function(s) {
  d <- if (s == "pooled") pos else pos %>% filter(survey == s)
  tibble(set = s, n_fish = sum(d$catch),
         !!!setNames(lapply(THR, function(x) round(pct_of(d, x), 1)), paste0("le_", THR)))
}) %>% as.data.frame())

TLAB <- sprintf("%s psu – %.0f%% of catch (TPS %.0f%%, EDSM %.0f%%)",
                c("0.1", "1", "6"),
                sapply(THR, function(x) pct_of(pos, x)),
                sapply(THR, function(x) pct_of(filter(pos, survey == "TPS"),  x)),
                sapply(THR, function(x) pct_of(filter(pos, survey == "EDSM"), x)))
TCOL <- c("#F9A825", "#E65100", "#7B1FA2")
names(TCOL) <- as.character(THR); names(TLAB) <- as.character(THR)

## ----------------------------- fields --------------------------------------
IDW_POWER <- 2.5
make_field <- function(w) {
  f   <- field %>% filter(wy == w)
  d   <- dstack[[f$stn]]
  wts <- 1 / (d + 50)^IDW_POWER
  num <- wts
  for (i in seq_len(nlyr(wts))) num[[i]] <- wts[[i]] * f$salinity[i]
  (sum(num) / sum(wts)) |> mask(main) |> project("EPSG:4326", method = "bilinear")
}
flds <- setNames(lapply(dates$wy, make_field), dates$wy)

## ----------------------------- constants -----------------------------------
hydro <- st_sf(geometry = readRDS("data/spatial/delta_hydro_tiger.rds")); st_crs(hydro) <- 4326
XLIM <- c(-122.510, -121.750)        # TPS sampling extent + margin
YLIM <- c(  38.000,   38.245)
WEST_LIM <- -122.22552               # westernmost salinity station (Carquinez Br)
LAND <- "#e6e6e6"
TOPO_LV <- c(0.06, 0.08, 0.12, 0.2, 0.3, 0.5, 0.8, 1.5, 2.5, 4, 8, 12, 16, 20)
BASE_PT <- 7

contour_sf <- function(r, levels) {
  mm <- minmax(r)
  levels <- levels[levels > mm[1] & levels < mm[2]]
  if (!length(levels)) return(NULL)
  cl <- suppressWarnings(as.contour(r, levels = levels))
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  st_as_sf(cl) %>% mutate(level = as.numeric(as.character(level)))
}

thm <- theme_void(base_size = BASE_PT, base_family = "Helvetica") +
  theme(plot.title = element_text(size = BASE_PT, face = "bold", hjust = 0,
                                  margin = margin(b = 1)),
        plot.margin = margin(1, 2, 1, 2),
        legend.position = "right",
        legend.title = element_text(size = BASE_PT),
        legend.text  = element_text(size = BASE_PT - 1),
        legend.key.height = unit(8, "pt"), legend.key.width = unit(8, "pt"),
        legend.margin = margin(0, 0, 0, 0))

panel <- function(w, version) {
  r   <- flds[[as.character(w)]]
  rdf <- as.data.frame(r, xy = TRUE, na.rm = TRUE) |> setNames(c("lon","lat","salinity"))
  d   <- dates %>% filter(wy == w)
  pz  <- zero %>% filter(wy == w); pp <- pos %>% filter(wy == w)
  thr <- contour_sf(r, THR)

  g <- ggplot() +
    annotate("rect", xmin = XLIM[1], xmax = XLIM[2], ymin = YLIM[1], ymax = YLIM[2],
             fill = LAND, colour = NA) +
    geom_sf(data = hydro, fill = "white", colour = NA) +
    geom_raster(data = rdf, aes(lon, lat, fill = salinity))

  if (version == "topo") {
    tp <- contour_sf(r, setdiff(TOPO_LV, THR))
    if (!is.null(tp)) g <- g +
      geom_sf(data = tp, colour = "white",  linewidth = 0.38, alpha = 0.60) +
      geom_sf(data = tp, colour = "grey10", linewidth = 0.15, alpha = 0.80)
  }
  if (!is.null(thr)) g <- g +
    geom_sf(data = thr, colour = "white", linewidth = 1.15, alpha = 0.80, show.legend = FALSE) +
    geom_sf(data = thr, aes(colour = factor(level)), linewidth = 0.62, show.legend = FALSE)

  g +
    # dummy layer so every panel contributes an identical contour legend
    geom_line(data = tibble(lon = NA_real_, lat = NA_real_, lvl = as.character(THR)),
              aes(lon, lat, colour = lvl), key_glyph = "path") +
    # west of the westernmost salinity station the field is extrapolated
    annotate("rect", xmin = XLIM[1], xmax = WEST_LIM, ymin = YLIM[1], ymax = YLIM[2],
             fill = "white", alpha = 0.42, colour = NA) +
    annotate("segment", x = WEST_LIM, xend = WEST_LIM, y = YLIM[1], yend = YLIM[2],
             linetype = "22", linewidth = 0.22, colour = "grey20") +
    # dummy layer so every panel contributes an identical symbol legend
    geom_point(data = tibble(lon = NA_real_, lat = NA_real_,
                             grp = factor(GRP, levels = GRP)),
               aes(lon, lat, shape = grp), fill = "#D81B60", colour = "black") +
    # zero-catch effort: + = TPS, x = EDSM, white halo so both read on any fill
    geom_point(data = pz, aes(lon, lat, shape = grp), size = 0.36, stroke = 0.36,
               colour = "white", alpha = 0.55, show.legend = FALSE) +
    geom_point(data = pz, aes(lon, lat, shape = grp), size = 0.32, stroke = 0.16,
               colour = "grey15", alpha = 0.62, show.legend = FALSE) +
    # positive catch: diamond = TPS, circle = EDSM
    geom_point(data = pp, aes(lon, lat, size = catch, shape = grp),
               fill = "#D81B60", colour = "black", stroke = 0.20, alpha = 0.95,
               show.legend = c(shape = FALSE, size = TRUE)) +
    geom_point(data = STN, aes(lon, lat), shape = 24, size = 0.7,
               fill = "white", colour = "black", stroke = 0.20,
               inherit.aes = FALSE) +
    annotate("segment", x = XLIM[1] + 0.025, xend = XLIM[1] + 0.025 + 10/87.6,
             y = YLIM[1] + 0.018, yend = YLIM[1] + 0.018, linewidth = 0.5, colour = "grey10") +
    annotate("text", x = XLIM[1] + 0.025 + 5/87.6, y = YLIM[1] + 0.033,
             label = "10 km", size = (BASE_PT - 1.5)/.pt, colour = "grey10") +
    scale_fill_gradientn(
      colours = c("#f7fcf5","#c7e9c0","#74c476","#31a354","#1f6ea8","#2b3f72","#141f3d"),
      trans = "log10", limits = c(0.05, 20),
      breaks = c(0.05, 0.2, 1, 5, 20), labels = c("0.05","0.2","1","5","20"),
      oob = scales::squish, name = "Salinity (psu)") +
    scale_colour_manual(values = TCOL, labels = TLAB, limits = as.character(THR),
                        drop = FALSE, name = NULL) +
    scale_shape_manual(values = GSH, limits = GRP, drop = FALSE, name = "Survey / tow") +
    scale_size_area(max_size = 2.6, breaks = c(1, 3, 10, 18),
                    limits = c(1, max(pos$catch)), name = "Fish per tow") +
    guides(shape  = guide_legend(order = 1, override.aes =
             list(size = c(1.9, 1.9, 1.4, 1.4), stroke = 0.3,
                  fill = c("#D81B60", "#D81B60", NA, NA),
                  colour = c("black", "black", "grey15", "grey15"), alpha = 1)),
           size   = guide_legend(order = 2, override.aes =
             list(shape = 21, fill = "#D81B60", colour = "black")),
           colour = guide_legend(order = 3, override.aes = list(linewidth = 0.8)),
           fill   = guide_colourbar(order = 4, barheight = unit(34, "pt"),
                                    barwidth = unit(7, "pt"))) +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    labs(title = sprintf("WY%d  —  %s   (max %.1f psu in view)",
                         w, format(d$date, "%d %b %Y"), max(rdf$salinity))) +
    thm
}


build <- function(version, stem) {
  ps  <- lapply(dates$wy, panel, version = version)
  fig <- (ps[[1]] | ps[[2]]) / (ps[[3]] | ps[[4]]) / (ps[[5]] | patchwork::guide_area()) +
    patchwork::plot_layout(guides = "collect")
  W <- 17.8/2.54; H <- 13.5/2.54
  base <- file.path(OUT_APPX, stem)
  open_pdf_device(paste0(base, ".pdf"), W, H); print(fig); dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = W, height = H, units = "in",
                 res = 600, compression = "lzw", background = "white"); print(fig); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = W, height = H, units = "in",
                res = 200, background = "white"); print(fig); dev.off()
  cat("wrote", stem, "\n")
}

build("clean", "sal_salinity_contour_catch_maps_bysurvey")
build("topo",  "sal_salinity_contour_catch_maps_bysurvey_topo")
cat("DONE\n")

# =============================================================================
# Part 5 - extent of the low-salinity envelope
# =============================================================================
# Numbers that accompany the salinity contour map.
suppressMessages({library(tidyverse); library(sf); library(terra)})

inp <- readRDS(file.path(DERIVED, "sal_field_inputs.rds"))
dstack <- rast(inp$dstack); main <- rast(inp$main)
STN <- inp$STN; dates <- inp$dates; field <- inp$field

IDW_POWER <- 2.5
make_field <- function(w) {
  f <- field %>% filter(wy == w); d <- dstack[[f$stn]]
  wts <- 1/(d + 50)^IDW_POWER; num <- wts
  for (i in seq_len(nlyr(wts))) num[[i]] <- wts[[i]] * f$salinity[i]
  (sum(num)/sum(wts)) |> mask(main)
}
flds <- setNames(lapply(dates$wy, make_field), dates$wy)

rd <- function(f, lab) read_csv(f, show_col_types = FALSE) %>%
  transmute(survey = lab, wy = water_year, lon = Station_Longitude,
            lat = Station_Latitude, catch = n_l_75)
tows <- bind_rows(rd("data/tow_12_6_tps_wq_nona.csv","TPS"),
                  rd("data/tow_12_6_edsm_wq_salfix.csv","EDSM"))
tpts <- st_as_sf(tows, coords = c("lon","lat"), crs = 4326) |> st_transform(inp$UTM)

CELL_KM2 <- (inp$RES/1000)^2
THR <- c(0.1, 1, 6)

res <- map_dfr(dates$wy, function(w) {
  r  <- flds[[as.character(w)]]
  tw <- tows %>% filter(wy == w)
  tp <- tpts[tpts$wy == w, ]
  sal_at_tow <- terra::extract(r, vect(tp))[, 2]
  map_dfr(THR, function(th) {
    inzone <- !is.na(r) & r <= th
    area   <- global(inzone, "sum", na.rm = TRUE)[1,1] * CELL_KM2
    # westernmost (most seaward) longitude of the <= th envelope
    pz <- as.points(ifel(inzone, 1, NA))
    wlon <- if (nrow(pz)) min(crds(project(pz, "EPSG:4326"))[,1]) else NA_real_
    inside <- !is.na(sal_at_tow) & sal_at_tow <= th
    tibble(wy = w, date = dates$date[dates$wy == w], threshold = th,
           area_km2 = round(area), west_lon = round(wlon, 3),
           pct_tows_inside = round(100*mean(inside), 1),
           pct_fish_inside = round(100*sum(tw$catch[inside])/max(sum(tw$catch),1), 1),
           fish_total = sum(tw$catch))
  })
})
cat("=== extent of the low-salinity envelope at the freshest moment of each water year ===\n")
print(as.data.frame(res))

cat("\n=== field salinity at the mapped date, summary ===\n")
print(map_dfr(dates$wy, function(w) {
  v <- values(flds[[as.character(w)]]); v <- v[!is.na(v)]
  tibble(wy = w, date = dates$date[dates$wy == w],
         min = round(min(v),2), median = round(median(v),2), max = round(max(v),2),
         pct_water_below_1 = round(100*mean(v <= 1), 1))
}) %>% as.data.frame())

saveRDS(res, file.path(DERIVED, "sal_map_stats.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 6 - share of fish at or below each threshold, by survey
# =============================================================================
suppressMessages(library(tidyverse))
rd <- function(f, lab) read_csv(f, show_col_types = FALSE) %>%
  transmute(survey = lab, wy = water_year, lon = Station_Longitude,
            lat = Station_Latitude, salinity = Salinity, catch = n_l_75)
tows <- bind_rows(rd("data/tow_12_6_tps_wq_nona.csv","TPS"),
                  rd("data/tow_12_6_edsm_wq_salfix.csv","EDSM"))
cat("=== tows / fish by survey ===\n")
print(tows %>% group_by(survey) %>%
  summarise(tows = n(), pos_tows = sum(catch > 0), fish = sum(catch),
            max_tow = max(catch), .groups="drop") %>% as.data.frame())

pos <- tows %>% filter(catch > 0)
pct <- function(d, x) 100 * mean(rep(d$salinity, d$catch) <= x)
THR <- c(0.1, 0.5, 1, 2, 6)
cat("\n=== % of FISH at or below each threshold ===\n")
tab <- map_dfr(c("TPS","EDSM","pooled"), function(s) {
  d <- if (s == "pooled") pos else pos %>% filter(survey == s)
  tibble(set = s, n_fish = sum(d$catch),
         !!!setNames(lapply(THR, function(x) round(pct(d, x), 1)), paste0("le_", THR)))
})
print(as.data.frame(tab))

cat("\n=== quantiles of salinity at catch (fish-weighted) ===\n")
print(map_dfr(c("TPS","EDSM","pooled"), function(s) {
  d <- if (s == "pooled") pos else pos %>% filter(survey == s)
  v <- rep(d$salinity, d$catch)
  tibble(set = s, n = length(v), med = round(median(v),3),
         p90 = round(quantile(v,.90),3), p95 = round(quantile(v,.95),3),
         p99 = round(quantile(v,.99),3), max = round(max(v),2))
}) %>% as.data.frame())

cat("\n=== positive tows west of Carquinez (-122.2255) ===\n")
print(tows %>% mutate(west = lon < -122.22552) %>% group_by(survey, west) %>%
  summarise(tows = n(), fish = sum(catch), .groups="drop") %>% as.data.frame())
