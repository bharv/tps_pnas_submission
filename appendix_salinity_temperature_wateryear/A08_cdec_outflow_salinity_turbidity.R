# =============================================================================
# A08_cdec_outflow_salinity_turbidity.R            (appendix Figure S1 + stats)
# -----------------------------------------------------------------------------
# Daily and 30-day Delta outflow (CDEC DTO), salinity at four stations from
# Collinsville to Martinez (PSS-78 from specific conductance at 25 C) and
# turbidity at three stations, Nov 1 - Jun 30 of WY2017-2021, one row per water
# year. Supports the appendix argument that salinity keeps a much longer memory of
# an outflow pulse than turbidity does, so the two respond to outflow on very
# different time scales.
#
# Parts
#   1. tidy the CDEC downloads           -> data/derived_appendix/cdec_tidy.rds
#   2. pulse-recovery statistics          -> data/derived_appendix/cdec_recovery.rds
#      (days for salinity / turbidity to recover half-way after each year's
#       largest outflow pulse)
#   3. memory statistics + the figure     -> data/derived_appendix/cdec_stats.rds,
#      output_figs_appendix/sal_cdec_outflow_salinity_turbidity_by_year.{pdf,tiff,png}
#      (best-fitting exponentially weighted outflow memory, between-year variance
#       share, e-folding decorrelation time). 17.8 cm wide, 7 pt type.
#
# Input : data/web_downloads/cdec_timeseries/*.csv  (A01 / 00_data/00_download_web_data.R)
# Cost  : < 1 min
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

# =============================================================================
# Part 1 - tidy the CDEC series
# =============================================================================
# Build the CDEC time-series table for the appendix figure:
# daily + 30-day Delta outflow, salinity (from specific conductance), turbidity,
# Nov 1 - Jun 30 of each water year 2017-2021.
suppressMessages({library(tidyverse); library(zoo); library(oce)})
CD <- file.path(WEB, "cdec_timeseries")

read_cdec <- function(file) {
  read_csv(file.path(CD, file), col_types = cols(.default = col_character())) %>%
    transmute(stn   = STATION_ID,
              dt    = as.POSIXct(`DATE TIME`, format = "%Y%m%d %H%M", tz = "Etc/GMT+8"),
              date  = as.Date(substr(`DATE TIME`, 1, 8), format = "%Y%m%d"),
              value = suppressWarnings(as.numeric(VALUE))) %>%
    filter(!is.na(date)) %>%
    # CDEC uses large negatives and absurd highs as missing flags
    mutate(value = ifelse(value < 0 | value > 1e6, NA, value))
}

CFS_TO_CMS <- 0.028316846592

## ---- Delta outflow (CDEC DTO, the Delta Outflow Index) ----------------------
out <- read_cdec("DTO_23_D.csv") %>%
  group_by(date) %>% summarise(out_cfs = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(date) %>%
  mutate(out_cms   = out_cfs * CFS_TO_CMS,
         out30_cms = zoo::rollmeanr(out_cms, 30, fill = NA, na.rm = TRUE))

## ---- Salinity from specific conductance -------------------------------------
# Sensor 324 is explicitly specific conductance @25 C; sensor 100 at these DWR Delta
# stations is likewise reported as specific conductance (uS/cm).  PSS-78 via oce::swSCTp
# at t = 25 C, p = 0, the same conversion used for the EDSM salinity field.
ec_files <- c(MAL = "MAL_100_D.csv", BDL = "BDL_324_D.csv",
              PCT = "PCT_100_D.csv", MRZ = "MRZ_100_D.csv", CLL = "CLL_100_D.csv",
              GZL = "GZL_100_E.csv")
sal <- map_dfr(names(ec_files), function(st) {
  read_cdec(ec_files[[st]]) %>%
    group_by(date) %>% summarise(ec = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(stn = st)
}) %>%
  filter(!is.na(ec), ec > 0) %>%
  mutate(salinity = oce::swSCTp(conductivity = ec / 1000, temperature = rep(25, length(ec)), pressure = rep(0, length(ec)),
                                conductivityUnit = "mS/cm", eos = "unesco"))

## ---- Turbidity ---------------------------------------------------------------
tb_files <- c(MAL = "MAL_27_H.csv", GZL = "GZL_27_E.csv", MRZ = "MRZ_27_E.csv")
turb <- map_dfr(names(tb_files), function(st) {
  read_cdec(tb_files[[st]]) %>%
    filter(value < 2000) %>%
    group_by(date) %>% summarise(turbidity = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(stn = st)
})

## ---- assemble, Nov 1 - Jun 30 of each water year -----------------------------
wy_of <- function(d) ifelse(as.integer(format(d, "%m")) >= 10,
                            as.integer(format(d, "%Y")) + 1L, as.integer(format(d, "%Y")))
in_window <- function(d) { m <- as.integer(format(d, "%m")); m >= 11 | m <= 6 }

tidy_all <- bind_rows(
  out  %>% transmute(date, stn = "DTO",  variable = "Outflow (daily)",   value = out_cms),
  out  %>% transmute(date, stn = "DTO",  variable = "Outflow (30-d mean)", value = out30_cms),
  sal  %>% transmute(date, stn, variable = "Salinity",  value = salinity),
  turb %>% transmute(date, stn, variable = "Turbidity", value = turbidity)
) %>%
  mutate(wy = wy_of(date)) %>%
  filter(wy %in% 2017:2021, in_window(date), !is.na(value)) %>%
  # common x axis: days since Nov 1 of the water year
  mutate(nov1 = as.Date(sprintf("%d-11-01", wy - 1L)),
         doy  = as.numeric(date - nov1))

saveRDS(tidy_all, file.path(DERIVED, "cdec_tidy.rds"))

cat("\n--- coverage: n days by variable/station/water year ---\n")
print(tidy_all %>% count(variable, stn, wy) %>%
        pivot_wider(names_from = wy, values_from = n, values_fill = 0) %>%
        arrange(variable, stn), n = 40)

cat("\n--- salinity range by station/wy (psu) ---\n")
print(tidy_all %>% filter(variable == "Salinity") %>% group_by(stn, wy) %>%
        summarise(min = round(min(value), 2), med = round(median(value), 2),
                  max = round(max(value), 2), .groups = "drop"), n = 40)

cat("\n--- turbidity range by station/wy (NTU) ---\n")
print(tidy_all %>% filter(variable == "Turbidity") %>% group_by(stn, wy) %>%
        summarise(min = round(min(value), 1), med = round(median(value), 1),
                  max = round(max(value), 1), .groups = "drop"), n = 40)

cat("\n--- peak 30-d outflow by wy (cms) ---\n")
print(tidy_all %>% filter(variable == "Outflow (30-d mean)") %>% group_by(wy) %>%
        summarise(mean = round(mean(value)), peak = round(max(value)), .groups = "drop"))
cat("\nDONE\n")

# =============================================================================
# Part 2 - recovery after the largest outflow pulse of each year
# =============================================================================
# Event-based recovery: after the largest outflow pulse of each water year, how many
# days does each variable take to come half-way back to its pre-pulse level?
suppressMessages({library(tidyverse)})
tidy_all <- readRDS(file.path(DERIVED, "cdec_tidy.rds"))

out_d <- tidy_all %>% filter(variable == "Outflow (daily)") %>% distinct(wy, doy, date, value)

peaks <- out_d %>% group_by(wy) %>% slice_max(value, n = 1, with_ties = FALSE) %>%
  ungroup() %>% transmute(wy, peak_doy = doy, peak_date = date, peak_cms = round(value))
cat("=== largest daily outflow pulse of each water year ===\n"); print(as.data.frame(peaks))

series <- tidy_all %>% filter(variable %in% c("Salinity", "Turbidity")) %>%
  select(variable, stn, wy, doy, value) %>% left_join(peaks, by = "wy")

recov <- series %>% group_by(variable, stn, wy) %>% group_modify(function(d, k) {
  pk   <- d$peak_doy[1]
  pre  <- d %>% filter(doy >= pk - 45, doy <= pk - 10) %>% pull(value)
  post <- d %>% filter(doy >= pk,      doy <= pk + 25)
  if (length(pre) < 15 || nrow(post) < 15) return(tibble())
  base <- median(pre, na.rm = TRUE)
  # salinity is pushed DOWN by the pulse, turbidity is pushed UP
  resp <- if (k$variable == "Salinity") min(post$value, na.rm = TRUE) else max(post$value, na.rm = TRUE)
  resp_doy <- post$doy[which(post$value == resp)[1]]
  half <- resp + 0.5 * (base - resp)
  tail <- d %>% filter(doy > resp_doy) %>% arrange(doy)
  hit  <- if (k$variable == "Salinity") which(tail$value >= half) else which(tail$value <= half)
  censored <- length(hit) == 0
  days <- if (censored) max(tail$doy, na.rm = TRUE) - resp_doy else tail$doy[hit[1]] - resp_doy
  tibble(baseline = round(base, 2), response = round(resp, 2),
         pct_change = round(100 * (resp - base) / base),
         days_to_half = days, censored = censored)
}) %>% ungroup()

cat("\n=== days to half-recovery after the peak (censored = never recovered in the season) ===\n")
print(as.data.frame(recov %>% mutate(days = ifelse(censored, paste0(">", days_to_half),
                                                   as.character(days_to_half))) %>%
                      select(variable, stn, wy, baseline, response, pct_change, days)))

cat("\n=== summary by variable (censored counted at their lower bound) ===\n")
print(as.data.frame(recov %>% group_by(variable) %>%
  summarise(n = n(), n_censored = sum(censored),
            median_days = median(days_to_half),
            mean_days = round(mean(days_to_half), 1),
            median_pct_change = median(pct_change))))

cat("\n=== response magnitude: % change from pre-pulse baseline ===\n")
print(as.data.frame(recov %>% group_by(variable) %>%
  summarise(median_pct = median(pct_change), min_pct = min(pct_change), max_pct = max(pct_change))))

## season medians by year, the between-year contrast seen in the figure
cat("\n=== season (Nov-Jun) median by station and water year ===\n")
print(tidy_all %>% filter(variable %in% c("Salinity","Turbidity")) %>%
  group_by(variable, stn, wy) %>% summarise(med = round(median(value), 1), .groups="drop") %>%
  pivot_wider(names_from = wy, values_from = med) %>%
  mutate(fold_range = round(pmax(`2017`,`2018`,`2019`,`2020`,`2021`) /
                            pmin(`2017`,`2018`,`2019`,`2020`,`2021`), 1)) %>%
  as.data.frame())

saveRDS(recov, file.path(DERIVED, "cdec_recovery.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 3 - outflow-memory statistics and Figure S1
# =============================================================================
# Appendix figure: Delta outflow, salinity and turbidity, Nov-Jun, WY2017-2021.
# Supports the claim that salinity retains a much longer memory of an outflow pulse
# than turbidity does.
suppressMessages({library(tidyverse); library(patchwork); library(zoo)})

tidy_all <- readRDS(file.path(DERIVED, "cdec_tidy.rds"))

SAL_STN  <- c("CLL", "BDL", "PCT", "MRZ")   # upstream -> downstream
TURB_STN <- c("MAL", "GZL", "MRZ")
STN_LAB <- c(CLL = "Collinsville (CLL)", BDL = "Montezuma Sl. (BDL)",
             PCT = "Port Chicago (PCT)", MRZ = "Martinez (MRZ)",
             MAL = "Mallard Is. (MAL)",  GZL = "Grizzly Bay (GZL)")

COL <- c(`Outflow (daily)` = "#9FC5E8", `Outflow (30-d mean)` = "#0B3D91",
         Salinity = "#C2452D", Turbidity = "#6B7D2E")

## ===========================================================================
## 1. Outflow-memory timescale: which exponentially-weighted outflow average
##    best predicts each series?  The best tau is the "memory" of the variable.
## ===========================================================================
daily_out <- tidy_all %>% filter(variable == "Outflow (daily)") %>%
  distinct(date, value) %>% arrange(date) %>% rename(out = value)
# a continuous daily outflow series (needed for the EWMA, incl. Jul-Oct gaps)
full_out <- tibble(date = seq(min(daily_out$date) - 150, max(daily_out$date), by = "day")) %>%
  left_join(daily_out, by = "date") %>%
  mutate(out = zoo::na.approx(out, na.rm = FALSE)) %>%
  fill(out, .direction = "downup")

ewma <- function(x, tau) {
  a <- 1 - exp(-1 / tau); y <- numeric(length(x)); y[1] <- x[1]
  for (i in 2:length(x)) y[i] <- a * x[i] + (1 - a) * y[i - 1]
  y
}
taus <- c(1, 2, 3, 5, 7, 10, 14, 20, 30, 45, 60, 90, 120, 180)
ew <- map_dfc(taus, ~ tibble(!!sprintf("tau%03d", .x) := ewma(full_out$out, .x)))
ew$date <- full_out$date

memory_tbl <- tidy_all %>%
  filter(variable %in% c("Salinity", "Turbidity")) %>%
  inner_join(ew, by = "date") %>%
  group_by(variable, stn) %>%
  group_modify(function(d, k) {
    r2 <- map_dbl(taus, function(t) {
      x <- log(d[[sprintf("tau%03d", t)]]); y <- log(pmax(d$value, 0.01))
      suppressWarnings(cor(x, y, use = "complete.obs"))^2
    })
    tibble(tau = taus, r2 = r2)
  }) %>% ungroup()

best_tau <- memory_tbl %>% group_by(variable, stn) %>%
  slice_max(r2, n = 1) %>% ungroup() %>% arrange(variable, stn)
cat("\n=== Best-fitting outflow memory timescale (tau, days) ===\n")
print(as.data.frame(best_tau %>% mutate(r2 = round(r2, 3))))

cat("\n=== R2 of log(variable) on log(EWMA outflow) at each tau ===\n")
print(memory_tbl %>% mutate(r2 = round(r2, 3)) %>%
        pivot_wider(names_from = tau, values_from = r2) %>% as.data.frame())

## between- vs within-year variance share, same two variables
vc <- tidy_all %>% filter(variable %in% c("Salinity", "Turbidity")) %>%
  group_by(variable, stn) %>%
  summarise(between_r2 = summary(lm(log(pmax(value, .01)) ~ factor(wy)))$r.squared,
            .groups = "drop") %>% mutate(between_r2 = round(between_r2, 3))
cat("\n=== Share of variance in log(variable) explained by water year ===\n")
print(as.data.frame(vc))

## ---------------------------------------------------------------------------
## 1b. Decorrelation timescale: how long does each variable remember ITS OWN
##     state?  e-folding time = first lag at which the within-year ACF < 1/e.
##     This is the direct measure of "salinity lingers, turbidity clears".
## ---------------------------------------------------------------------------
efold <- function(x, maxlag = 90) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 60) return(NA_real_)
  a <- acf(x, lag.max = maxlag, plot = FALSE, na.action = na.pass)$acf[-1]
  w <- which(a < exp(-1))
  if (!length(w)) return(NA_real_) else return(as.numeric(w[1]))
}

dec_tbl <- tidy_all %>%
  filter(variable %in% c("Salinity", "Turbidity")) %>%
  arrange(variable, stn, wy, doy) %>%
  group_by(variable, stn, wy) %>%
  summarise(efold_days = efold(value), .groups = "drop")

cat("\n=== e-folding decorrelation time of the daily series (days) ===\n")
print(dec_tbl %>% pivot_wider(names_from = wy, values_from = efold_days) %>% as.data.frame())
cat("\n=== median e-folding time across stations x years ===\n")
print(dec_tbl %>% group_by(variable) %>%
        summarise(median_days = median(efold_days, na.rm = TRUE),
                  min = min(efold_days, na.rm = TRUE),
                  max = max(efold_days, na.rm = TRUE)) %>% as.data.frame())

# same thing at the two stations that carry BOTH variables' nearest pairing
cat("\n=== paired at Martinez/Carquinez (MRZ), both variables, by year ===\n")
print(dec_tbl %>% filter(stn == "MRZ") %>%
        pivot_wider(names_from = variable, values_from = efold_days) %>% as.data.frame())

## ===========================================================================
## 2. Figure
## ===========================================================================
BASE <- 7
thm <- theme_bw(base_size = BASE, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.15, colour = "grey88"),
        panel.border      = element_rect(linewidth = 0.25, colour = "grey25"),
        axis.ticks        = element_line(linewidth = 0.25),
        axis.text         = element_text(size = BASE - 1),
        axis.title        = element_text(size = BASE),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        strip.text        = element_text(size = BASE - 0.5, margin = margin(1.2, 1.2, 1.2, 1.2)),
        plot.title        = element_text(size = BASE + 0.5, face = "bold",
                                         margin = margin(b = 1.5)),
        legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = BASE - 1),
        legend.key.height = unit(6, "pt"),
        legend.key.width  = unit(16, "pt"),
        legend.margin     = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin       = margin(2, 3, 1, 2))

xbreaks <- as.numeric(as.Date(c("2016-11-01","2016-12-01","2017-01-01","2017-02-01",
                                "2017-03-01","2017-04-01","2017-05-01","2017-06-01")) -
                        as.Date("2016-11-01"))
xlabs <- c("N","D","J","F","M","A","M","J")
xsc <- scale_x_continuous(breaks = xbreaks, labels = xlabs, limits = c(0, 242),
                          expand = expansion(mult = 0.01))

wy_lab <- function(x) paste("WY", x)

## -- panel A: outflow -------------------------------------------------------
pa_dat <- tidy_all %>% filter(str_starts(variable, "Outflow")) %>%
  mutate(variable = factor(variable, levels = c("Outflow (daily)", "Outflow (30-d mean)")))

pA <- ggplot(pa_dat, aes(doy, value, colour = variable, linewidth = variable)) +
  geom_line() +
  facet_wrap(~ wy, ncol = 1, strip.position = "right", labeller = labeller(wy = wy_lab)) +
  scale_colour_manual(values = COL[1:2]) +
  scale_linewidth_manual(values = c(0.22, 0.55), guide = "none") +
  scale_y_log10(breaks = c(100, 300, 1000, 3000),
                labels = c("100","300","1,000","3,000")) +
  xsc + thm +
  labs(title = "a  Delta outflow", x = "Month", y = expression(Outflow~(m^3~s^-1))) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 0.6)))

## -- panel B: salinity ------------------------------------------------------
pb_dat <- tidy_all %>% filter(variable == "Salinity", stn %in% SAL_STN) %>%
  mutate(stn = factor(stn, levels = SAL_STN, labels = STN_LAB[SAL_STN]))

pB <- ggplot(pb_dat, aes(doy, value, linetype = stn)) +
  geom_line(colour = COL[["Salinity"]], linewidth = 0.3) +
  facet_wrap(~ wy, ncol = 1, strip.position = "right", labeller = labeller(wy = wy_lab)) +
  scale_linetype_manual(values = c("solid", "22", "4121", "11")) +
  scale_y_continuous(limits = c(0, 21), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xsc + thm +
  labs(title = "b  Salinity", x = "Month", y = "Salinity (psu)") +
  guides(linetype = guide_legend(ncol = 1,
                                 override.aes = list(colour = COL[["Salinity"]], linewidth = 0.4)))

## -- panel C: turbidity -----------------------------------------------------
pc_dat <- tidy_all %>% filter(variable == "Turbidity", stn %in% TURB_STN) %>%
  mutate(stn = factor(stn, levels = TURB_STN, labels = STN_LAB[TURB_STN]))

pC <- ggplot(pc_dat, aes(doy, value, linetype = stn)) +
  geom_line(colour = COL[["Turbidity"]], linewidth = 0.3) +
  facet_wrap(~ wy, ncol = 1, strip.position = "right", labeller = labeller(wy = wy_lab)) +
  scale_linetype_manual(values = c("solid", "22", "4121")) +
  scale_y_continuous(limits = c(0, 260), breaks = seq(0, 250, 50),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xsc + thm +
  labs(title = "c  Turbidity", x = "Month", y = "Turbidity (NTU)") +
  guides(linetype = guide_legend(ncol = 1,
                                 override.aes = list(colour = COL[["Turbidity"]], linewidth = 0.4)))

fig <- (pA | pB | pC) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

W <- 17.8 / 2.54; H <- 19.0 / 2.54
base <- file.path(OUT_APPX, "sal_cdec_outflow_salinity_turbidity_by_year")

open_pdf_device(paste0(base, ".pdf"), W, H)
print(fig); dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = W, height = H, units = "in",
               res = 600, compression = "lzw", background = "white")
print(fig); dev.off()

ragg::agg_png(paste0(base, ".png"), width = W, height = H, units = "in",
              res = 170, background = "white")
print(fig); dev.off()

cat("\nwrote:\n"); print(list.files(OUT_APPX, pattern = "sal_cdec", full.names = TRUE))
saveRDS(list(best_tau = best_tau, memory_tbl = memory_tbl, vc = vc, dec_tbl = dec_tbl),
        file.path(DERIVED, "cdec_stats.rds"))
cat("\nDONE\n")
