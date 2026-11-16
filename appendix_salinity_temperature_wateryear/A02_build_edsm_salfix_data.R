# =============================================================================
# A02_build_edsm_salfix_data.R
# -----------------------------------------------------------------------------
# Builds data/tow_12_6_edsm_wq_salfix.csv: the EDSM model frame with a CORRECTED
# salinity column and a temperature z-score, for appendix EDSM fits 22-27.
#
# Why: EDSM reports specific conductivity (uS/cm). The original frame
# (00_data/01_build_model_datasets.Rmd) converted it to salinity as if it were
# mS/cm, so its Salinity / salinity_z are wrong. No main-text model uses salinity,
# so the published frame is left untouched; this script writes a separate file.
#
# What changes vs data/tow_12_6_edsm_wq_nona.csv (and nothing else):
#   Salinity       recomputed with PSS-78 (oce::swSCTp) from raw conductivity
#                  rejoined on SampleDate|StationCode|TowNumber, uS/cm -> mS/cm;
#                  3 physically impossible readings (> 35 psu: sensor failures at
#                  61,108 / 121,319 / 189,013 uS/cm) winsorized to the highest
#                  plausible value. Rows are KEPT so the positional K-fold vector
#                  still lines up.
#   salinity_z     re-standardised from the corrected Salinity
#   temperature_z  new (standardised Temperature)
#   + audit columns Salinity_broken, Salinity_raw, salinity_winsorized,
#                   salinity_z_broken, Specific_Conductivity, Temperature_src
# Every other modelling column is asserted byte-identical, and the row order is
# asserted to match the cached station fold vector, so fits 22-27 stay
# ELPD-comparable with fits 0-21.
#
# Inputs : data/tow_12_6_edsm_wq_nona.csv, data/All_Surveys_Combined_v7.3.csv,
#          fits_edsm_Vfix/kfold5_station_folds.rds
# Output : data/tow_12_6_edsm_wq_salfix.csv   (refuses to overwrite)
# Cost   : seconds
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")
suppressMessages({library(dplyr); library(readr); library(oce)})

OUT <- Sys.getenv("SALFIX_OUT", EDSM_SALFIX)   # env var only used to test-build elsewhere
if (file.exists(OUT)) stop(OUT, " already exists; delete it first to rebuild.")

e0 <- read_csv("data/tow_12_6_edsm_wq_nona.csv", show_col_types = FALSE)
e  <- e0 %>% select(-any_of("...1"))

# --- 1. rejoin raw specific conductivity on a real key ----------------------------
src <- read_csv("data/All_Surveys_Combined_v7.3.csv", show_col_types = FALSE) %>%
  filter(Survey == "EDSM") %>%
  select(SampleDate, StationCode, TowNumber, Specific_Conductivity,
         Temperature_src = Temperature)
k <- function(d) paste(d$SampleDate, d$StationCode, d$TowNumber, sep = "|")
stopifnot(!anyDuplicated(k(src)), !anyDuplicated(k(e)), all(k(e) %in% k(src)))
e <- e %>% left_join(src, by = c("SampleDate", "StationCode", "TowNumber"))
stopifnot(nrow(e) == 4241,
          isTRUE(all.equal(e$Temperature, e$Temperature_src)),   # join landed on the right rows
          !any(is.na(e$Specific_Conductivity)))

# --- 2. recompute salinity (uS/cm -> mS/cm, PSS-78, surface pressure) ---------------
e <- e %>% mutate(
  Salinity_broken = Salinity,
  Salinity_raw    = swSCTp(c = Specific_Conductivity / 1000, t = Temperature,
                           p = rep(0, n()), conductivityUnit = "mS/cm", eos = "unesco"))

# --- 3. winsorize the physically impossible values -----------------------------------
BAD <- e$Salinity_raw > 35
CAP <- max(e$Salinity_raw[!BAD])
message("winsorizing ", sum(BAD), " values > 35 psu to ", round(CAP, 4))
e <- e %>% mutate(Salinity = pmin(Salinity_raw, CAP), salinity_winsorized = BAD)

# --- 4. standardise ONLY salinity (replaced) and temperature (new) -------------------
e <- e %>% mutate(salinity_z_broken = salinity_z,
                  salinity_z    = as.numeric(scale(Salinity)),
                  temperature_z = as.numeric(scale(Temperature)))

# --- 5. nothing else moved; row order still matches the fold vector ------------------
keep <- c("out30_z", "distance_m_z", "distance_along_channel_z", "secchi_z", "dox_z",
          "water_year_day", "Volume", "n_l_75", "StationCode")
for (cc in keep) stopifnot(identical(e[[cc]], e0[[cc]]))
folds <- readRDS(file.path(FIT_EDSM, "kfold5_station_folds.rds"))
stopifnot(length(folds) == nrow(e),
          all(tapply(folds, e$StationCode, function(x) length(unique(x))) == 1))

write_csv(e, OUT)
message("wrote ", OUT, " (", nrow(e), " rows, ", ncol(e), " cols)")
