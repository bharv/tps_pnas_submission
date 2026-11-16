# =============================================================================
# A10_appendix_tables_and_stats.R
# -----------------------------------------------------------------------------
# Numbers quoted in the appendix tables and text.
#   Part 1
#     (a) paired ELPD contrasts on the EDSM grouped 5-fold CV: what salinity,
#         Secchi and temperature each add to the fit-9 backbone; the "mirror"
#         (salinity added to a Secchi model vs Secchi added to a salinity model);
#         what distance along channel adds in each covariate context; ranking of 27
#     (b) outflow by water year in m3/s: tow-weighted 30-day mean vs plain
#         Nov-Jun seasonal mean (needs data/derived_appendix/cdec_tidy.rds from A08)
#     (c) stacking weight carried by models CONTAINING each covariate
#   Part 2  model-averaged ("inclusion": absent models contribute 0) and
#     conditional (only models containing the term, weights renormalised) stacked
#     coefficient summaries for EDSM with salinity (27 models), EDSM without
#     salinity (22) and TPS (17), from 40,000 mixture draws
#   Part 3  sensitivity of the EDSM stacked coefficients to which salinity /
#     temperature models are in the set
#
# Inputs : fits_edsm_Vfix/kfold5_station_criteria_bundle_27.rds (A05),
#          cache/sal_fast_figure_cache_kfold.rds (A07), data/derived_appendix/cdec_tidy.rds (A08),
#          data/tow_12_6_tps_wq_nona.csv, data/tow_12_6_edsm_wq_salfix.csv
# Outputs: data/derived_appendix/appendix_stats.rds, table_s10.rds (+ printed tables)
# Cost   : ~1 min
# =============================================================================

source("appendix_salinity_temperature_wateryear/A00_appendix_setup.R")

# =============================================================================
# Part 1 - ELPD contrasts, outflow by year, weight on models containing each covariate
# =============================================================================
# Recomputations for the appendix revision:
#  (a) extra paired ELPD contrasts for Table S6 (Secchi analogue, the mirror, and dac)
#  (b) unit conversions + unambiguous seasonal outflow for Table S1
#  (c) revised Table S7 with weight-on-including-models and conditional intervals
suppressMessages({library(tidyverse); library(loo)})
CFS <- 0.028316846592

## =========================================================================
## (a) paired ELPD contrasts, EDSM grouped 5-fold CV
## =========================================================================
b <- readRDS("fits_edsm_Vfix/kfold5_station_criteria_bundle_27.rds")
pw <- function(n) b[[paste0("fit_edsm_", n)]]$pointwise[, "elpd_kfold"]
elpd <- function(n) sum(pw(n))

pair <- function(a, bnum, label) {
  d <- pw(a) - pw(bnum)
  se <- sqrt(length(d)) * sd(d)
  tibble(comparison = label, A = a, B = bnum,
         elpd_A = round(elpd(a), 2), elpd_B = round(elpd(bnum), 2),
         dELPD = round(sum(d), 2), SE = round(se, 2), ratio = round(sum(d)/se, 1))
}

cat("=== EDSM paired ELPD contrasts (grouped 5-fold CV by StationCode) ===\n\n")

cat("-- What each water-quality covariate buys on the SAME fit-9 backbone --\n")
t1 <- bind_rows(
  pair(22,  9, "fit 22 vs 9: + salinity only"),
  pair(10,  9, "fit 10 vs 9: + Secchi only"),
  pair(23,  9, "fit 23 vs 9: + temperature only"),
  pair(26,  9, "fit 26 vs 9: + Secchi AND salinity")
); print(as.data.frame(t1))

cat("\n-- THE MIRROR: salinity added to a Secchi model vs Secchi added to a salinity model --\n")
t2 <- bind_rows(
  pair(26, 10, "fit 26 vs 10: + salinity, Secchi already present"),
  pair(26, 22, "fit 26 vs 22: + Secchi, salinity already present"),
  pair(25, 21, "fit 25 vs 21: + salinity, Secchi already present"),
  pair(24, 19, "fit 24 vs 19: + salinity, Secchi already present")
); print(as.data.frame(t2))

cat("\n-- What distance downstream (dac) buys, in each covariate context --\n")
t3 <- bind_rows(
  pair( 2,  1, "fit  2 vs 1 : + dac, outflow only"),
  pair(19, 17, "fit 19 vs 17: + dac, outflow + Secchi"),
  pair( 9,  6, "fit  9 vs 6 : + dac (interacted), outflow x dist offshore"),
  pair(10, 20, "fit 10 vs 20: + dac (interacted), outflow x dist offshore + Secchi")
); print(as.data.frame(t3))

cat("\n-- model ranking, all 27 --\n")
rk <- tibble(fit = 0:26, elpd = map_dbl(0:26, elpd)) %>% arrange(desc(elpd)) %>%
  mutate(d_from_best = round(elpd - max(elpd), 2), elpd = round(elpd, 2))
print(as.data.frame(head(rk, 12)))

## =========================================================================
## (b) outflow in m3/s, and an unambiguous seasonal mean
## =========================================================================
cat("\n\n=== Outflow: tow-weighted 30-d mean vs plain seasonal mean, by water year ===\n")
tidy_all <- readRDS(file.path(DERIVED, "cdec_tidy.rds"))
seasonal <- tidy_all %>% filter(variable == "Outflow (daily)") %>% distinct(wy, date, value) %>%
  group_by(wy) %>% summarise(season_mean_cms = round(mean(value)),
                             season_peak_cms = round(max(value)), .groups = "drop")

tow_sum <- function(path, lab) {
  read_csv(path, show_col_types = FALSE) %>%
    group_by(wy = water_year) %>%
    summarise(tows = n(), fish = sum(n_l_75), pct_pos = round(100*mean(n_l_75 > 0), 2),
              tow30_cfs = round(mean(OUT_Moving_Mean_30d)),
              tow30_cms = round(mean(OUT_Moving_Mean_30d) * CFS),
              mean_sal  = round(mean(Salinity), 2),
              sd_sal    = round(sd(Salinity), 2),
              iqr_sal   = round(IQR(Salinity), 2),
              mean_secchi = round(mean(Secchi), 1), .groups = "drop") %>%
    left_join(seasonal, by = "wy") %>% mutate(survey = lab)
}
tps_tbl  <- tow_sum("data/tow_12_6_tps_wq_nona.csv", "TPS")
edsm_tbl <- tow_sum("data/tow_12_6_edsm_wq_salfix.csv", "EDSM")
cat("\n-- EDSM --\n"); print(as.data.frame(edsm_tbl %>% select(-survey)))
cat("\n-- TPS --\n");  print(as.data.frame(tps_tbl  %>% select(-survey)))

cat("\n=== key cfs -> m3/s conversions used in the text ===\n")
key <- c(183153, 12713, 50237, 11948, 8040, 21340, 80055, 12462, 8702,
         212416, 35051, 98046, 129675)
print(data.frame(cfs = key, cms = round(key * CFS)))

## =========================================================================
## (c) Table S7: model-averaged vs conditional stacked coefficients
## =========================================================================
cat("\n\n=== Stacking weights over the 27-model EDSM set ===\n")
lpd <- sapply(0:26, pw)                       # n x 27 pointwise matrix
colnames(lpd) <- paste0("fit", 0:26)
w <- loo::stacking_weights(lpd)
wtab <- tibble(fit = 0:26, weight = as.numeric(w)) %>% arrange(desc(weight))
print(as.data.frame(wtab %>% filter(weight > 1e-4) %>% mutate(weight = round(weight, 4))))

fixef_by_fit <- c(
  "0"="", "1"="out30_z","2"="out30_z + distance_along_channel_z","3"="out30_z + distance_m_z",
  "4"="out30_z + distance_along_channel_z + distance_m_z","5"="out30_z * distance_along_channel_z",
  "6"="out30_z * distance_m_z","7"="out30_z * distance_along_channel_z + distance_m_z",
  "8"="out30_z * distance_m_z + distance_along_channel_z",
  "9"="out30_z * (distance_m_z + distance_along_channel_z)",
  "10"="out30_z * (distance_m_z + distance_along_channel_z) + secchi_z",
  "11"="out30_z * (distance_m_z + distance_along_channel_z) + habitat2020 + secchi_z",
  "12"="out30_z * (distance_m_z + distance_along_channel_z) + habitat2020 + dox_z",
  "13"="out30_z * (distance_m_z + distance_along_channel_z) + habitat2020 + secchi_z + dox_z",
  "14"="out30_z * (distance_m_z + distance_along_channel_z) + habitat2020",
  "15"="out30_z * (distance_m_z + distance_along_channel_z) + dox_z",
  "16"="out30_z * (distance_m_z + distance_along_channel_z) + secchi_z + dox_z",
  "17"="out30_z + secchi_z","18"="out30_z + distance_m_z + secchi_z",
  "19"="out30_z + distance_along_channel_z + secchi_z","20"="out30_z * distance_m_z + secchi_z",
  "21"="out30_z * distance_along_channel_z + secchi_z",
  "22"="out30_z * (distance_m_z + distance_along_channel_z) + salinity_z",
  "23"="out30_z * (distance_m_z + distance_along_channel_z) + temperature_z",
  "24"="out30_z + distance_along_channel_z + secchi_z + salinity_z",
  "25"="out30_z * distance_along_channel_z + secchi_z + salinity_z",
  "26"="out30_z * (distance_m_z + distance_along_channel_z) + secchi_z + salinity_z")

has <- function(fit, term) grepl(term, fixef_by_fit[[as.character(fit)]], fixed = TRUE)
cat("\n=== total stacking weight on models CONTAINING each covariate ===\n")
for (tm in c("out30_z","distance_along_channel_z","distance_m_z","secchi_z","salinity_z","temperature_z")) {
  keep <- sapply(0:26, has, term = tm)
  cat(sprintf("  %-26s weight = %.4f  (%d of 27 models)\n", tm, sum(as.numeric(w)[keep]), sum(keep)))
}
saveRDS(list(t1=t1,t2=t2,t3=t3,rk=rk,w=wtab,edsm_tbl=edsm_tbl,tps_tbl=tps_tbl,seasonal=seasonal),
        file.path(DERIVED, "appendix_stats.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 2 - stacked coefficient table (inclusion vs conditional mixtures)
# =============================================================================
# Revised stacked-coefficient table:
#   (i) model-averaged ("inclusion") mixture - models lacking the covariate contribute 0
#  (ii) conditional mixture - only models containing the covariate, weights renormalised
# (iii) the stacking weight carried by models that contain the covariate at all
suppressMessages({library(tidyverse); library(loo)})
set.seed(123)

cc <- readRDS("cache/sal_fast_figure_cache_kfold.rds")
bd <- cc$beta_draws_cache

pwmat <- function(lst, nms) {
  M <- sapply(nms, function(nm) {
    cr <- lst[[nm]]
    cr$pointwise[, grep("^elpd", colnames(cr$pointwise))[1]]
  })
  storage.mode(M) <- "double"
  as.matrix(M)
}

## ---------------- EDSM weights from the authoritative bundle -----------------
bun <- readRDS("fits_edsm_Vfix/kfold5_station_criteria_bundle_27.rds")
Me <- pwmat(bun, paste0("fit_edsm_", 0:26))
we <- as.numeric(loo::stacking_weights(Me)); names(we) <- as.character(0:26)

## salinity-free subset (drop 22, 23, 24, 25, 26 -> 22 models)
nosal <- setdiff(0:26, c(22, 23, 24, 25, 26))
we_ns <- as.numeric(loo::stacking_weights(pwmat(bun, paste0("fit_edsm_", nosal))))
names(we_ns) <- as.character(nosal)

cat("=== EDSM stacking weights, 27-model set (>1e-4) ===\n")
print(round(we[we > 1e-4], 4))
cat("\n=== EDSM stacking weights, 22-model salinity-free set (>1e-4) ===\n")
print(round(we_ns[we_ns > 1e-4], 4))

## ---------------- TPS weights ------------------------------------------------
tnames <- names(cc$tps_criteria_cache)
cat("\ntps criteria names:", paste(tnames, collapse=", "), "\n")
Mt <- pwmat(cc$tps_criteria_cache, tnames)
wt <- as.numeric(loo::stacking_weights(Mt))
names(wt) <- gsub("[^0-9]", "", tnames)
cat("\n=== TPS stacking weights (>1e-4) ===\n"); print(round(wt[wt > 1e-4], 4))

## ---------------- mixture machinery ------------------------------------------
NDRAW <- 40000
mix <- function(weights, draws_for, param) {
  nm <- names(weights)
  present <- sapply(nm, function(k) {
    d <- draws_for(k); !is.null(d) && param %in% names(d)
  })
  w  <- weights / sum(weights)
  wp <- sum(w[present])
  if (wp == 0) return(NULL)

  # (i) inclusion mixture: absent models contribute an exact zero
  pick <- sample(nm, NDRAW, replace = TRUE, prob = w)
  inc  <- vapply(pick, function(k) {
    if (!present[[k]]) return(0)
    d <- draws_for(k); d[[param]][sample.int(nrow(d), 1)]
  }, numeric(1))

  # (ii) conditional mixture: only models containing the term
  wc   <- w[present] / sum(w[present])
  pick2 <- sample(names(wc), NDRAW, replace = TRUE, prob = wc)
  cond <- vapply(pick2, function(k) {
    d <- draws_for(k); d[[param]][sample.int(nrow(d), 1)]
  }, numeric(1))

  tibble(param = param,
         w_in_models = round(wp, 4),
         n_models    = sum(present),
         inc_mean = round(mean(inc), 2),
         inc_lo   = round(unname(quantile(inc, .025)), 2),
         inc_hi   = round(unname(quantile(inc, .975)), 2),
         cond_mean = round(mean(cond), 2),
         cond_lo   = round(unname(quantile(cond, .025)), 2),
         cond_hi   = round(unname(quantile(cond, .975)), 2),
         p_neg     = round(mean(cond < 0), 3))
}

edsm_draws <- function(k) bd[[paste0("fit_edsm_catch", k)]]
tps_draws  <- function(k) bd[[paste0("fit_tps_catch",  k)]]

PARAMS <- c("b_out30_z", "b_distance_along_channel_z", "b_distance_m_z",
            "b_secchi_z", "b_salinity_z", "b_temperature_z")

cat("\n\n================ EDSM, 27-model set (salinity INCLUDED) ================\n")
e_all <- map_dfr(PARAMS, ~ mix(we, edsm_draws, .x)); print(as.data.frame(e_all))

cat("\n================ EDSM, 22-model set (salinity EXCLUDED) ================\n")
e_ns <- map_dfr(PARAMS, ~ mix(we_ns, edsm_draws, .x)); print(as.data.frame(e_ns))

cat("\n================ TPS, 17-model set ================\n")
t_all <- map_dfr(PARAMS, ~ mix(wt, tps_draws, .x)); print(as.data.frame(t_all))

saveRDS(list(we = we, we_ns = we_ns, wt = wt,
             e_all = e_all, e_ns = e_ns, t_all = t_all),
        file.path(DERIVED, "table_s10.rds"))
cat("\nDONE\n")

# =============================================================================
# Part 3 - sensitivity of EDSM stacked coefficients to the model set
# =============================================================================
suppressMessages({library(tidyverse); library(loo)})
set.seed(123)
bun <- readRDS("fits_edsm_Vfix/kfold5_station_criteria_bundle_27.rds")
cc  <- readRDS("cache/sal_fast_figure_cache_kfold.rds"); bd <- cc$beta_draws_cache
pwm <- function(nms) { M <- sapply(nms, function(nm) bun[[nm]]$pointwise[,"elpd_kfold"]); as.matrix(M) }
mixq <- function(nums, param, NDRAW=40000) {
  w <- as.numeric(loo::stacking_weights(pwm(paste0("fit_edsm_", nums)))); names(w) <- as.character(nums)
  w <- w/sum(w)
  pres <- sapply(names(w), function(k) param %in% names(bd[[paste0("fit_edsm_catch",k)]]))
  if (!any(pres)) return(NULL)
  pick <- sample(names(w), NDRAW, TRUE, prob=w)
  inc <- vapply(pick, function(k) if(!pres[[k]]) 0 else { d<-bd[[paste0("fit_edsm_catch",k)]]; d[[param]][sample.int(nrow(d),1)] }, numeric(1))
  wc <- w[pres]/sum(w[pres]); p2 <- sample(names(wc), NDRAW, TRUE, prob=wc)
  cond <- vapply(p2, function(k){ d<-bd[[paste0("fit_edsm_catch",k)]]; d[[param]][sample.int(nrow(d),1)] }, numeric(1))
  tibble(param, n=length(nums), w_in=round(sum(w[pres]),4),
         inc=sprintf("%.2f (%.2f, %.2f)", mean(inc), quantile(inc,.025), quantile(inc,.975)),
         cond=sprintf("%.2f (%.2f, %.2f)", mean(cond), quantile(cond,.025), quantile(cond,.975)))
}
P <- c("b_out30_z","b_distance_along_channel_z","b_distance_m_z","b_secchi_z","b_salinity_z")
cat("=== 27-model set ===\n"); print(as.data.frame(map_dfr(P, ~mixq(0:26, .x))))
cat("\n=== 23-model set (drop 22,24,25,26; keep temperature fit 23) ===\n")
print(as.data.frame(map_dfr(P, ~mixq(setdiff(0:26, c(22,24,25,26)), .x))))
cat("\n=== 22-model set (drop 22,23,24,25,26) ===\n")
print(as.data.frame(map_dfr(P, ~mixq(setdiff(0:26, 22:26), .x))))
cat("\n=== original 22-model set (fits 0-21 only) ===\n")
print(as.data.frame(map_dfr(P, ~mixq(0:21, .x))))
cat("\nDONE\n")
