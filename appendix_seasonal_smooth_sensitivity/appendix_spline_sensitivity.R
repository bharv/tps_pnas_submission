## Regenerate the s(water_year_day) sensitivity analysis and save its outputs:
##   output_figs_appendix/appendix_spline_k_comparison.png  (figure)
##   cache/appendix_spline_sensitivity.rds                  (smooths + stats + loo table)
## Run from the REPOSITORY ROOT (the folder holding the .Rproj), e.g.
##   Rscript appendix_seasonal_smooth_sensitivity/appendix_spline_sensitivity.R
## Fits 5 brms models in memory (none are saved); takes a while.
suppressMessages({ library(tidyverse); library(brms) })
options(mc.cores = parallel::detectCores())
if (!length(list.files(".", pattern = "\\.Rproj$")))
  stop("Run from the repository root (the folder with the .Rproj file).")
for (dd in c("cache", "output_figs_appendix")) dir.create(dd, showWarnings = FALSE)

d <- read_csv("data/tow_12_6_tps_nona.csv", show_col_types = FALSE)
wyd_mean <- mean(d$water_year_day); wyd_sd <- sd(d$water_year_day)
d$water_year_day_z <- (d$water_year_day - wyd_mean) / wyd_sd
priors <- c(prior(normal(0, 2.5), class = "b"), prior(gamma(2, 0.1), class = "shape"))

fit_k <- function(k, zscore = FALSE) {
  v  <- if (zscore) "water_year_day_z" else "water_year_day"
  kt <- if (is.na(k)) "" else paste0(", k = ", k)
  form <- as.formula(sprintf("n_l_75 ~ out30_z + s(%s%s) + offset(log(Volume))", v, kt))
  brm(formula = form, family = negbinomial("log"), data = d, prior = priors,
      iter = 3000, warmup = 2000, chains = 4,
      control = list(max_treedepth = 12, adapt_delta = 0.99), seed = 1, refresh = 0)
}

specs <- tribble(~tag, ~k, ~z,
  "k=10 (default)",   10L, FALSE,
  "k=6",               6L, FALSE,
  "k=4",               4L, FALSE,
  "k=3",               3L, FALSE,
  "k=10 standardized",10L, TRUE)

smooths <- list(); loos <- list(); stat <- list()
for (i in seq_len(nrow(specs))) {
  tg <- specs$tag[i]; message("=== ", tg, " ===")
  fit <- fit_k(specs$k[i], specs$z[i])
  cs  <- conditional_smooths(fit, spaghetti = FALSE)[[1]]
  xcol <- if (specs$z[i]) "water_year_day_z" else "water_year_day"
  smooths[[tg]] <- tibble(
    x_raw = if (specs$z[i]) cs[[xcol]] * wyd_sd + wyd_mean else cs[[xcol]],
    est = cs$estimate__, lo = cs$lower__, hi = cs$upper__, tag = tg)
  lo <- loo(fit); loos[[tg]] <- lo
  np <- nuts_params(fit)
  stat[[tg]] <- tibble(
    tag = tg, k = specs$k[i], standardized = specs$z[i],
    sds = as.data.frame(summary(fit)$splines)[1, "Estimate"],
    divergences = sum(subset(np, Parameter == "divergent__")$Value),
    elpd_loo = lo$estimates["elpd_loo", "Estimate"],
    se_elpd  = lo$estimates["elpd_loo", "SE"],
    max_pareto_k = max(lo$diagnostics$pareto_k))
  rm(fit); gc()
}

stat_tbl <- bind_rows(stat)
cmp <- as.data.frame(loo_compare(loos)); cmp$tag <- rownames(cmp)
sm  <- bind_rows(smooths); sm$tag <- factor(sm$tag, levels = specs$tag)

p <- ggplot(sm, aes(x_raw, est)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, fill = "#2c7fb8") +
  geom_line(colour = "#08519c", linewidth = .7) +
  facet_wrap(~tag, nrow = 1) +
  labs(x = "Water-year day", y = "Partial effect of s(water-year day) (log scale)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave("output_figs_appendix/appendix_spline_k_comparison.png", p, width = 17.8, height = 5,
       units = "cm", dpi = 600)

saveRDS(list(stats = stat_tbl, loo_compare = cmp, smooths = sm,
             wyd_mean = wyd_mean, wyd_sd = wyd_sd, n = nrow(d),
             prop_zero = mean(d$n_l_75 == 0), max_catch = max(d$n_l_75)),
        "cache/appendix_spline_sensitivity.rds")

cat("\nSaved figure -> output_figs_appendix/appendix_spline_k_comparison.png\n")
cat("Saved cache  -> cache/appendix_spline_sensitivity.rds\n")
print(stat_tbl); cat("\n"); print(cmp)
cat("\nDONE\n")
