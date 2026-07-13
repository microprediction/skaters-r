# Parity against the Python-generated vectors, same 1e-6 discipline as the
# JS twin. Run: Rscript tests/parity.R
# Runs in two contexts: repo root (Rscript tests/parity.R) sources R/;
# R CMD check uses the installed package and system.file.
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
  .vec_path <- system.file("parity", "vectors.json", package = "skaters")
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) source(.f)
  .vec_path <- "inst/parity/vectors.json"
}
source("R/runstats.R"); source("R/ema.R"); source("R/ensemble.R"); source("R/bayesian.R")
source("R/multiscale.R"); source("R/sticky.R"); source("R/tails.R"); source("R/terminal.R")
source("R/parade.R"); source("R/api.R")
v <- jsonlite::fromJSON(.vec_path, simplifyVector = FALSE)
series <- unlist(v$series)
ATOL <- 1e-6; RTOL <- 1e-6
probe <- function(d, p, qlo, qhi) c(dist_mean(d), dist_std(d), dist_logpdf(d, p),
                                    dist_cdf(d, p), dist_quantile(d, qlo),
                                    dist_quantile(d, qhi), dist_crps(d, p))
scenarios <- list()
for (k in c(1L, 3L)) {
  suf <- if (k == 1L) "" else sprintf("_k%d", k)
  add <- function(name, sk) scenarios[[paste0(name, suf)]] <<- list(k = k, sk = sk)
  add("leaf", leaf(k))
  add("diff", conjugate(leaf(k), difference(), k))
  add("ema_t", conjugate(leaf(k), ema_transform(0.1), k))
  add("standardize", conjugate(leaf(k), standardize(), k))
  add("theta", conjugate(leaf(k), theta(0.1), k))
  add("drift", conjugate(leaf(k), drift(0.05, 0.01), k))
  add("holt", conjugate(leaf(k), holt_linear(0.1, 0.05), k))
  add("garch", conjugate(leaf(k), garch(), k))
  add("seasonal", conjugate(leaf(k), seasonal_difference(7L), k))
  add("power", conjugate(leaf(k), power_transform(0.5), k))
  add("ar1", conjugate(leaf(k), ar(1L), k))
  add("ar2", conjugate(leaf(k), ar(2L, decay = 1), k))
  add("frac", conjugate(leaf(k), fractional_difference(0.4, 30L), k))
  add("grouped_ar", conjugate(leaf(k), grouped_ar(8L), k))
  add("yeojohnson_log", conjugate(leaf(k), yeo_johnson(0.0), k))
  add("yeojohnson_half", conjugate(leaf(k), yeo_johnson(0.5), k))
  add("ou", conjugate(leaf(k), ou_transform(0.1), k))
  add("ou_sqrt", conjugate(conjugate(leaf(k), ou_transform(0.1), k), yeo_johnson(0.5), k))
  add("ema_skater", ema(0.05, k))
  add("pw_ensemble", precision_weighted_ensemble(list(ema(0.05, k), ema(0.2, k)), k))
  add("multiscale", multiscale(function(kk) conjugate(leaf(kk), ema_transform(0.1), kk), k))
  add("bayes_ensemble", bayesian_ensemble(
    list(ema(0.05, k), conjugate(leaf(k), difference(), k)),
    k = k, learning_rate = 0.5, complexity_penalty = 0.02, depths = c(1, 1)))
}
scenarios[["scale_mixture_leaf"]] <- list(k = 1L, sk = scale_mixture_leaf(1L))
scenarios[["crps_leaf"]] <- list(k = 1L, sk = crps_leaf(1L))
scenarios[["garch_leaf"]] <- list(k = 1L, sk = garch_leaf(1L))
scenarios[["scalemix_ema"]] <- list(
  k = 1L, sk = conjugate(scale_mixture_leaf(1L), ema_transform(0.1), 1L))
scenarios[["gpd_tails"]] <- list(
  k = 1L, sk = gpdtails(conjugate(leaf(1L), ema_transform(0.1), 1L),
                        k = 1L, level = 0.9, nexc = 50L, warmup = 100L))
scenarios[["pol_laplace"]] <- list(k = 1L, sk = laplace(k = 1L))
scenarios[["pol_laplace_k3"]] <- list(k = 3L, sk = laplace(k = 3L))

fails <- 0L; checked <- 0L
check_block <- function(scenarios, series, expected_block) {
  for (name in names(scenarios)) {
    sc <- scenarios[[name]]
    expected <- expected_block[[name]]$out
    st <- NULL; row <- 0L
    for (i in seq_along(series)) {
      r <- sc$sk(series[i], st); st <- r$state
      if (i - 1 >= v$burn) {
        row <- row + 1L
        for (h in seq_len(sc$k)) {
          got <- probe(r$dists[[h]], v$probe, v$q_lo, v$q_hi)
          exp_ <- unlist(expected[[row]][[h]])
          exp_ <- suppressWarnings(as.numeric(exp_))
          for (j in seq_along(got)) {
            checked <<- checked + 1L
            if (is.na(exp_[j])) next
            if (abs(got[j] - exp_[j]) > ATOL + RTOL * abs(exp_[j])) {
              fails <<- fails + 1L
              if (fails < 8) cat(sprintf("FAIL %s row %d h %d probe %d: got %.9g want %.9g\n",
                                          name, row, h, j, got[j], exp_[j]))
            }
          }
        }
      }
    }
    cat(sprintf("ok   %-9s\n", name))
  }
}
check_block(scenarios, series, v$scenarios)

# Sticky/dirac on the repeat-heavy series (exact repeats + 0.25-grid jumps).
repeat_series <- unlist(v$repeat_series)
repeat_scenarios <- list(
  sticky_ema = list(k = 1L, sk = sticky(conjugate(leaf(1L), ema_transform(0.1), 1L), k = 1L))
)
check_block(repeat_scenarios, repeat_series, v$repeat_scenarios)

cat(sprintf("%d values checked\n", checked))
if (fails > 0) { cat(sprintf("PARITY FAILED: %d\n", fails)); quit(status = 1) }
cat("PARITY OK\n")
