# Parity against the Python-generated vectors, same 1e-6 discipline as the
# JS twin. Run: Rscript tests/parity.R
source("R/dist.R"); source("R/leaf.R"); source("R/transform.R"); source("R/conjugate.R")
v <- jsonlite::fromJSON("inst/parity/vectors.json", simplifyVector = FALSE)
series <- unlist(v$series)
ATOL <- 1e-6; RTOL <- 1e-6
probe <- function(d, p, qlo, qhi) c(dist_mean(d), dist_std(d), dist_logpdf(d, p),
                                    dist_cdf(d, p), dist_quantile(d, qlo),
                                    dist_quantile(d, qhi), dist_crps(d, p))
scenarios <- list(
  leaf      = list(k = 1L, sk = leaf(1L)),
  leaf_k3   = list(k = 3L, sk = leaf(3L)),
  diff      = list(k = 1L, sk = conjugate(leaf(1L), difference(), 1L)),
  diff_k3   = list(k = 3L, sk = conjugate(leaf(3L), difference(), 3L)),
  ema_t     = list(k = 1L, sk = conjugate(leaf(1L), ema_transform(0.1), 1L)),
  ema_t_k3  = list(k = 3L, sk = conjugate(leaf(3L), ema_transform(0.1), 3L))
)
fails <- 0L; checked <- 0L
for (name in names(scenarios)) {
  sc <- scenarios[[name]]
  expected <- v$scenarios[[name]]$out
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
          checked <- checked + 1L
          if (is.na(exp_[j])) next
          if (abs(got[j] - exp_[j]) > ATOL + RTOL * abs(exp_[j])) {
            fails <- fails + 1L
            if (fails < 8) cat(sprintf("FAIL %s row %d h %d probe %d: got %.9g want %.9g\n",
                                        name, row, h, j, got[j], exp_[j]))
          }
        }
      }
    }
  }
  cat(sprintf("ok   %-9s\n", name))
}
cat(sprintf("%d values checked\n", checked))
if (fails > 0) { cat(sprintf("PARITY FAILED: %d\n", fails)); quit(status = 1) }
cat("PARITY OK\n")
