# granger R parity against the Python reference (skaters/parity/granger_vectors.json).
#   Rscript inst/parity/check_granger.R ../skaters/parity/granger_vectors.json
args <- commandArgs(trailingOnly = TRUE)
vpath <- if (length(args)) args[1] else "../skaters/parity/granger_vectors.json"
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
V <- jsonlite::fromJSON(vpath, simplifyVector = FALSE)
ATOL <- 1e-6; RTOL <- 1e-6
worst <- 0; bad <- 0
for (key in names(V)) {
  freq <- strsplit(key, "/")[[1]][1]
  g <- granger(granger_params(freq))
  ys <- unlist(V[[key]]$ys); rows <- V[[key]]$rows
  s <- NULL
  for (t in seq_len(length(ys) - 1)) {
    r <- g(ys[t], s); s <- r$state; m <- r$dists[[1]]
    got <- c(dist_logpdf(m, ys[t + 1]), dist_mean(m), dist_std(m), m$w[1], dist_quantile(m, 0.1), dist_quantile(m, 0.9))
    want <- unlist(rows[[t]])
    dev <- abs(got - want)
    worst <- max(worst, max(dev / pmax(1, abs(want))))
    if (any(dev > ATOL + RTOL * abs(want))) {
      if (bad < 5) cat(sprintf("%s step %d: r %s python %s\n", key, t - 1, paste(signif(got, 10), collapse = " "), paste(signif(want, 10), collapse = " ")))
      bad <- bad + 1
    }
  }
}
cat(sprintf("[granger R parity] worst relative deviation %.2e, failing steps %d\n", worst, bad))
quit(status = if (bad) 1 else 0)
