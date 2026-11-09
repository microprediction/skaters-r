# Convert granger parameter files (JSON, from the Python package's
# skaters/granger_params) into the RDS files this package ships, so granger
# needs no JSON reader at run time.
#   Rscript tools/granger_params.R ../skaters/src/skaters/granger_params
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "../skaters/src/skaters/granger_params"
for (f in list.files(src, pattern = "\\.json$", full.names = TRUE)) {
  J <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  mat <- function(rows) matrix(unlist(rows), nrow = length(rows), byrow = TRUE)
  P <- list(
    frequency = J$frequency,
    pi_lap = J$pi_lap,
    trunk = list(learning_rate = J$trunk$learning_rate, complexity_penalty = J$trunk$complexity_penalty,
                 forget = J$trunk$forget, prior_log_weights = unlist(J$trunk$prior_log_weights)),
    tleaf = list(eta = J$tleaf$eta, scale_alpha = J$tleaf$scale_alpha),
    candidates = lapply(J$candidates, function(c) list(
      chain = lapply(c$chain, function(e) list(name = e[[1]], kw = e[[2]])),
      depth = as.numeric(c$depth))),
    gate = list(l1_weight = mat(J$gate$l1_weight), l1_bias = unlist(J$gate$l1_bias),
                l2_weight = mat(J$gate$l2_weight), l2_bias = unlist(J$gate$l2_bias)),
    source = J$source
  )
  out <- file.path("inst", "extdata", "granger", sub("\\.json$", ".rds", basename(f)))
  saveRDS(P, out)
  cat(sprintf("[granger_params] %s -> %s\n", f, out))
}
