# Parity against the Python-generated vectors, 1e-6, both harnesses.
# Engine lives in R/parity_check.R (run_parity, run_transform_parity); this
# script is one entry point, tests/testthat/test-parity.R and
# test-transform-parity.R are the other. Run: Rscript tests/parity.R
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
  .vec_path <- system.file("parity", "vectors.json", package = "skaters")
  .tvec_path <- system.file("parity", "transform_vectors.json", package = "skaters")
  .rp <- skaters:::run_parity
  .rtp <- skaters:::run_transform_parity
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) source(.f)
  .vec_path <- "inst/parity/vectors.json"
  .tvec_path <- "inst/parity/transform_vectors.json"
  .rp <- run_parity
  .rtp <- run_transform_parity
}
.res <- .rp(.vec_path, verbose = TRUE)
.tres <- .rtp(.tvec_path, verbose = TRUE)
if (.res$fails > 0L || .tres$fails > 0L) quit(status = 1)
