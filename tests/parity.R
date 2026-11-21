# Parity against the Python-generated vectors, 1e-6, both harnesses.
# Engine lives in R/parity_check.R (run_parity, run_transform_parity,
# run_wave2_parity, run_state_parity); this script is one entry point, the testthat files
# test-parity.R, test-transform-parity.R, test-wave2-parity.R and
# test-state.R are the
# other. Run: Rscript tests/parity.R
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
  .path <- function(f) system.file("parity", f, package = "skaters")
  .rp <- skaters:::run_parity
  .rtp <- skaters:::run_transform_parity
  .rwp <- skaters:::run_wave2_parity
  .rsp <- skaters:::run_state_parity
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) source(.f)
  .path <- function(f) file.path("inst", "parity", f)
  .rp <- run_parity
  .rtp <- run_transform_parity
  .rwp <- run_wave2_parity
  .rsp <- run_state_parity
}
.res <- .rp(.path("vectors.json"), verbose = TRUE)
.tres <- .rtp(.path("transform_vectors.json"), verbose = TRUE)
.wres <- .rwp(.path("wave2_vectors.json"), verbose = TRUE)
.sres <- .rsp(.path("state_vectors.json.gz"), verbose = TRUE)
if (.res$fails > 0L || .tres$fails > 0L || .wres$fails > 0L || .sres$fails > 0L) quit(status = 1)
