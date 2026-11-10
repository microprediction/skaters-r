# Parity against the Python-generated vectors, 1e-6, both harnesses.
# Engine lives in R/parity_check.R (run_parity); this script is one entry point,
# tests/testthat/test-parity.R is the other. Run: Rscript tests/parity.R
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
  .vec_path <- system.file("parity", "vectors.json", package = "skaters")
  .rp <- skaters:::run_parity
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) source(.f)
  .vec_path <- "inst/parity/vectors.json"
  .rp <- run_parity
}
.res <- .rp(.vec_path, verbose = TRUE)
if (.res$fails > 0L) quit(status = 1)
