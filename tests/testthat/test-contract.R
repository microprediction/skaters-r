# Runs the hand-rolled contract script unchanged inside testthat. One assertion
# engine, two harnesses: the script stays runnable directly (Rscript
# tests/contract.R) and keeps its own failure counter; here we assert on it.
test_that("contract script passes", {
  # also runs standalone under R CMD check (tests/contract.R); skipping the
  # duplicate here keeps CRAN check time down without losing coverage
  skip_on_cran()
  script <- file.path("..", "contract.R")
  skip_if_not(file.exists(script))
  env <- new.env(parent = globalenv())
  # the script quits on failure when run standalone; under testthat we capture
  # the counter instead
  env$quit <- function(status = 0) stop(sprintf("script exited with status %d", status))
  expect_no_error(sys.source(script, envir = env))
})
