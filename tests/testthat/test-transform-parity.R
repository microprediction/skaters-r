# Transform parity against inst/parity/transform_vectors.json, generated from
# the reference by tools/gen_transform_vectors.py: seasonal_scale,
# gaussianize, log_transform, flesaker, the exact flags, pushforward
# composition and the spliced paths. Engine in R/parity_check.R.
.transform_vectors_path <- function() {
  p <- system.file("parity", "transform_vectors.json", package = "skaters")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("..", "..", "inst", "parity", "transform_vectors.json"),
                 file.path("inst", "parity", "transform_vectors.json"))) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  ""
}

test_that("new transforms match the reference to 1e-6", {
  skip_if_not(nzchar(.transform_vectors_path()), "transform_vectors.json not found")
  skip_if_not_installed("jsonlite")
  rtp <- if (exists("run_transform_parity")) run_transform_parity else skaters:::run_transform_parity
  res <- rtp(.transform_vectors_path(), verbose = FALSE)
  for (nm in names(res$scenarios)) {
    expect_equal(res$scenarios[[nm]]$fails, 0L, info = paste("scenario", nm, "drifted from the reference"))
  }
  expect_length(res$missing, 0)
  expect_equal(res$fails, 0L)
  expect_gt(res$values_checked, 100000)
})

test_that("log_transform rejects observations on or past the bound", {
  expect_error(log_transform()$forward(0.0, NULL))
  expect_error(log_transform(upper = 1.0)$forward(2.0, NULL))
  expect_error(log_transform(lower = 0, upper = 1))
  expect_true(is.nan(log_transform()$forward(NaN, NULL)$y))
})
