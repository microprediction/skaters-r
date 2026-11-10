# The parity contract under testthat. Same engine as tests/parity.R
# (R/parity_check.R), so the assertions are ported, never rewritten.
test_that("parity vectors are present", {
  skip_if_not(nzchar(skaters_vectors_path()),
              "vectors.json not found; regenerate with parity/gen_vectors.py")
  expect_true(file.exists(skaters_vectors_path()))
})

test_that("all 55 scenarios match the reference to 1e-6", {
  skip_if_not(nzchar(skaters_vectors_path()))
  skip_if_not_installed("jsonlite")
  res <- run_parity(skaters_vectors_path(), verbose = FALSE)
  for (nm in names(res$scenarios)) {
    expect_equal(res$scenarios[[nm]]$fails, 0L,
                 info = paste("scenario", nm, "drifted from the reference"))
  }
  expect_equal(res$fails, 0L)
  expect_gt(res$values_checked, 100000)
  expect_true(res$structure_ok)
})
