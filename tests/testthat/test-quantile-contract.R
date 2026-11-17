# Cross-port inverse-CDF contract: x accuracy, spread, atoms and reviewed
# regressions. The cases are inst/parity/quantile_contract.json, copied
# unchanged from the Python reference (parity/quantile_contract.json there,
# described in parity/QUANTILE_CONTRACT.md) and shared with the Python,
# JavaScript and Rust runners. The assertions mirror
# tests/test_quantile_inverse_contract.py in the reference.
#
# Tolerance key: the fixture's `atol`, the cross-port key the JavaScript and
# Rust runners use. The Python runner prefers a `python_atol` override when a
# case carries one; that key is Python-specific and is ignored here.

quantile_contract_path <- function() {
  p <- system.file("parity", "quantile_contract.json", package = "skaters")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("..", "..", "inst", "parity", "quantile_contract.json"),
                 file.path("inst", "parity", "quantile_contract.json"))) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  ""
}

contract_dist <- function(components) {
  col <- function(i) vapply(components, function(x) as.numeric(x[[i]]), 0.0)
  dist_new(col(1), col(2), col(3))
}

test_that("quantile contract cases are present", {
  expect_true(nzchar(quantile_contract_path()))
})

test_that("all 64 quantile contract cases hold", {
  skip_if_not(nzchar(quantile_contract_path()))
  skip_if_not_installed("jsonlite")
  cases <- jsonlite::fromJSON(quantile_contract_path(), simplifyVector = FALSE)
  expect_equal(length(cases), 64L)
  for (case in cases) {
    d <- contract_dist(case$components)
    q <- dist_quantile(d, case$p)
    info <- sprintf("%s: q=%.17g", case$name, q)
    if (!is.null(case$expected)) {
      expect_true(abs(q - case$expected) <= case$atol, info = info)
    } else {
      expect_true(case$lower < q && q < case$upper, info = info)
    }
    if (!is.null(case$cdf_tol)) {
      cq <- dist_cdf(d, q)
      expect_true(abs(cq - case$p) <= case$cdf_tol,
                  info = sprintf("%s cdf=%.17g", info, cq))
    }
  }
})

test_that("tiny spread is preserved", {
  d <- dist_new(1, 0, 1e-11)
  spread <- dist_quantile(d, 0.99) - dist_quantile(d, 0.01)
  expect_lt(abs(spread - 4.6526957480816816e-11), 2e-20)
})

test_that("explicit budget and resolution stop", {
  d <- dist_new(1, 0, 1)
  expect_identical(dist_quantile(d, 0.5, max_iter = 0L), 0)
  expect_identical(dist_quantile(d, 0.5, max_iter = 1L), -4)
  for (tol in c(0, -1)) {
    q <- dist_quantile(d, 0.5, tol = tol)
    expect_lt(abs(q), 1e-15)
    expect_gte(dist_cdf(d, q), 0.5)
  }
})

test_that("unit-scale CDF call count is unchanged", {
  # two bracket checks, 34 historical bisections
  where <- environment(dist_quantile)
  calls <- 0L
  suppressMessages(trace(".cdf_terms_sum", tracer = function() calls <<- calls + 1L,
                         where = where, print = FALSE))
  on.exit(suppressMessages(untrace(".cdf_terms_sum", where = where)), add = TRUE)
  dist_quantile(dist_new(1, 0, 1), 0.9)
  expect_equal(calls, 36L)
})

test_that("variance overflow is explicit", {
  d <- dist_new(c(0.5, 0.5), c(-1e155, 1e155), c(1, 1))
  expect_error(dist_quantile(d, 0.5), "finite mixture variance")
})

test_that("p outside (0, 1) is rejected", {
  d <- dist_new(1, 0, 1)
  for (p in c(0, 1, -0.1, 1.1, NaN)) {
    expect_error(dist_quantile(d, p), "quantile needs 0 < p < 1")
  }
})

test_that("the spliced wrapper forwards the default budget to the body", {
  expect_identical(formals(dist_quantile)$max_iter, 2100L)
  expect_identical(formals(spliced_quantile)$max_iter, 2100L)
})
