# dist_affine, dist_shift and dist_scale on the spliced predictive. Mirrors
# tests/test_spliced_affine.py in the reference. Conjugating laplace with a
# transform calls these on the predictive, which is spliced once the GPD tail
# splice is live. The splice is defined in the body's PIT space, so the map
# moves only the body; cdf, quantile and logpdf of the mapped law must equal
# the change of variables.

.sa_series <- function(n, seed = 3L) {
  set.seed(seed)
  x <- 0.0
  out <- numeric(n)
  for (i in seq_len(n)) {
    x <- 0.9 * x + rnorm(1) * (1.0 + 0.5 * sin((i - 1) / 7.0))
    out[i] <- x
  }
  out
}

.sa_spliced <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      f <- laplace()
      st <- NULL
      for (y in .sa_series(1500L)) {
        r <- f(y, st)
        st <- r$state
      }
      cached <<- r$dists[[1]]
    }
    cached
  }
})

.sa_check_map <- function(d, e, a, b) {
  expect_true(isTRUE(e$spliced))
  for (p in c(0.001, 0.05, 0.5, 0.95, 0.999)) {
    want <- a * dist_quantile(d, if (a > 0) p else 1 - p) + b
    expect_equal(dist_quantile(e, p), want, tolerance = 1e-6, scale = max(abs(want), 1e-3))
  }
  for (x in c(-2.0, 0.3, 4.0)) {
    u <- dist_cdf(d, x)
    expect_equal(dist_cdf(e, a * x + b), if (a > 0) u else 1 - u, tolerance = 1e-9, scale = 1)
    expect_equal(dist_logpdf(e, a * x + b), dist_logpdf(d, x) - log(abs(a)), tolerance = 1e-6, scale = 1)
  }
}

test_that("dist_affine on a spliced dist is the change of variables", {
  d <- .sa_spliced()
  expect_true(isTRUE(d$spliced))
  for (ab in list(c(2.5, -1.0), c(0.1, 3.0), c(-1.5, 0.7))) {
    .sa_check_map(d, dist_affine(d, ab[1], ab[2]), ab[1], ab[2])
  }
})

test_that("dist_shift and dist_scale on a spliced dist are the change of variables", {
  d <- .sa_spliced()
  .sa_check_map(d, dist_shift(d, 0.7), 1.0, 0.7)
  for (g in c(2.5, 0.1, -1.5)) {
    .sa_check_map(d, dist_scale(d, g), g, 0.0)
  }
})

test_that("spliced affine and scale reject a zero multiplier", {
  d <- .sa_spliced()
  expect_error(dist_affine(d, 0.0, 1.0))
  expect_error(dist_scale(d, 0.0))
  expect_error(spliced_affine(d, 0.0, 1.0))
  expect_error(spliced_scale(d, 0.0))
})

test_that("transforms wrap laplace past the splice", {
  # Before the fix these failed in dist_new once the GPD splice was live.
  # seasonal_scale is not yet ported to R.
  transforms <- list(
    standardize = function() standardize(0.05),
    garch = function() garch(),
    ema_transform = function() ema_transform(0.05),
    difference = function() difference(),
    seasonal_difference = function() seasonal_difference(7L),
    seasonal_anchor = function() seasonal_anchor(7L),
    ou_transform = function() ou_transform()
  )
  ys <- .sa_series(1500L)
  for (nm in names(transforms)) {
    f <- conjugate(laplace(), transforms[[nm]](), k = 1L)
    st <- NULL
    for (y in ys) {
      r <- f(y, st)
      st <- r$state
    }
    d <- r$dists[[1]]
    expect_true(is.finite(dist_logpdf(d, 0.0)), info = nm)
    expect_lt(dist_quantile(d, 0.1), dist_quantile(d, 0.9))
  }
})
