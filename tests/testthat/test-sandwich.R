# The sandwich inverse. The assertions that matter are the identity (a
# no-op inner model must reproduce laplace alone) and the density integral
# (a wrong Jacobian still gives plausible quantiles, but not unit mass).

test_that("a no-op inner model is the identity on the body", {
  body <- dist_combine(list(dist_gaussian(0.0, 1.0), dist_gaussian(1.5, 0.4)),
                       c(0.7, 0.3))
  d <- sandwich_new(body, z_mu = 0.0, z_sigma = 1.0)
  for (x in c(-2.5, -0.7, 0.0, 0.9, 2.2)) {
    expect_equal(sandwich_cdf(d, x), dist_cdf(body, x), tolerance = 1e-7)
    # the change-of-variables terms cancel exactly here, so this is not approximate
    expect_equal(sandwich_logpdf(d, x), dist_logpdf(body, x), tolerance = 1e-12)
  }
  for (p in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    expect_equal(sandwich_quantile(d, p), dist_quantile(body, p), tolerance = 1e-6)
  }
})

test_that("the density integrates to one", {
  body <- dist_combine(list(dist_gaussian(0.0, 1.0), dist_gaussian(2.0, 0.5)),
                       c(0.6, 0.4))
  for (cfg in list(c(0.0, 1.0), c(0.5, 0.8), c(-0.9, 1.6))) {
    d <- sandwich_new(body, z_mu = cfg[1], z_sigma = cfg[2])
    lo <- sandwich_quantile(d, 1e-5)
    hi <- sandwich_quantile(d, 1 - 1e-5)
    grid <- seq(lo, hi, length.out = 4001L)
    dens <- vapply(grid, function(x) exp(sandwich_logpdf(d, x)), 0.0)
    mass <- sum((dens[-1] + dens[-length(dens)]) / 2) * (grid[2] - grid[1])
    expect_equal(mass, 1.0, tolerance = 1e-3,
                 info = paste("z_mu", cfg[1], "z_sigma", cfg[2]))
  }
})

test_that("cdf and quantile invert each other", {
  d <- sandwich_new(dist_gaussian(0.3, 1.2), z_mu = 0.4, z_sigma = 0.9)
  for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
    expect_equal(sandwich_cdf(d, sandwich_quantile(d, p)), p, tolerance = 1e-6)
  }
})

test_that("quantiles are monotone and shift with the inner mean", {
  d <- sandwich_new(dist_gaussian(0.0, 1.0), z_mu = 0.0, z_sigma = 1.0)
  q <- vapply(c(0.05, 0.2, 0.5, 0.8, 0.95), function(p) sandwich_quantile(d, p), 0.0)
  expect_true(all(diff(q) > 0))
  med <- function(mu) {
    sandwich_quantile(sandwich_new(dist_gaussian(0.0, 1.0), z_mu = mu), 0.5)
  }
  expect_true(med(-1.0) < med(0.0))
  expect_true(med(0.0) < med(1.0))
  expect_equal(med(0.0), 0.0, tolerance = 1e-6)
})

test_that("dist_* dispatch reaches the sandwich", {
  body <- dist_gaussian(1.0, 2.0)
  d <- sandwich_new(body, z_mu = 0.2, z_sigma = 1.1)
  expect_equal(dist_cdf(d, 1.0), sandwich_cdf(d, 1.0))
  expect_equal(dist_logpdf(d, 1.0), sandwich_logpdf(d, 1.0))
  expect_equal(dist_quantile(d, 0.3), sandwich_quantile(d, 0.3))
  expect_equal(dist_mean(d), sandwich_mean(d))
  expect_equal(dist_var(d), sandwich_var(d))
  expect_equal(dist_crps(d, 1.0), sandwich_crps(d, 1.0))
  expect_true(is.finite(dist_std(d)))
})

test_that("a sandwich composes with a spliced body", {
  set.seed(11)
  # warmup defaults to 500, so the splice needs a series well past that
  y <- rt(1500, df = 3)
  f <- gpdtails(laplace(k = 1), k = 1L)
  st <- NULL
  pend <- NULL
  for (v in y) {
    r <- f(v, st)
    pend <- r$dists
    st <- r$state
  }
  body <- pend[[1]]
  skip_if_not(isTRUE(body$spliced), "splice not active on this sample")
  d <- sandwich_new(body, z_mu = 0.3, z_sigma = 1.2)
  expect_true(is.finite(sandwich_logpdf(d, 0.5)))
  expect_true(sandwich_quantile(d, 0.2) < sandwich_quantile(d, 0.8))
  expect_equal(sandwich_cdf(d, sandwich_quantile(d, 0.7)), 0.7, tolerance = 1e-6)
})

test_that("rosenblatt returns a z stream and a usable transport", {
  set.seed(3)
  y <- cumsum(rnorm(300)) * 0.1
  r <- rosenblatt(y, k = 3L)
  expect_equal(length(r$z), length(y))
  # no predictive had been issued when the first observation arrived
  expect_true(is.na(r$z[1]))
  expect_true(is.na(r$pit[1]))
  expect_false(any(is.na(r$z[-1])))
  expect_true(all(r$pit[-1] > 0 & r$pit[-1] < 1))
  expect_equal(r$k, 3L)
  expect_true(is.finite(dist_mean(sandwich_transport(r, 1L))))
  # steps past k reuse the k-step transport, a disclosed approximation
  expect_identical(sandwich_transport(r, 9L), sandwich_transport(r, 3L))
})

test_that("the z stream of a calibrated forecaster is roughly standard normal", {
  set.seed(5)
  y <- rnorm(3000)
  z <- rosenblatt(y, k = 1L)$z
  z <- z[is.finite(z)]
  expect_gt(length(z), 2500)
  expect_lt(abs(mean(z)), 0.15)
  expect_lt(abs(stats::sd(z) - 1.0), 0.2)
})
