# Transforms that read mixture components, fed laplace's tail-spliced
# predictive. Mirrors tests/test_spliced_nonaffine.py in the reference. The
# nonlinear inverses (power_transform, yeo_johnson, gaussianize) emit the
# exact pushforward of the spliced law; the k>1 widening of ema_transform,
# seasonal_difference and seasonal_anchor widens the body and lets the splice
# follow it.

.sn_series <- function(n, seed = 5L) {
  set.seed(seed)
  x <- 0.0
  out <- numeric(n)
  for (i in seq_len(n)) {
    x <- 0.9 * x + rnorm(1) * (1.0 + 0.5 * sin((i - 1) / 7.0))
    out[i] <- x
  }
  out
}

.sn_run <- function(f, n = 1500L) {
  st <- NULL
  for (y in .sn_series(n)) {
    r <- f(y, st)
    st <- r$state
  }
  r$dists
}

.sn_spliced <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- .sn_run(laplace())[[1]]
    }
    cached
  }
})

.sn_sane <- function(d, info) {
  expect_true(is.finite(dist_logpdf(d, dist_quantile(d, 0.5))), info = info)
  qs <- vapply(c(0.001, 0.1, 0.5, 0.9, 0.999), function(p) dist_quantile(d, p), 0.0)
  expect_true(all(diff(qs) > 0), info = info)
  expect_true(is.finite(dist_mean(d)), info = info)
  expect_gt(dist_std(d), 0)
}

test_that("nonlinear inverses of a spliced dist are the exact pushforward", {
  d <- .sn_spliced()
  expect_true(isTRUE(d$spliced))
  cases <- list(
    list(tr = power_transform(0.5), inv = function(z) sign(z) * abs(z)^2.0),
    list(tr = yeo_johnson(0.0), inv = function(z) if (z >= 0) expm1(z) else 1.0 - sqrt(1.0 - 2.0 * z))
  )
  for (cs in cases) {
    e <- cs$tr$inverse_k(list(d), list())[[1]]
    expect_true(isTRUE(e$pushforward))
    for (p in c(0.001, 0.05, 0.5, 0.95, 0.999)) {
      want <- cs$inv(dist_quantile(d, p))
      expect_equal(dist_quantile(e, p), want, tolerance = 1e-4, scale = max(abs(want), 1e-2))
    }
    for (z in c(-1.0, 0.2, 1.5)) {
      # table resolution
      expect_equal(dist_cdf(e, cs$inv(z)), dist_cdf(d, z), tolerance = 1e-3, scale = 1)
    }
    expect_true(is.finite(dist_crps(e, 0.3)))
    expect_gt(dist_crps(e, 0.3), 0)
  }
})

test_that("widening a spliced dist moves the body and keeps the splice", {
  d <- .sn_spliced()
  e <- .widen(d, 2.0, 0.5)
  expect_true(isTRUE(e$spliced))
  expect_identical(
    c(e$t_lo, e$t_up, e$zeta_lo, e$zeta_up, e$g_lo, e$s_lo, e$g_up, e$s_up),
    c(d$t_lo, d$t_up, d$zeta_lo, d$zeta_up, d$g_lo, d$s_lo, d$g_up, d$s_up)
  )
  expect_equal(dist_mean(e$body), dist_mean(d$body) + 2.0, tolerance = 1e-9)
  expect_equal(dist_var(e$body), dist_var(d$body) + 0.5, tolerance = 1e-9)
})

test_that("widening a mixture is unchanged", {
  d <- dist_new(c(0.3, 0.7), c(-1.0, 2.0), c(0.5, 1.5))
  e <- .widen(d, 0.4, 0.25)
  expect_identical(e, dist_new(d$w, d$m + 0.4, sqrt(d$s * d$s + 0.25)))
})

test_that("conjugated laplace runs past the splice", {
  # The splice is live from about step 510.
  cases <- list(
    power_transform = list(make = function() power_transform(0.5), k = 1L),
    yeo_johnson = list(make = function() yeo_johnson(0.5), k = 1L),
    gaussianize = list(make = function() gaussianize(), k = 1L),
    ema_transform = list(make = function() ema_transform(0.05), k = 3L),
    seasonal_difference = list(make = function() seasonal_difference(2L), k = 3L),
    seasonal_anchor = list(make = function() seasonal_anchor(2L), k = 3L)
  )
  for (nm in names(cases)) {
    k <- cases[[nm]]$k
    ds <- .sn_run(conjugate(laplace(k), cases[[nm]]$make(), k), n = 600L)
    expect_equal(length(ds), k, info = nm)
    for (d in ds) {
      .sn_sane(d, nm)
    }
  }
})
