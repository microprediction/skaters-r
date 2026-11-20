# Wave 2 parity against inst/parity/wave2_vectors.json, generated from the
# reference by tools/gen_wave2_vectors.py: homogenize, residual_transform and
# its diagnostics, dantzig and search, pushforward affine maps, spec JSON.
.wave2_vectors_path <- function() {
  p <- system.file("parity", "wave2_vectors.json", package = "skaters")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("..", "..", "inst", "parity", "wave2_vectors.json"),
                 file.path("inst", "parity", "wave2_vectors.json"))) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  ""
}

test_that("wave 2 surface matches the reference to 1e-6", {
  skip_if_not(nzchar(.wave2_vectors_path()), "wave2_vectors.json not found")
  skip_if_not_installed("jsonlite")
  rwp <- if (exists("run_wave2_parity")) run_wave2_parity else skaters:::run_wave2_parity
  res <- rwp(.wave2_vectors_path(), verbose = FALSE)
  for (nm in names(res$scenarios)) {
    expect_equal(res$scenarios[[nm]]$fails, 0L, info = paste("scenario", nm))
  }
  for (nm in names(res$sections)) {
    expect_equal(res$sections[[nm]], 0L, info = paste("section", nm))
  }
  expect_length(res$missing, 0)
  expect_equal(res$fails, 0L)
})

test_that("to_json writes Python's compact json.dumps", {
  expect_identical(to_json(ema_spec(1e-5, 1L)), '{"op":"ema","alpha":1e-05,"k":1}')
  expect_identical(to_json(ema_spec(1, 1L)), '{"op":"ema","alpha":1.0,"k":1}')
  expect_identical(to_json(frac_spec(0.4, 30L)), '{"op":"frac","d":0.4,"window":30}')
  expect_identical(to_json(list(a = "é\n", b = list(), c = NULL, d = TRUE)),
                   '{"a":"\\u00e9\\n","b":[],"c":null,"d":true}')
  expect_identical(spec_name(ema_spec(1e-5, 1L)), "ema(1e-05)")
  expect_identical(spec_name(ema_spec(1e16, 1L)), "ema(10000000000000000)")
  sp <- from_json('{"op":"ema","alpha":1e-05,"k":1}')
  expect_true(is.integer(sp$k))
  expect_true(is.double(sp$alpha))
})

test_that("pushforward affine rejects a zero or non-finite multiplier", {
  d <- yeo_johnson(0.0, exact = TRUE)$inverse_k(list(dist_gaussian(0.2, 0.7)), list())[[1]]
  expect_error(dist_scale(d, 0.0))
  expect_error(pushforward_affine(d, 0.0, 1.0))
  expect_error(pushforward_affine(d, Inf, 1.0))
  e <- dist_scale(d, -2.0)
  expect_equal(dist_quantile(e, 0.3), -2.0 * dist_quantile(d, 0.7), tolerance = 1e-9)
})

test_that("the deep-tail spliced quantile resolves the body tail mass", {
  # Mirrors test_deep_tail_quantile_is_resolved in the reference.
  set.seed(5L)
  x <- 0.0
  ys <- numeric(1500)
  for (i in seq_along(ys)) {
    x <- 0.9 * x + rnorm(1) * (1.0 + 0.5 * sin((i - 1) / 7.0))
    ys[i] <- x
  }
  f <- laplace()
  st <- NULL
  for (y in ys) {
    r <- f(y, st)
    st <- r$state
  }
  d <- r$dists[[1]]
  expect_true(isTRUE(d$spliced))
  hit <- 0L
  for (p in c(1e-9, 1e-11, 1 - 1e-9, 1 - 1e-11)) {
    xq <- dist_quantile(d, p)
    upper <- p > 0.5
    expect_identical(xq < dist_quantile(d, 0.5), !upper)
    b <- d$body
    mass <- sum(b$w * pnorm((xq - b$m) / b$s, lower.tail = !upper))
    zt <- if (upper) {
      d$t_up + .gpd_isf((1 - p) / d$zeta_up, d$g_up, d$s_up)
    } else {
      d$t_lo - .gpd_isf(p / d$zeta_lo, d$g_lo, d$s_lo)
    }
    t <- max(pnorm(-abs(zt)), 1e-12)
    if (t < 1e-8) {
      hit <- hit + 1L
      expect_equal(mass, t, tolerance = 1e-6)
    }
  }
  expect_gte(hit, 2L)
})
