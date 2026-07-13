# Adversarial release gate, mirroring the reference repo's
# parity/adversarial.mjs and tests/test_tails_robustness.py.
#
# Runs the deployed default (laplace, GPD tails) over the pathological
# streams a production deployment will eventually meet: constant series,
# lattice/repeats, a monster spike, an extreme finite tick, scale collapse,
# vol whiplash on a trend. Asserts the contract at checkpoints: forecasts
# stay finite and well-formed, the cdf stays a cdf, quantiles stay ordered,
# and the detector recovers after the insult.
#
# Run: Rscript tests/robustness.R   (also runs under R CMD check)
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) source(.f)
}

# Deterministic RNG (same LCG as the JS gate; no set.seed dependence on
# R's generator versioning).
lcg <- function(seed) {
  s <- seed %% 2^32
  function() {
    s <<- (1664525 * s + 1013904223) %% 2^32
    s / 2^32
  }
}
gauss <- function(rand) {
  u <- 0.0; v <- 0.0
  while (u == 0) u <- rand()
  while (v == 0) v <- rand()
  sqrt(-2.0 * log(u)) * cos(2.0 * pi * v)
}

failures <- 0L
check <- function(cond, label) {
  if (!isTRUE(cond)) {
    failures <<- failures + 1L
    cat(sprintf("FAIL %s\n", label))
  }
}

assert_wellformed <- function(d, y_near, label) {
  lp <- dist_logpdf(d, y_near)
  check(!is.na(lp), sprintf("%s: logpdf NaN", label))   # +/-Inf tolerated
  cv <- dist_cdf(d, y_near)
  check(cv >= 0.0 && cv <= 1.0, sprintf("%s: cdf out of [0,1]", label))
  qs <- vapply(c(0.001, 0.25, 0.5, 0.75, 0.999),
               function(p) dist_quantile(d, p), 0.0)
  check(all(is.finite(qs)), sprintf("%s: non-finite quantile", label))
  check(all(diff(qs) >= -1e-9), sprintf("%s: quantiles unordered", label))
  probes <- c(qs[1] - 1.0, qs, qs[length(qs)] + 1.0)
  cs <- vapply(probes, function(x) dist_cdf(d, x), 0.0)
  check(all(diff(cs) >= -1e-9), sprintf("%s: cdf not monotone", label))
}

soak <- function(ys, label) {
  f <- laplace(1L)
  st <- NULL; dists <- NULL
  for (t in seq_along(ys)) {
    r <- f(ys[t], st); st <- r$state; dists <- r$dists
    if (t > 1 && (t - 1) %% 997 == 0) assert_wellformed(dists[[1]], ys[t], label)
  }
  assert_wellformed(dists[[1]], ys[length(ys)], label)
  invisible(st)
}

# 1. constant series
soak(rep(3.7, 4000), "constant")
cat("ok   constant\n")

# 2. lattice / repeats (exact repeats + 0.25-grid jumps)
{
  rand <- lcg(17)
  v <- 1.0
  ys <- numeric(6000)
  for (i in 1:6000) {
    if (rand() >= 0.7) v <- v + c(-0.25, 0.25, 0.5)[floor(rand() * 3) + 1]
    ys[i] <- v
  }
  soak(ys, "lattice")
  cat("ok   lattice\n")
}

# 3. monster spike, then recovery (no deafness, no permanent alarm)
{
  rand <- lcg(23)
  f <- laplace(1L)
  st <- NULL
  for (i in 1:3000) { r <- f(gauss(rand), st); st <- r$state }
  r <- f(1e9, st); st <- r$state                     # the insult
  assert_wellformed(r$dists[[1]], 0.0, "spike:after")
  alarms <- 0L; n <- 0L
  for (i in 1:3000) {
    y <- gauss(rand)
    r <- f(y, st); st <- r$state
    z <- st$z[1]
    if (i > 1001 && !is.na(z)) {                     # mjs: i > 1000, 0-based
      n <- n + 1L
      if (abs(z) > 2.5758) alarms <- alarms + 1L     # ~1e-2 two-sided
    }
    if ((i - 1) %% 500 == 0) assert_wellformed(r$dists[[1]], y, "spike:recovery")
  }
  check(n > 1500, "spike: too few matured ticks")
  check(alarms / n < 0.06,
        sprintf("spike: alarm rate %.4f after recovery", alarms / n))
  cat(sprintf("ok   spike (alarm rate %.4f on %d matured ticks)\n", alarms / n, n))
}

# 3b. extreme finite tick: the input gate must keep the tree alive
{
  rand <- lcg(29)
  f <- laplace(1L)
  st <- NULL
  for (i in 1:1500) { r <- f(gauss(rand), st); st <- r$state }
  r <- f(1e300, st); st <- r$state                   # near the double limit
  assert_wellformed(r$dists[[1]], 0.0, "extreme:after")
  for (i in 1:1500) {
    r <- f(gauss(rand), st); st <- r$state
    if ((i - 1) %% 500 == 0) assert_wellformed(r$dists[[1]], 0.0, "extreme:recovery")
  }
  cat("ok   extreme finite tick\n")
}

# 4. scale collapse and recovery
{
  rand <- lcg(31)
  ys <- c(vapply(1:2000, function(i) gauss(rand), 0.0),
          rep(0.0, 2000),
          vapply(1:2000, function(i) gauss(rand), 0.0))
  soak(ys, "collapse")
  cat("ok   collapse\n")
}

# 5. vol whiplash on a trend
{
  rand <- lcg(41)
  ys <- numeric(8000)
  lvl <- 0.0
  for (t in 0:7999) {
    vol <- if ((t %/% 700) %% 2 == 1) 10.0 else 1.0
    lvl <- lvl + 0.05 + vol * gauss(rand)
    ys[t + 1] <- lvl
  }
  soak(ys, "whiplash")
  cat("ok   whiplash\n")
}

# 6. long-series soak: 3000 gaussian ticks so the default 500-tick tail
# warm-up activates and the SplicedDist path is exercised in-package.
{
  rand <- lcg(53)
  f <- laplace(1L)
  st <- NULL; dists <- NULL
  spliced_seen <- FALSE
  for (i in 1:3000) {
    y <- gauss(rand)
    r <- f(y, st); st <- r$state; dists <- r$dists
    if (isTRUE(dists[[1]]$spliced)) spliced_seen <- TRUE
    if ((i - 1) %% 500 == 0 && i > 1) assert_wellformed(dists[[1]], y, "soak3000")
  }
  d <- dists[[1]]
  check(spliced_seen, "soak3000: splice never activated")
  check(isTRUE(d$spliced), "soak3000: final predictive not spliced")
  check(is.finite(d$t_lo) && is.finite(d$t_up) && d$t_lo < d$t_up,
        "soak3000: splice thresholds malformed")
  check(d$zeta_lo > 0 && d$zeta_lo < 1 && d$zeta_up > 0 && d$zeta_up < 1,
        "soak3000: splice tail masses out of (0,1)")
  check(is.finite(d$g_lo) && d$s_lo > 0 && is.finite(d$g_up) && d$s_up > 0,
        "soak3000: GPD parameters malformed")
  assert_wellformed(d, 0.0, "soak3000:final")
  # Tail state must round-trip through serialization (the Python test
  # json-dumps it; RDS is the R equivalent).
  tf <- tempfile(fileext = ".rds")
  saveRDS(st, tf)
  check(identical(readRDS(tf), st), "soak3000: state not RDS-stable")
  unlink(tf)
  cat("ok   soak3000 (splice active)\n")
}

if (failures > 0) {
  cat(sprintf("ROBUSTNESS FAILED: %d violation(s)\n", failures))
  quit(status = 1)
}
cat("ROBUSTNESS OK (constant, lattice, spike, extreme, collapse, whiplash, soak)\n")
