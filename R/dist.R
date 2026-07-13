# Gaussian mixture distribution: the distributional prediction type.
# Faithful port of skaters/dist.py; parity-checked against vectors.json.

dist_new <- function(w, m, s) {
  stopifnot(length(w) > 0, sum(w) > 0)
  list(w = w / sum(w), m = m, s = s)
}

dist_gaussian <- function(mean = 0.0, std = 1.0) dist_new(1.0, mean, std)

dist_combine <- function(dists, weights = NULL) {
  n <- length(dists)
  if (is.null(weights)) weights <- rep(1.0 / n, n)
  wt <- sum(weights)
  w <- c(); m <- c(); s <- c()
  for (i in seq_len(n)) {
    d <- dists[[i]]
    w <- c(w, weights[i] / wt * d$w); m <- c(m, d$m); s <- c(s, d$s)
  }
  dist_new(w, m, s)
}

dist_mean <- function(d) sum(d$w * d$m)

dist_var <- function(d) {
  mu <- dist_mean(d)
  sum(d$w * (d$s^2 + (d$m - mu)^2))
}

dist_std <- function(d) { v <- dist_var(d); if (v > 0) sqrt(v) else 0.0 }

dist_logpdf <- function(d, x) {
  ok <- d$w > 0 & d$s > 0
  if (any(d$s <= 0 & x == d$m & d$w > 0)) return(Inf)
  w <- d$w[ok]; m <- d$m[ok]; s <- d$s[ok]
  if (length(w) == 0) return(-Inf)
  z <- (x - m) / s
  t <- log(w) - 0.5 * z * z - log(s) - 0.5 * log(2 * pi)
  b <- max(t)
  if (!is.finite(b)) return(-Inf)
  b + log(sum(exp(t - b)))
}

dist_cdf <- function(d, x) sum(d$w * pnorm(x, d$m, d$s))

.abs_expectation <- function(m, s) {
  if (s <= 0) return(abs(m))
  m * (2 * pnorm(m / s) - 1) + 2 * s * dnorm(m / s)
}

dist_crps <- function(d, x) {
  t1 <- sum(vapply(seq_along(d$w),
                   function(i) d$w[i] * .abs_expectation(d$m[i] - x, d$s[i]),
                   0.0))
  t2 <- 0.0
  for (i in seq_along(d$w)) for (j in seq_along(d$w)) {
    t2 <- t2 + d$w[i] * d$w[j] *
      .abs_expectation(d$m[i] - d$m[j], sqrt(d$s[i]^2 + d$s[j]^2))
  }
  t1 - 0.5 * t2
}

dist_quantile <- function(d, p, tol = 1e-9, max_iter = 100L) {
  stopifnot(p > 0, p < 1)
  mu <- dist_mean(d); sigma <- sqrt(dist_var(d))
  lo <- mu - 8 * sigma; hi <- mu + 8 * sigma
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (dist_cdf(d, mid) < p) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  0.5 * (lo + hi)
}

dist_shift <- function(d, delta) dist_new(d$w, d$m + delta, d$s)
