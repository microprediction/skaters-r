# Gaussian mixture distribution: the distributional prediction type.
# Faithful port of skaters/dist.py; parity-checked against vectors.json.

#' Gaussian mixture distribution: the distributional prediction type
#'
#' Constructors and probes for the Gaussian mixture that every skater emits:
#' `dist_new(w, m, s)` builds one, `dist_gaussian` a single component,
#' `dist_combine` a weighted mixture of mixtures. Probes: `dist_mean`,
#' `dist_var`, `dist_std`, `dist_logpdf`, `dist_cdf`, `dist_crps`,
#' `dist_quantile`. Affine maps: `dist_shift`, `dist_scale`, `dist_affine`.
#' `dist_prune` reduces the component count by ulp-tolerant closest-pair
#' merging. Port of the Python reference `skaters/dist.py`, parity-checked
#' at 1e-6.
#'
#' The probes also accept the tail-spliced distributions produced by
#' [gpdtails()] and dispatch to the `spliced_*` functions for them.
#'
#' @param w numeric vector of mixture weights (normalized internally).
#' @param m numeric vector of component means.
#' @param s numeric vector of component standard deviations.
#' @return Constructors and affine maps return a dist: `list(w, m, s)`.
#'   `dist_mean`, `dist_var`, `dist_std`, `dist_logpdf`, `dist_cdf`,
#'   `dist_crps`, and `dist_quantile` return a scalar.
#' @examples
#' d <- dist_new(c(0.5, 0.5), c(0, 1), c(1, 2))
#' dist_mean(d)
#' dist_quantile(d, 0.9)
#' dist_crps(d, 0.3)
#' dist_mean(dist_shift(d, 10))
#' @rdname dist
#' @export
dist_new <- function(w, m, s) {
  stopifnot(length(w) > 0, sum(w) > 0)
  list(w = w / sum(w), m = m, s = s)
}

#' @param mean mean of the single Gaussian component.
#' @param std standard deviation of the single Gaussian component.
#' @rdname dist
#' @export
dist_gaussian <- function(mean = 0.0, std = 1.0) dist_new(1.0, mean, std)

#' @param dists list of dists to combine into one mixture.
#' @param weights weights over `dists`; equal by default.
#' @rdname dist
#' @export
dist_combine <- function(dists, weights = NULL) {
  n <- length(dists)
  if (is.null(weights)) {
    weights <- rep(1.0 / n, n)
  }
  wt <- sum(weights)
  w <- c()
  m <- c()
  s <- c()
  for (i in seq_len(n)) {
    d <- dists[[i]]
    w <- c(w, weights[i] / wt * d$w)
    m <- c(m, d$m)
    s <- c(s, d$s)
  }
  dist_new(w, m, s)
}

#' @param d a dist.
#' @rdname dist
#' @export
dist_mean <- function(d) {
  if (isTRUE(d$spliced)) {
    return(spliced_mean(d))
  }
  sum(d$w * d$m)
}

#' @rdname dist
#' @export
dist_var <- function(d) {
  if (isTRUE(d$spliced)) {
    return(spliced_var(d))
  }
  mu <- dist_mean(d)
  sum(d$w * (d$s^2 + (d$m - mu)^2))
}

#' @rdname dist
#' @export
dist_std <- function(d) {
  v <- dist_var(d)
  if (v > 0) sqrt(v) else 0.0
}

#' @param x point at which to evaluate the density, cdf, or CRPS.
#' @rdname dist
#' @export
dist_logpdf <- function(d, x) {
  if (isTRUE(d$spliced)) {
    return(spliced_logpdf(d, x))
  }
  ok <- d$w > 0 & d$s > 0
  if (any(d$s <= 0 & x == d$m & d$w > 0)) {
    return(Inf)
  }
  w <- d$w[ok]
  m <- d$m[ok]
  s <- d$s[ok]
  if (length(w) == 0) {
    return(-Inf)
  }
  z <- (x - m) / s
  t <- log(w) - 0.5 * z * z - log(s) - 0.5 * log(2 * pi)
  b <- max(t)
  if (!is.finite(b)) {
    return(-Inf)
  }
  b + log(sum(exp(t - b)))
}

#' @rdname dist
#' @export
dist_cdf <- function(d, x) {
  if (isTRUE(d$spliced)) {
    return(spliced_cdf(d, x))
  }
  sum(d$w * pnorm(x, d$m, d$s))
}

.abs_expectation <- function(m, s) {
  if (s <= 0) {
    return(abs(m))
  }
  m * (2 * pnorm(m / s) - 1) + 2 * s * dnorm(m / s)
}

#' @rdname dist
#' @export
dist_crps <- function(d, x) {
  if (isTRUE(d$spliced)) {
    return(spliced_crps(d, x))
  }
  t1 <- sum(vapply(seq_along(d$w), function(i) d$w[i] * .abs_expectation(d$m[i] - x, d$s[i]), 0.0))
  t2 <- 0.0
  for (i in seq_along(d$w)) {
    for (j in seq_along(d$w)) {
      t2 <- t2 + d$w[i] * d$w[j] * .abs_expectation(d$m[i] - d$m[j], sqrt(d$s[i]^2 + d$s[j]^2))
    }
  }
  t1 - 0.5 * t2
}

#' @param p probability in (0, 1).
#' @param tol bisection tolerance for the quantile search.
#' @param max_iter bisection iteration cap.
#' @rdname dist
#' @export
dist_quantile <- function(d, p, tol = 1e-9, max_iter = 100L) {
  if (isTRUE(d$spliced)) {
    return(spliced_quantile(d, p, tol = tol, max_iter = max_iter))
  }
  stopifnot(p > 0, p < 1)
  mu <- dist_mean(d)
  sigma <- sqrt(dist_var(d))
  lo <- mu - 8 * sigma
  hi <- mu + 8 * sigma
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (dist_cdf(d, mid) < p) {
      lo <- mid
    } else {
      hi <- mid
    }
    if (hi - lo < tol) break
  }
  0.5 * (lo + hi)
}

#' @param delta shift added to every component mean.
#' @rdname dist
#' @export
dist_shift <- function(d, delta) dist_new(d$w, d$m + delta, d$s)

#' @param factor scale applied to means and (in absolute value) to sds.
#' @rdname dist
#' @export
dist_scale <- function(d, factor) {
  stopifnot(factor != 0)
  dist_new(d$w, d$m * factor, d$s * abs(factor))
}

#' @param a multiplier of the affine map `a * x + b`; must be nonzero.
#' @param b offset of the affine map.
#' @rdname dist
#' @export
dist_affine <- function(d, a, b) {
  stopifnot(a != 0)
  dist_new(d$w, a * d$m + b, abs(a) * d$s)
}

# Reduce component count by merging closest pairs. Exact port of
# Dist.prune: sorted components, ulp-tolerant first-pair-within-threshold
# selection (so platforms that disagree at the last ulp merge the same
# pairs in the same order), moment-matched merges.
#' @param max_components component budget; closest pairs are moment-matched
#'   and merged until the mixture fits.
#' @rdname dist
#' @export
dist_prune <- function(d, max_components = 20L) {
  max_components <- max(1L, max_components)
  if (length(d$w) <= max_components) {
    return(d)
  }
  ord <- order(d$m, d$s, d$w)
  w <- d$w[ord]
  m <- d$m[ord]
  s <- d$s[ord]
  scale <- abs(m[1]) + abs(m[length(m)]) + 1e-12
  while (length(w) > max_components) {
    n <- length(w)
    best_dist <- Inf
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        dd <- abs(m[i] - m[j])
        if (dd < best_dist) best_dist <- dd
      }
    }
    thresh <- best_dist + 1e-9 * scale
    best_i <- NA_integer_
    best_j <- NA_integer_
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        if (abs(m[i] - m[j]) <= thresh) {
          best_i <- i
          best_j <- j
          break
        }
      }
      if (!is.na(best_i)) break
    }
    if (is.na(best_i)) {
      best_i <- 1L
      best_j <- 2L
    } # NaN means: still terminate
    wi <- w[best_i]
    mi <- m[best_i]
    si <- s[best_i]
    wj <- w[best_j]
    mj <- m[best_j]
    sj <- s[best_j]
    w_new <- wi + wj
    if (w_new < 1e-300) {
      m_new <- 0.5 * (mi + mj)
      s_new <- max(si, sj, 1e-12)
    } else {
      m_new <- (wi * mi + wj * mj) / w_new
      v_new <- (wi * (si * si + (mi - m_new)^2) + wj * (sj * sj + (mj - m_new)^2)) / w_new
      s_new <- sqrt(max(v_new, 0.0))
    }
    w[best_i] <- w_new
    m[best_i] <- m_new
    s[best_i] <- s_new
    w <- w[-best_j]
    m <- m[-best_j]
    s <- s[-best_j]
  }
  dist_new(w, m, s)
}
