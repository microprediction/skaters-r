# Gaussian mixture distribution: the distributional prediction type.
# Faithful port of skaters/dist.py; parity-checked against vectors.json.

dist_new <- function(w, m, s) {
  stopifnot(length(w) > 0, sum(w) > 0)
  list(w = w / sum(w), m = m, s = s)
}

dist_gaussian <- function(mean = 0.0, std = 1.0) dist_new(1.0, mean, std)

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

dist_mean <- function(d) {
  if (isTRUE(d$mixture)) {
    return(mixture_mean(d))
  }
  if (isTRUE(d$spliced)) {
    return(spliced_mean(d))
  }
  sum(d$w * d$m)
}

dist_var <- function(d) {
  if (isTRUE(d$mixture)) {
    return(mixture_var(d))
  }
  if (isTRUE(d$spliced)) {
    return(spliced_var(d))
  }
  mu <- dist_mean(d)
  sum(d$w * (d$s^2 + (d$m - mu)^2))
}

dist_std <- function(d) {
  v <- dist_var(d)
  if (v > 0) sqrt(v) else 0.0
}

dist_logpdf <- function(d, x) {
  if (isTRUE(d$mixture)) {
    return(mixture_logpdf(d, x))
  }
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

dist_cdf <- function(d, x) {
  if (isTRUE(d$mixture)) {
    return(mixture_cdf(d, x))
  }
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

dist_crps <- function(d, x) {
  if (isTRUE(d$mixture)) {
    return(mixture_crps(d, x))
  }
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

# One mixture CDF value accumulated component by component, in component
# order. The reference adds the terms in a plain left-to-right loop; base
# sum() accumulates in long double where the platform has one, which can
# differ in the last bit. The quantile's bracket checks, bisection, atom snap
# and float-resolution comparison all use this form.
.cdf_terms_sum <- function(w, m, s, x) {
  total <- 0.0
  for (v in w * pnorm(x, m, s)) {
    total <- total + v
  }
  total
}

# Inverse CDF by bisection, a mirror of Dist.quantile in skaters/dist.py
# (contract: parity/QUANTILE_CONTRACT.md there, cases in
# inst/parity/quantile_contract.json here). The x stopping tolerance is
# tol * min(1, smallest positive active component std). An ordinary tolerance
# exit returns the midpoint, except that a probability inside an atom's jump
# returns the atom exactly. At float resolution the rounded midpoint is kept
# when 0 < tol <= 2 * gap and moving the unresolved components across the
# adjacent pair leaves both endpoint CDFs unchanged; otherwise hi is returned.
# The default budget spans the binary64 exponent range; a smaller explicit
# max_iter can still return an unconverged midpoint.
dist_quantile <- function(d, p, tol = 1e-9, max_iter = 2100L) {
  if (isTRUE(d$mixture)) {
    return(mixture_quantile(d, p))
  }
  if (isTRUE(d$spliced)) {
    return(spliced_quantile(d, p, tol = tol, max_iter = max_iter))
  }
  if (!isTRUE(p > 0 && p < 1)) {
    stop(sprintf("quantile needs 0 < p < 1, got %s", format(p)))
  }
  w <- d$w
  m <- d$m
  s <- d$s
  cdf <- function(x) .cdf_terms_sum(w, m, s, x)
  mu <- dist_mean(d)
  sigma <- sqrt(dist_var(d))
  if (!is.finite(sigma)) {
    stop("quantile needs a finite mixture variance")
  }
  active <- s[w > 0 & s > 0]
  scale <- if (length(active) > 0) min(active) else 1.0
  x_tol <- tol * min(1.0, scale)
  lo <- mu - 8 * sigma
  hi <- mu + 8 * sigma
  # Geometric expansion of the moment bracket. The location-scaled seed also
  # moves a collapsed bracket at large-magnitude atoms. Strictly below p
  # avoids excluding an earlier solution on a CDF plateau.
  seed <- max(16 * sigma, 1e-12 * max(1.0, abs(mu)))
  width <- seed
  ok <- FALSE
  for (i in seq_len(120)) {
    if (cdf(lo) < p) {
      ok <- TRUE
      break
    }
    lo <- lo - width
    width <- width * 2.0
  }
  if (!ok) stop(sprintf("could not bracket quantile p=%s from below", format(p)))
  width <- seed
  ok <- FALSE
  for (i in seq_len(120)) {
    if (cdf(hi) >= p) {
      ok <- TRUE
      break
    }
    hi <- hi + width
    width <- width * 2.0
  }
  if (!ok) stop(sprintf("could not bracket quantile p=%s from above", format(p)))
  if (!(is.finite(lo) && is.finite(hi))) {
    stop("quantile needs a finite search bracket")
  }
  n <- length(w)
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (!is.finite(mid)) {
      mid <- 0.5 * lo + 0.5 * hi
    }
    if (!(lo < mid && mid < hi)) {
      # The reference before this contract may have stopped on this adjacent
      # pair when its width first fell below tol. Keep that rounded midpoint
      # only if unresolved components do not change either endpoint CDF when
      # moved across the pair. No weight cutoff is used.
      gap <- hi - lo
      if (0 < tol && tol <= 2 * gap) {
        lower <- 0.0
        upper <- 0.0
        swapped_lower <- 0.0
        swapped_upper <- 0.0
        for (j in seq_len(n)) {
          cl <- w[j] * pnorm(lo, m[j], s[j])
          cu <- w[j] * pnorm(hi, m[j], s[j])
          lower <- lower + cl
          upper <- upper + cu
          unresolved <- w[j] > 0 && s[j] <= gap
          swapped_lower <- swapped_lower + (if (unresolved) cu else cl)
          swapped_upper <- swapped_upper + (if (unresolved) cl else cu)
        }
        if (swapped_lower == lower && swapped_upper == upper) {
          return(mid)
        }
      }
      return(hi)
    }
    if (cdf(mid) < p) {
      lo <- mid
    } else {
      hi <- mid
    }
    if (hi - lo < x_tol) {
      # A midpoint just below an atom has the wrong CDF even if its x error
      # is tiny. Snap only when p is inside the atom's jump.
      for (j in seq_len(n)) {
        if (w[j] > 0 && s[j] == 0 && lo <= m[j] && m[j] <= hi) {
          left <- 0.0
          for (k in seq_len(n)) {
            left <- left + w[k] * (if (s[k] == 0) as.numeric(m[k] < m[j]) else pnorm(m[j], m[k], s[k]))
          }
          if (left < p && p <= cdf(m[j])) {
            return(m[j])
          }
        }
      }
      break
    }
  }
  mid <- 0.5 * (lo + hi)
  if (is.finite(mid)) mid else 0.5 * lo + 0.5 * hi
}

dist_shift <- function(d, delta) dist_new(d$w, d$m + delta, d$s)

dist_scale <- function(d, factor) {
  stopifnot(factor != 0)
  dist_new(d$w, d$m * factor, d$s * abs(factor))
}

dist_affine <- function(d, a, b) {
  stopifnot(a != 0)
  dist_new(d$w, a * d$m + b, abs(a) * d$s)
}

# Reduce component count by merging closest pairs. Exact port of
# Dist.prune: sorted components, ulp-tolerant first-pair-within-threshold
# selection (so platforms that disagree at the last ulp merge the same
# pairs in the same order), moment-matched merges.
dist_prune <- function(d, max_components = 20L) {
  max_components <- max(1L, max_components)
  if (length(d$w) <= max_components) {
    return(d)
  }
  ord <- order(d$m, d$s, d$w)
  w <- d$w[ord]
  m <- d$m[ord]
  s <- d$s[ord]
  # translation-invariant tie scale: the span of the sorted means (skaters#226)
  scale <- (m[length(m)] - m[1]) + 1e-12
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
