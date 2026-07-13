# GPD tail splice: censored-ML generalized-Pareto tails spliced into the
# body's predictive beyond frozen z-thresholds. Port of skaters/tails.py.

.TAILS_EPS <- 1e-12
.LOG_SQRT2PI <- 0.5 * log(2 * pi)

# Acklam's rational approximation to the standard normal quantile, polished
# with one Halley step. Ported verbatim (deterministic, cross-platform).
.ACK_A <- c(-3.969683028665376e+01, 2.209460984245205e+02,
            -2.759285104469687e+02, 1.383577518672690e+02,
            -3.066479806614716e+01, 2.506628277459239e+00)
.ACK_B <- c(-5.447609879822406e+01, 1.615858368580409e+02,
            -1.556989798598866e+02, 6.680131188771972e+01,
            -1.328068155288572e+01)
.ACK_C <- c(-7.784894002430293e-03, -3.223964580411365e-01,
            -2.400758277161838e+00, -2.549732539343734e+00,
            4.374664141464968e+00, 2.938163982698783e+00)
.ACK_D <- c(7.784695709041462e-03, 3.224671290700398e-01,
            2.445134137142996e+00, 3.754408661907416e+00)

.phi_inv <- function(p) {
  p <- min(max(p, .TAILS_EPS), 1.0 - .TAILS_EPS)
  A <- .ACK_A; B <- .ACK_B; C <- .ACK_C; D <- .ACK_D
  if (p < 0.02425) {
    q <- sqrt(-2.0 * log(p))
    x <- ((((((C[1] * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) * q + C[6]) /
      ((((D[1] * q + D[2]) * q + D[3]) * q + D[4]) * q + 1.0))
  } else if (p <= 0.97575) {
    q <- p - 0.5
    r <- q * q
    x <- ((((((A[1] * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * r + A[6]) * q /
      (((((B[1] * r + B[2]) * r + B[3]) * r + B[4]) * r + B[5]) * r + 1.0))
  } else {
    q <- sqrt(-2.0 * log(1.0 - p))
    x <- -((((((C[1] * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5]) * q + C[6]) /
      ((((D[1] * q + D[2]) * q + D[3]) * q + D[4]) * q + 1.0))
  }
  e <- pnorm(x) - p
  u <- e * sqrt(2.0 * pi) * exp(0.5 * x * x)
  x - u / (1.0 + 0.5 * x * u)
}

# --- GPD helpers ---

.gpd_logpdf <- function(e, gamma, sigma) {
  if (abs(gamma) < 1e-9) return(-log(sigma) - e / sigma)
  arg <- 1.0 + gamma * e / sigma
  if (arg <= 0.0) return(-745.0)
  -log(sigma) - (1.0 / gamma + 1.0) * log(arg)
}

.gpd_sf <- function(e, gamma, sigma) {
  if (e <= 0.0) return(1.0)
  if (abs(gamma) < 1e-9) return(exp(-e / sigma))
  arg <- 1.0 + gamma * e / sigma
  if (arg <= 0.0) return(0.0)
  arg^(-1.0 / gamma)
}

.gpd_isf <- function(p, gamma, sigma) {
  p <- min(max(p, 1e-300), 1.0)
  if (abs(gamma) < 1e-9) return(-sigma * log(p))
  sigma / gamma * (p^(-gamma) - 1.0)
}

.TAU_GRID <- c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.0, 8.0)

# Censored-ML GPD fit (Grimshaw profile over a fixed tau grid).
.fit_ml <- function(exc, s1) {
  n <- length(exc)
  emean <- s1 / n
  if (n < 20 || emean <= 0.0) return(c(0.0, max(emean, 1e-12)))
  emax <- max(exc)
  best_g <- 0.0; best_s <- max(emean, 1e-12); best_ll <- -1e300
  taus <- c(.TAU_GRID / emean, -0.5 / emax, -0.25 / emax, -0.1 / emax)
  for (tau in taus) {
    if (tau <= -1.0 / emax || abs(tau) < 1e-12) next
    g <- 0.0
    for (e in exc) g <- g + log1p(tau * e)
    g <- g / n
    if (g <= 1e-9) next
    sigma <- g / tau
    if (sigma <= 0.0) next
    ll <- -n * log(sigma) - (1.0 + 1.0 / g) * n * g
    if (ll > best_ll) { best_ll <- ll; best_g <- g; best_s <- sigma }
  }
  c(best_g, best_s)
}

# --- the spliced predictive ---

spliced_new <- function(body, t_lo, t_up, zeta_lo, zeta_up, g_lo, s_lo, g_up, s_up) {
  plo <- pnorm(t_lo); pup <- pnorm(t_up)
  interior <- max(pup - plo, 1e-12)
  cc <- max(1.0 - zeta_lo - zeta_up, 1e-12) / interior
  list(spliced = TRUE, body = body, t_lo = t_lo, t_up = t_up,
       zeta_lo = zeta_lo, zeta_up = zeta_up,
       g_lo = g_lo, s_lo = s_lo, g_up = g_up, s_up = s_up,
       plo = plo, pup = pup, c = cc)
}

.spliced_z <- function(d, x) {
  u <- min(max(dist_cdf(d$body, x), .TAILS_EPS), 1.0 - .TAILS_EPS)
  .phi_inv(u)
}

spliced_cdf <- function(d, x) {
  z <- .spliced_z(d, x)
  if (z < d$t_lo) return(d$zeta_lo * .gpd_sf(d$t_lo - z, d$g_lo, d$s_lo))
  if (z > d$t_up) return(1.0 - d$zeta_up * .gpd_sf(z - d$t_up, d$g_up, d$s_up))
  d$zeta_lo + d$c * (pnorm(z) - d$plo)
}

spliced_logpdf <- function(d, x) {
  base <- dist_logpdf(d$body, x)
  if (!is.finite(base)) return(base)
  z <- .spliced_z(d, x)
  corr <- if (z < d$t_lo) {
    log(max(d$zeta_lo, 1e-300)) + .gpd_logpdf(d$t_lo - z, d$g_lo, d$s_lo) -
      (-0.5 * z * z - .LOG_SQRT2PI)
  } else if (z > d$t_up) {
    log(max(d$zeta_up, 1e-300)) + .gpd_logpdf(z - d$t_up, d$g_up, d$s_up) -
      (-0.5 * z * z - .LOG_SQRT2PI)
  } else {
    log(d$c)
  }
  base + corr
}

spliced_quantile <- function(d, p, tol = 1e-9, max_iter = 100L) {
  stopifnot(p > 0, p < 1)
  if (p < d$zeta_lo) {
    z <- d$t_lo - .gpd_isf(p / d$zeta_lo, d$g_lo, d$s_lo)
  } else if (p > 1.0 - d$zeta_up) {
    z <- d$t_up + .gpd_isf((1.0 - p) / d$zeta_up, d$g_up, d$s_up)
  } else {
    u <- d$plo + (p - d$zeta_lo) / d$c
    z <- .phi_inv(min(max(u, .TAILS_EPS), 1.0 - .TAILS_EPS))
  }
  ub <- min(max(pnorm(z), .TAILS_EPS), 1.0 - .TAILS_EPS)
  dist_quantile(d$body, ub, tol = tol, max_iter = max_iter)
}

# Numeric moments and CRPS over a fixed 65-node quantile grid (midpoint rule).
.SPLICED_GRID_N <- 65L

spliced_qgrid <- function(d) {
  n <- .SPLICED_GRID_N
  vapply(seq_len(n), function(i) spliced_quantile(d, (i - 0.5) / n), 0.0)
}

spliced_mean <- function(d) { q <- spliced_qgrid(d); sum(q) / length(q) }

spliced_var <- function(d) {
  q <- spliced_qgrid(d)
  m <- sum(q) / length(q)
  sum((q - m) * (q - m)) / length(q)
}

spliced_crps <- function(d, x) {
  q <- spliced_qgrid(d)
  n <- length(q)
  t1 <- sum(abs(q - x)) / n
  t2 <- 2.0 * sum(q * (2.0 * (seq_len(n) - 0.5) / n - 1.0)) / n
  t1 - t2 * 0.5
}

# --- the wrapper ---

.GUARD_SF <- 1e-3
.ADAPT_AFTER <- 10L
.REFIT_EVERY <- 25L

.tail_new <- function() {
  list(t = NULL, exc = numeric(0), s1 = 0.0, nx = 0L, r = 0.0,
       g = 0.0, s = 1.0, since = 0L, run = 0L)
}

.tail_add <- function(tail, e, nexc) {
  # Contamination guard: winsorize intake beyond the fitted 1-in-1000 excess,
  # but a RUN of guarded ticks is a changepoint: escape after .ADAPT_AFTER.
  if (length(tail$exc) >= 20) {
    cap <- .gpd_isf(.GUARD_SF, tail$g, tail$s)
    if (e > cap) {
      if (tail$run < .ADAPT_AFTER) e <- cap
      tail$run <- tail$run + 1L
    } else {
      tail$run <- 0L
    }
  }
  tail$exc <- c(tail$exc, e)
  tail$s1 <- tail$s1 + e
  tail$nx <- tail$nx + 1L
  if (length(tail$exc) > nexc) {
    tail$s1 <- tail$s1 - tail$exc[1]
    tail$exc <- tail$exc[-1]
  }
  tail$since <- tail$since + 1L
  if (tail$since >= .REFIT_EVERY || length(tail$exc) <= 25) {
    fit <- .fit_ml(tail$exc, tail$s1)
    tail$g <- fit[1]; tail$s <- fit[2]
    tail$since <- 0L
  }
  tail
}

gpdtails <- function(base, k, level = 0.98, nexc = 500L, warmup = 500L,
                     rate_alpha = 0.002) {
  stopifnot(k >= 1, level > 0.5, level < 1.0, nexc >= 50, warmup >= 100,
            rate_alpha > 0.0, rate_alpha < 0.1)
  force(base)   # bind now: callers rebind f <- gpdtails(f, ...)
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(base = NULL, pending = list(),
                    tails = lapply(seq_len(k), function(m)
                      list(up = .tail_new(), lo = .tail_new(),
                           warm = numeric(0), n = 0L)))
    }
    # Resolve arrivals against the BODY predictions made for them.
    pend <- state$pending
    n <- length(pend)
    for (m in seq_len(k)) {
      if (m > n) next
      d <- pend[[n - m + 1]][[m]]
      u <- dist_cdf(d, y)
      if (!is.finite(u)) next
      z <- .phi_inv(min(max(u, .TAILS_EPS), 1.0 - .TAILS_EPS))
      th <- state$tails[[m]]
      up <- th$up; lo <- th$lo
      if (is.null(up$t)) {
        th$warm <- c(th$warm, z)
        if (length(th$warm) >= warmup) {
          w <- sort(th$warm)
          iu <- min(floor(level * length(w)), length(w) - 1)   # 0-based
          up$t <- w[iu + 1]
          lo$t <- w[length(w) - iu]
          for (x in w) {
            if (x > up$t) up <- .tail_add(up, x - up$t, nexc)
            else if (x < lo$t) lo <- .tail_add(lo, lo$t - x, nexc)
          }
          th$n <- length(w)
          up$r <- up$nx / length(w)
          lo$r <- lo$nx / length(w)
          th$warm <- numeric(0)
        }
      } else {
        th$n <- th$n + 1L
        up$r <- up$r + rate_alpha * ((if (z > up$t) 1.0 else 0.0) - up$r)
        lo$r <- lo$r + rate_alpha * ((if (z < lo$t) 1.0 else 0.0) - lo$r)
        if (z > up$t) up <- .tail_add(up, z - up$t, nexc)
        else if (z < lo$t) lo <- .tail_add(lo, lo$t - z, nexc)
      }
      th$up <- up; th$lo <- lo
      state$tails[[m]] <- th
    }

    r <- base(y, state$base)
    state$base <- r$state
    dists <- r$dists
    pend[[length(pend) + 1L]] <- dists
    if (length(pend) > k) pend <- pend[-1]
    state$pending <- pend

    out <- vector("list", k)
    for (m in seq_len(k)) {
      d <- dists[[m]]
      th <- state$tails[[m]]
      up <- th$up; lo <- th$lo
      if (is.null(up$t) || length(up$exc) < 8 || length(lo$exc) < 8 || th$n <= 0) {
        out[[m]] <- d
        next
      }
      out[[m]] <- spliced_new(d, lo$t, up$t, lo$r, up$r,
                              lo$g, lo$s, up$g, up$s)
    }
    list(dists = out, state = state)
  }
}
