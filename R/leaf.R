# Centered Gaussian residual leaf (Welford online variance).
# Port of skaters/leaf.py::leaf.

leaf <- function(k = 1L) {
  function(y, state = NULL) {
    if (is.null(state)) state <- list(n = 0L, mean = 0.0, m2 = 0.0)
    n <- state$n + 1L
    delta <- y - state$mean
    mean <- state$mean + delta / n
    m2 <- state$m2 + delta * (y - mean)
    state <- list(n = n, mean = mean, m2 = m2)
    v <- if (n < 2L) Inf else m2 / (n - 1)
    std <- if (is.finite(v) && v > 0) sqrt(v) else max(abs(y), 1e-8)
    d <- dist_gaussian(0.0, std)
    list(dists = rep(list(d), k), state = state)
  }
}

.SCALE_BASIS <- c(0.7, 1.0, 1.6, 3.0, 6.0)

# Residual model as a fixed Gaussian scale mixture, weights fit online
# by recency-weighted likelihood EM. Port of leaf.py::scale_mixture_leaf.
scale_mixture_leaf <- function(k = 1L, gamma = 0.02, scale_alpha = 0.01,
                               scales = .SCALE_BASIS) {
  C <- scales
  K <- length(C)
  one_idx <- which.min(abs(C - 1.0))
  function(y, state = NULL) {
    if (is.null(state)) {
      w <- rep(1e-6, K)
      w[one_idx] <- 1.0
      state <- list(v = 0.0, w = w, n = 0L)
    }
    state$n <- state$n + 1L
    a <- if (scale_alpha > 1.0 / state$n) scale_alpha else 1.0 / state$n
    state$v <- (1 - a) * state$v + a * y * y
    v <- state$v
    sigma <- if (is.finite(v) && v > 0) sqrt(v) else max(abs(y), 1e-8)
    z <- y / sigma
    w <- state$w
    dens <- w * exp(-0.5 * z * z / (C * C)) / C
    total <- sum(dens)
    if (total > 0) {
      g <- if (gamma > 1.0 / state$n) gamma else 1.0 / state$n
      state$w <- (1 - g) * w + g * dens / total
    }
    d <- dist_new(state$w, rep(0.0, K), C * sigma)
    list(dists = rep(list(d), k), state = state)
  }
}

# --- CRPS leaf ---

.A0 <- 2.0 / sqrt(2.0 * pi)   # A(0, s) = 2 phi(0) s

.abs_normal <- function(m, s) {
  if (s <= 0) return(abs(m))
  z <- m / s
  m * (2.0 * pnorm(z) - 1.0) + 2.0 * s * (exp(-0.5 * z * z) / sqrt(2.0 * pi))
}

# round(0.4 * 1.28^i, 4) for i in 0..14, precomputed (matches Python's FINE).
.FINE <- c(0.4, 0.512, 0.6554, 0.8389, 1.0737, 1.3744, 1.7592, 2.2518,
           2.8823, 3.6893, 4.7224, 6.0446, 7.7371, 9.9035, 12.6765)

# Scale-mixture leaf with weights fit by online CRPS exponentiated gradient.
# Port of leaf.py::crps_leaf.
crps_leaf <- function(k = 1L, eta = 1.0, scale_alpha = 0.01, scales = .FINE) {
  C <- scales
  K <- length(C)
  B <- outer(C, C, function(a, b) sqrt(a * a + b * b) * .A0)
  one_idx <- which.min(abs(C - 1.0))
  function(y, state = NULL) {
    if (is.null(state)) {
      w <- rep(1e-6, K)
      w[one_idx] <- 1.0
      state <- list(v = 0.0, w = w, n = 0L)
    }
    state$n <- state$n + 1L
    a <- if (scale_alpha > 1.0 / state$n) scale_alpha else 1.0 / state$n
    state$v <- (1 - a) * state$v + a * y * y
    sig <- if (is.finite(state$v) && state$v > 0) sqrt(state$v) else max(abs(y), 1e-8)
    z <- y / sig
    w <- state$w
    g <- vapply(seq_len(K), function(cc) .abs_normal(-z, C[cc]) - sum(w * B[cc, ]), 0.0)
    gm <- sum(g) / K
    e <- -eta * (g - gm)
    emax <- max(e)
    nw <- w * exp(e - emax)
    Z <- sum(nw)
    if (!(isTRUE(Z > 0.0) && is.finite(Z))) {
      nw <- w
      Z <- sum(nw)
    }
    state$w <- nw / Z
    d <- dist_new(state$w, rep(0.0, K), C * sig)
    list(dists = rep(list(d), k), state = state)
  }
}

# --- GARCH(1,1)-t leaf ---

.GARCH_A <- c(0.02, 0.04, 0.06, 0.09, 0.12, 0.16, 0.20)
.GARCH_B <- c(0.72, 0.78, 0.84, 0.88, 0.92, 0.95, 0.97)
.GARCH_OMEGA_MULT <- c(0.5, 0.7, 1.0, 1.4, 2.0)

# Terminal leaf with GARCH(1,1) conditional variance and scale-mixture tails.
# Port of leaf.py::garch_leaf (variance-targeted QMLE grid refit).
garch_leaf <- function(k = 1L, gamma = 0.02, refit_every = 40L, min_obs = 80L,
                       window = 400L, scales = .SCALE_BASIS) {
  C <- scales
  K <- length(C)
  one_idx <- which.min(abs(C - 1.0))
  function(y, state = NULL) {
    if (is.null(state)) {
      w <- rep(1e-6, K)
      w[one_idx] <- 1.0
      state <- list(h = 0.0, s2 = 0.0, n = 0L, omega = 0.0, alpha = 0.05,
                    beta = 0.90, buf = numeric(0), w = w, last_r2 = 0.0)
    }
    s <- state
    s$n <- s$n + 1L
    a0 <- if (0.02 > 1.0 / s$n) 0.02 else 1.0 / s$n
    s$s2 <- (1 - a0) * s$s2 + a0 * y * y
    if (s$s2 <= 0) s$s2 <- max(y * y, 1e-12)

    h <- if (s$n == 1L) s$s2 else s$omega + s$alpha * s$last_r2 + s$beta * s$h
    if (h <= 1e-300) h <- s$s2
    s$h <- h
    s$last_r2 <- y * y
    s$buf <- c(s$buf, y)
    if (length(s$buf) > window) s$buf <- s$buf[-1]

    if (s$n >= min_obs && s$n %% refit_every == 0L && length(s$buf) >= min_obs) {
      resid <- s$buf
      s2 <- sum(resid * resid) / length(resid)
      if (s2 > 0) {
        best_om <- NA_real_; best_al <- NA_real_; best_be <- NA_real_
        best_v <- Inf
        for (al in .GARCH_A) for (be in .GARCH_B) {
          if (al + be >= 0.999) next
          base <- (1.0 - al - be) * s2
          for (cc in .GARCH_OMEGA_MULT) {
            om <- if (base * cc > 1e-12) base * cc else 1e-12
            hh <- om / (1.0 - al - be)
            v <- 0.0
            for (r in resid) {
              hh <- om + al * (r * r) + be * hh
              if (hh <= 1e-300) hh <- 1e-300
              v <- v + log(hh) + (r * r) / hh
            }
            if (v < best_v) {   # first-wins on ties
              best_v <- v; best_om <- om; best_al <- al; best_be <- be
            }
          }
        }
        s$omega <- best_om; s$alpha <- best_al; s$beta <- best_be
      }
    }

    sigma <- if (is.finite(h) && h > 0) sqrt(h) else max(abs(y), 1e-8)
    z <- y / sigma
    w <- s$w
    dens <- w * exp(-0.5 * z * z / (C * C)) / C
    total <- sum(dens)
    if (total > 0) {
      g <- if (gamma > 1.0 / s$n) gamma else 1.0 / s$n
      s$w <- (1 - g) * w + g * dens / total
    }
    d <- dist_new(s$w, rep(0.0, K), C * sigma)
    list(dists = rep(list(d), k), state = s)
  }
}
