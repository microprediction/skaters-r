# Invertible transforms: (forward, inverse_k) pairs.
# Ports of skaters/transform.py.

difference <- function() {
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(last = y)))
    }
    list(y = y - tstate$last, state = list(last = y))
  }
  inverse_k <- function(dists, tstate) {
    anchor <- tstate$last
    out <- vector("list", length(dists))
    cm <- 0.0
    cv <- 0.0
    for (i in seq_along(dists)) {
      d <- dists[[i]]
      cm <- cm + dist_mean(d)
      cv <- cv + dist_var(d)
      std <- if (cv > 0) sqrt(cv) else max(dist_std(d), 1e-12)
      out[[i]] <- dist_gaussian(anchor + cm, std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

ema_transform <- function(alpha = 0.05) {
  stopifnot(alpha > 0, alpha < 1)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(level = y)))
    }
    residual <- y - tstate$level
    list(y = residual, state = list(level = tstate$level + alpha * residual))
  }
  inverse_k <- function(dists, tstate) {
    lapply(dists, dist_shift, delta = tstate$level)
  }
  list(forward = forward, inverse_k = inverse_k)
}

standardize <- function(alpha = 0.05, eps = 1e-8) {
  force(alpha)
  force(eps)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(mu = y, var = 0.0)))
    }
    mu <- tstate$mu
    v <- tstate$var
    diff <- y - mu
    # Emit against the PRIOR state, so the forward map is the affine change of
    # coordinates z = (y - mu) / sigma that the inverse applies. Scaling by the
    # post-update std makes the emission self-normalized (bounded by
    # 1/sqrt(alpha)) and the affine inverse is then not its inverse.
    # Cold start: before the variance is informative, scale by the first
    # nonzero residual, so the first informative emission is +-1.
    sigma <- if (v > eps * eps) {
      sqrt(v)
    } else if (abs(diff) > eps) {
      abs(diff)
    } else {
      eps
    }
    y_prime <- diff / sigma
    mu_new <- mu + alpha * diff
    v_new <- (1 - alpha) * v + alpha * diff * diff
    list(y = y_prime, state = list(mu = mu_new, var = v_new))
  }
  inverse_k <- function(dists, tstate) {
    sigma <- if (tstate$var > 1e-16) sqrt(tstate$var) else 1e-8
    lapply(dists, dist_affine, a = sigma, b = tstate$mu)
  }
  list(forward = forward, inverse_k = inverse_k)
}

ou_transform <- function(kappa = 0.1, alpha = 0.02) {
  stopifnot(kappa > 0, kappa <= 1, alpha > 0, alpha < 1)
  phi <- 1.0 - kappa
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate) || !is.finite(y)) {
      y0 <- if (is.finite(y)) y else 0.0
      return(list(y = 0.0, state = list(m = y0, fc = y0, y = y0)))
    }
    resid <- y - tstate$fc
    if (!is.finite(resid)) {
      resid <- 0.0
    }
    m <- tstate$m + alpha * (y - tstate$m)
    fc <- m + phi * (y - m)
    list(y = resid, state = list(m = m, fc = fc, y = y))
  }
  inverse_k <- function(dists, tstate) {
    m <- tstate$m
    ylast <- tstate$y
    out <- vector("list", length(dists))
    for (h in seq_along(dists)) {
      center <- m + (phi^h) * (ylast - m)
      g <- if (phi < 1.0 - 1e-9) sqrt((1.0 - phi^(2 * h)) / (1.0 - phi * phi)) else sqrt(h)
      out[[h]] <- dist_shift(dist_scale(dists[[h]], g), center)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

theta <- function(alpha = 0.1) {
  stopifnot(alpha > 0, alpha < 1)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(ses = y, t = 1, sum_t = 1.0, sum_t2 = 1.0, sum_y = y, sum_ty = y)))
    }
    s <- tstate
    s$t <- s$t + 1
    t <- s$t
    slope_prev <- if (is.null(s$slope)) 0.0 else s$slope
    forecast <- s$ses + slope_prev / 2
    residual <- y - forecast
    s$ses <- alpha * y + (1 - alpha) * s$ses
    s$sum_t <- s$sum_t + t
    s$sum_t2 <- s$sum_t2 + t * t
    s$sum_y <- s$sum_y + y
    s$sum_ty <- s$sum_ty + t * y
    n <- t
    denom <- n * s$sum_t2 - s$sum_t^2
    s$slope <- if (abs(denom) > 1e-12) (n * s$sum_ty - s$sum_t * s$sum_y) / denom else 0.0
    list(y = residual, state = s)
  }
  inverse_k <- function(dists, tstate) {
    ses <- tstate$ses
    slope <- if (is.null(tstate$slope)) 0.0 else tstate$slope
    out <- vector("list", length(dists))
    cv <- 0.0
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      cv <- cv + dist_var(d)
      forecast <- ses + h * slope / 2 + dist_mean(d)
      std <- if (cv > 0) sqrt(cv) else max(dist_std(d), 1e-12)
      out[[h]] <- dist_gaussian(forecast, std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

drift <- function(alpha = 0.002, shrinkage = 0.001) {
  stopifnot(alpha > 0, alpha < 1, shrinkage >= 0, shrinkage < 1)
  decay <- 1 - alpha - shrinkage
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(last = y, mu = 0.0)))
    }
    dy <- y - tstate$last
    residual <- dy - tstate$mu
    mu <- decay * tstate$mu + alpha * dy
    list(y = residual, state = list(last = y, mu = mu))
  }
  inverse_k <- function(dists, tstate) {
    anchor <- tstate$last
    mu <- tstate$mu
    out <- vector("list", length(dists))
    cm <- 0.0
    cv <- 0.0
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      cm <- cm + dist_mean(d)
      cv <- cv + dist_var(d)
      total_mean <- anchor + h * mu + cm
      std <- if (cv > 0) sqrt(cv) else max(dist_std(d), 1e-12)
      out[[h]] <- dist_gaussian(total_mean, std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

holt_linear <- function(alpha = 0.1, beta = 0.05) {
  stopifnot(alpha > 0, alpha < 1, beta > 0, beta < 1)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(level = y, trend = 0.0)))
    }
    l_prev <- tstate$level
    b_prev <- tstate$trend
    l_new <- alpha * y + (1 - alpha) * (l_prev + b_prev)
    b_new <- beta * (l_new - l_prev) + (1 - beta) * b_prev
    residual <- y - (l_prev + b_prev)
    list(y = residual, state = list(level = l_new, trend = b_new))
  }
  inverse_k <- function(dists, tstate) {
    level <- tstate$level
    trend <- tstate$trend
    out <- vector("list", length(dists))
    cv <- 0.0
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      cv <- cv + dist_var(d)
      forecast <- level + h * trend + dist_mean(d)
      std <- if (cv > 0) sqrt(cv) else max(dist_std(d), 1e-12)
      out[[h]] <- dist_gaussian(forecast, std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

garch <- function(omega = 0.01, alpha = 0.1, beta = 0.85, mean_alpha = 0.05,
                  eps = 1e-8) {
  stopifnot(omega > 0, alpha >= 0, beta >= 0)
  force(eps)
  force(mean_alpha)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      persist <- alpha + beta
      var0 <- if (persist < 1) omega / (1 - persist) else omega / eps
      return(list(y = 0.0, state = list(var = var0, last_dev = 0.0, mu = y)))
    }
    # GARCH conditional variance of the DEVIATION from a running mean, not of
    # the raw value. alpha * y^2 treats y as a mean-zero return; on a level
    # series that makes the "volatility" of order |y|, and the inverse then
    # re-inflates it -- unstable, and it widened real level-series forecasts.
    # Tracking a mean (like `standardize`) makes garch shift-invariant. On the
    # mean-zero returns garch is meant for, mu stays ~0 and this reduces to the
    # original y/sigma exactly.
    mu <- tstate$mu
    dev <- y - mu
    v <- omega + alpha * tstate$last_dev^2 + beta * tstate$var
    sigma <- if (v > eps * eps) sqrt(v) else eps
    y_prime <- dev / sigma
    mu_new <- mu + mean_alpha * dev
    list(y = y_prime, state = list(var = v, last_dev = dev, mu = mu_new))
  }
  inverse_k <- function(dists, tstate) {
    sigma <- if (tstate$var > 1e-16) sqrt(tstate$var) else 1e-8
    lapply(dists, dist_affine, a = sigma, b = tstate$mu)
  }
  list(forward = forward, inverse_k = inverse_k)
}

seasonal_difference <- function(period = 12L) {
  stopifnot(period >= 1)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(buffer = y)))
    }
    buf <- tstate$buffer
    y_prime <- if (length(buf) >= period) y - buf[length(buf) - period + 1] else 0.0
    buf <- c(buf, y)
    if (length(buf) > 2 * period) {
      buf <- buf[-1]
    }
    list(y = y_prime, state = list(buffer = buf))
  }
  inverse_k <- function(dists, tstate) {
    buf <- tstate$buffer
    k <- length(dists)
    recovered_means <- numeric(k)
    recovered_vars <- numeric(k)
    out <- vector("list", k)
    for (h in seq_len(k)) {
      lag_idx <- (h - 1) - period
      if (lag_idx < 0) {
        buf_idx <- length(buf) - period + h # 1-based
        anchor_mean <- if (buf_idx >= 1 && buf_idx <= length(buf)) buf[buf_idx] else 0.0
        anchor_var <- 0.0
      } else {
        anchor_mean <- recovered_means[lag_idx + 1]
        anchor_var <- recovered_vars[lag_idx + 1]
      }
      d <- dists[[h]]
      recovered_means[h] <- dist_mean(d) + anchor_mean
      recovered_vars[h] <- dist_var(d) + anchor_var
      out[[h]] <- if (anchor_var > 0.0) {
        dist_new(d$w, d$m + anchor_mean, sqrt(d$s * d$s + anchor_var))
      } else {
        dist_shift(d, anchor_mean)
      }
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

# Residual from a hedged seasonal anchor.
#
# Forward:   y'_t = y_t - a_t,  a_t = weight * phaseEMA_{p(t)} + (1-weight) * y_{t-s}
# Inverse:   shift each horizon's Dist by its anchor; for h >= s the naive
#            component is a value recovered earlier in the same call and its
#            variance is convolved in (as in seasonal_difference).
#
# The phase-EMA is a recency-weighted mean of same-phase values (memory ~1/alpha
# cycles); the seasonal-naive is the single value one period ago. The naive alone
# adapts instantly but is one noisy draw; the phase-EMA averages the noise but
# lags level shifts. The blend beats either alone on seasonal series.
# weight = 0 recovers seasonal_difference.
#
# Forecasting from same-phase component series follows Viole's NNS package
# (NNS.ARMA, CRAN, since 2017).
seasonal_anchor <- function(period, alpha = 0.2, weight = 0.5) {
  stopifnot(period >= 1, alpha > 0, alpha < 1, weight >= 0, weight <= 1)
  force(alpha)
  force(weight)
  anchor_of <- function(ema_p, snaive) {
    if (is.null(ema_p)) snaive else weight * ema_p + (1.0 - weight) * snaive
  }
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(ema = vector("list", period),
                                        buffer = c(y), n = 1)))
    }
    st <- tstate
    buf <- st$buffer
    p <- st$n %% period # 0-based phase, as in Python
    snaive <- if (length(buf) >= period) buf[length(buf) - period + 1] else buf[length(buf)]
    y_prime <- y - anchor_of(st$ema[[p + 1]], snaive)
    e <- st$ema[[p + 1]]
    st$ema[[p + 1]] <- if (is.null(e)) y else e + alpha * (y - e)
    buf <- c(buf, y)
    if (length(buf) > 2 * period) {
      buf <- buf[-1]
    }
    st$buffer <- buf
    st$n <- st$n + 1
    list(y = y_prime, state = st)
  }
  inverse_k <- function(dists, tstate) {
    buf <- tstate$buffer
    k <- length(dists)
    recovered_means <- numeric(k)
    recovered_vars <- numeric(k)
    out <- vector("list", k)
    for (h in seq_len(k)) {
      hh <- h - 1 # 0-based horizon
      p <- (tstate$n + hh) %% period
      lag_idx <- hh - period
      if (lag_idx < 0) {
        buf_idx <- length(buf) - period + hh + 1 # 1-based
        snaive <- if (buf_idx >= 1 && buf_idx <= length(buf)) {
          buf[buf_idx]
        } else {
          buf[length(buf)]
        }
        snaive_var <- 0.0
      } else {
        snaive <- recovered_means[lag_idx + 1]
        snaive_var <- recovered_vars[lag_idx + 1]
      }
      a_mean <- anchor_of(tstate$ema[[p + 1]], snaive)
      a_var <- ((1.0 - weight)^2) * snaive_var
      d <- dists[[h]]
      recovered_means[h] <- dist_mean(d) + a_mean
      recovered_vars[h] <- dist_var(d) + a_var
      out[[h]] <- if (a_var > 0.0) {
        dist_new(d$w, d$m + a_mean, sqrt(d$s * d$s + a_var))
      } else {
        dist_shift(d, a_mean)
      }
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

power_transform <- function(p = 0.5) {
  stopifnot(p > 0, p < 1)
  inv_p <- 1.0 / p
  pw_fwd <- function(y) if (y < 0) -abs(y)^p else abs(y)^p
  forward <- function(y, tstate = NULL) {
    list(y = pw_fwd(y), state = if (is.null(tstate)) list() else tstate)
  }
  inverse_k <- function(dists, tstate) {
    out <- vector("list", length(dists))
    for (i in seq_along(dists)) {
      d <- dists[[i]]
      n <- length(d$w)
      m_out <- numeric(n)
      s_out <- numeric(n)
      for (j in seq_len(n)) {
        mu <- d$m[j]
        sigma <- d$s[j]
        orig_mean <- if (mu < 0) -abs(mu)^inv_p else abs(mu)^inv_p
        abs_mu <- abs(mu)
        deriv <- if (abs_mu > 1e-12) inv_p * abs_mu^(inv_p - 1) else inv_p
        m_out[j] <- orig_mean
        s_out[j] <- max(sigma * deriv, 1e-12)
      }
      out[[i]] <- dist_new(d$w, m_out, s_out)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

yeo_johnson <- function(lmbda = 0.0) {
  L <- as.numeric(lmbda)
  yj_fwd <- function(y) {
    if (y >= 0.0) {
      if (L == 0.0) {
        return(log1p(y))
      }
      return(((y + 1.0)^L - 1.0) / L)
    }
    if (L == 2.0) {
      return(-log1p(-y))
    }
    -(((-y + 1.0)^(2.0 - L) - 1.0) / (2.0 - L))
  }
  yj_inv <- function(yp) {
    if (yp >= 0.0) {
      if (L == 0.0) {
        return(expm1(min(yp, 350.0)))
      }
      base <- max(L * yp + 1.0, 1e-12)
      return(base^(1.0 / L) - 1.0)
    }
    if (L == 2.0) {
      return(1.0 - exp(min(-yp, 350.0)))
    }
    base <- max(-(2.0 - L) * yp + 1.0, 1e-12)
    1.0 - base^(1.0 / (2.0 - L))
  }
  yj_dinv <- function(yp) {
    if (yp >= 0.0) {
      if (L == 0.0) {
        return(exp(min(yp, 350.0)))
      }
      base <- max(L * yp + 1.0, 1e-12)
      return(base^(1.0 / L - 1.0))
    }
    if (L == 2.0) {
      return(exp(min(-yp, 350.0)))
    }
    base <- max(-(2.0 - L) * yp + 1.0, 1e-12)
    base^(1.0 / (2.0 - L) - 1.0)
  }
  forward <- function(y, tstate = NULL) {
    list(y = yj_fwd(y), state = if (is.null(tstate)) list() else tstate)
  }
  inverse_k <- function(dists, tstate) {
    out <- vector("list", length(dists))
    for (i in seq_along(dists)) {
      d <- dists[[i]]
      m_out <- vapply(d$m, yj_inv, 0.0)
      s_out <- pmax(d$s * vapply(d$m, yj_dinv, 0.0), 1e-12)
      out[[i]] <- dist_new(d$w, m_out, s_out)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

# --- Fractional differencing ---

.frac_diff_weights <- function(d, window) {
  w <- numeric(window)
  w[1] <- 1.0
  if (window > 1) {
    for (i in seq_len(window - 1)) {
      w[i + 1] <- -w[i] * (d - i + 1) / i
    }
  }
  w
}

fractional_difference <- function(d = 0.4, window = 50L) {
  w_fwd <- .frac_diff_weights(d, window)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      tstate <- list(buffer = numeric(0))
    }
    buf <- c(tstate$buffer, y)
    if (length(buf) > window) {
      buf <- buf[-1]
    }
    n <- length(buf)
    y_prime <- sum(w_fwd[seq_len(n)] * buf[n:1])
    list(y = y_prime, state = list(buffer = buf))
  }
  inverse_k <- function(dists, tstate) {
    buf <- tstate$buffer
    out <- vector("list", length(dists))
    for (i in seq_along(dists)) {
      d_in <- dists[[i]]
      buf <- c(buf, 0.0)
      n <- length(buf)
      shift <- 0.0
      jmax <- min(n, window) - 1
      if (jmax >= 1) {
        for (j in seq_len(jmax)) {
          shift <- shift - w_fwd[j + 1] * buf[n - j]
        }
      }
      recovered <- dist_mean(d_in) + shift
      buf[n] <- recovered
      out[[i]] <- dist_gaussian(recovered, dist_std(d_in))
      if (length(buf) > window) buf <- buf[-1]
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

# --- AR(p) with online recursive least squares ---

.rls_mat_vec <- function(M, v, n) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- sum(M[i, ] * v)
  }
  out
}

# Spectral radius of the AR companion matrix (largest |characteristic root|).
# The h-step AR forecast is governed by powers of this matrix: it converges
# (a well-posed forecast) iff the radius is < 1, and diverges geometrically
# otherwise. Closed form for p <= 2; power iteration for higher orders.
# Port of _ar_spectral_radius in transform.py -- the 60 fixed iterations and
# the "or 1.0" zero-guard are load-bearing for 1e-6 parity.
.ar_spectral_radius <- function(phi) {
  p <- length(phi)
  if (p == 1) {
    return(abs(phi[1]))
  }
  if (p == 2) {
    a <- phi[1]
    b <- phi[2]
    disc <- a * a + 4.0 * b
    if (disc >= 0.0) {
      r <- sqrt(disc)
      return(max(abs((a + r) / 2.0), abs((a - r) / 2.0)))
    }
    return(sqrt((a / 2.0)^2 + (sqrt(-disc) / 2.0)^2)) # |complex root|
  }
  v <- rep(1.0, p)
  rho <- 0.0
  for (.i in seq_len(60)) {
    nv <- c(sum(phi * v), v[-p])
    m <- max(abs(nv))
    if (m == 0) m <- 1.0
    v <- nv / m
    rho <- m
  }
  rho
}

# Damp AR coefficients into the stationary region for forecasting.
#
# An online least-squares fit -- especially from a handful of warm-up points --
# can land outside the stationary region (a near-unit or explosive root). Its
# multi-step forecast then diverges (an AR(2) fit to 3 points reached a 13-step
# mean of ~1e22 on real GIFT-Eval series). Scaling phi_j by gamma^(j+1) scales
# every companion eigenvalue by gamma, so gamma = margin / rho brings the radius
# to margin and leaves an already-stationary fit untouched. Constrained
# forecasting, not a magnitude clip.
.ar_stationary <- function(phi, margin = 0.999) {
  rho <- .ar_spectral_radius(phi)
  if (rho <= margin) {
    return(phi)
  }
  g <- margin / rho
  phi * g^(seq_along(phi))
}

ar <- function(order = 2L, lam = 0.99, ridge = 1.0, decay = 0.0) {
  stopifnot(order >= 1, lam > 0, lam <= 1, decay >= 0)
  force(ridge)
  p <- order
  init_P <- function() {
    P <- matrix(0.0, p, p)
    for (j in seq_len(p)) {
      P[j, j] <- if (decay > 0) ridge / (j^decay) else ridge
    }
    P
  }
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      tstate <- list(buffer = numeric(0), phi = rep(0.0, p), P = init_P(), n = 0)
    }
    st <- tstate
    st$n <- st$n + 1
    buf <- st$buffer
    if (length(buf) >= p) {
      x <- buf[seq(length(buf), by = -1L, length.out = p)]
      prediction <- sum(st$phi * x)
      residual <- y - prediction
      P <- st$P
      Px <- .rls_mat_vec(P, x, p)
      denom <- lam + sum(x * Px)
      if (abs(denom) > 1e-15) {
        K <- Px / denom
        st$phi <- st$phi + K * residual
        P <- (P - K %o% Px) / lam
        if (!all(is.finite(P)) || max(abs(P)) > 1e10) {
          P <- init_P()
        }
        st$P <- P
      }
    } else {
      residual <- y
    }
    buf <- c(buf, y)
    if (length(buf) > 2 * p + 10) {
      buf <- buf[-1]
    }
    st$buffer <- buf
    list(y = residual, state = st)
  }
  inverse_k <- function(dists, tstate) {
    buf <- tstate$buffer
    # Constrain to stationarity: an online fit from few warm-up points can land
    # on an explosive root whose multi-step forecast diverges.
    phi <- .ar_stationary(tstate$phi)
    k <- length(dists)
    # Impulse responses psi_0..psi_{k-1}: psi_0 = 1, psi_i = sum_j phi_j psi_{i-j}.
    # The h-step variance is sigma^2 * sum_{i<=h} psi_i^2, NOT a recursion on the
    # recovered variances -- successive horizons share the same innovation, so
    # summing phi_j^2 * var_{h-j} both mis-weights it and double counts.
    psi <- numeric(k)
    psi[1] <- 1.0
    for (i in seq_len(k - 1)) {
      s <- 0.0
      for (j in seq_len(p)) {
        idx <- i - j # 0-based psi index
        if (idx >= 0) s <- s + phi[j] * psi[idx + 1]
      }
      psi[i + 1] <- s
    }
    recovered_means <- numeric(k)
    out <- vector("list", k)
    cum_psi2 <- 0.0
    for (h in seq_len(k)) {
      ar_mean <- 0.0
      for (j in seq_len(p)) {
        lag_h <- (h - 1) - j # 0-based previous horizon index
        if (lag_h < 0) {
          buf_idx <- length(buf) + lag_h + 1
          if (buf_idx >= 1 && buf_idx <= length(buf)) {
            ar_mean <- ar_mean + phi[j] * buf[buf_idx]
          }
        } else if (lag_h < h - 1) {
          ar_mean <- ar_mean + phi[j] * recovered_means[lag_h + 1]
        }
      }
      d <- dists[[h]]
      total_mean <- dist_mean(d) + ar_mean
      cum_psi2 <- cum_psi2 + psi[h] * psi[h]
      total_var <- cum_psi2 * dist_var(d)
      total_std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
      recovered_means[h] <- total_mean
      out[[h]] <- dist_gaussian(total_mean, total_std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

.build_groups <- function(max_lag) {
  groups <- integer(0)
  g <- 1L
  size <- 1L
  assigned <- 0L
  while (assigned < max_lag) {
    for (i in seq_len(size)) {
      if (assigned >= max_lag) {
        break
      }
      groups <- c(groups, g)
      assigned <- assigned + 1L
    }
    g <- g + 1L
    size <- size * 2L
  }
  groups
}

grouped_ar <- function(max_lag = 16L, lam = 0.99, ridge = 1.0) {
  stopifnot(max_lag >= 1, lam > 0, lam <= 1, ridge > 0)
  groups <- .build_groups(max_lag) # 1-based group id per lag j = 1..max_lag
  n_groups <- max(groups)
  init_P <- function() diag(ridge, n_groups, n_groups)
  group_regressor <- function(buf) {
    x <- numeric(n_groups)
    for (j in seq_len(max_lag)) {
      g <- groups[j]
      x[g] <- x[g] + buf[length(buf) - j + 1]
    }
    x
  }
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      tstate <- list(buffer = numeric(0), theta = rep(0.0, n_groups), P = init_P(), n = 0)
    }
    st <- tstate
    st$n <- st$n + 1
    buf <- st$buffer
    if (length(buf) >= max_lag) {
      x <- group_regressor(buf)
      prediction <- sum(st$theta * x)
      residual <- y - prediction
      P <- st$P
      Px <- .rls_mat_vec(P, x, n_groups)
      denom <- lam + sum(x * Px)
      if (abs(denom) > 1e-15) {
        K <- Px / denom
        st$theta <- st$theta + K * residual
        P <- (P - K %o% Px) / lam
        if (!all(is.finite(P)) || max(abs(P)) > 1e10) {
          P <- init_P()
        }
        st$P <- P
      }
    } else {
      residual <- y
    }
    buf <- c(buf, y)
    if (length(buf) > max_lag + 10) {
      buf <- buf[-1]
    }
    st$buffer <- buf
    list(y = residual, state = st)
  }
  inverse_k <- function(dists, tstate) {
    buf <- tstate$buffer
    # Expand group coefficients to per-lag, then constrain to stationarity for a
    # well-posed forecast (same as ar(); grouped AR is AR-family).
    phi <- .ar_stationary(tstate$theta[groups])
    k <- length(dists)
    # Impulse responses, as in ar(): the h-step variance is
    # sigma^2 * sum_{i<=h} psi_i^2, not a recursion on recovered variances.
    psi <- numeric(k)
    psi[1] <- 1.0
    for (i in seq_len(k - 1)) {
      s <- 0.0
      for (j in seq_len(max_lag)) {
        idx <- i - j
        if (idx >= 0) s <- s + phi[j] * psi[idx + 1]
      }
      psi[i + 1] <- s
    }
    recovered_means <- numeric(k)
    out <- vector("list", k)
    cum_psi2 <- 0.0
    for (h in seq_len(k)) {
      ar_mean <- 0.0
      for (j in seq_len(max_lag)) {
        lag_h <- (h - 1) - j
        if (lag_h < 0) {
          buf_idx <- length(buf) + lag_h + 1
          if (buf_idx >= 1 && buf_idx <= length(buf)) {
            ar_mean <- ar_mean + phi[j] * buf[buf_idx]
          }
        } else if (lag_h < h - 1) {
          ar_mean <- ar_mean + phi[j] * recovered_means[lag_h + 1]
        }
      }
      d <- dists[[h]]
      total_mean <- dist_mean(d) + ar_mean
      cum_psi2 <- cum_psi2 + psi[h] * psi[h]
      total_var <- cum_psi2 * dist_var(d)
      total_std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
      recovered_means[h] <- total_mean
      out[[h]] <- dist_gaussian(total_mean, total_std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}
