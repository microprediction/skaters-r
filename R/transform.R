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
    # Shift by the level and convolve in the level's own innovation response:
    # Var(y_{t+h+1}) = var_h + alpha^2 sum_{i<h} var_i (ETS(A,N,N); skaters#254).
    # Component-wise, so the leaf's shape is kept.
    level <- tstate$level
    out <- vector("list", length(dists))
    add <- 0.0
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      out[[h]] <- if (add > 0.0) {
        .widen(d, level, add)
      } else {
        dist_shift(d, level)
      }
      add <- add + alpha * alpha * dist_var(d)
    }
    out
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
      # Var(y_{t+h}) = sum_{i<h} phi^{2(h-1-i)} var_i: each future innovation
      # decays from its own time (skaters#256), as a scale on this horizon's Dist.
      v_h <- dist_var(dists[[h]])
      g <- if (v_h > 0.0) {
        tot <- 0.0
        for (i in 0:(h - 1)) {
          tot <- tot + (phi^(2 * (h - 1 - i))) * dist_var(dists[[i + 1]])
        }
        sqrt(tot / v_h)
      } else if (phi < 1.0 - 1e-9) {
        sqrt((1.0 - phi^(2 * h)) / (1.0 - phi * phi))
      } else {
        sqrt(h)
      }
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
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      forecast <- ses + h * slope / 2 + dist_mean(d)
      # SES response: c_0 = 1, c_j = alpha; Var = sum_j c_j^2 var_{h-j}
      # (skaters#253). The running-OLS slope's O(1/t) response is dropped.
      total_var <- 0.0
      for (j in 0:(h - 1)) {
        cc <- if (j == 0) 1.0 else alpha
        total_var <- total_var + cc * cc * dist_var(dists[[h - j]])
      }
      std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
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
    rho <- decay + alpha # 1 - shrinkage: persistence of the drift estimate
    out <- vector("list", length(dists))
    cm <- 0.0
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      cm <- cm + dist_mean(d)
      total_mean <- anchor + h * mu + cm
      # The innovation at t+i+1 reaches y_{t+h+1} directly and through the
      # drift estimate on each later step: c = 1 + alpha * sum_{m<h-i} rho^m
      # (skaters#257); i is 0-based, h here is 1-based.
      total_var <- 0.0
      for (i in 0:(h - 1)) {
        acc <- 0.0
        f <- 1.0
        for (.m in seq_len(h - 1 - i)) {
          acc <- acc + f
          f <- f * rho
        }
        cc <- 1.0 + alpha * acc
        total_var <- total_var + cc * cc * dist_var(dists[[i + 1]])
      }
      std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
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
    for (h in seq_along(dists)) {
      d <- dists[[h]]
      forecast <- level + h * trend + dist_mean(d)
      # ETS(A,A,N): c_0 = 1, c_j = alpha + alpha*beta*j;
      # Var = sum_j c_j^2 var_{h-j} (skaters#242). j is 0-based as in Python.
      total_var <- 0.0
      for (j in 0:(h - 1)) {
        cc <- if (j == 0) 1.0 else alpha + alpha * beta * j
        total_var <- total_var + cc * cc * dist_var(dists[[h - j]])
      }
      std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
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

# Shift by delta and add independent Gaussian variance add. Exact for a
# Gaussian mixture (each component widens). A tail-spliced predictive has no
# components to widen, so its body widens and the GPD splice, which lives in
# the body's PIT space, follows it. Port of _widen in transform.py.
.widen <- function(d, delta, add) {
  if (isTRUE(d$spliced)) {
    return(spliced_new(.widen(d$body, delta, add), d$t_lo, d$t_up, d$zeta_lo, d$zeta_up,
                       d$g_lo, d$s_lo, d$g_up, d$s_up))
  }
  cm <- dist_components(d)
  dist_new(cm$w, cm$m + delta, sqrt(cm$s * cm$s + add))
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
        .widen(d, anchor_mean, anchor_var)
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
    # Innovation-response rows for each recovered value and each phase's EMA:
    # the forward updates a phase's EMA from recovered values too, so an anchor
    # at h >= period responds to the recovered value one period back through
    # both the naive term and the EMA term (skaters#258).
    rec_resp <- vector("list", k)
    ema_resp <- rep(list(numeric(k)), period)
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
        sn_resp <- NULL
      } else {
        snaive <- recovered_means[lag_idx + 1]
        sn_resp <- rec_resp[[lag_idx + 1]]
      }
      e <- tstate$ema[[p + 1]]
      a_mean <- anchor_of(e, snaive)
      a_resp <- numeric(k)
      for (i in seq_len(hh)) {
        sn <- if (!is.null(sn_resp)) sn_resp[i] else 0.0
        a_resp[i] <- if (!is.null(e)) weight * ema_resp[[p + 1]][i] + (1.0 - weight) * sn else sn
      }
      a_var <- 0.0
      for (i in seq_len(hh)) {
        a_var <- a_var + a_resp[i] * a_resp[i] * dist_var(dists[[i]])
      }
      r <- a_resp
      r[h] <- 1.0
      d <- dists[[h]]
      recovered_means[h] <- dist_mean(d) + a_mean
      rec_resp[[h]] <- r
      ema_resp[[p + 1]] <- if (is.null(e)) r else (1.0 - alpha) * ema_resp[[p + 1]] + alpha * r
      out[[h]] <- if (a_var > 0.0) {
        .widen(d, a_mean, a_var)
      } else {
        dist_shift(d, a_mean)
      }
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

# Phase-indexed running standardization (seasonal heteroskedasticity).
#
# Forward:   y'_t = (y_t - mu_{p(t)}) / sigma_{p(t)},   p(t) = t mod period
# Inverse:   affine by the stats of the phase being forecast: dists[h] is the
#            h-step-ahead predictive, so its phase is (n + h - 1) mod period
#            where n counts observations seen.
#
# Each phase keeps its own EWMA mean and variance. center = FALSE divides by
# the phase sigma without subtracting the phase mean; compose it over
# seasonal_difference(period). A phase's variance is shrunk toward the global
# EWMA variance in log space by visit count, w = visits / (visits + 3). The
# global stats run at alpha / period.
seasonal_scale <- function(period, alpha = 0.05, center = FALSE, eps = 1e-8) {
  stopifnot(period >= 1, alpha > 0, alpha < 1)
  force(center)
  force(eps)
  a_g <- alpha / period
  K0 <- 3.0 # shrinkage prior strength (visits)

  # Bias-corrected global variance (the slow EWMA starts at zero).
  vg_of <- function(st) {
    updates <- st$n - 1
    if (updates <= 0) {
      return(0.0)
    }
    st$var_g / (1.0 - (1.0 - a_g)^updates)
  }
  # Shrunk (mu, sigma) for 0-based phase p.
  phase_stats <- function(st, p) {
    w <- st$vis[p + 1] / (st$vis[p + 1] + K0)
    vp <- max(st$var[p + 1], eps * eps)
    vg <- max(vg_of(st), eps * eps)
    var_eff <- exp(w * log(vp) + (1.0 - w) * log(vg))
    mu_eff <- if (center) (w * st$mu[p + 1] + (1.0 - w) * st$mu_g) else 0.0
    c(mu_eff, sqrt(var_eff))
  }
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) {
      return(list(y = 0.0, state = list(
        mu_g = y, var_g = 0.0, mu = numeric(period), var = numeric(period),
        vis = integer(period), n = 1L
      )))
    }
    st <- tstate
    p <- st$n %% period
    # Center against the prior mean, scale by the updated variance.
    mu_eff <- phase_stats(st, p)[1]
    dg <- y - st$mu_g
    st$mu_g <- st$mu_g + a_g * dg
    st$var_g <- (1 - a_g) * st$var_g + a_g * dg * dg
    if (st$vis[p + 1] == 0) { # first visit seeds from global stats
      st$mu[p + 1] <- y
      st$var[p + 1] <- vg_of(st)
    } else {
      dp <- y - st$mu[p + 1]
      st$mu[p + 1] <- st$mu[p + 1] + alpha * dp
      st$var[p + 1] <- (1 - alpha) * st$var[p + 1] + alpha * dp * dp
    }
    st$vis[p + 1] <- st$vis[p + 1] + 1L
    st$n <- st$n + 1L
    sigma <- phase_stats(st, p)[2]
    list(y = (y - mu_eff) / sigma, state = st)
  }
  inverse_k <- function(dists, tstate) {
    out <- vector("list", length(dists))
    for (h in seq_along(dists)) {
      ms <- phase_stats(tstate, (tstate$n + h - 1) %% period)
      out[[h]] <- dist_affine(dists[[h]], ms[2], ms[1])
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

power_transform <- function(p = 0.5, exact = FALSE) {
  stopifnot(p > 0, p < 1)
  inv_p <- 1.0 / p
  pw_fwd <- function(y) if (y < 0) -abs(y)^p else abs(y)^p
  pw_inv <- function(z) if (z < 0) -abs(z)^inv_p else abs(z)^inv_p
  forward <- function(y, tstate = NULL) {
    list(y = pw_fwd(y), state = if (is.null(tstate)) list() else tstate)
  }
  inverse_k <- function(dists, tstate) {
    out <- vector("list", length(dists))
    for (i in seq_along(dists)) {
      d <- dist_components(dists[[i]])
      if (is.null(d)) { # tail-spliced: exact pushforward
        out[[i]] <- .pushforward(dists[[i]], pw_inv)
        next
      }
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
  list(forward = forward, inverse_k = if (exact) .exact_inverse_k(pw_inv) else inverse_k)
}

yeo_johnson <- function(lmbda = 0.0, exact = FALSE) {
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
      d <- dist_components(dists[[i]])
      if (is.null(d)) { # tail-spliced: exact pushforward
        out[[i]] <- .pushforward(dists[[i]], yj_inv)
        next
      }
      m_out <- vapply(d$m, yj_inv, 0.0)
      s_out <- pmax(d$s * vapply(d$m, yj_dinv, 0.0), 1e-12)
      out[[i]] <- dist_new(d$w, m_out, s_out)
    }
    out
  }
  list(forward = forward, inverse_k = if (exact) .exact_inverse_k(yj_inv) else inverse_k)
}

# Log coordinate for a series with a one-sided bound.
#
# Forward:   y' = log(y - lower)        (lower bound; the default, lower = 0)
#            y' = -log(upper - y)       (upper bound)
# Inverse:   the exact pushforward through the inverse map, so the predictive
#            in y keeps the skew of the coordinate change, its logpdf includes
#            the Jacobian, and every quantile lies on the right side of the
#            bound.
#
# A finite observation on or past the bound is an error. Missing observations
# pass through. Quantiles that round to the bound are rounded inward by one
# floating-point step. Experimental, alongside flesaker.
log_transform <- function(lower = NULL, upper = NULL) {
  if (!is.null(lower) && !is.null(upper)) {
    stop("log_transform takes one bound: lower or upper")
  }
  b <- as.numeric(if (!is.null(upper)) upper else if (is.null(lower)) 0.0 else lower)
  if (!is.finite(b)) {
    stop("log_transform needs a finite bound")
  }
  direction <- if (!is.null(upper)) -1.0 else 1.0
  bound_label <- if (!is.null(upper)) "upper" else "lower"
  log_two <- log(2.0)

  fwd_query <- function(y) {
    gap <- direction * (y - b)
    if (gap <= 0.0) {
      return(-direction * Inf)
    }
    if (is.infinite(gap)) {
      # A gap between finite floats can overflow even when its log fits.
      half_gap <- direction * (0.5 * y - 0.5 * b)
      return(direction * (log(half_gap) + log_two))
    }
    direction * log(gap)
  }
  inverse_at_scale <- function(z, weight) {
    bound <- weight * b
    log_gap <- direction * z + log(weight)
    gap <- exp(log_gap)
    if (is.infinite(gap)) {
      # The shift by the bound can bring an overflowing gap into range.
      return(2.0 * (0.5 * bound + direction * exp(log_gap - log_two)))
    }
    bound + direction * gap
  }
  lt_inv <- function(z) {
    y <- inverse_at_scale(z, 1.0)
    if (isTRUE(y == b)) .pf_nextafter(b, direction) else y
  }
  log_jacobian <- function(z) -direction * z
  moments <- function(m, s) {
    variance <- s * s
    log_mean_gap <- direction * m + 0.5 * variance
    mean <- inverse_at_scale(direction * log_mean_gap, 1.0)
    if (s == 0.0) {
      return(c(mean, 0.0))
    }
    # log(sqrt(expm1(s^2))) stays stable as the spread approaches zero or
    # grows large, where the correction to s^2 / 2 rounds to zero.
    log_relative_std <- if (variance < 1e-16) {
      log(s)
    } else if (variance > 40.0) {
      0.5 * variance
    } else {
      0.5 * log(expm1(variance))
    }
    c(mean, exp(log_mean_gap + log_relative_std))
  }
  forward <- function(y, tstate = NULL) {
    st <- if (is.null(tstate)) list() else tstate
    if (!is.finite(y)) {
      return(list(y = y, state = st))
    }
    z <- fwd_query(y)
    if (z == -direction * Inf) {
      stop(sprintf("log_transform: observation %s must lie inside %s=%s", format(y), bound_label, format(b)))
    }
    list(y = z, state = st)
  }
  inverse_k <- function(dists, tstate) {
    lapply(dists, function(d) {
      pushforward_new(d,
        inv_fn = lt_inv, fwd_fn = fwd_query, log_jacobian_fn = log_jacobian,
        component_moments_fn = moments, scaled_inv_fn = inverse_at_scale
      )
    })
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
    # Recovered earlier horizons are uncertain and feed back through the same
    # weights (skaters#255): each recovered value carries a row of responses to
    # the future innovations, built by the same recursion as its mean, and
    # Var(y_{t+h+1}) = sum_i c_{h,i}^2 var_i. NULL rows mark known values.
    H <- length(dists)
    buf <- tstate$buffer
    cbuf <- vector("list", length(buf))
    out <- vector("list", H)
    for (i in seq_along(dists)) {
      d_in <- dists[[i]]
      buf <- c(buf, 0.0)
      cbuf <- c(cbuf, list(NULL))
      n <- length(buf)
      shift <- 0.0
      coef <- numeric(H)
      coef[i] <- 1.0
      jmax <- min(n, window) - 1
      if (jmax >= 1) {
        for (j in seq_len(jmax)) {
          shift <- shift - w_fwd[j + 1] * buf[n - j]
          cj <- cbuf[[n - j]]
          if (!is.null(cj)) {
            for (ii in seq_len(i - 1)) {
              coef[ii] <- coef[ii] - w_fwd[j + 1] * cj[ii]
            }
          }
        }
      }
      recovered <- dist_mean(d_in) + shift
      buf[n] <- recovered
      cbuf[[n]] <- coef
      total_var <- 0.0
      for (ii in seq_len(i)) {
        total_var <- total_var + coef[ii] * coef[ii] * dist_var(dists[[ii]])
      }
      std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d_in), 1e-12)
      out[[i]] <- dist_gaussian(recovered, std)
      if (length(buf) > window) {
        buf <- buf[-1]
        cbuf <- cbuf[-1]
      }
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
  # Two seeds (all-ones and alternating +-1), max of the two: a single seed can
  # be an eigenvector of a non-dominant root and miss the true radius
  # (phi = c(3, -2, 0, 0, 0) read 1 instead of 2; skaters#232/#235).
  best <- 0.0
  for (seed in list(rep(1.0, p), (-1.0)^(seq_len(p) - 1L))) {
    v <- seed
    rho <- 0.0
    for (.i in seq_len(60)) {
      nv <- c(sum(phi * v), v[-p])
      m <- max(abs(nv))
      if (m == 0) m <- 1.0
      v <- nv / m
      rho <- m
    }
    best <- max(best, rho)
  }
  best
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
      tstate <- list(buffer = numeric(0), phi = numeric(p), phi_used = numeric(p), P = init_P(), n = 0)
    }
    st <- tstate
    if (is.null(st$phi_used)) {
      # a state from before skaters#232/#235 carries only the raw coefficients
      st$phi_used <- .ar_stationary(st$phi)
    }
    st$n <- st$n + 1
    buf <- st$buffer
    if (length(buf) >= p) {
      x <- buf[seq(length(buf), by = -1L, length.out = p)]
      # the residual EMITTED (and inverse_k's forecast) uses the stationarity-
      # projected coefficients from the end of the previous tick (skaters#232);
      # the RLS keeps fitting the raw problem with its own prediction error
      residual <- y - sum(st$phi_used * x)
      raw_residual <- y - sum(st$phi * x)
      P <- st$P
      Px <- .rls_mat_vec(P, x, p)
      denom <- lam + sum(x * Px)
      if (abs(denom) > 1e-15) {
        K <- Px / denom
        st$phi <- st$phi + K * raw_residual
        P <- (P - K %o% Px) / lam
        if (!all(is.finite(P)) || max(abs(P)) > 1e10) {
          P <- init_P()
        }
        st$P <- P
      }
      # re-project after every update so the next emission and inverse_k agree
      st$phi_used <- .ar_stationary(st$phi)
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
    # The same stationarity-projected coefficients forward() emitted residuals
    # under (skaters#232); a pre-#235 state has no cache yet, so derive it.
    phi <- if (!is.null(tstate$phi_used)) tstate$phi_used else .ar_stationary(tstate$phi)
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
      # Var(y_{t+h+1}) = sum_i psi_{h-i}^2 var_i: each future innovation enters
      # with its own variance (skaters#245); plain accumulation in this order.
      total_var <- 0.0
      for (i in seq_len(h)) {
        total_var <- total_var + psi[h - i + 1] * psi[h - i + 1] * dist_var(dists[[i]])
      }
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
      tstate <- list(buffer = numeric(0), theta = numeric(n_groups), phi_used = numeric(max_lag),
                     P = init_P(), n = 0)
    }
    st <- tstate
    if (is.null(st$phi_used)) {
      st$phi_used <- .ar_stationary(st$theta[groups])
    }
    st$n <- st$n + 1
    buf <- st$buffer
    if (length(buf) >= max_lag) {
      x <- group_regressor(buf)
      # emitted residual under the previous tick's projected per-lag coefficients
      # (skaters#232); the RLS fits the raw grouped problem with its own error
      lag_x <- buf[seq(length(buf), by = -1L, length.out = max_lag)]
      residual <- y - sum(st$phi_used * lag_x)
      raw_residual <- y - sum(st$theta * x)
      P <- st$P
      Px <- .rls_mat_vec(P, x, n_groups)
      denom <- lam + sum(x * Px)
      if (abs(denom) > 1e-15) {
        K <- Px / denom
        st$theta <- st$theta + K * raw_residual
        P <- (P - K %o% Px) / lam
        if (!all(is.finite(P)) || max(abs(P)) > 1e10) {
          P <- init_P()
        }
        st$P <- P
      }
      st$phi_used <- .ar_stationary(st$theta[groups])
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
    phi <- if (!is.null(tstate$phi_used)) tstate$phi_used else .ar_stationary(tstate$theta[groups])
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
      # Var(y_{t+h+1}) = sum_i psi_{h-i}^2 var_i: each future innovation enters
      # with its own variance (skaters#245); plain accumulation in this order.
      total_var <- 0.0
      for (i in seq_len(h)) {
        total_var <- total_var + psi[h - i + 1] * psi[h - i + 1] * dist_var(dists[[i]])
      }
      total_std <- if (total_var > 0) sqrt(total_var) else max(dist_std(d), 1e-12)
      recovered_means[h] <- total_mean
      out[[h]] <- dist_gaussian(total_mean, total_std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

# --- Gaussianize: probit of the smoothed fence-post empirical CDF ---
#
# z = Phi^{-1}(F_hat(y)) with fence-post weights, F_hat built from the sorted
# past inputs at most max_knots equal-mass centers. smooth = FALSE
# interpolates the fence-post empirical CDF linearly; smooth = TRUE uses a
# Gaussian-kernel CDF (Silverman bandwidth) through a monotone C^1 cubic
# (Fritsch-Butland tangents). Knot-segment slopes are guarded at slope_cap
# times the interquartile slope (guard "move" caps, "delete" prunes, "emin"
# projects into the smallest rank corridor) and the forward output is clamped
# to |z| <= z_max. The inverse maps each component through an equal-probability
# particle approximation (n_particles local Gaussians) and prunes to 12; a
# tail-spliced predictive gets the exact pushforward. Before warmup
# observations the transform standardizes by expanding moments.
gaussianize <- function(max_knots = 39L, warmup = 30L, n_particles = 9L, smooth = FALSE,
                        slope_cap = 20.0, z_max = 8.0, guard = "move") {
  force(max_knots)
  force(warmup)
  force(smooth)
  force(guard)
  SLOPE_CAP <- slope_cap
  Z_MAX <- z_max
  SQ2 <- sqrt(2.0)
  J <- as.integer(n_particles)
  mids <- qnorm(((0:(J - 1)) + 0.5) / J)
  edges <- if (J > 1) qnorm((1:(J - 1)) / J) else numeric(0)
  .erf <- function(x) 2.0 * pnorm(x * SQ2) - 1.0

  # Cap reference: SLOPE_CAP times the interquartile inverse slope, or NULL.
  slope_ref <- function(ys, zs) {
    K <- length(ys)
    i1 <- K %/% 4
    i3 <- (3 * K) %/% 4
    if (i3 <= i1 || zs[i3 + 1] <= zs[i1 + 1]) {
      return(NULL)
    }
    cap <- SLOPE_CAP * (ys[i3 + 1] - ys[i1 + 1]) / (zs[i3 + 1] - zs[i1 + 1])
    if (cap <= 0) NULL else cap
  }

  # Delete knots whose segment slope exceeds the cap, walking outward from the
  # median knot.
  prune_slopes <- function(ys, zs) {
    cap <- slope_ref(ys, zs)
    if (is.null(cap)) {
      return(list(ys, zs))
    }
    K <- length(ys)
    c0 <- K %/% 2 # 0-based
    keep <- c0
    if (c0 + 1 <= K - 1) {
      for (i in (c0 + 1):(K - 1)) {
        j <- keep[length(keep)]
        if ((ys[i + 1] - ys[j + 1]) / (zs[i + 1] - zs[j + 1]) <= cap) keep <- c(keep, i)
      }
    }
    if (c0 - 1 >= 0) {
      for (i in (c0 - 1):0) {
        j <- min(keep)
        if ((ys[j + 1] - ys[i + 1]) / (zs[j + 1] - zs[i + 1]) <= cap) keep <- c(keep, i)
      }
    }
    keep <- sort(keep)
    if (length(keep) < 2) {
      return(list(ys, zs))
    }
    list(ys[keep + 1], zs[keep + 1])
  }

  # Smallest-rank-displacement guard.
  emin_knots <- function(buf, ys, zs, idxs) {
    cap <- slope_ref(ys, zs)
    if (is.null(cap)) {
      return(list(ys, zs))
    }
    nb <- length(buf)
    K <- length(ys)
    dz <- diff(zs)
    corridor <- function(e) {
      lo <- buf[pmax(0, idxs - e) + 1]
      hi <- buf[pmin(nb - 1, idxs + e) + 1]
      U <- numeric(K)
      U[K] <- hi[K]
      Db <- numeric(K)
      Db[K] <- lo[K]
      if (K >= 2) {
        for (j in (K - 1):1) {
          U[j] <- min(hi[j], U[j + 1])
          Db[j] <- max(lo[j], Db[j + 1] - cap * dz[j])
        }
      }
      Lw <- numeric(K)
      Lw[1] <- Db[1]
      if (K >= 2) {
        for (j in 2:K) Lw[j] <- max(Db[j], Lw[j - 1])
      }
      if (any(Lw > U)) {
        return(NULL)
      }
      list(L = Lw, U = U)
    }
    if (!is.null(corridor(0))) {
      return(list(ys, zs))
    }
    lo_e <- 0
    hi_e <- nb
    while (hi_e - lo_e > 1) {
      mid <- (lo_e + hi_e) %/% 2
      if (!is.null(corridor(mid))) hi_e <- mid else lo_e <- mid
    }
    res <- corridor(hi_e)
    if (is.null(res)) {
      return(list(ys, zs))
    }
    Lw <- res$L
    U <- res$U
    w <- numeric(K)
    w[1] <- max(Lw[1], min(ys[1], U[1]))
    if (K >= 2) {
      for (j in 2:K) {
        a_j <- max(Lw[j], w[j - 1])
        b_j <- min(U[j], w[j - 1] + cap * dz[j - 1])
        w[j] <- max(a_j, min(ys[j], b_j))
      }
    }
    # drop any knot that failed to stay strictly increasing
    ks <- w[1]
    kz <- zs[1]
    if (K >= 2) {
      for (j in 2:K) {
        if (w[j] > ks[length(ks)]) {
          ks <- c(ks, w[j])
          kz <- c(kz, zs[j])
        }
      }
    }
    list(ks, kz)
  }

  # Cap segment slopes, walking outward from the median knot.
  cap_slopes <- function(ys, zs) {
    cap <- slope_ref(ys, zs)
    if (is.null(cap)) {
      return(ys)
    }
    K <- length(ys)
    c0 <- K %/% 2 + 1 # 1-based median knot
    if (c0 + 1 <= K) {
      for (i in (c0 + 1):K) ys[i] <- min(ys[i], ys[i - 1] + cap * (zs[i] - zs[i - 1]))
    }
    if (c0 - 1 >= 1) {
      for (i in (c0 - 1):1) ys[i] <- max(ys[i], ys[i + 1] - cap * (zs[i + 1] - zs[i]))
    }
    ys
  }

  # Fritsch-Butland harmonic-mean tangents: monotone C^1 cubic.
  tangents <- function(xs, vs) {
    d <- diff(vs) / diff(xs)
    n <- length(xs)
    t <- numeric(n)
    t[1] <- d[1]
    if (n > 2) {
      for (i in 2:(n - 1)) {
        a <- d[i - 1]
        b <- d[i]
        t[i] <- if (a * b > 0) 2.0 * a * b / (a + b) else 0.0
      }
    }
    t[n] <- d[n - 1]
    t
  }
  herm <- function(s, dx, v0, v1, t0, t1) {
    s2 <- s * s
    s3 <- s2 * s
    (v0 * (2 * s3 - 3 * s2 + 1) + dx * t0 * (s3 - 2 * s2 + s) +
      v1 * (-2 * s3 + 3 * s2) + dx * t1 * (s3 - s2))
  }
  herm_ds <- function(s, dx, v0, v1, t0, t1) {
    s2 <- s * s
    (v0 * (6 * s2 - 6 * s) + dx * t0 * (3 * s2 - 4 * s + 1) +
      v1 * (-6 * s2 + 6 * s) + dx * t1 * (3 * s2 - 2 * s))
  }

  build <- function(buf) {
    n <- length(buf)
    m <- min(max_knots, n)
    cs <- numeric(0)
    idxs <- integer(0)
    for (i in 0:(m - 1)) {
      idx <- as.integer((i + 0.5) * n / m) # 0-based
      y <- buf[idx + 1]
      if (length(cs) && y <= cs[length(cs)]) next
      cs <- c(cs, y)
      idxs <- c(idxs, idx)
    }
    if (length(cs) < 2) {
      return(NULL)
    }
    if (!smooth) {
      zs <- qnorm((idxs + 1) / (n + 1))
      if (guard == "delete") {
        r <- prune_slopes(cs, zs)
        cs <- r[[1]]
        zs <- r[[2]]
      } else if (guard == "emin") {
        r <- emin_knots(buf, cs, zs, idxs)
        cs <- r[[1]]
        zs <- r[[2]]
      } else {
        cs <- cap_slopes(cs, zs)
      }
      return(list(cs = cs, zs = zs, ts = NULL))
    }
    # kernel-smoothed CDF on the sketch, Silverman bandwidth
    mm <- length(cs)
    i1 <- mm %/% 4
    i3 <- (3 * mm) %/% 4
    iqr <- cs[i3 + 1] - cs[i1 + 1]
    mean <- sum(cs) / mm
    sd <- sqrt(sum((cs - mean)^2) / mm)
    scale <- if (iqr > 0) min(sd, iqr / 1.34) else sd
    if (scale <= 0) {
      scale <- max(abs(cs[mm] - cs[1]), 1e-8)
    }
    h <- max(0.9 * scale * n^-0.2, 1e-12)
    R <- 8.0 * h
    zs <- numeric(mm)
    for (ci in seq_len(mm)) {
      cc <- cs[ci]
      lo <- sum(cs < cc - R)
      hi <- sum(cs <= cc + R)
      inner <- if (hi > lo) sum(0.5 * (1.0 + .erf((cc - cs[(lo + 1):hi]) / (SQ2 * h)))) else 0.0
      f <- (lo + inner) / mm
      zs[ci] <- qnorm((1.0 + n * f) / (n + 2.0))
    }
    cs2 <- cs[1]
    zs2 <- zs[1]
    if (mm >= 2) {
      for (ci in 2:mm) {
        if (zs[ci] > zs2[length(zs2)] + 1e-12) {
          cs2 <- c(cs2, cs[ci])
          zs2 <- c(zs2, zs[ci])
        }
      }
    }
    if (length(cs2) < 2) {
      return(NULL)
    }
    cs2 <- cap_slopes(cs2, zs2)
    list(cs = cs2, zs = zs2, ts = tangents(cs2, zs2))
  }

  g_fwd <- function(y, knots) {
    cs <- knots$cs
    zs <- knots$zs
    ts <- knots$ts
    K <- length(cs)
    if (y <= cs[1]) {
      s0 <- if (!is.null(ts)) ts[1] else (zs[2] - zs[1]) / (cs[2] - cs[1])
      return(zs[1] + s0 * (y - cs[1]))
    }
    if (y >= cs[K]) {
      s1 <- if (!is.null(ts)) ts[K] else (zs[K] - zs[K - 1]) / (cs[K] - cs[K - 1])
      return(zs[K] + s1 * (y - cs[K]))
    }
    i <- sum(cs <= y) # bisect_right - 1, as a 1-based index
    if (is.null(ts)) {
      return(zs[i] + (zs[i + 1] - zs[i]) * (y - cs[i]) / (cs[i + 1] - cs[i]))
    }
    dx <- cs[i + 1] - cs[i]
    herm((y - cs[i]) / dx, dx, zs[i], zs[i + 1], ts[i], ts[i + 1])
  }

  g_inv <- function(z, knots) {
    cs <- knots$cs
    zs <- knots$zs
    ts <- knots$ts
    K <- length(cs)
    if (z <= zs[1]) {
      s0 <- if (!is.null(ts)) ts[1] else (zs[2] - zs[1]) / (cs[2] - cs[1])
      return(cs[1] + (z - zs[1]) / max(s0, 1e-12))
    }
    if (z >= zs[K]) {
      s1 <- if (!is.null(ts)) ts[K] else (zs[K] - zs[K - 1]) / (cs[K] - cs[K - 1])
      return(cs[K] + (z - zs[K]) / max(s1, 1e-12))
    }
    i <- sum(zs <= z)
    if (is.null(ts)) {
      return(cs[i] + (cs[i + 1] - cs[i]) * (z - zs[i]) / (zs[i + 1] - zs[i]))
    }
    dx <- cs[i + 1] - cs[i]
    s <- (z - zs[i]) / (zs[i + 1] - zs[i])
    for (it in 1:6) {
      dd <- herm_ds(s, dx, zs[i], zs[i + 1], ts[i], ts[i + 1])
      if (dd <= 1e-14) break
      s <- s - (herm(s, dx, zs[i], zs[i + 1], ts[i], ts[i + 1]) - z) / dd
      s <- if (s < 0.0) 0.0 else if (s > 1.0) 1.0 else s
    }
    cs[i] + s * dx
  }

  g_moments <- function(st) {
    var <- if (st$n > 2) st$m2 / (st$n - 1) else 0.0
    c(st$mu, if (var > 0) sqrt(var) else 1e-8)
  }

  forward <- function(y, tstate = NULL) {
    st <- if (is.null(tstate)) list(buf = numeric(0), n = 0L, mu = 0.0, m2 = 0.0, knots = NULL) else tstate
    knots <- st$knots
    if (!is.null(knots)) {
      z <- g_fwd(y, knots)
    } else {
      ms <- g_moments(st)
      mu <- ms[1]
      sigma <- ms[2]
      if (sigma <= 1e-8) {
        # cold start: scale by the first observed deviation
        sigma <- max(abs(y - mu), 1e-8)
      }
      z <- if (st$n > 0) (y - mu) / sigma else 0.0
    }
    z <- max(-Z_MAX, min(Z_MAX, z))
    st$buf <- append(st$buf, y, after = sum(st$buf <= y)) # bisect.insort
    n <- st$n + 1L
    d0 <- y - st$mu
    mu <- st$mu + d0 / n
    st$m2 <- st$m2 + d0 * (y - mu)
    st$n <- n
    st$mu <- mu
    st["knots"] <- list(if (n >= warmup) build(st$buf) else NULL)
    list(y = z, state = st)
  }

  inverse_k <- function(dists, tstate) {
    knots <- tstate$knots
    if (is.null(knots)) {
      ms <- g_moments(tstate)
      return(lapply(dists, dist_affine, a = ms[2], b = ms[1]))
    }
    out <- vector("list", length(dists))
    for (di in seq_along(dists)) {
      d <- dist_components(dists[[di]])
      if (is.null(d)) { # tail-spliced: exact pushforward
        out[[di]] <- .pushforward(dists[[di]], function(z) g_inv(z, knots))
        next
      }
      w_out <- numeric(0)
      m_out <- numeric(0)
      s_out <- numeric(0)
      for (ci in seq_along(d$w)) {
        w <- d$w[ci]
        mu <- d$m[ci]
        sig <- d$s[ci]
        pts <- vapply(mids, function(q) g_inv(mu + sig * q, knots), 0.0)
        eds <- vapply(edges, function(e) g_inv(mu + sig * e, knots), 0.0)
        for (j in seq_len(J)) {
          lo <- if (j > 1) eds[j - 1] else NULL
          hi <- if (j < J) eds[j] else NULL
          std <- if (!is.null(lo) && !is.null(hi)) {
            (hi - lo) / 2.0
          } else if (!is.null(hi)) {
            hi - pts[j]
          } else {
            pts[j] - lo
          }
          w_out <- c(w_out, w / J)
          m_out <- c(m_out, pts[j])
          s_out <- c(s_out, max(std, 1e-9))
        }
      }
      out[[di]] <- dist_prune(dist_new(w_out, m_out, s_out), 12L)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}
