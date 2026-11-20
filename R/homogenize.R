# Online homogenized scale-mixture correction: a small Bayesian-hedged pool
# of two-point Gaussian scale-mixture hypotheses for a skater's own PIT
# residual stream. Port of skaters/homogenize.py.
#
# Each candidate (rho, gain, delta) carries a Huber-capped fixed-gain filter q
# on H_2(z) = z^2 - 1, giving the scale mixture
#   g(z) = 1/2 N(0, v(1-delta)) + 1/2 N(0, v(1+delta)),  v = clip(1 + q)
# (one N(0, v) when delta = 0), and a discounted Bayesian log-weight. An
# identity candidate frozen at v = 1 holds most of the prior weight. Pooling
# is exact: the pooled density is the union of the candidates' (scale,
# weight) pairs. The pooled correction H composes onto the base predictive:
#   F~(y) = H(Phi^-1(F(y))),  Q~(p) = F.quantile(Phi(H.quantile(p)))

.HG_EPS <- 1e-12
.HG_LOG_SQRT2PI <- log(sqrt(2.0 * pi))
.HG_V_MIN <- 0.25
.HG_V_MAX <- 9.0
.HG_DELTA_MAX <- 0.8
.HG_INNOVATION_CAP <- 8.0
.HG_ADAPT_AFTER <- 10L
.HG_WEIGHT_DECAY <- 0.9999
.HG_LEARNING_RATE <- 0.1
.HG_SCORE_CAP <- 20.0

.hg_clip <- function(x, lo, hi) if (x < lo) lo else if (x > hi) hi else x
.hg_clamp01 <- function(u) min(max(u, .HG_EPS), 1.0 - .HG_EPS)
.hg_phi_inv <- function(p) dist_quantile(dist_gaussian(0.0, 1.0), p)
.hg_phi_logpdf <- function(z) -0.5 * z * z - .HG_LOG_SQRT2PI

# Log-density of a scale mixture list(weights, scales) at z (log-sum-exp).
.hg_logg <- function(z, correction) {
  terms <- numeric(0)
  best <- -Inf
  for (i in seq_along(correction$weights)) {
    w <- correction$weights[i]
    if (w <= 0.0) next
    s <- correction$scales[i]
    t <- log(w) - log(s) + .hg_phi_logpdf(z / s)
    terms <- c(terms, t)
    if (t > best) best <- t
  }
  if (best == -Inf) {
    return(-Inf)
  }
  best + log(sum(exp(terms - best)))
}

# Candidate config grid plus a distinguished identity candidate. The scale
# mixture is a density only for 0 <= delta < 1. Python's make_candidates,
# exported there as homogenize_candidates.
homogenize_candidates <- function(rho_grid = c(0.0, 0.8, 0.95, 0.99), gain_grid = c(0.03, 0.12),
                                  delta_grid = c(0.0, 0.25, 0.50), delta_max = .HG_DELTA_MAX) {
  if (!(0.0 <= delta_max && delta_max < 1.0)) {
    stop(sprintf("delta_max must satisfy 0 <= delta_max < 1, got %s", format(delta_max)))
  }
  for (delta in delta_grid) {
    if (!(0.0 <= delta && delta < 1.0)) {
      stop(sprintf("delta_grid entries must satisfy 0 <= delta < 1, got %s", format(delta)))
    }
  }
  cands <- list(list(rho = 0.0, gain = 0.0, delta = 0.0, is_identity = TRUE))
  for (rho in rho_grid) {
    for (gain in gain_grid) {
      for (delta in delta_grid) {
        cands[[length(cands) + 1L]] <- list(
          rho = as.numeric(rho), gain = as.numeric(gain),
          delta = as.numeric(.hg_clip(delta, 0.0, delta_max)), is_identity = FALSE
        )
      }
    }
  }
  cands
}

.hg_cand_scales <- function(v, delta) {
  if (delta <= 0.0) sqrt(v) else c(sqrt(v * (1.0 - delta)), sqrt(v * (1.0 + delta)))
}

.hg_cand_weights <- function(delta) if (delta <= 0.0) 1.0 else c(0.5, 0.5)

.hg_cand_logpdf <- function(z, v, delta) {
  .hg_logg(z, list(weights = .hg_cand_weights(delta), scales = .hg_cand_scales(v, delta)))
}

.hg_init_state <- function(candidates) {
  n_other <- max(sum(!vapply(candidates, function(c) isTRUE(c$is_identity), FALSE)), 1)
  dynamic <- lapply(candidates, function(c) {
    list(q = 0.0, log_weight = if (isTRUE(c$is_identity)) log(n_other) else 0.0, guard_run = 0L)
  })
  list(dynamic = dynamic, issued = NULL, n = 0L)
}

# Fixed-gain filter on x = H_2(z), Huber-capped, with a changepoint escape.
.hg_update_q <- function(d, cfg, x) {
  q_minus <- cfg$rho * d$q
  resid <- x - q_minus
  capped <- abs(resid) > .HG_INNOVATION_CAP
  d$guard_run <- if (capped) d$guard_run + 1L else 0L
  used <- if (capped && d$guard_run <= .HG_ADAPT_AFTER) {
    .hg_clip(resid, -.HG_INNOVATION_CAP, .HG_INNOVATION_CAP)
  } else {
    resid
  }
  d$q <- q_minus + cfg$gain * used
  d
}

.hg_omega <- function(dynamic) {
  lw <- vapply(dynamic, function(d) d$log_weight, 0.0)
  raw <- exp(lw - max(lw))
  raw / sum(raw)
}

# Exact Bayesian pooling: the union of every candidate's (scale, weight).
.hg_pool <- function(candidates, dynamic) {
  omega <- .hg_omega(dynamic)
  weights <- numeric(0)
  scales <- numeric(0)
  for (i in seq_along(candidates)) {
    cfg <- candidates[[i]]
    v <- .hg_clip(1.0 + dynamic[[i]]$q, .HG_V_MIN, .HG_V_MAX)
    scales <- c(scales, .hg_cand_scales(v, cfg$delta))
    weights <- c(weights, omega[i] * .hg_cand_weights(cfg$delta))
  }
  list(weights = weights, scales = scales)
}

# Score each candidate's pre-update correction against z, update weights,
# update each filter, then pool the correction issued for the next tick.
# Returns list(correction, state).
.hg_cell_step <- function(z, state, candidates) {
  if (is.null(state)) {
    state <- .hg_init_state(candidates)
  }
  dyn <- state$dynamic
  x <- z * z - 1.0
  for (i in seq_along(candidates)) {
    v <- .hg_clip(1.0 + dyn[[i]]$q, .HG_V_MIN, .HG_V_MAX)
    r <- .hg_cand_logpdf(z, v, candidates[[i]]$delta) - .hg_phi_logpdf(z)
    dyn[[i]]$log_weight <- .HG_WEIGHT_DECAY * dyn[[i]]$log_weight +
      .HG_LEARNING_RATE * .hg_clip(r, -.HG_SCORE_CAP, .HG_SCORE_CAP)
  }
  for (i in seq_along(candidates)) {
    dyn[[i]] <- .hg_update_q(dyn[[i]], candidates[[i]], x)
  }
  correction <- .hg_pool(candidates, dyn)
  state$dynamic <- dyn
  state$issued <- correction
  state$n <- state$n + 1L
  list(correction = correction, state = state)
}

# --- the homogenized predictive ---

# The base predictive reshaped by a pooled scale-mixture correction
# list(weights, scales). Answers the dist_* probes. Port of HomogenizedDist.
homogenized_dist <- function(base, correction) {
  list(homogenized = TRUE, base = base, correction = correction, cache = new.env(parent = emptyenv()))
}

.hd_z <- function(d, y) .hg_phi_inv(.hg_clamp01(dist_cdf(d$base, y)))

.hd_z_cdf <- function(d, z) {
  total <- 0.0
  for (i in seq_along(d$correction$weights)) {
    total <- total + d$correction$weights[i] * pnorm(z / d$correction$scales[i])
  }
  .hg_clamp01(total)
}

.hd_z_quantile <- function(d, p, tol = 1e-9, max_iter = 100L) {
  lo <- -40.0
  hi <- 40.0
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (.hd_z_cdf(d, mid) < p) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  0.5 * (lo + hi)
}

homogenized_cdf <- function(d, y) .hd_z_cdf(d, .hd_z(d, y))

homogenized_logpdf <- function(d, y) {
  z <- .hd_z(d, y)
  dist_logpdf(d$base, y) + .hg_logg(z, d$correction) - .hg_phi_logpdf(z)
}

homogenized_quantile <- function(d, p, tol = 1e-9, max_iter = 2100L) {
  stopifnot(p > 0, p < 1)
  hp <- .hd_z_quantile(d, p, tol = tol, max_iter = max_iter)
  u <- .hg_clamp01(pnorm(hp))
  dist_quantile(d$base, u, tol = tol, max_iter = max_iter)
}

homogenized_mean <- function(d) .qgrid_mean(.cd_qgrid(d, homogenized_quantile))
homogenized_var <- function(d) .qgrid_var(.cd_qgrid(d, homogenized_quantile))
homogenized_std <- function(d) .qgrid_std(.cd_qgrid(d, homogenized_quantile))
homogenized_crps <- function(d, x) .qgrid_crps(.cd_qgrid(d, homogenized_quantile), x)

# Borrow strength across horizons, in place on the cell list (returned).
.hg_couple <- function(cells, mode, lam) {
  live_idx <- which(!vapply(cells, is.null, FALSE))
  if (length(live_idx) < 2 || lam <= 0.0) {
    return(cells)
  }
  n <- length(cells[[live_idx[1]]]$dynamic)
  if (mode == "weights") {
    for (i in seq_len(n)) {
      mean_w <- sum(vapply(live_idx, function(j) cells[[j]]$dynamic[[i]]$log_weight, 0.0)) / length(live_idx)
      for (j in live_idx) {
        lw <- cells[[j]]$dynamic[[i]]$log_weight
        cells[[j]]$dynamic[[i]]$log_weight <- (1.0 - lam) * lw + lam * mean_w
      }
    }
  } else if (mode == "state") {
    for (i in seq_len(n)) {
      prev <- NULL
      for (j in live_idx) {
        if (!is.null(prev)) {
          q <- cells[[j]]$dynamic[[i]]$q
          cells[[j]]$dynamic[[i]]$q <- (1.0 - lam) * q + lam * prev
        }
        prev <- cells[[j]]$dynamic[[i]]$q
      }
    }
  }
  cells
}

# Wrap a skater with an online homogenized scale-mixture correction, one
# candidate pool per horizon (couple = "none"), or coupled across horizons by
# "weights", "state" or one "shared" pool. The state carries pending (the
# last k predictive lists), cell (per-horizon pool state), correction and raw.
homogenize <- function(base, candidates = NULL, k = 1L, couple = "none", lam = 0.25) {
  force(base)
  force(k)
  force(lam)
  stopifnot(k >= 1, couple %in% c("none", "weights", "state", "shared"))
  candidates <- if (is.null(candidates)) homogenize_candidates() else candidates
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(base = NULL, raw = NULL, pending = list(),
                    cell = vector("list", k), correction = vector("list", k))
    }
    s <- state
    pend <- s$pending
    n <- length(pend)
    for (m in seq_len(k)) {
      if (m > n) next
      d_m <- pend[[n - m + 1]][[m]]
      z <- .hg_phi_inv(.hg_clamp01(dist_cdf(d_m, y)))
      if (is.nan(z)) next
      slot <- if (couple == "shared") 1L else m
      st <- .hg_cell_step(z, s$cell[[slot]], candidates)
      s$cell[slot] <- list(st$state)
      if (couple == "shared") {
        for (j in seq_len(k)) s$correction[j] <- list(st$correction)
      } else {
        s$correction[m] <- list(st$correction)
      }
    }
    if (couple %in% c("weights", "state")) {
      s$cell <- .hg_couple(s$cell, couple, lam)
      for (m in seq_len(k)) {
        if (!is.null(s$cell[[m]])) {
          s$correction[m] <- list(.hg_pool(candidates, s$cell[[m]]$dynamic))
        }
      }
    }
    r <- base(y, s$base)
    s["base"] <- list(r$state)
    dists <- r$dists
    s$raw <- dists[[1]]
    pend[[length(pend) + 1L]] <- dists
    if (length(pend) > k) pend <- pend[-1]
    s$pending <- pend
    out <- vector("list", k)
    for (m in seq_len(k)) {
      d <- if (m <= length(dists)) dists[[m]] else dists[[length(dists)]]
      corr <- s$correction[[m]]
      out[[m]] <- if (is.null(corr)) d else homogenized_dist(d, corr)
    }
    list(dists = out, state = s)
  }
}
