# granger: laplace with trained parameters, mixed with laplace so it cannot do
# much worse. Port of the Python reference skaters/granger.py (horizon 1),
# parity-checked against it (inst/parity/check_granger.R).
#
# After Clive Granger, whose work with Bates (1969) founded forecast
# combination. The issued predictive is the Aggregating-Algorithm mixture of
# laplace and its trained twin under log loss, so its total log score is never
# more than log(1 / pi_lap) behind laplace's (log 4 at the shipped 0.25), up to
# floating-point rounding.

#' Shipped granger parameters ("monthly", "weekly" or "hourly").
#' @export
granger_params <- function(frequency = c("monthly", "weekly", "hourly")) {
  frequency <- match.arg(frequency)
  path <- system.file("extdata", "granger", paste0(frequency, ".rds"), package = "skaters")
  if (!nzchar(path)) {
    path <- file.path("inst", "extdata", "granger", paste0(frequency, ".rds"))
  }
  if (!file.exists(path)) stop("no shipped granger parameters for ", frequency)
  readRDS(path)
}

.granger_transform <- function(name, kw) {
  no_exact <- function() if (isTRUE(kw$exact)) stop("granger: ", name, "(exact=TRUE) is not ported")
  switch(name,
    ema_transform = ema_transform(kw$alpha),
    difference = difference(),
    drift = drift(kw$alpha, kw$shrinkage),
    theta = theta(kw$alpha),
    ar = ar(as.integer(kw$order), kw$lam, kw$ridge, kw$decay),
    holt_linear = holt_linear(kw$alpha, kw$beta),
    seasonal_difference = seasonal_difference(as.integer(kw$period)),
    seasonal_anchor = seasonal_anchor(as.integer(kw$period), kw$alpha, kw$weight),
    standardize = standardize(kw$alpha, kw$eps),
    fractional_difference = fractional_difference(kw$d, as.integer(kw$window)),
    garch = garch(kw$omega, kw$alpha, kw$beta, kw$mean_alpha, kw$eps),
    power_transform = { no_exact(); power_transform(kw$p) },
    yeo_johnson = { no_exact(); yeo_johnson(kw$lmbda) },
    stop("granger: unknown transform ", name)
  )
}

.granger_candidate <- function(chain, k) {
  f <- leaf(k)
  for (j in rev(seq_along(chain))) {
    f <- conjugate(f, .granger_transform(chain[[j]]$name, chain[[j]]$kw), k)
  }
  f
}

.granger_expert <- function(P, k, gate_hook) {
  cands <- lapply(P$candidates, function(c) .granger_candidate(c$chain, k))
  depths <- vapply(P$candidates, function(c) c$depth, 0.0)
  eta <- P$tleaf$eta
  scale_alpha <- P$tleaf$scale_alpha
  f <- terminal_leaf_ensemble(
    cands,
    k = k,
    leaf_fn = function(k = 1L) crps_leaf(k = k, eta = eta, scale_alpha = scale_alpha),
    learning_rate = P$trunk$learning_rate,
    complexity_penalty = P$trunk$complexity_penalty,
    depths = depths,
    prior_log_weights = P$trunk$prior_log_weights,
    max_components = 20L,
    forget = P$trunk$forget,
    gate_hook = gate_hook
  )
  f <- sticky(f, k = k)
  f <- gpdtails(f, k)
  parade(f, k)
}

.granger_trunk_state <- function(s) {
  for (i in 1:3) {
    if (is.null(s)) return(NULL)
    s <- s$base
  }
  s
}

.granger_gate <- function(G, feats) {
  hidden <- tanh(G$l1_bias + as.vector(G$l1_weight %*% feats))
  G$l2_bias + as.vector(G$l2_weight %*% hidden)
}

.granger_features <- function(state, n) {
  z <- state$z
  z <- if (is.null(z) || is.na(z)) 0.0 else max(-8.0, min(8.0, z))
  ts <- .granger_trunk_state(state$learned)
  if (!is.null(ts) && length(ts$log_w) == n) {
    lw <- pmax(ts$log_w - max(ts$log_w), -20.0) / 20.0
  } else {
    lw <- rep(0.0, n)
  }
  c(z, state$repeat_frac, state$abs_z_mean, lw)
}

.granger_mixture <- function(dists, log_w) {
  m <- max(log_w)
  lw <- log_w - (log(sum(exp(log_w - m))) + m)
  list(mixture = TRUE, dists = dists, log_w = lw, w = exp(lw))
}

mixture_logpdf <- function(d, x) {
  v <- vapply(d$dists, function(q) dist_logpdf(q, x), 0.0)
  t <- ifelse(is.na(v), -Inf, d$log_w + v)   # a NaN component contributes nothing
  m <- max(t)
  if (m == Inf) return(Inf)
  if (m == -Inf) return(-Inf)
  m + log(sum(exp(t - m)))
}

mixture_cdf <- function(d, x) sum(d$w * vapply(d$dists, function(q) dist_cdf(q, x), 0.0))

mixture_mean <- function(d) sum(d$w * vapply(d$dists, dist_mean, 0.0))

mixture_var <- function(d) {
  ms <- vapply(d$dists, dist_mean, 0.0)
  vs <- vapply(d$dists, dist_var, 0.0)
  mu <- sum(d$w * ms)
  sum(d$w * (vs + (ms - mu)^2))
}

mixture_quantile <- function(d, p, tol = 1e-9, max_iter = 200L) {
  stopifnot(p > 0, p < 1)
  qs <- vapply(d$dists, function(q) dist_quantile(q, p), 0.0)
  lo <- min(qs)
  hi <- max(qs)
  # absolute tol at scale >= 1, relative below it so small-scale series keep their digits
  tol <- tol * min(1.0, max(abs(lo), abs(hi)))
  if (hi - lo <= tol) return(0.5 * (lo + hi))
  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    if (mixture_cdf(d, mid) < p) lo <- mid else hi <- mid
    if (hi - lo <= tol) break
  }
  0.5 * (lo + hi)
}

mixture_crps <- function(d, x, n = 65L) {
  q <- vapply(seq_len(n), function(i) mixture_quantile(d, (i - 0.5) / n), 0.0)
  t1 <- sum(abs(q - x)) / n
  t2 <- 2.0 * sum(q * (2.0 * ((seq_len(n) - 1) + 0.5) / n - 1.0)) / n
  t1 - 0.5 * t2
}

#' granger: laplace's trained twin mixed with laplace (horizon 1).
#'
#' @param params a parameter list from [granger_params()].
#' @return a skater `f(y, state)` returning `list(dists, state)`.
#' @export
granger <- function(params = granger_params("monthly"), k = 1L) {
  if (k != 1L) stop("granger is trained at horizon 1 only")
  P <- params
  if (!(P$pi_lap > 0 && P$pi_lap < 1)) stop("pi_lap must lie in (0, 1)")
  log_prior <- c(log(P$pi_lap), log(1.0 - P$pi_lap))
  ref <- laplace(1L)
  cur <- new.env(parent = emptyenv())
  cur$offset <- NULL
  ex <- .granger_expert(P, 1L, function() cur$offset)
  n <- length(P$candidates)

  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(lap = NULL, learned = NULL, pend = NULL, cum = c(0.0, 0.0), t = 0L,
                    z = NA_real_, repeat_frac = 0.0, abs_z_mean = 0.0, n = 0L, last_y = NA_real_)
    }
    yv <- if (is.null(y)) NaN else as.numeric(y)
    finite <- is.finite(yv)
    if (finite && !is.null(state$pend)) {
      for (j in 1:2) {
        lp <- dist_logpdf(state$pend[[j]], yv)
        # +Inf (an exact hit on a lattice atom) counts as +1e6; NaN or -Inf as -1e6
        state$cum[j] <- state$cum[j] + (if (is.na(lp)) -1e6 else if (lp == Inf) 1e6 else if (is.finite(lp)) lp else -1e6)
      }
      # re-centre: the weights depend only on the difference
      state$cum <- state$cum - max(state$cum)
    }
    cur$offset <- .granger_gate(P$gate, .granger_features(state, n))
    r <- tryCatch(ex(yv, state$learned), finally = cur$offset <- NULL)
    state$learned <- r$state
    lrn <- r$dists[[1]]
    rl <- ref(yv, state$lap)
    state$lap <- rl$state
    lap <- rl$dists[[1]]
    state$pend <- list(lap, lrn)
    pz <- state$learned$z
    state$z <- if (is.null(pz)) NA_real_ else pz[1]
    if (finite) {
      nn <- state$n + 1L
      rep_ <- if (!is.na(state$last_y) && yv == state$last_y) 1.0 else 0.0
      state$repeat_frac <- state$repeat_frac + (rep_ - state$repeat_frac) / nn
      z0 <- state$z
      if (!is.na(z0) && is.finite(z0)) {
        state$abs_z_mean <- state$abs_z_mean + (min(abs(z0), 8.0) - state$abs_z_mean) / nn
      }
      state$n <- nn
      state$last_y <- yv
    }
    out <- .granger_mixture(list(lap, lrn), log_prior + state$cum)
    state$t <- state$t + 1L
    list(dists = list(out), state = state)
  }
}
