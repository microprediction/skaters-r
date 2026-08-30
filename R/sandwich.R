# The sandwich: run any model in laplace's z-coordinates and map its forecasts
# back with no approximation in the accounting.
#
# The front half is the Rosenblatt transform z_t = Phi^-1(F_t(y_t)), where F_t
# is the predictive laplace issued for y_t. parade already emits exactly that
# stream in state$z; rosenblatt() below is the batch form, and it also keeps
# the transport needed to map forecasts home.
#
# This file is the other half: given the inner model's z-space Gaussian and the
# h-step transport, assemble the y-space predictive. Nothing here knows what the
# inner model is, so the package takes on no dependency to support one.

.SANDWICH_EPS <- 1e-12

# Rosenblatt front-end over a whole series.
#
# z[1] is NA because no predictive had been issued when the first observation
# arrived. The Python mirror (prophet-laplace) substitutes 0 there instead,
# because Prophet requires a gap-free series; R models take NA, so this reports
# the absence rather than inventing a median draw.
rosenblatt <- function(y, k = 1L, fn = NULL) {
  y <- as.numeric(y)
  stopifnot(length(y) > 0)
  if (is.null(fn)) {
    fn <- laplace(k = k)
  }
  z <- rep(NA_real_, length(y))
  pit <- rep(NA_real_, length(y))
  st <- NULL
  pend <- NULL
  for (i in seq_along(y)) {
    if (!is.null(pend)) {
      u <- min(max(dist_cdf(pend[[1]], y[i]), .SANDWICH_EPS), 1.0 - .SANDWICH_EPS)
      pit[i] <- u
      z[i] <- .phi_inv(u)
    }
    r <- fn(y[i], st)
    pend <- r$dists
    st <- r$state
  }
  list(z = z, pit = pit, transport = pend, state = st, k = length(pend))
}

# The transport for a given step ahead. Steps past k reuse the k-step
# transport, which is an approximation and is disclosed as one.
sandwich_transport <- function(r, step = 1L) {
  stopifnot(step >= 1L)
  r$transport[[min(as.integer(step), r$k)]]
}

# A y-space predictive built from the inner model's Gaussian in z and the
# transport laplace issued for that step.
#
# With z_mu = 0 and z_sigma = 1 this is the identity on body: the sandwich
# reduces to laplace alone when the inner model says nothing. test-sandwich.R
# asserts that, since it is the invariant that catches a wrong Jacobian.
sandwich_new <- function(body, z_mu = 0.0, z_sigma = 1.0) {
  stopifnot(is.list(body), is.finite(z_mu), is.finite(z_sigma))
  list(
    sandwiched = TRUE,
    body = body,
    z_mu = z_mu,
    z_sigma = max(z_sigma, 1e-9)
  )
}

.sandwich_z <- function(d, x) {
  u <- min(max(dist_cdf(d$body, x), .SANDWICH_EPS), 1.0 - .SANDWICH_EPS)
  .phi_inv(u)
}

sandwich_cdf <- function(d, x) {
  pnorm((.sandwich_z(d, x) - d$z_mu) / d$z_sigma)
}

# log f_Y(y) = log f_Z(z) + log f_t(y) - log phi(z)
sandwich_logpdf <- function(d, x) {
  base <- dist_logpdf(d$body, x)
  if (!is.finite(base)) {
    return(base)
  }
  z <- .sandwich_z(d, x)
  r <- (z - d$z_mu) / d$z_sigma
  log_fz <- -0.5 * r * r - log(d$z_sigma) - .LOG_SQRT2PI
  log_phi <- -0.5 * z * z - .LOG_SQRT2PI
  base + log_fz - log_phi
}

sandwich_quantile <- function(d, p, tol = 1e-9, max_iter = 100L) {
  stopifnot(p > 0, p < 1)
  z <- d$z_mu + d$z_sigma * .phi_inv(p)
  u <- min(max(pnorm(z), .SANDWICH_EPS), 1.0 - .SANDWICH_EPS)
  dist_quantile(d$body, u, tol = tol, max_iter = max_iter)
}

# Numeric moments and CRPS over a fixed 65-node quantile grid (midpoint rule),
# the same treatment spliced bodies get.
.SANDWICH_GRID_N <- 65L

sandwich_qgrid <- function(d) {
  n <- .SANDWICH_GRID_N
  vapply(seq_len(n), function(i) sandwich_quantile(d, (i - 0.5) / n), 0.0)
}

sandwich_mean <- function(d) {
  q <- sandwich_qgrid(d)
  sum(q) / length(q)
}

sandwich_var <- function(d) {
  q <- sandwich_qgrid(d)
  m <- sum(q) / length(q)
  sum((q - m) * (q - m)) / length(q)
}

sandwich_crps <- function(d, x) {
  q <- sandwich_qgrid(d)
  n <- length(q)
  t1 <- sum(abs(q - x)) / n
  t2 <- 2.0 * sum(q * (2.0 * (seq_len(n) - 0.5) / n - 1.0)) / n
  t1 - t2 * 0.5
}
