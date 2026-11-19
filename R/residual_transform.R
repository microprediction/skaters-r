# Online residual-transform wrapper: learn to reshape a skater's own errors.
# Port of skaters/residual_transform.py.
#
# Given the base forecast F_t and a residual model H_t (a predictive for the
# next z = Phi^-1(F(y)), learned from past z), the corrected forecast is
#   F~_t(y) = H_t(Phi^-1(F_t(y)))
#   Q~_t(p) = F_t.quantile(Phi(H_t.quantile(p)))
# Nested models: m0 identity, m1 dynamic scale, m2 scale plus leverage, with
# (a, b, c) learned online by an RMSProp-normalized RTRL-lite step.

.RT_EPS <- 1e-12
.RT_LOG_SQRT2PI <- log(sqrt(2.0 * pi))
.RT_A_MIN <- 0.0
.RT_A_MAX <- 0.995
.RT_BC_MAX <- 3.0
.RT_ELL_MIN <- -20.0
.RT_ELL_MAX <- 20.0
.RT_SENS_MAX <- 50.0
.RT_GRAD_MAX <- 3.0
.RT_RMS_RHO <- 0.01
.RT_RMS_EPS <- 1e-6

# Quantile-grid quadrature for mean, var and crps (midpoint rule).
.RT_QGRID_N <- 32L
.RT_QGRID_P <- ((0:(32L - 1L)) + 0.5) / 32L
.RT_QGRID_W <- 1.0 / 32L

.rt_clamp01 <- function(u) min(max(u, .RT_EPS), 1.0 - .RT_EPS)
.rt_clip <- function(x, lo, hi) if (x < lo) lo else if (x > hi) hi else x
# Standard normal quantile by the Dist bisection, as the reference does.
.rt_phi_inv <- function(p) dist_quantile(dist_gaussian(0.0, 1.0), p)

# --- transformation families: H_t as cdf, quantile, logpdf of (z or p, theta) ---

.rt_identity_family <- list(
  name = "identity",
  cdf = function(z, theta) pnorm(z),
  quantile = function(p, theta) .rt_phi_inv(p),
  logpdf = function(z, theta) -0.5 * z * z - .RT_LOG_SQRT2PI
)

.rt_gss_family <- list(
  name = "gaussian_shift_scale",
  cdf = function(z, theta) pnorm((z - theta$m) / exp(0.5 * theta$ell)),
  quantile = function(p, theta) theta$m + exp(0.5 * theta$ell) * .rt_phi_inv(p),
  logpdf = function(z, theta) {
    sigma <- exp(0.5 * theta$ell)
    zz <- (z - theta$m) / sigma
    -0.5 * zz * zz - log(sigma) - .RT_LOG_SQRT2PI
  }
)

# --- the corrected predictive ---

# The base predictive reshaped by the residual model family(theta). Answers
# the dist_* probes. family is a list(cdf, quantile, logpdf); the identity
# and Gaussian shift-scale families are .rt_identity_family and
# .rt_gss_family. Port of CorrectedDist.
corrected_dist <- function(base, family, theta) {
  list(corrected = TRUE, base = base, family = family, theta = theta,
       cache = new.env(parent = emptyenv()))
}

.cd_z <- function(d, y) .rt_phi_inv(.rt_clamp01(dist_cdf(d$base, y)))

corrected_cdf <- function(d, y) d$family$cdf(.cd_z(d, y), d$theta)

corrected_logpdf <- function(d, y) {
  z <- .cd_z(d, y)
  # chain rule: log f~(y) = log f(y) + log h(z) - log phi(z)
  dist_logpdf(d$base, y) + d$family$logpdf(z, d$theta) + 0.5 * z * z + .RT_LOG_SQRT2PI
}

corrected_quantile <- function(d, p, tol = 1e-9, max_iter = 2100L) {
  stopifnot(p > 0, p < 1)
  hp <- d$family$quantile(p, d$theta)
  u <- .rt_clamp01(pnorm(hp))
  dist_quantile(d$base, u, tol = tol, max_iter = max_iter)
}

.cd_qgrid <- function(d, qfun) {
  if (is.null(d$cache$qgrid)) {
    d$cache$qgrid <- vapply(.RT_QGRID_P, function(p) qfun(d, p), 0.0)
  }
  d$cache$qgrid
}

.qgrid_mean <- function(q) sum(q) * .RT_QGRID_W
.qgrid_var <- function(q) {
  mu <- sum(q) * .RT_QGRID_W
  sum((q - mu)^2) * .RT_QGRID_W
}
.qgrid_std <- function(q) {
  v <- .qgrid_var(q)
  if (v > 0) sqrt(v) else 0.0
}
# CRPS via the quantile-form identity (Gneiting and Raftery 2007).
.qgrid_crps <- function(q, x) {
  total <- 0.0
  for (i in seq_along(q)) {
    ind <- if (x <= q[i]) 1.0 else 0.0
    total <- total + (q[i] - x) * (ind - .RT_QGRID_P[i])
  }
  2.0 * total * .RT_QGRID_W
}

corrected_mean <- function(d) .qgrid_mean(.cd_qgrid(d, corrected_quantile))
corrected_var <- function(d) .qgrid_var(.cd_qgrid(d, corrected_quantile))
corrected_std <- function(d) .qgrid_std(.cd_qgrid(d, corrected_quantile))
corrected_crps <- function(d, x) .qgrid_crps(.cd_qgrid(d, corrected_quantile), x)

# --- RTRL-lite online learners for m1 (scale) and m2 (scale + leverage) ---

.rt_init_scale_learner <- function(learn_b) {
  st <- list(ell = 0.0, a = 0.5, c = 0.05, ea = 0.0, ec = 0.0, g2a = 0.0, g2c = 0.0)
  if (learn_b) {
    st$b <- 0.0
    st$eb <- 0.0
    st$g2b <- 0.0
  }
  st
}

.rt_rms_step <- function(lr, grad, g2_prev) {
  g2_new <- (1.0 - .RT_RMS_RHO) * g2_prev + .RT_RMS_RHO * grad * grad
  c(lr * grad / sqrt(g2_new + .RT_RMS_EPS), g2_new)
}

.rt_update_scale_learner <- function(st, z, lr, learn_b) {
  ell_prev <- st$ell
  a <- st$a
  cc <- st$c
  b <- if (is.null(st$b)) 0.0 else st$b
  grad_ell <- .rt_clip(0.5 - 0.5 * z * z * exp(-ell_prev), -.RT_GRAD_MAX, .RT_GRAD_MAX)
  sa <- .rt_rms_step(lr, grad_ell * st$ea, st$g2a)
  sc <- .rt_rms_step(lr, grad_ell * st$ec, st$g2c)
  a_new <- .rt_clip(a - sa[1], .RT_A_MIN, .RT_A_MAX)
  c_new <- .rt_clip(cc - sc[1], -.RT_BC_MAX, .RT_BC_MAX)
  if (learn_b) {
    sb <- .rt_rms_step(lr, grad_ell * st$eb, st$g2b)
    b_new <- .rt_clip(b - sb[1], -.RT_BC_MAX, .RT_BC_MAX)
  } else {
    b_new <- 0.0
  }
  ea_new <- .rt_clip(ell_prev + a_new * st$ea, -.RT_SENS_MAX, .RT_SENS_MAX)
  ec_new <- .rt_clip((z * z - 1.0) + a_new * st$ec, -.RT_SENS_MAX, .RT_SENS_MAX)
  ell_new <- .rt_clip(a_new * ell_prev + b_new * z + c_new * (z * z - 1.0), .RT_ELL_MIN, .RT_ELL_MAX)
  out <- list(ell = ell_new, a = a_new, c = c_new, ea = ea_new, ec = ec_new, g2a = sa[2], g2c = sc[2])
  if (learn_b) {
    out$b <- b_new
    out$eb <- .rt_clip(z + a_new * st$eb, -.RT_SENS_MAX, .RT_SENS_MAX)
    out$g2b <- sb[2]
  }
  out
}

.rt_theta_of <- function(model, learner) {
  if (model == "m0") {
    return(list(m = 0.0, ell = 0.0))
  }
  theta <- list(m = 0.0, ell = learner$ell, a = learner$a, c = learner$c)
  if (model == "m2") {
    theta$b <- learner$b
  }
  theta
}

# Wrap a k = 1 skater with an online residual-transform correction. The
# returned state carries z (the resolved residual, NULL on the first tick),
# raw (the base forecast) and theta (the H_t parameters), as in Python.
residual_transform <- function(base, model = "m2", lr = 0.01, k = 1L) {
  force(base)
  force(lr)
  stopifnot(k == 1)
  if (!(model %in% c("m0", "m1", "m2"))) {
    stop(sprintf("model must be one of m0, m1, m2, got %s", model))
  }
  learn_b <- model == "m2"
  family <- if (model == "m0") .rt_identity_family else .rt_gss_family
  function(y, state = NULL) {
    if (is.null(state)) {
      learner0 <- if (model == "m0") NULL else .rt_init_scale_learner(learn_b)
      state <- list(base = NULL, raw = NULL, z = NULL, learner = learner0,
                    theta = .rt_theta_of(model, learner0))
    }
    z <- NULL
    if (!is.null(state$raw)) {
      z <- .rt_phi_inv(.rt_clamp01(dist_cdf(state$raw, y)))
      if (model != "m0") {
        state$learner <- .rt_update_scale_learner(state$learner, z, lr, learn_b)
      }
    }
    state["z"] <- list(z)
    r <- base(y, state$base)
    state["base"] <- list(r$state)
    f_t <- r$dists[[1]]
    state$raw <- f_t
    theta <- .rt_theta_of(model, state$learner)
    state$theta <- theta
    list(dists = list(corrected_dist(f_t, family, theta)), state = state)
  }
}

# Online EWMA diagnostics on a z-stream: moments and the leverage statistic
# L_k = E[z_t (z_{t+k}^2 - 1)] for k = 1..K. Port of ResidualDiagnostics as
# a pure value: residual_diagnostics(K, halflife) makes one,
# residual_diagnostics_update(rd, z) returns the updated one, and
# residual_diagnostics_summary(rd) gives the reference's summary dict.
residual_diagnostics <- function(K = 5L, halflife = 50.0) {
  list(K = as.integer(K), alpha = 1.0 - 0.5^(1.0 / halflife), n = 0L,
       m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0, buf = numeric(0),
       L = numeric(K), cross_sq = numeric(K))
}

residual_diagnostics_update <- function(rd, z) {
  rd$n <- rd$n + 1L
  a <- if (rd$alpha > 1.0 / rd$n) rd$alpha else 1.0 / rd$n
  z2 <- z * z
  rd$m1 <- (1 - a) * rd$m1 + a * z
  rd$m2 <- (1 - a) * rd$m2 + a * z2
  rd$m3 <- (1 - a) * rd$m3 + a * z * z2
  rd$m4 <- (1 - a) * rd$m4 + a * z2 * z2
  nb <- length(rd$buf)
  for (lag in seq_len(rd$K)) {
    if (nb >= lag) {
      z_lag <- rd$buf[nb - lag + 1]
      rd$L[lag] <- (1 - a) * rd$L[lag] + a * z_lag * (z2 - 1.0)
      rd$cross_sq[lag] <- (1 - a) * rd$cross_sq[lag] + a * z_lag * z_lag * z2
    }
  }
  rd$buf <- c(rd$buf, z)
  if (length(rd$buf) > rd$K) {
    rd$buf <- rd$buf[(length(rd$buf) - rd$K + 1):length(rd$buf)]
  }
  rd
}

residual_diagnostics_summary <- function(rd) {
  list(n = rd$n, `E[z]` = rd$m1, `E[z2]` = rd$m2, `E[z3]` = rd$m3, `E[z4]` = rd$m4,
       L = rd$L, cross_sq = rd$cross_sq)
}
