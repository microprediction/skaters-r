# Exact pushforward of a predictive through a monotone coordinate map.
# Port of the table-plus-inverse path of skaters/pushforward.py
# (PushforwardDist with ys, zs and inv_fn; no fwd_fn or analytic Jacobian).
#
# The inner predictive stays in the transformed z-coordinate and queries go
# through the change of variables:
#   F_Y(y)     = F_Z(z(y))
#   log p_Y(y) = log p_Z(z(y)) + log dz/dy
#   Q_Y(p)     = y(Q_Z(p))
# The transforms use it for a tail-spliced predictive, which has no mixture
# components for the component-wise delta method to map.

.PF_SQRT2 <- sqrt(2.0)
.PF_SQRTPI <- sqrt(pi)

# 7-point Gauss-Hermite nodes and weights (physicists').
.PF_GH7_X <- c(
  -2.6519613568352334, -1.6735516287674714, -0.8162878828589647, 0.0,
  0.8162878828589647, 1.6735516287674714, 2.6519613568352334
)
.PF_GH7_W <- c(
  9.717812450995192e-04, 5.451558281912703e-02, 4.256072526101278e-01,
  8.102646175568073e-01, 4.256072526101278e-01, 5.451558281912703e-02,
  9.717812450995192e-04
)

# Fourth-power odds give the tails finer resolution under nonlinear maps.
.pf_quantile_rule <- function(n) {
  i <- 0:n
  edges <- i^4 / (i^4 + (n - i)^4)
  lo <- edges[-length(edges)]
  hi <- edges[-1]
  list(p = (lo + hi) * 0.5, w = hi - lo)
}

.PF_QUANTILE_RULE <- .pf_quantile_rule(257L)

# Piecewise-linear with linear extrapolation from the end segments. The
# segment search matches Python's bisect_left.
interp_table <- function(x, xs, fs) {
  n <- length(xs)
  i <- sum(xs < x) # bisect_left, 0-based insertion point
  if (i <= 0) {
    i <- 1L
  } else if (i >= n) {
    i <- n - 1L
  }
  x0 <- xs[i]
  x1 <- xs[i + 1]
  fs[i] + (fs[i + 1] - fs[i]) * (x - x0) / (x1 - x0)
}

# (ys, zs) knots of a monotone map y = inv_fn(z) on [z_lo, z_hi]. Knots where
# inv_fn fails to strictly increase are dropped so every segment has a
# positive, finite slope.
table_from_map <- function(inv_fn, z_lo, z_hi, nodes = 400L) {
  z_hi <- max(z_hi, z_lo + 1e-9)
  step <- (z_hi - z_lo) / nodes
  ys <- numeric(0)
  zs <- numeric(0)
  for (j in 0:nodes) {
    z <- z_lo + step * j
    y <- inv_fn(z)
    if (length(ys) == 0 || y > ys[length(ys)]) {
      ys <- c(ys, y)
      zs <- c(zs, z)
    }
  }
  if (length(ys) < 2) {
    y0 <- inv_fn(z_lo)
    ys <- c(y0, y0 + 1e-9)
    zs <- c(z_lo, z_lo + 1e-9)
  }
  list(ys = ys, zs = zs)
}

pushforward_new <- function(inner, ys, zs, inv_fn) {
  # The environment caches the mapped quantile grid, built lazily as in Python.
  list(pushforward = TRUE, inner = inner, ys = ys, zs = zs, inv_fn = inv_fn, cache = new.env(parent = emptyenv()))
}

.pf_qgrid <- function(d) {
  if (is.null(d$cache$qgrid)) {
    d$cache$qgrid <- vapply(.PF_QUANTILE_RULE$p, function(p) pushforward_quantile(d, p), 0.0)
  }
  d$cache$qgrid
}

.pf_z <- function(d, y) interp_table(y, d$ys, d$zs)

.pf_y <- function(d, z) d$inv_fn(z)

# dz/dy at y: central difference of the exact inverse, segment secant otherwise.
.pf_slope <- function(d, y) {
  z <- .pf_z(d, y)
  h <- 1e-5 * (1.0 + abs(z))
  dy <- d$inv_fn(z + h) - d$inv_fn(z - h)
  if (dy > 0.0) {
    return(2.0 * h / dy)
  }
  ys <- d$ys
  zs <- d$zs
  i <- sum(ys < y)
  i <- min(max(i, 1L), length(ys) - 1L)
  (zs[i + 1] - zs[i]) / (ys[i + 1] - ys[i])
}

pushforward_logpdf <- function(d, x) {
  z <- .pf_z(d, x)
  if (is.nan(z)) {
    return(NaN)
  }
  if (is.infinite(z)) {
    return(-Inf)
  }
  dist_logpdf(d$inner, z) + log(.pf_slope(d, x))
}

pushforward_cdf <- function(d, x) {
  z <- .pf_z(d, x)
  if (is.nan(z)) {
    return(NaN)
  }
  if (is.infinite(z)) {
    return(if (z < 0) 0.0 else 1.0)
  }
  dist_cdf(d$inner, z)
}

pushforward_quantile <- function(d, p, tol = 1e-9, max_iter = 2100L) {
  .pf_y(d, dist_quantile(d$inner, p, tol = tol, max_iter = max_iter))
}

# CRPS by quadrature over the table in z-space.
pushforward_crps <- function(d, x) {
  if (!is.finite(x)) {
    return(abs(x))
  }
  zstar <- .pf_z(d, x)
  zs <- d$zs
  nz <- length(zs)
  s_lo <- 1.0 / .pf_slope(d, d$ys[1])
  s_hi <- 1.0 / .pf_slope(d, d$ys[length(d$ys)])
  base <- c(
    zs[1] - c(8.0, 6.0, 4.0, 2.0, 1.0), zs, zs[nz] + c(1.0, 2.0, 4.0, 6.0, 8.0),
    zstar - 1.0, zstar, zstar + 1.0
  )
  base <- sort(unique(base))
  nodes <- base[1]
  for (z in base[-1]) {
    lo <- nodes[length(nodes)]
    gap <- z - lo
    if (gap > 0.1) {
      pieces <- min(as.integer(gap / 0.1) + 1L, 40L)
      step <- gap / pieces
      if (pieces > 1) {
        nodes <- c(nodes, lo + step * seq_len(pieces - 1))
      }
    }
    nodes <- c(nodes, z)
  }
  total <- 0.0
  prev_z <- nodes[1]
  c_prev <- dist_cdf(d$inner, prev_z)
  for (z in nodes[-1]) {
    cz <- dist_cdf(d$inner, z)
    ind <- if (prev_z >= zstar) 1.0 else 0.0
    g0 <- (c_prev - ind)^2
    g1 <- (cz - ind)^2
    yp <- if (z <= zs[1]) {
      s_lo
    } else if (prev_z >= zs[nz]) {
      s_hi
    } else {
      (.pf_y(d, z) - .pf_y(d, prev_z)) / (z - prev_z)
    }
    total <- total + 0.5 * (g0 + g1) * (z - prev_z) * yp
    prev_z <- z
    c_prev <- cz
  }
  total
}

# Moment-matched y-space components. A spliced or already transformed inner
# gets one Gaussian fitted to the mapped quantile grid; a plain mixture maps
# each component by centered Gauss-Hermite quadrature.
pushforward_components <- function(d) {
  inner <- d$inner
  if (isTRUE(inner$spliced) || isTRUE(inner$pushforward) || is.null(inner$w) || isTRUE(inner$mixture)) {
    q <- .pf_qgrid(d)
    w <- .PF_QUANTILE_RULE$w
    mu <- sum(w * q)
    sd <- if (!is.finite(mu)) Inf else sqrt(sum(w * (q - mu)^2))
    return(list(w = 1.0, m = mu, s = sd))
  }
  keep <- inner$w > 0.0
  wi <- inner$w[keep]
  mi <- inner$m[keep]
  si <- inner$s[keep]
  m_out <- numeric(length(wi))
  s_out <- numeric(length(wi))
  gw <- .PF_GH7_W / .PF_SQRTPI
  for (j in seq_along(wi)) {
    node_ys <- vapply(mi[j] + si[j] * .PF_SQRT2 * .PF_GH7_X, function(z) .pf_y(d, z), 0.0)
    mu <- sum(gw * node_ys)
    m_out[j] <- mu
    s_out[j] <- if (!is.finite(mu)) Inf else max(sqrt(sum(gw * (node_ys - mu)^2)), 1e-12)
  }
  list(w = wi, m = m_out, s = s_out)
}

pushforward_mean <- function(d) {
  cm <- pushforward_components(d)
  sum(cm$w * cm$m)
}

pushforward_var <- function(d) {
  cm <- pushforward_components(d)
  if (length(cm$w) == 1) {
    return(cm$s[1] * cm$s[1])
  }
  mu <- sum(cm$w * cm$m)
  if (!is.finite(mu)) {
    return(Inf)
  }
  sum(cm$w * (cm$s^2 + (cm$m - mu)^2))
}

# The exact pushforward of d through the monotone map inv_fn. The table window
# adapts to the predictive: [m - 10s, m + 10s] over a mixture, or, with no
# mixture to read (laplace's tail-spliced predictive), its extreme quantiles
# padded by their own span on each side.
.pushforward <- function(d, inv_fn) {
  if (!is.null(d$w) && !isTRUE(d$mixture) && !isTRUE(d$pushforward)) {
    lo <- min(d$m - 10.0 * d$s)
    hi <- max(d$m + 10.0 * d$s)
  } else {
    q_lo <- dist_quantile(d, 1e-6)
    q_hi <- dist_quantile(d, 1.0 - 1e-6)
    lo <- 2.0 * q_lo - q_hi
    hi <- 2.0 * q_hi - q_lo
  }
  tab <- table_from_map(inv_fn, lo, hi)
  pushforward_new(d, tab$ys, tab$zs, inv_fn)
}
