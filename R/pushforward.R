# Exact pushforward of a predictive through a monotone coordinate map.
# Port of skaters/pushforward.py (PushforwardDist).
#
# The inner predictive stays in the transformed z-coordinate and queries go
# through the change of variables:
#   F_Y(y)     = F_Z(z(y))
#   log p_Y(y) = log p_Z(z(y)) + log dz/dy
#   Q_Y(p)     = y(Q_Z(p))
# The map is a strictly increasing table (ys, zs), an analytic inverse inv_fn,
# or both; fwd_fn, log_jacobian_fn (log dz/dy at z), component_moments_fn
# (mapped mean and std of a Gaussian component) and scaled_inv_fn(z, weight)
# (the inverse times a positive weight, for an inverse that overflows) are
# optional. Consumers that need mixture components read the moment-matched
# compatibility view through dist_components.

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

# math.hypot over a vector: the Euclidean norm without overflow or underflow.
.pf_hypot <- function(v) {
  if (length(v) == 0) {
    return(0.0)
  }
  a <- abs(v)
  if (any(is.infinite(a))) {
    return(Inf)
  }
  if (anyNA(a)) {
    return(NaN)
  }
  m <- max(a)
  if (m == 0) {
    return(0.0)
  }
  m * sqrt(sum((a / m)^2))
}

# math.nextafter(x, direction * Inf) for finite x.
.pf_nextafter <- function(x, direction) {
  if (x == 0) {
    return(direction * 2^-1074)
  }
  a <- abs(x)
  e <- floor(log2(a))
  if (2^e > a) e <- e - 1
  if (2^(e + 1) <= a) e <- e + 1
  ulp <- 2^max(e - 52, -1074)
  away <- (x > 0) == (direction > 0)
  if (!away && a == 2^e && e > -1022) {
    ulp <- ulp / 2 # stepping toward zero from a power of two
  }
  x + direction * ulp
}

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

pushforward_new <- function(inner, ys = NULL, zs = NULL, inv_fn = NULL, fwd_fn = NULL,
                            log_jacobian_fn = NULL, component_moments_fn = NULL,
                            scaled_inv_fn = NULL) {
  # The environment caches the components and the quantile grid, built lazily
  # as in Python.
  list(
    pushforward = TRUE, inner = inner, ys = ys, zs = zs, inv_fn = inv_fn, fwd_fn = fwd_fn,
    log_jacobian_fn = log_jacobian_fn, component_moments_fn = component_moments_fn,
    scaled_inv_fn = scaled_inv_fn, cache = new.env(parent = emptyenv())
  )
}

# --- the map ---

# No clamp (skaters#221): past the table ends interp_table extrapolates
# linearly, which is the map's own behavior there.
.pf_z <- function(d, y) {
  if (!is.null(d$fwd_fn)) {
    return(d$fwd_fn(y))
  }
  interp_table(y, d$ys, d$zs)
}

.pf_y <- function(d, z) {
  if (!is.null(d$inv_fn)) {
    return(d$inv_fn(z))
  }
  interp_table(z, d$zs, d$ys)
}

# dz/dy at y: central difference of the exact map when available, segment
# secant otherwise (end segments extend beyond the table).
.pf_slope <- function(d, y) {
  if (!is.null(d$fwd_fn)) {
    h <- 1e-5 * (1.0 + abs(y))
    za <- d$fwd_fn(y + h)
    zb <- d$fwd_fn(y - h)
    if (is.finite(za) && is.finite(zb) && za > zb) {
      return((za - zb) / (2.0 * h))
    }
  }
  if (!is.null(d$inv_fn)) {
    z <- .pf_z(d, y)
    h <- 1e-5 * (1.0 + abs(z))
    dy <- d$inv_fn(z + h) - d$inv_fn(z - h)
    if (isTRUE(dy > 0.0)) {
      return(2.0 * h / dy)
    }
  }
  ys <- d$ys
  zs <- d$zs
  i <- sum(ys < y)
  i <- min(max(i, 1L), length(ys) - 1L)
  (zs[i + 1] - zs[i]) / (ys[i + 1] - ys[i])
}

# --- queries (exact) ---

pushforward_logpdf <- function(d, x) {
  z <- .pf_z(d, x)
  if (is.nan(z)) {
    return(NaN)
  }
  if (is.infinite(z)) { # x on or past a bound of the map
    return(-Inf)
  }
  if (!is.null(d$log_jacobian_fn)) {
    return(dist_logpdf(d$inner, z) + d$log_jacobian_fn(z))
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

# CRPS by quadrature over quantiles or over the table in z-space.
pushforward_crps <- function(d, x) {
  if (!is.finite(x)) {
    return(abs(x))
  }
  if (!is.null(d$fwd_fn) || is.null(d$zs)) {
    return(.pf_crps_from_quantiles(d, x))
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

# Twice the pinball loss integrated over the quantile rule, including
# observations outside the support. Weights go in first if the unweighted
# difference overflows.
.pf_crps_from_quantiles <- function(d, x) {
  q <- .pf_qgrid(d)
  rule <- .PF_QUANTILE_RULE
  total <- 0.0
  for (i in seq_along(q)) {
    p <- rule$p[i]
    weight <- 2.0 * rule$w[i] * (if (q[i] >= x) (1.0 - p) else p)
    gap <- abs(q[i] - x)
    if (is.finite(gap)) {
      total <- total + weight * gap
    } else {
      weighted_q <- weight * q[i]
      if (is.infinite(q[i]) && !is.null(d$scaled_inv_fn)) {
        weighted_q <- d$scaled_inv_fn(dist_quantile(d$inner, p), weight)
      }
      total <- total + abs(weighted_q - weight * x)
    }
  }
  total
}

.pf_qgrid <- function(d) {
  if (is.null(d$cache$qgrid)) {
    d$cache$qgrid <- vapply(.PF_QUANTILE_RULE$p, function(p) pushforward_quantile(d, p), 0.0)
  }
  d$cache$qgrid
}

# --- affine maps ---

# The law of a*Y + b, still a pushforward of (a reflection of) the same inner.
# a > 0 composes the map; a < 0 reflects the inner with dist_scale(inner, -1)
# so the composed map stays increasing. Port of PushforwardDist.affine.
pushforward_affine <- function(d, a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (a == 0.0 || !is.finite(a) || !is.finite(b)) {
    stop("affine needs a finite nonzero a and a finite b")
  }
  r <- if (a > 0) 1.0 else -1.0 # z' = r * z
  inner <- if (r > 0) d$inner else dist_scale(d$inner, -1.0)
  inv <- d$inv_fn
  fwd <- d$fwd_fn
  lj <- d$log_jacobian_fn
  cmf <- d$component_moments_fn
  sinv <- d$scaled_inv_fn
  log_a <- log(abs(a))
  ys <- NULL
  zs <- NULL
  if (!is.null(d$ys)) {
    ys <- a * d$ys + b
    zs <- r * d$zs
    if (r < 0) {
      ys <- rev(ys)
      zs <- rev(zs)
    }
  }
  pushforward_new(
    inner, ys, zs,
    inv_fn = if (is.null(inv)) NULL else function(z) a * inv(r * z) + b,
    fwd_fn = if (is.null(fwd)) NULL else function(y) r * fwd((y - b) / a),
    log_jacobian_fn = if (is.null(lj)) NULL else function(z) lj(r * z) - log_a,
    component_moments_fn = if (is.null(cmf)) NULL else function(m, s) {
      ms <- cmf(r * m, s)
      c(a * ms[1] + b, abs(a) * ms[2])
    },
    scaled_inv_fn = if (is.null(sinv)) NULL else function(z, w) a * sinv(r * z, w) + w * b
  )
}

pushforward_shift <- function(d, delta) pushforward_affine(d, 1.0, delta)

pushforward_scale <- function(d, factor) {
  if (factor == 0) {
    stop("scale factor must be nonzero")
  }
  pushforward_affine(d, factor, 0.0)
}

# --- Dist-compatibility view ---

# Centered Gauss-Hermite quadrature (skaters#225).
.pf_gaussian_moments <- function(d, m, s) {
  gw <- .PF_GH7_W / .PF_SQRTPI
  node_ys <- vapply(m + s * .PF_SQRT2 * .PF_GH7_X, function(z) .pf_y(d, z), 0.0)
  mu <- sum(gw * node_ys)
  if (!is.finite(mu)) {
    return(c(mu, Inf))
  }
  sd <- .pf_hypot(sqrt(gw) * (node_ys - mu))
  c(mu, max(sd, 1e-12))
}

# Moment-matched y-space components (lazy, cached), as list(w, m, s). A
# spliced, granger-mixture or already transformed inner gets one Gaussian
# fitted to the mapped quantile grid; a plain mixture maps each component.
pushforward_components <- function(d) {
  if (!is.null(d$cache$comps)) {
    return(d$cache$comps)
  }
  inner_comps <- if (isTRUE(d$inner$pushforward)) NULL else dist_components(d$inner)
  if (is.null(inner_comps)) {
    q <- .pf_qgrid(d)
    w <- .PF_QUANTILE_RULE$w
    mu <- sum(w * q)
    sd <- if (!is.finite(mu)) Inf else .pf_hypot(sqrt(w) * (q - mu))
    comps <- list(w = 1.0, m = mu, s = sd)
  } else {
    keep <- inner_comps$w > 0.0
    wi <- inner_comps$w[keep]
    mi <- inner_comps$m[keep]
    si <- inner_comps$s[keep]
    m_out <- numeric(length(wi))
    s_out <- numeric(length(wi))
    for (j in seq_along(wi)) {
      ms <- if (!is.null(d$component_moments_fn)) {
        d$component_moments_fn(mi[j], si[j])
      } else {
        .pf_gaussian_moments(d, mi[j], si[j])
      }
      m_out[j] <- ms[1]
      s_out[j] <- ms[2]
    }
    comps <- list(w = wi, m = m_out, s = s_out)
  }
  d$cache$comps <- comps
  comps
}

pushforward_mean <- function(d) {
  cm <- pushforward_components(d)
  sum(cm$w * cm$m)
}

pushforward_std <- function(d) {
  cm <- pushforward_components(d)
  if (length(cm$w) == 1) {
    return(cm$s[1])
  }
  mu <- sum(cm$w * cm$m)
  if (!is.finite(mu)) {
    return(Inf)
  }
  .pf_hypot(sqrt(cm$w) * vapply(seq_along(cm$w), function(j) .pf_hypot(c(cm$s[j], cm$m[j] - mu)), 0.0))
}

pushforward_var <- function(d) {
  sd <- pushforward_std(d)
  sd * sd
}

# The exact pushforward of d through the monotone map inv_fn. The table window
# adapts to the predictive: [m - 10s, m + 10s] over a mixture, or, with no
# mixture to read (laplace's tail-spliced predictive), its extreme quantiles
# padded by their own span on each side.
.pushforward <- function(d, inv_fn) {
  if (!isTRUE(d$pushforward) && !is.null(dist_components(d))) {
    lo <- min(d$m - 10.0 * d$s)
    hi <- max(d$m + 10.0 * d$s)
  } else {
    q_lo <- dist_quantile(d, 1e-6)
    q_hi <- dist_quantile(d, 1.0 - 1e-6)
    lo <- 2.0 * q_lo - q_hi
    hi <- 2.0 * q_hi - q_lo
  }
  tab <- table_from_map(inv_fn, lo, hi)
  pushforward_new(d, tab$ys, tab$zs, inv_fn = inv_fn)
}

# An inverse_k that emits the exact pushforward.
.exact_inverse_k <- function(inv_fn) {
  force(inv_fn)
  function(dists, tstate) lapply(dists, .pushforward, inv_fn = inv_fn)
}
