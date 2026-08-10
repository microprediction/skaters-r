# User-facing API: laplace, the general forecaster.
# Port of skaters/api.py (the candidate population and the composition).

.objective_leaf <- function(objective, scale_alpha = 0.03) {
  force(scale_alpha)
  if (objective == "crps") {
    return(function(k = 1L) crps_leaf(k = k, scale_alpha = scale_alpha))
  }
  if (objective == "likelihood") {
    return(function(k = 1L) scale_mixture_leaf(k = k, scale_alpha = scale_alpha))
  }
  stop(sprintf("objective must be 'crps' or 'likelihood', got %s", objective))
}

# Generate the full candidate population (shared by all policies).
# Returns list(candidates=..., depths=...). Order mirrors api.py exactly.
build_candidates <- function(k, leaf_fn = leaf) {
  force(k)
  force(leaf_fn)
  candidates <- list()
  depths <- numeric(0)
  add <- function(cand, depth) {
    candidates[[length(candidates) + 1L]] <<- cand
    depths <<- c(depths, depth)
  }

  # Depth 0: just noise (baseline)
  add(leaf_fn(k = k), 0)

  # Depth 1: single EMA at various speeds
  for (alpha in c(0.01, 0.05, 0.1, 0.3)) {
    add(conjugate(leaf_fn(k = k), ema_transform(alpha), k = k), 1)
  }

  # Depth 1: differencing + leaf (pure random walk)
  add(conjugate(leaf_fn(k = k), difference(), k = k), 1)

  # Depth 1: drift + leaf, multiple speeds
  for (as in list(c(0.05, 0.01), c(0.01, 0.002), c(0.002, 0.001), c(0.0005, 0.0002))) {
    add(conjugate(leaf_fn(k = k), drift(alpha = as[1], shrinkage = as[2]), k = k), 1)
  }

  # Depth 1: Theta method
  for (a in c(0.05, 0.1, 0.3)) {
    add(conjugate(leaf_fn(k = k), theta(alpha = a), k = k), 1)
  }

  # Depth 1: AR
  add(conjugate(leaf_fn(k = k), ar(1L), k = k), 1)
  add(conjugate(leaf_fn(k = k), ar(2L, decay = 1), k = k), 1)

  # Depth 1: Holt linear
  for (ab in list(c(0.1, 0.02), c(0.1, 0.05), c(0.3, 0.1))) {
    add(conjugate(leaf_fn(k = k), holt_linear(alpha = ab[1], beta = ab[2]), k = k), 1)
  }

  # Depth 1: seasonal differencing (common periods)
  for (period in c(7L, 12L, 24L)) {
    add(conjugate(leaf_fn(k = k), seasonal_difference(period), k = k), 1)
  }

  # Depth 1: hedged seasonal anchor -- phase-EMA blended 50/50 with the
  # seasonal-naive. The naive adapts instantly but is one noisy draw; the
  # phase-EMA averages same-phase noise but lags shifts; the hedge beats either
  # alone on genuinely seasonal series (median -8% CRPS on M4-Hourly) at
  # unmeasurable cost elsewhere. MUST stay in this position: the ensemble aligns
  # depths and weights by candidate INDEX, so inserting elsewhere silently
  # repermutes the pool relative to the reference.
  for (period in c(7L, 12L, 24L)) {
    add(conjugate(leaf_fn(k = k), seasonal_anchor(period), k = k), 1)
  }

  # Depth 2: seasonal differencing + EMA
  for (period in c(7L, 12L, 24L)) {
    for (alpha in c(0.05, 0.1)) {
      add(conjugate(conjugate(leaf_fn(k = k), ema_transform(alpha), k = k), seasonal_difference(period), k = k), 2)
    }
  }

  # Depth 2: differencing + EMA
  for (alpha in c(0.05, 0.1, 0.3)) {
    add(conjugate(conjugate(leaf_fn(k = k), ema_transform(alpha), k = k), difference(), k = k), 2)
  }

  # Depth 2: standardize + EMA
  for (alpha in c(0.05, 0.1)) {
    add(conjugate(conjugate(leaf_fn(k = k), ema_transform(alpha), k = k), standardize(), k = k), 2)
  }

  # Depth 2: fractional diff + EMA
  for (d in c(0.2, 0.4)) {
    add(
      conjugate(
        conjugate(leaf_fn(k = k), ema_transform(0.1), k = k),
        fractional_difference(d = d, window = 30L),
        k = k
      ),
      2
    )
  }

  # Depth 2: drift + EMA
  for (ds in list(c(0.002, 0.001), c(0.0005, 0.0002))) {
    for (a_ema in c(0.05, 0.1)) {
      add(
        conjugate(
          conjugate(leaf_fn(k = k), ema_transform(a_ema), k = k),
          drift(alpha = ds[1], shrinkage = ds[2]),
          k = k
        ),
        2
      )
    }
  }

  # Depth 2: drift + Holt linear
  add(
    conjugate(
      conjugate(leaf_fn(k = k), holt_linear(0.1, 0.05), k = k),
      drift(alpha = 0.001, shrinkage = 0.0005),
      k = k
    ),
    2
  )

  # Depth 2: GARCH + EMA
  add(conjugate(conjugate(leaf_fn(k = k), ema_transform(0.1), k = k), garch(), k = k), 2)

  # Depth 2: power transform + EMA
  add(conjugate(conjugate(leaf_fn(k = k), ema_transform(0.1), k = k), power_transform(0.5), k = k), 2)

  # Depth 2: thinking fast and slow — fast tracker outside, slow scale inside
  fast_trackers <- function() {
    list(
      ema_transform(0.3),
      ema_transform(0.5),
      holt_linear(alpha = 0.4, beta = 0.2),
      ar(1L),
      drift(alpha = 0.05, shrinkage = 0.01),
      difference()
    )
  }
  for (scale_alpha in c(0.02, 0.05)) {
    for (tracker in fast_trackers()) {
      slow_scale <- standardize(alpha = scale_alpha)
      add(conjugate(conjugate(leaf_fn(k = k), slow_scale, k = k), tracker, k = k), 2)
    }
  }

  # Coordinate prior (Yeo-Johnson)
  for (L in c(0.0, 0.5)) {
    for (inner_tx in list(difference(), ema_transform(0.1))) {
      add(conjugate(conjugate(leaf_fn(k = k), inner_tx, k = k), yeo_johnson(L), k = k), 2)
    }
  }

  # Mean-reversion prior (Ornstein-Uhlenbeck), multi-step only (gated on k > 1)
  if (k > 1) {
    for (L in c(0.0, 0.5)) {
      for (kappa in c(0.03, 0.1, 0.3)) {
        add(conjugate(conjugate(leaf_fn(k = k), ou_transform(kappa, 0.02), k = k), yeo_johnson(L), k = k), 2)
      }
    }
  }

  list(candidates = candidates, depths = depths)
}

# One laplace instance on one clock: the likelihood-weighted trunk with a
# terminal leaf, plus the lattice projection.
.laplace_single_scale <- function(k, objective, use_sticky, leaf_arg, scale_alpha) {
  cd <- build_candidates(k)
  f <- terminal_leaf_ensemble(
    cd$candidates,
    k = k,
    leaf_fn = if (!is.null(leaf_arg)) leaf_arg else .objective_leaf(objective, scale_alpha),
    learning_rate = 0.8,
    complexity_penalty = 0.005,
    depths = cd$depths,
    max_components = 20L,
    forget = 0.99
  )
  if (use_sticky) {
    f <- sticky(f, k = k)
  }
  f
}

# The general forecaster: likelihood-weighted trunk, CRPS terminal leaf,
# lattice projection, multi-scale at k > 1, GPD tail splice, parade wrapper.
laplace <- function(
  k = 1L,
  objective = c("crps", "likelihood"),
  sticky = TRUE,
  leaf = NULL,
  scales = NULL,
  scale_alpha = 0.03,
  tails = c("gpd", "gaussian")
) {
  objective <- match.arg(objective)
  tails <- match.arg(tails)
  use_sticky <- sticky
  leaf_arg <- leaf
  f <- multiscale(
    function(kk) .laplace_single_scale(kk, objective, use_sticky, leaf_arg, scale_alpha),
    k = k,
    scales = scales
  )
  if (tails == "gpd") {
    f <- gpdtails(f, k = k)
  } # conditional tail fit
  parade(f, k = k) # PIT/z against the spliced predictive
}
