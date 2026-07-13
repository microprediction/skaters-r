# Bayesian model averaging ensemble with shrinkage and complexity penalty.
# Port of skaters/bayesian.py.

bayesian_ensemble <- function(
  skaters,
  k = 1L,
  learning_rate = 0.5,
  complexity_penalty = 0.0,
  depths = NULL,
  prior_log_weights = NULL,
  max_components = 20L
) {
  force(k)
  force(max_components)
  n <- length(skaters)
  stopifnot(n > 0, learning_rate > 0, learning_rate <= 1, complexity_penalty >= 0)
  if (is.null(depths)) {
    depths <- rep(0.0, n)
  }
  if (is.null(prior_log_weights)) {
    prior_log_weights <- rep(0.0, n)
  }
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(
        sub = vector("list", n),
        queues = lapply(seq_len(n), function(i) lapply(seq_len(k), function(h) list())),
        log_w = lapply(seq_len(n), function(i) rep(prior_log_weights[i], k)),
        n_obs = 0L
      )
    }
    state$n_obs <- state$n_obs + 1L
    all_dists <- vector("list", n)
    for (i in seq_len(n)) {
      r <- skaters[[i]](y, state$sub[[i]])
      state$sub[[i]] <- r$state
      all_dists[[i]] <- r$dists
    }
    # Enqueue current predictions, resolve matured ones (h-step lag at 1-based h).
    for (i in seq_len(n)) {
      for (h in seq_len(k)) {
        q <- state$queues[[i]][[h]]
        q[[length(q) + 1L]] <- all_dists[[i]][[h]]
        if (length(q) > h) {
          past_dist <- q[[1]]
          q <- q[-1]
          lp <- dist_logpdf(past_dist, y)
          # Bounded loss: clamp to [-20, 20]; the isTRUE arm also catches NaN.
          if (isTRUE(lp > 20.0)) {
            lp <- 20.0
          } else if (!isTRUE(lp >= -20.0)) {
            lp <- -20.0
          }
          state$log_w[[i]][h] <- state$log_w[[i]][h] +
            learning_rate * lp -
            complexity_penalty * depths[i]
        }
        state$queues[[i]][[h]] <- q
      }
    }
    combined <- vector("list", k)
    for (h in seq_len(k)) {
      log_ws <- vapply(seq_len(n), function(i) state$log_w[[i]][h], 0.0)
      max_lw <- max(log_ws)
      weights <- if (is.finite(max_lw)) exp(log_ws - max_lw) else rep(1.0, n)
      horizon_dists <- lapply(seq_len(n), function(i) all_dists[[i]][[h]])
      d <- dist_combine(horizon_dists, weights)
      if (length(d$w) > max_components) {
        d <- dist_prune(d, max_components)
      }
      combined[[h]] <- d
    }
    list(dists = combined, state = state)
  }
}
