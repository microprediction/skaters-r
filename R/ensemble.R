# Precision-weighted ensemble of skaters.
# Port of skaters/ensemble.py.

precision_weighted_ensemble <- function(skaters, k = 1L, floor = 1e-6) {
  force(k)
  force(floor)
  n <- length(skaters)
  stopifnot(n > 0)
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(
        sub = vector("list", n),
        queues = lapply(seq_len(n), function(i) lapply(seq_len(k), function(h) numeric(0))),
        stats = lapply(seq_len(n), function(i) lapply(seq_len(k), function(h) running_var_init()))
      )
    }
    all_dists <- vector("list", n)
    for (i in seq_len(n)) {
      r <- skaters[[i]](y, state$sub[[i]])
      state$sub[[i]] <- r$state
      all_dists[[i]] <- r$dists
    }
    # Resolve pending predictions (Dist means) and update error stats.
    # Horizon h (1-based here) is h-step-ahead: buffer h predictions, then resolve.
    for (i in seq_len(n)) {
      for (h in seq_len(k)) {
        q <- c(state$queues[[i]][[h]], dist_mean(all_dists[[i]][[h]]))
        if (length(q) > h) {
          pred_mean <- q[1]
          q <- q[-1]
          error <- y - pred_mean
          state$stats[[i]][[h]] <- running_var_update(state$stats[[i]][[h]], error)
        }
        state$queues[[i]][[h]] <- q
      }
    }
    combined <- vector("list", k)
    for (h in seq_len(k)) {
      weights <- numeric(n)
      for (i in seq_len(n)) {
        mse <- running_mse_get(state$stats[[i]][[h]])
        w <- if (is.finite(mse) && mse > 0) 1.0 / mse else floor
        weights[i] <- max(w, floor)
      }
      horizon_dists <- lapply(seq_len(n), function(i) all_dists[[i]][[h]])
      combined[[h]] <- dist_combine(horizon_dists, weights)
    }
    list(dists = combined, state = state)
  }
}
