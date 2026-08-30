# Terminal-leaf ensemble: mix for the mean, model the residual once.
# Port of skaters/terminal.py.

#' @param leaf_fn leaf factory fitted to the ensemble's residual stream,
#'   e.g. [crps_leaf()].
#' @param forget per-step decay applied to accumulated member weights.
#' @rdname ensembles
#' @export
terminal_leaf_ensemble <- function(
  skaters,
  leaf_fn = crps_leaf,
  k = 1L,
  learning_rate = 0.5,
  complexity_penalty = 0.0,
  depths = NULL,
  prior_log_weights = NULL,
  max_components = 20L,
  forget = 1.0
) {
  force(k)
  force(learning_rate)
  force(complexity_penalty)
  force(max_components)
  force(forget)
  n <- length(skaters)
  stopifnot(n > 0)
  if (is.null(depths)) {
    depths <- numeric(n)
  }
  prior <- if (is.null(prior_log_weights)) numeric(n) else prior_log_weights
  # One terminal leaf per horizon; closures live here, never in state.
  tleafs <- lapply(seq_len(k), function(h) leaf_fn(k = 1L))

  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(
        sub = vector("list", n),
        qdist = lapply(seq_len(n), function(i) list()), # h=1 Dist queue for weighting
        log_w = prior,
        leaf_state = vector("list", k),
        leaf_pred = vector("list", k),
        mean_q = lapply(seq_len(k), function(h) numeric(0))
      )
    }

    all_dists <- vector("list", n)
    for (i in seq_len(n)) {
      r <- skaters[[i]](y, state$sub[[i]])
      state$sub[[i]] <- r$state
      all_dists[[i]] <- r$dists
    }

    # Update model weights from the resolved one-step prediction.
    for (i in seq_len(n)) {
      q <- state$qdist[[i]]
      if (length(q) > 0) {
        lp <- dist_logpdf(q[[1]], y)
        q <- q[-1]
        # Bounded loss: clamp to [-20, 20]; the isTRUE arm also catches NaN.
        if (isTRUE(lp > 20.0)) {
          lp <- 20.0
        } else if (!isTRUE(lp >= -20.0)) {
          lp <- -20.0
        }
        state$log_w[i] <- forget * state$log_w[i] + learning_rate * lp - complexity_penalty * depths[i]
      }
      q[[length(q) + 1L]] <- all_dists[[i]][[1]]
      state$qdist[[i]] <- q
    }

    log_w <- state$log_w
    max_lw <- max(log_w)
    w <- exp(log_w - max_lw)
    tot <- sum(w)

    combined <- vector("list", k)
    for (h in seq_len(k)) {
      mu_h <- sum(vapply(seq_len(n), function(i) w[i] * dist_mean(all_dists[[i]][[h]]), 0.0)) / tot

      # Resolve the h-step-ahead combined-mean prediction made h steps ago
      # into a residual, and update this horizon's terminal leaf.
      mq <- state$mean_q[[h]]
      if (length(mq) >= h) {
        resid <- y - mq[1]
        mq <- mq[-1]
        lr <- tleafs[[h]](resid, state$leaf_state[[h]])
        state$leaf_state[[h]] <- lr$state
        state$leaf_pred[[h]] <- lr$dists[[1]]
      }

      pred <- if (!is.null(state$leaf_pred[[h]])) {
        dist_shift(state$leaf_pred[[h]], mu_h)
      } else {
        # Warm-up: fall back to the candidate mixture until the leaf has data.
        p <- dist_combine(lapply(seq_len(n), function(i) all_dists[[i]][[h]]), w)
        if (length(p$w) > max_components) {
          p <- dist_prune(p, max_components)
        }
        p
      }
      combined[[h]] <- pred
      state$mean_q[[h]] <- c(mq, mu_h)
    }

    list(dists = combined, state = state)
  }
}
