# Terminal-leaf ensemble: mix for the mean, model the residual once.
# Port of skaters/terminal.py.

.PIT_DISPERSION_TARGET <- 1.0 / 12.0 # E[(u - 0.5)^2] for u ~ Uniform(0, 1)
.PIT_EWMA_DECAY <- 0.98
.LOSS_FAST_DECAY <- 0.90
.LOSS_SLOW_DECAY <- 0.995
.TEMP_AMBIG_EPS <- 0.15 * .PIT_DISPERSION_TARGET # below this, the loss slope decides

terminal_leaf_ensemble <- function(
  skaters,
  leaf_fn = crps_leaf,
  k = 1L,
  learning_rate = 0.5,
  complexity_penalty = 0.0,
  depths = NULL,
  prior_log_weights = NULL,
  max_components = 20L,
  forget = 1.0,
  adaptive_temperature = FALSE,
  temp_rho = 0.03,
  temp_band = 3.0,
  gate_hook = NULL
) {
  force(adaptive_temperature)
  force(temp_rho)
  force(temp_band)
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
        mean_q = lapply(seq_len(k), function(h) numeric(0)),
        # adaptive selection temperature (skaters#213/#215), or NULL
        adapt = if (adaptive_temperature) {
          list(log_eta = log(learning_rate), pit_ewma = .PIT_DISPERSION_TARGET,
               loss_fast = NULL, loss_slow = NULL, pred_q = list())
        } else {
          NULL
        }
      )
    }

    all_dists <- vector("list", n)
    for (i in seq_len(n)) {
      r <- skaters[[i]](y, state$sub[[i]])
      state$sub[[i]] <- r$state
      all_dists[[i]] <- r$dists
    }

    # Adaptive temperature: resolve the issued predictive from one step ago
    # against y; its PIT dispersion (or, in the dead band, its own logpdf
    # trend) heats or cools eta for this step's weight update.
    eta_t <- learning_rate
    if (adaptive_temperature) {
      ad <- state$adapt
      if (length(ad$pred_q) > 0) {
        prev_pred <- ad$pred_q[[1]]
        ad$pred_q <- ad$pred_q[-1]
        lp_c <- dist_logpdf(prev_pred, y)
        if (!isTRUE(lp_c >= -20.0)) {
          lp_c <- -20.0
        } else if (lp_c > 20.0) {
          lp_c <- 20.0
        }
        u <- dist_cdf(prev_pred, y)
        if (!isTRUE(u >= 0.0 && u <= 1.0)) {
          u <- 0.5 # degenerate or NaN: no signal
        }
        ad$pit_ewma <- .PIT_EWMA_DECAY * ad$pit_ewma + (1.0 - .PIT_EWMA_DECAY) * (u - 0.5)^2
        if (is.null(ad$loss_fast)) {
          ad$loss_fast <- lp_c
          ad$loss_slow <- lp_c
        } else {
          ad$loss_fast <- .LOSS_FAST_DECAY * ad$loss_fast + (1.0 - .LOSS_FAST_DECAY) * lp_c
          ad$loss_slow <- .LOSS_SLOW_DECAY * ad$loss_slow + (1.0 - .LOSS_SLOW_DECAY) * lp_c
        }
        e_t <- ad$pit_ewma - .PIT_DISPERSION_TARGET
        dloss_t <- ad$loss_fast - ad$loss_slow
        direction <- if (abs(e_t) > .TEMP_AMBIG_EPS) {
          if (e_t > 0.0) -1.0 else 1.0
        } else if (dloss_t > 0.0) {
          1.0
        } else if (dloss_t < 0.0) {
          -1.0
        } else {
          0.0
        }
        ad$log_eta <- ad$log_eta + temp_rho * direction
        ad$log_eta <- min(max(ad$log_eta, log(learning_rate / temp_band)), log(learning_rate * temp_band))
      }
      eta_t <- exp(ad$log_eta)
      state$adapt <- ad
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
        state$log_w[i] <- forget * state$log_w[i] + eta_t * lp - complexity_penalty * depths[i]
      }
      q[[length(q) + 1L]] <- all_dists[[i]][[1]]
      state$qdist[[i]] <- q
    }

    log_w <- state$log_w
    # granger's gate shifts the combining weights for this step only; the
    # recursion above is untouched. laplace passes no hook.
    off <- if (is.null(gate_hook)) NULL else gate_hook()
    if (!is.null(off)) log_w <- log_w + off
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
    if (adaptive_temperature) {
      state$adapt$pred_q[[length(state$adapt$pred_q) + 1L]] <- combined[[1]]
    }

    list(dists = combined, state = state)
  }
}
