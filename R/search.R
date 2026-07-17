# Adaptive search over the transform tree. Port of skaters/search.py.
#
# Beam search over the transform grammar: score candidates by cumulative
# clamped log-likelihood, expand top performers with new transforms
# (replaying recent history so children join warm), prune losers.
# Named adaptive_search here: `search` would mask base::search.
#
# One deliberate departure from the Python reference, which stores live
# callables in the state: here the state holds only recipes and plain
# data (the package convention), and candidate skaters are rebuilt on
# demand from their recipes through a per-instance memo. Rebuilding is
# deterministic, so values match the reference exactly and the state
# survives saveRDS/readRDS.

# The grammar: list of (name, factory, cost) triples.
.SEARCH_TRANSFORMS <- list(
  list("ema_t(0.05)", function() ema_transform(0.05), 1),
  list("ema_t(0.1)", function() ema_transform(0.1), 1),
  list("ema_t(0.3)", function() ema_transform(0.3), 1),
  list("diff", function() difference(), 1),
  list("std(0.05)", function() standardize(0.05), 1),
  list("frac(0.3)", function() fractional_difference(0.3, 30L), 3),
  list("garch", function() garch(), 1),
  list("pow(0.5)", function() power_transform(0.5), 1),
  list("ar(2)", function() ar(2L), 2),
  list("ar(5)", function() ar(5L, decay = 1), 3),
  list("gar(16)", function() grouped_ar(16L), 2),
  list("theta(0.1)", function() theta(0.1), 1),
  list("theta(0.3)", function() theta(0.3), 1),
  list("drift", function() drift(), 1),
  list("drift(0.01)", function() drift(alpha = 0.01, shrinkage = 0.002), 1),
  list("holt(0.1,0.05)", function() holt_linear(0.1, 0.05), 1),
  list("holt(0.3,0.1)", function() holt_linear(0.3, 0.1), 1),
  list("seas(7)", function() seasonal_difference(7L), 1),
  list("seas(12)", function() seasonal_difference(12L), 1),
  list("seas(24)", function() seasonal_difference(24L), 1)
)

# The grammar in play: base transforms plus one seasonal per detected
# period, in detection order (derived from state, never stored in it).
.search_transforms <- function(detected_periods) {
  out <- .SEARCH_TRANSFORMS
  for (p in detected_periods) {
    pp <- p
    out[[length(out) + 1L]] <- list(
      sprintf("seas(%d)", pp),
      local({
        q <- pp
        function() seasonal_difference(q)
      }),
      2
    )
  }
  out
}

.search_entry <- function(depth, recipe, k, cost = 0.0) {
  list(
    s = NULL,
    depth = depth,
    recipe = recipe,
    cost = cost,
    age = 0L,
    warmed = FALSE,
    log_w = rep(0.0, k),
    queues = rep(list(list()), k),
    dists = NULL
  )
}

.search_init_pool <- function(k, cost_budget = Inf) {
  pool <- list()
  e <- .search_entry(0L, character(0), k, cost = 1.0)
  e$warmed <- TRUE
  pool[[1L]] <- e
  for (tr in .SEARCH_TRANSFORMS) {
    cand_cost <- 1.0 + tr[[3]]
    if (cand_cost > cost_budget) {
      next
    }
    e <- .search_entry(1L, tr[[1]], k, cost = cand_cost)
    e$warmed <- TRUE
    pool[[length(pool) + 1L]] <- e
  }
  pool
}

.search_build_from_recipe <- function(recipe, k, transforms) {
  lookup <- new.env(parent = emptyenv())
  for (tr in transforms) {
    assign(tr[[1]], tr[[2]], envir = lookup)
  } # last wins
  f <- leaf(k = k)
  for (t_name in recipe) {
    f <- conjugate(f, get(t_name, envir = lookup)(), k = k)
  }
  f
}

.search_expand <- function(pool, k, top_n, max_depth, transforms, cost_budget) {
  scores <- vapply(pool, function(e) sum(e$log_w) / k, 0.0)
  # Python sorts (score, index) tuples reverse=True: score desc, index desc.
  ord <- order(-scores, -seq_along(pool), method = "radix")
  existing <- vapply(pool, function(e) paste(e$recipe, collapse = "|"), "")
  children <- list()
  for (pi in ord[seq_len(min(top_n, length(ord)))]) {
    parent <- pool[[pi]]
    if (parent$depth >= max_depth) {
      next
    }
    for (tr in transforms) {
      t_name <- tr[[1]]
      t_cost <- tr[[3]]
      child_cost <- parent$cost + t_cost
      if (child_cost > cost_budget) {
        next
      }
      if (
        length(parent$recipe) > 0 &&
          parent$recipe[length(parent$recipe)] == t_name
      ) {
        next
      }
      new_recipe <- c(parent$recipe, t_name)
      key <- paste(new_recipe, collapse = "|")
      if (key %in% existing) {
        next
      }
      existing <- c(existing, key)
      children[[length(children) + 1L]] <- .search_entry(
        length(new_recipe),
        new_recipe,
        k,
        cost = child_cost
      )
    }
  }
  children
}

.search_prune <- function(pool, threshold, max_pool, k) {
  if (length(pool) <= 1) {
    return(pool)
  }
  avg <- function(e) sum(e$log_w) / k
  best <- max(vapply(pool, avg, 0.0))
  i <- 1L
  while (i <= length(pool)) {
    if (avg(pool[[i]]) < best + threshold && length(pool) > 1) {
      pool[[i]] <- NULL
    } else {
      i <- i + 1L
    }
  }
  while (length(pool) > max_pool) {
    worst <- which.min(vapply(pool, avg, 0.0)) # first argmin, as in Python
    pool[[worst]] <- NULL
  }
  pool
}

adaptive_search <- function(
  k = 1L,
  learning_rate = 0.5,
  complexity_penalty = 0.02,
  max_pool = 30L,
  expand_interval = 100L,
  expand_top_n = 3L,
  max_depth = 3L,
  replay_buffer = 500L,
  prune_threshold = -50.0,
  max_components = 20L,
  cost_budget = Inf
) {
  force(k)
  force(learning_rate)
  force(complexity_penalty)
  force(max_pool)
  force(expand_interval)
  force(expand_top_n)
  force(max_depth)
  force(replay_buffer)
  force(prune_threshold)
  force(max_components)
  force(cost_budget)
  pd_func <- period_detector()

  # Per-instance memo: recipe key -> live skater. Skaters are pure
  # functions of the recipe (all mutable state lives in the entry), so
  # a fresh instance resuming a saved state rebuilds identical ones.
  memo <- new.env(parent = emptyenv())
  get_skater <- function(recipe, detected_periods) {
    key <- paste0("r:", paste(recipe, collapse = "|"))
    if (is.null(memo[[key]])) {
      memo[[key]] <- .search_build_from_recipe(
        recipe,
        k,
        .search_transforms(detected_periods)
      )
    }
    memo[[key]]
  }

  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(
        pool = .search_init_pool(k, cost_budget = cost_budget),
        n_obs = 0L,
        buffer = numeric(0),
        pd_state = NULL,
        detected_periods = integer(0)
      )
    }
    state$n_obs <- state$n_obs + 1L
    state$buffer <- c(state$buffer, y)
    if (length(state$buffer) > replay_buffer) {
      state$buffer <- state$buffer[-seq_len(length(state$buffer) - replay_buffer)]
    }
    pool <- state$pool

    # 1. Run all active candidates.
    for (i in seq_along(pool)) {
      f <- get_skater(pool[[i]]$recipe, state$detected_periods)
      r <- f(y, pool[[i]]$s)
      pool[[i]]$s <- r$state
      pool[[i]]$dists <- r$dists
      pool[[i]]$age <- pool[[i]]$age + 1L
    }

    # 2+3. Queue current predictions, resolve matured ones. The Dist issued
    # h+1 steps ago targeted the current y. Bounded loss: clamp logpdf to
    # [-20, 20]; the lower clamp also catches NaN.
    for (i in seq_along(pool)) {
      e <- pool[[i]]
      for (h in seq_len(k)) {
        q <- e$queues[[h]]
        q[[length(q) + 1L]] <- e$dists[[h]]
        if (length(q) > h) {
          past <- q[[1L]]
          q <- q[-1L]
          if (e$warmed) {
            lp <- dist_logpdf(past, y)
            if (is.na(lp) || lp < -20.0) {
              lp <- -20.0
            }
            if (lp > 20.0) {
              lp <- 20.0
            }
            e$log_w[h] <- e$log_w[h] +
              learning_rate * lp -
              complexity_penalty * e$depth
          }
        }
        e$queues[[h]] <- q
      }
      pool[[i]] <- e
    }

    # 3b. Run the period detector.
    rp <- pd_func(y, state$pd_state)
    state$pd_state <- rp$state

    # 4. Periodically expand and prune.
    if (state$n_obs %% expand_interval == 0 && state$n_obs > 10) {
      detected <- top_periods(rp$scores, threshold = 0.3, max_periods = 3L)
      for (period in detected) {
        if (!(period %in% state$detected_periods)) {
          state$detected_periods <- c(state$detected_periods, as.integer(period))
        }
      }
      transforms <- .search_transforms(state$detected_periods)
      children <- .search_expand(pool, k, expand_top_n, max_depth, transforms = transforms, cost_budget = cost_budget)
      for (child in children) {
        # Replay recent history through the child so it joins warm.
        f <- get_skater(child$recipe, state$detected_periods)
        for (yb in state$buffer) {
          r <- f(yb, child$s)
          child$s <- r$state
          child$dists <- r$dists
          child$age <- child$age + 1L
        }
        if (!is.null(child$dists)) {
          for (h in seq_len(k)) {
            child$queues[[h]] <- list(child$dists[[h]])
          }
        }
        child$warmed <- TRUE
        pool[[length(pool) + 1L]] <- child
      }
      pool <- .search_prune(pool, prune_threshold, max_pool, k)
    }

    # 5. Combine predictions via softmax weights.
    combined <- vector("list", k)
    for (h in seq_len(k)) {
      log_ws <- vapply(pool, function(e) e$log_w[h], 0.0)
      max_lw <- max(log_ws)
      weights <- if (is.finite(max_lw)) {
        exp(log_ws - max_lw)
      } else {
        rep(1.0, length(pool))
      }
      d <- dist_combine(lapply(pool, function(e) e$dists[[h]]), weights)
      if (length(d$w) > max_components) {
        d <- dist_prune(d, max_components)
      }
      combined[[h]] <- d
    }

    state$pool <- pool
    list(dists = combined, state = state)
  }
}
