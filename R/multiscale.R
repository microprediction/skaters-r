# Multi-scale ensemble: combine forecasters running on decimated clocks.
# Port of skaters/multiscale.py.

multiscale <- function(base, k, scales = NULL, forget = 0.99, max_components = 20L) {
  force(base)
  force(forget)
  force(max_components)
  if (is.null(scales)) {
    scales <- sort(unique(c(1L, as.integer(ceiling(sqrt(k))), as.integer(k))))
  }
  scales <- sort(unique(as.integer(scales[scales >= 1 & scales <= k])))
  stopifnot(length(scales) > 0, scales[1] == 1L)
  if (length(scales) == 1L) {
    return(base(k))
  }
  ns <- length(scales)
  subs <- lapply(scales, function(s) base(max(1L, as.integer(ceiling(k / s)))))

  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(
        t = 0L,
        phase = lapply(scales, function(s) vector("list", s)),
        pending = lapply(scales, function(s) vector("list", s)),
        latest = vector("list", ns),
        score = rep(NA_real_, ns) # NA = no score yet (Python None)
      )
    }
    t <- state$t
    for (i in seq_len(ns)) {
      s <- scales[i]
      ph <- (t %% s) + 1L
      prev <- state$pending[[i]][[ph]]
      if (!is.null(prev)) {
        lp <- max(dist_logpdf(prev, y), -20.0)
        m <- state$score[i]
        state$score[i] <- if (is.na(m)) lp else forget * m + (1.0 - forget) * lp
      }
      r <- subs[[i]](y, state$phase[[i]][[ph]])
      state$phase[[i]][[ph]] <- r$state
      state$pending[[i]][[ph]] <- r$dists[[1]]
      state$latest[[i]] <- r$dists
    }
    state$t <- t + 1L

    scores <- state$score
    top <- if (any(!is.na(scores))) max(scores[!is.na(scores)]) else 0.0
    out <- vector("list", k)
    for (h in seq_len(k)) {
      fcs <- list()
      wts <- numeric(0)
      for (i in seq_len(ns)) {
        s <- scales[i]
        if (s > h || is.null(state$latest[[i]])) {
          next
        }
        j <- max(1L, as.integer(floor(h / s + 0.5))) # half-up, portable
        dists <- state$latest[[i]]
        if (j > length(dists)) {
          next
        }
        fcs[[length(fcs) + 1L]] <- dists[[j]]
        m <- scores[i]
        wts <- c(wts, exp((if (is.na(m)) top else m) - top))
      }
      # One eligible scale: pass its Dist through untouched.
      out[[h]] <- if (length(fcs) == 1L) {
        fcs[[1]]
      } else {
        dist_prune(dist_combine(fcs, wts), max_components)
      }
    }
    list(dists = out, state = state)
  }
}
