# The prediction parade: online PIT/z calibration diagnostics in the state.
# Port of skaters/parade.py. Pass-through for the forecasts themselves.

.PARADE_EPS <- 1e-12
.STD_NORMAL <- dist_gaussian(0.0, 1.0)

parade <- function(base, k) {
  force(base); force(k)
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(base = NULL, pending = list(),
                    pit = rep(NA_real_, k), z = rep(NA_real_, k))
    }
    pend <- state$pending
    n <- length(pend)
    pit <- rep(NA_real_, k)
    z <- rep(NA_real_, k)
    for (m in seq_len(k)) {
      if (m <= n) {
        d <- pend[[n - m + 1]][[m]]       # issued m steps ago, horizon m
        u <- dist_cdf(d, y)
        if (!is.finite(u)) next           # leave this horizon's entry NA
        u <- min(max(u, .PARADE_EPS), 1.0 - .PARADE_EPS)
        pit[m] <- u
        z[m] <- dist_quantile(.STD_NORMAL, u)
      }
    }
    # Gate the observation before the tree consumes it (magnitude-relative,
    # not sigma-relative): identity on any comfortably representable stream.
    y_fed <- y
    if (is.numeric(y_fed) && is.finite(y_fed)) {
      y_fed <- min(max(y_fed, -1e60), 1e60)
      if (n > 0) {
        d1 <- pend[[n]][[1]]              # the 1-step predictive for y
        if (isTRUE(d1$spliced)) d1 <- d1$body
        mp <- dist_mean(d1); sp <- dist_std(d1)
        if (is.finite(mp) && is.finite(sp)) {
          w <- 1e12 * (1.0 + abs(mp) + sp)
          y_fed <- min(max(y_fed, mp - w), mp + w)
        }
      }
    }
    r <- base(y_fed, state$base)
    state$base <- r$state
    pend[[length(pend) + 1L]] <- r$dists
    if (length(pend) > k) pend <- pend[-1]   # deque maxlen = k
    state$pending <- pend
    state$pit <- pit
    state$z <- z
    list(dists = r$dists, state = state)
  }
}
