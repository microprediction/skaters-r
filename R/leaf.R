# Centered Gaussian residual leaf (Welford online variance).
# Port of skaters/leaf.py::leaf.

leaf <- function(k = 1L) {
  function(y, state = NULL) {
    if (is.null(state)) state <- list(n = 0L, mean = 0.0, m2 = 0.0)
    n <- state$n + 1L
    delta <- y - state$mean
    mean <- state$mean + delta / n
    m2 <- state$m2 + delta * (y - mean)
    state <- list(n = n, mean = mean, m2 = m2)
    v <- if (n < 2L) Inf else m2 / (n - 1)
    std <- if (is.finite(v) && v > 0) sqrt(v) else max(abs(y), 1e-8)
    d <- dist_gaussian(0.0, std)
    list(dists = rep(list(d), k), state = state)
  }
}
