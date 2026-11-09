# Online covariance estimation: running (Welford), EMA, and Ledoit-Wolf
# shrinkage toward identity. Port of skaters/cov/.
#
# API pattern (mirrors the Python tuple):
#   r <- f(y, state);  r$mean, r$cov (flat row-major n*n), r$state

running_cov <- function(y, state = NULL) {
  n <- length(y)
  if (is.null(state)) {
    state <- list(n = 0L, mean = numeric(n), C = numeric(n * n))
  }
  state$n <- state$n + 1L
  k <- state$n
  delta <- y - state$mean
  state$mean <- state$mean + delta / k
  delta2 <- y - state$mean
  # C += delta %o% delta2, flat row-major: index (i,j) -> (i-1)*n + j
  state$C <- state$C + as.vector(t(outer(delta, delta2)))
  cov <- if (k < 2) numeric(n * n) else state$C / (k - 1)
  list(mean = state$mean, cov = cov, state = state)
}

ema_cov <- function(y, state = NULL, alpha = 0.05) {
  n <- length(y)
  if (is.null(state)) {
    state <- list(mean = as.numeric(y), cov = numeric(n * n), n = 1L)
    return(list(mean = as.numeric(y), cov = numeric(n * n), state = state))
  }
  state$n <- state$n + 1L
  delta <- y - state$mean
  state$mean <- state$mean + alpha * delta
  state$cov <- (1 - alpha) * (state$cov + alpha * as.vector(t(outer(delta, delta))))
  list(mean = state$mean, cov = state$cov, state = state)
}

# skaters#224: delegate to ema_cov (a rank-one EMA update, PSD by
# construction) and shrink THAT covariance's off-diagonal:
# (1 - rho) * C + rho * diag(C) stays PSD whenever C is. The prior version
# updated each pairwise correlation independently, which is not PSD-safe.
ledoit_wolf_cov <- function(y, state = NULL, alpha = 0.05, shrinkage = 0.5) {
  n <- length(y)
  r <- ema_cov(y, if (is.null(state)) NULL else state$inner, alpha = alpha)
  cov <- r$cov
  shrunk <- numeric(n * n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      idx <- (i - 1L) * n + j
      shrunk[idx] <- if (i == j) cov[idx] else (1 - shrinkage) * cov[idx]
    }
  }
  list(mean = r$mean, cov = shrunk, state = list(inner = r$state))
}
