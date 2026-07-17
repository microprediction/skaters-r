# Online covariance estimation: running (Welford), EMA, and Ledoit-Wolf
# shrinkage toward identity. Port of skaters/cov/.
#
# API pattern (mirrors the Python tuple):
#   r <- f(y, state);  r$mean, r$cov (flat row-major n*n), r$state

running_cov <- function(y, state = NULL) {
  n <- length(y)
  if (is.null(state)) {
    state <- list(n = 0L, mean = rep(0.0, n), C = rep(0.0, n * n))
  }
  state$n <- state$n + 1L
  k <- state$n
  delta <- y - state$mean
  state$mean <- state$mean + delta / k
  delta2 <- y - state$mean
  # C += delta %o% delta2, flat row-major: index (i,j) -> (i-1)*n + j
  state$C <- state$C + as.vector(t(outer(delta, delta2)))
  cov <- if (k < 2) rep(0.0, n * n) else state$C / (k - 1)
  list(mean = state$mean, cov = cov, state = state)
}

ema_cov <- function(y, state = NULL, alpha = 0.05) {
  n <- length(y)
  if (is.null(state)) {
    state <- list(mean = as.numeric(y), cov = rep(0.0, n * n), n = 1L)
    return(list(mean = as.numeric(y), cov = rep(0.0, n * n), state = state))
  }
  state$n <- state$n + 1L
  delta <- y - state$mean
  state$mean <- state$mean + alpha * delta
  state$cov <- (1 - alpha) * (state$cov + alpha * as.vector(t(outer(delta, delta))))
  list(mean = state$mean, cov = state$cov, state = state)
}

ledoit_wolf_cov <- function(y, state = NULL, alpha = 0.05, shrinkage = 0.5) {
  n <- length(y)
  if (is.null(state)) {
    corr <- rep(0.0, n * n)
    corr[seq(1L, n * n, by = n + 1L)] <- 1.0
    state <- list(mean = as.numeric(y), var = rep(0.0, n), corr = corr, n = 1L)
    return(list(mean = as.numeric(y), cov = rep(0.0, n * n), state = state))
  }
  state$n <- state$n + 1L
  delta <- y - state$mean
  state$mean <- state$mean + alpha * delta
  delta2 <- y - state$mean
  state$var <- (1 - alpha) * state$var + alpha * delta * delta2

  # Update correlations in standardized space, clamped to [-1, 1].
  s <- ifelse(state$var > 1e-16, sqrt(state$var), 1e-8)
  corr <- state$corr
  for (i in seq_len(n)) {
    if (i < n) {
      for (j in (i + 1L):n) {
        z_cross <- (delta[i] / s[i]) * (delta[j] / s[j])
        idx <- (i - 1L) * n + j
        r <- (1 - alpha) * corr[idx] + alpha * z_cross
        r <- max(-1.0, min(1.0, r))
        corr[idx] <- r
        corr[(j - 1L) * n + i] <- r
      }
    }
  }
  state$corr <- corr

  # Shrink correlation toward identity and reconstitute the covariance.
  shrunk <- rep(0.0, n * n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      idx <- (i - 1L) * n + j
      shrunk[idx] <- if (i == j) {
        state$var[i]
      } else {
        (1 - shrinkage) * corr[idx] * s[i] * s[j]
      }
    }
  }
  list(mean = state$mean, cov = shrunk, state = state)
}
