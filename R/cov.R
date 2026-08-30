# Online covariance estimation: running (Welford), EMA, and Ledoit-Wolf
# shrinkage toward identity. Port of skaters/cov/.
#
# API pattern (mirrors the Python tuple):
#   r <- f(y, state);  r$mean, r$cov (flat row-major n*n), r$state

#' Online covariance estimation
#'
#' Streaming estimators processing one observation vector at a time:
#' `running_cov` (Welford), `ema_cov` (exponentially weighted), and
#' `ledoit_wolf_cov` (EMA correlations shrunk toward identity).
#' Port of `skaters/cov/`.
#'
#' @param y numeric vector, one observation of the multivariate series.
#' @param state the list returned by the previous call, or `NULL` to start.
#' @return `list(mean, cov, state)`: the running mean vector, the covariance
#'   estimate as a flat row-major vector, and the state to pass back in.
#' @examples
#' st <- NULL
#' for (i in 1:10) {
#'   r <- running_cov(c(sin(i), cos(i)), st)
#'   st <- r$state
#' }
#' matrix(r$cov, 2, 2)
#' @rdname cov-estimators
#' @export
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

#' @param alpha EMA rate: weight given to the newest observation.
#' @rdname cov-estimators
#' @export
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

#' @param shrinkage weight of the identity target in the Ledoit-Wolf
#'   shrinkage of the EMA correlation matrix.
#' @rdname cov-estimators
#' @export
ledoit_wolf_cov <- function(y, state = NULL, alpha = 0.05, shrinkage = 0.5) {
  n <- length(y)
  if (is.null(state)) {
    corr <- numeric(n * n)
    corr[seq(1L, n * n, by = n + 1L)] <- 1.0
    state <- list(mean = as.numeric(y), var = numeric(n), corr = corr, n = 1L)
    return(list(mean = as.numeric(y), cov = numeric(n * n), state = state))
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
  shrunk <- numeric(n * n)
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
