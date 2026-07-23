# Online periodicity detection via running autocorrelation.
# Port of skaters/periodicity.py. O(n_lags) per observation.

.DEFAULT_LAGS <- c(2L, 3L, 4L, 5L, 6L, 7L, 12L, 14L, 24L, 28L, 30L, 52L, 60L, 90L, 168L, 365L)

# Returns a callable: r <- f(y, state); r$scores is a list of c(lag, acf)
# pairs sorted by |acf| descending (stable, so ties keep lag order).
period_detector <- function(lags = NULL, alpha = 0.01, min_observations = 50L) {
  if (is.null(lags)) {
    lags <- .DEFAULT_LAGS
  }
  force(alpha)
  force(min_observations)
  max_lag <- max(lags)

  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(buffer = numeric(0), n = 0L, mean = 0.0, var = 0.0, cross = numeric(length(lags)))
    }
    state$buffer <- c(state$buffer, y)
    state$n <- state$n + 1L

    diff_ <- y - state$mean
    state$mean <- state$mean + alpha * diff_
    state$var <- (1 - alpha) * (state$var + alpha * diff_ * diff_)

    mu <- state$mean
    buf <- state$buffer
    nb <- length(buf)
    for (i in seq_along(lags)) {
      L <- lags[i]
      if (nb > L) {
        cross <- (y - mu) * (buf[nb - L] - mu)
        state$cross[i] <- (1 - alpha) * state$cross[i] + alpha * cross
      }
    }
    if (nb > max_lag + 1) {
      state$buffer <- buf[-1]
    }

    if (state$n < min_observations || state$var < 1e-12) {
      return(list(scores = list(), state = state))
    }
    keep <- which(state$n > lags)
    acf <- if (state$var > 0) state$cross[keep] / state$var else numeric(length(keep))
    ord <- order(-abs(acf), method = "radix") # stable: ties keep lag order
    scores <- lapply(ord, function(j) c(lags[keep[j]], acf[j]))
    list(scores = scores, state = state)
  }
}

# Extract the best periods from detector scores.
top_periods <- function(scores, threshold = 0.3, max_periods = 3L) {
  out <- integer(0)
  for (s in scores[seq_len(min(length(scores), max_periods))]) {
    if (abs(s[2]) >= threshold) out <- c(out, as.integer(s[1]))
  }
  out
}
