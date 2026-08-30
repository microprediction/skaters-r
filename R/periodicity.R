# Online periodicity detection via running autocorrelation.
# Port of skaters/periodicity.py. O(n_lags) per observation.

.DEFAULT_LAGS <- c(2L, 3L, 4L, 5L, 6L, 7L, 12L, 14L, 24L, 28L, 30L, 52L, 60L, 90L, 168L, 365L)

# Returns a callable: r <- f(y, state); r$scores is a list of c(lag, acf)
# pairs sorted by |acf| descending (stable, so ties keep lag order).
#' Online periodicity detection
#'
#' `period_detector` maintains exponentially weighted autocorrelation
#' estimates at candidate lags, O(n_lags) per observation, returning ranked
#' `(lag, acf)` scores; `top_periods` extracts the significant ones. Port of
#' `skaters/periodicity.py`.
#'
#' @param lags integer candidate lags; a calendar-minded default if `NULL`.
#' @param alpha EMA rate of the running autocorrelation estimates.
#' @param min_observations observations before any score is emitted.
#' @return `period_detector` returns a function `f(y, state)` returning
#'   `list(scores, state)`, where `scores` is a list of `c(lag, acf)` pairs
#'   sorted by absolute autocorrelation. `top_periods` returns an integer
#'   vector of periods.
#' @examples
#' f <- period_detector(lags = c(2L, 3L), min_observations = 10L)
#' st <- NULL
#' for (y in rep(c(1, -1), 20)) {
#'   r <- f(y, st)
#'   st <- r$state
#' }
#' top_periods(r$scores)
#' @rdname periodicity
#' @export
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
#' @param scores the `scores` list from a `period_detector` call.
#' @param threshold minimum absolute autocorrelation to count as a period.
#' @param max_periods number of leading scores considered.
#' @rdname periodicity
#' @export
top_periods <- function(scores, threshold = 0.3, max_periods = 3L) {
  out <- integer(0)
  for (s in scores[seq_len(min(length(scores), max_periods))]) {
    if (abs(s[2]) >= threshold) out <- c(out, as.integer(s[1]))
  }
  out
}
