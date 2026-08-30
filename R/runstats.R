# Lightweight online running statistics (Welford).
# Port of skaters/runstats.py.

#' Running statistics (Welford)
#'
#' Lightweight online mean/variance state used by the leaves and ensembles.
#' Port of `skaters/runstats.py`.
#'
#' @return `running_var_init` returns a fresh state; `running_var_update`
#'   the updated state; `running_var_get` `c(mean, variance)` (variance
#'   `Inf` until two observations); `running_mse_get` the running mean
#'   squared error `mean^2 + variance`.
#' @examples
#' st <- running_var_init()
#' for (x in c(1, 2, 3)) st <- running_var_update(st, x)
#' running_var_get(st)
#' running_mse_get(st)
#' @rdname runstats
#' @export
running_var_init <- function() list(n = 0L, mean = 0.0, m2 = 0.0)

#' @param state a state from `running_var_init` or `running_var_update`.
#' @param x the arriving observation.
#' @rdname runstats
#' @export
running_var_update <- function(state, x) {
  n <- state$n + 1L
  delta <- x - state$mean
  mean <- state$mean + delta / n
  delta2 <- x - mean
  m2 <- state$m2 + delta * delta2
  list(n = n, mean = mean, m2 = m2)
}

#' @rdname runstats
#' @export
running_var_get <- function(state) {
  if (state$n < 2L) {
    return(c(state$mean, Inf))
  }
  c(state$mean, state$m2 / (state$n - 1))
}

#' @rdname runstats
#' @export
running_mse_get <- function(state) {
  if (state$n < 1L) {
    return(Inf)
  }
  mv <- running_var_get(state)
  if (!is.finite(mv[2])) {
    return(Inf)
  }
  mv[1] * mv[1] + mv[2]
}
