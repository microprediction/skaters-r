# Lightweight online running statistics (Welford).
# Port of skaters/runstats.py.

running_var_init <- function() list(n = 0L, mean = 0.0, m2 = 0.0)

running_var_update <- function(state, x) {
  n <- state$n + 1L
  delta <- x - state$mean
  mean <- state$mean + delta / n
  delta2 <- x - mean
  m2 <- state$m2 + delta * delta2
  list(n = n, mean = mean, m2 = m2)
}

running_var_get <- function(state) {
  if (state$n < 2L) return(c(state$mean, Inf))
  c(state$mean, state$m2 / (state$n - 1))
}

running_mse_get <- function(state) {
  if (state$n < 1L) return(Inf)
  mv <- running_var_get(state)
  if (!is.finite(mv[2])) return(Inf)
  mv[1] * mv[1] + mv[2]
}
