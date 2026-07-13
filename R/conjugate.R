# Conjugation: transform in front, skater inside, exact inverse behind.
# Port of skaters/conjugate.py.

conjugate <- function(skater, transform, k = 1L) {
  function(y, state = NULL) {
    if (is.null(state)) state <- list(t_state = NULL, s_state = NULL)
    f <- transform$forward(y, state$t_state)
    state$t_state <- f$state
    r <- skater(f$y, state$s_state)
    state$s_state <- r$state
    stopifnot(length(r$dists) == k)
    list(dists = transform$inverse_k(r$dists, state$t_state), state = state)
  }
}
