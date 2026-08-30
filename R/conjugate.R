# Conjugation: transform in front, skater inside, exact inverse behind.
# Port of skaters/conjugate.py.

#' Conjugate a skater by an invertible transform
#'
#' Runs the transform forward on each observation, lets the inner skater
#' forecast the transformed series, and pushes the predictive distributions
#' back through the inverse. The composition operator of the whole package.
#' Port of `skaters/conjugate.py`.
#'
#' @param skater a skater: a function `f(y, state)` returning
#'   `list(dists, state)`, such as [leaf()].
#' @param transform a transform pair `list(forward, inverse_k)`, such as
#'   [difference()]; see the transforms page.
#' @param k forecast horizon in steps; must match the inner skater's.
#' @return a skater `function(y, state)` returning `list(dists, state)` in
#'   the original coordinates.
#' @examples
#' f <- conjugate(leaf(k = 1), difference(), k = 1)
#' st <- NULL
#' for (y in c(10, 10.5, 11.2, 10.9)) {
#'   r <- f(y, st)
#'   st <- r$state
#' }
#' dist_mean(r$dists[[1]])
#' @export
conjugate <- function(skater, transform, k = 1L) {
  force(skater)
  force(transform)
  force(k) # bind eagerly: callers construct in loops
  function(y, state = NULL) {
    if (is.null(state)) {
      state <- list(t_state = NULL, s_state = NULL)
    }
    f <- transform$forward(y, state$t_state)
    state$t_state <- f$state
    r <- skater(f$y, state$s_state)
    state$s_state <- r$state
    stopifnot(length(r$dists) == k)
    list(dists = transform$inverse_k(r$dists, state$t_state), state = state)
  }
}
