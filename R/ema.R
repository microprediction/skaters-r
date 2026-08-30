# EMA skater: conjugation of a leaf with the EMA transform.
# Port of skaters/ema.py.

#' EMA skater
#'
#' Exponential moving average point forecast with running-variance bands, as
#' a skater `function(y, state)`. Port of `skaters/ema.py`.
#'
#' @param alpha EMA rate: weight given to the newest observation.
#' @param k forecast horizon in steps.
#' @return a skater: a function `f(y, state)` returning `list(dists, state)`
#'   with one predictive mixture per horizon (see [dist_new()]).
#' @examples
#' f <- ema(alpha = 0.1, k = 2)
#' st <- NULL
#' for (y in c(1, 1.2, 0.9, 1.1, 1.0)) {
#'   r <- f(y, st)
#'   st <- r$state
#' }
#' vapply(r$dists, dist_mean, numeric(1))
#' @export
ema <- function(alpha = 0.05, k = 1L) {
  conjugate(leaf(k), ema_transform(alpha), k)
}
