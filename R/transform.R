# Invertible transforms: (forward, inverse_k) pairs.
# Ports of skaters/transform.py::difference, ema_transform.

difference <- function() {
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) return(list(y = 0.0, state = list(last = y)))
    list(y = y - tstate$last, state = list(last = y))
  }
  inverse_k <- function(dists, tstate) {
    anchor <- tstate$last
    out <- vector("list", length(dists))
    cm <- 0.0; cv <- 0.0
    for (i in seq_along(dists)) {
      d <- dists[[i]]
      cm <- cm + dist_mean(d)
      cv <- cv + dist_var(d)
      std <- if (cv > 0) sqrt(cv) else max(dist_std(d), 1e-12)
      out[[i]] <- dist_gaussian(anchor + cm, std)
    }
    out
  }
  list(forward = forward, inverse_k = inverse_k)
}

ema_transform <- function(alpha = 0.05) {
  stopifnot(alpha > 0, alpha < 1)
  forward <- function(y, tstate = NULL) {
    if (is.null(tstate)) return(list(y = 0.0, state = list(level = y)))
    residual <- y - tstate$level
    list(y = residual, state = list(level = tstate$level + alpha * residual))
  }
  inverse_k <- function(dists, tstate) {
    lapply(dists, dist_shift, delta = tstate$level)
  }
  list(forward = forward, inverse_k = inverse_k)
}
