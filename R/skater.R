#' Wrap a skater closure in an S3 object
#'
#' A skater in this package is a closure `f(y, state)` returning
#' `list(dists, state)`. That closure is the parity surface, verified against
#' the reference implementation's vectors to 1e-6, and it never changes shape.
#' `skater()` wraps any such closure in a small S3 object so the usual R
#' generics work: [observe()] to feed data, [predict.skater()] for interval
#' forecasts, [quantile.skater()] for quantile paths, and `print()`.
#'
#' The wrapper is deliberately the only integration seam. A faster backend
#' producing the same closure contract (for example a future Rust-backed
#' `laplace_fast()`) plugs in here with no change to this file or to callers:
#' `skater(laplace_fast(k))` behaves like `skater(laplace(k))`.
#'
#' @param fn a skater closure, such as the value of [laplace()].
#' @param name optional label used by `print()`.
#' @return an object of class `"skater"`.
#' @examples
#' sk <- skater(laplace(k = 3))
#' sk <- observe(sk, c(0.4, -0.2, 0.1, 0.6, -0.3))
#' predict(sk)
#' quantile(sk, probs = c(0.1, 0.5, 0.9))
#' @export
skater <- function(fn, name = NULL) {
  stopifnot(is.function(fn))
  structure(
    list(fn = fn, state = NULL, dists = NULL, n = 0L,
         name = if (is.null(name)) deparse(substitute(fn)) else name),
    class = "skater"
  )
}

#' Feed observations to a skater
#'
#' Steps the wrapped closure once per element of `y`, in order. Returns the
#' updated object; assign it back, as with any functional update in R.
#'
#' @param object a `"skater"` object.
#' @param y numeric vector of observations, fed in order.
#' @param ... unused.
#' @return the updated `"skater"` object, whose latest predictive
#'   distributions are available through [predict.skater()].
#' @export
observe <- function(object, y, ...) UseMethod("observe")

#' @rdname observe
#' @export
observe.skater <- function(object, y, ...) {
  fn <- object$fn
  st <- object$state
  dists <- object$dists
  for (v in as.numeric(y)) {
    r <- fn(v, st)
    st <- r$state
    dists <- r$dists
  }
  object$state <- st
  object$dists <- dists
  object$n <- object$n + length(y)
  object
}

#' Forecast intervals from a skater
#'
#' Returns one row per horizon with the predictive mean and central intervals.
#' Quantiles are exact inversions of the predictive mixture, so an 80\%
#' interval is the 0.1 and 0.9 quantiles, with no normality assumption.
#'
#' @param object a `"skater"` object that has seen at least one observation.
#' @param level vector of central interval levels in percent.
#' @param ... unused.
#' @return a `data.frame` with columns `h`, `mean`, and `lo`/`hi` pairs per
#'   level.
#' @export
predict.skater <- function(object, level = c(80, 95), ...) {
  dists <- .skater_dists(object)
  out <- data.frame(h = seq_along(dists),
                    mean = vapply(dists, dist_mean, numeric(1)))
  for (lv in level) {
    a <- (1 - lv / 100) / 2
    out[[paste0("lo", lv)]] <-
      vapply(dists, function(d) dist_quantile(d, a), numeric(1))
    out[[paste0("hi", lv)]] <-
      vapply(dists, function(d) dist_quantile(d, 1 - a), numeric(1))
  }
  out
}

#' Quantile paths from a skater
#'
#' @param x a `"skater"` object that has seen at least one observation.
#' @param probs probabilities.
#' @param ... unused.
#' @return a matrix with one row per horizon and one column per probability.
#' @export
quantile.skater <- function(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), ...) {
  dists <- .skater_dists(x)
  m <- vapply(probs, function(p)
    vapply(dists, function(d) dist_quantile(d, p), numeric(1)),
    numeric(length(dists)))
  m <- matrix(m, nrow = length(dists))
  dimnames(m) <- list(h = seq_along(dists), probs = format(probs))
  m
}

#' @export
print.skater <- function(x, ...) {
  cat("<skater>", x$name, "\n")
  cat("  observations seen:", x$n, "\n")
  if (is.null(x$dists)) {
    cat("  no predictive yet; use observe()\n")
  } else {
    cat("  horizons:", length(x$dists), "\n")
    m <- vapply(x$dists, dist_mean, numeric(1))
    cat("  predictive means:", paste(signif(m, 4), collapse = " "), "\n")
  }
  invisible(x)
}

.skater_dists <- function(object) {
  if (is.null(object$dists)) {
    stop("this skater has seen no observations; call observe() first",
         call. = FALSE)
  }
  object$dists
}
