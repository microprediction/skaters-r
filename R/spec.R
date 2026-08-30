# Symbolic specification for skater pipelines. Port of skaters/spec.py.
#
# A spec is a plain list that fully describes how to build a skater. It
# serializes to JSON, compares, and materializes via spec_build(spec).
# Grammar (op key required):
#   leaf(k) | ema(alpha, k) | ensemble(k, skaters) | conjugate(skater, transform)
#   transforms: diff | frac(d, window) | std(alpha) | ema_t(alpha)

#' Symbolic specification of skater pipelines
#'
#' A spec is a plain list describing how to build a skater; it can be
#' serialized, compared, and materialized with `spec_build`. `spec_name`
#' derives the canonical name (transform chains read left to right, e.g.
#' `"diff|ema_t(0.1)|leaf"`). The `*_spec` functions are constructors for
#' the grammar. Port of `skaters/spec.py`.
#'
#' @param spec a spec list, from the `*_spec` constructors.
#' @return `spec_build` returns the materialized skater; `spec_name` a
#'   string; the constructors return spec lists.
#' @examples
#' sp <- conjugate_spec(leaf_spec(k = 1), diff_spec())
#' spec_name(sp)
#' f <- spec_build(sp)
#' r <- f(1.5, NULL)
#' @rdname spec
#' @export
spec_build <- function(spec) {
  op <- spec$op
  f <- if (op == "leaf") {
    leaf(k = spec$k)
  } else if (op == "ema") {
    ema(alpha = spec$alpha, k = spec$k)
  } else if (op == "ensemble") {
    subs <- lapply(spec$skaters, spec_build)
    fl <- if (is.null(spec$floor)) 1e-6 else spec$floor
    precision_weighted_ensemble(subs, k = spec$k, floor = fl)
  } else if (op == "conjugate") {
    inner <- spec_build(spec$skater)
    tr <- .spec_build_transform(spec$transform)
    conjugate(inner, tr, k = .spec_infer_k(spec$skater))
  } else {
    stop(sprintf("Unknown op: %s", op))
  }
  attr(f, "spec_name") <- spec_name(spec)
  f
}

.spec_build_transform <- function(spec) {
  op <- spec$op
  if (op == "diff") {
    return(difference())
  }
  if (op == "frac") {
    w <- if (is.null(spec$window)) 50L else spec$window
    return(fractional_difference(d = spec$d, window = w))
  }
  if (op == "std") {
    a <- if (is.null(spec$alpha)) 0.05 else spec$alpha
    return(standardize(alpha = a))
  }
  if (op == "ema_t") {
    return(ema_transform(alpha = spec$alpha))
  }
  stop(sprintf("Unknown transform op: %s", op))
}

.spec_infer_k <- function(spec) {
  if (!is.null(spec$k)) {
    return(spec$k)
  }
  if (!is.null(spec$skater)) {
    return(.spec_infer_k(spec$skater))
  }
  if (!is.null(spec$skaters)) {
    return(.spec_infer_k(spec$skaters[[1]]))
  }
  stop("Cannot infer k from spec")
}

.spec_fmt <- function(x) {
  if (x == as.integer(x)) {
    return(sprintf("%d", as.integer(x)))
  }
  sprintf("%.6g", x)
}

#' @rdname spec
#' @export
spec_name <- function(spec) {
  op <- spec$op
  if (op == "leaf") {
    return("leaf")
  }
  if (op == "ema") {
    return(sprintf("ema(%s)", .spec_fmt(spec$alpha)))
  }
  if (op == "ensemble") {
    inner <- paste(vapply(spec$skaters, spec_name, ""), collapse = ",")
    return(sprintf("ensemble(%s)", inner))
  }
  if (op == "conjugate") {
    return(sprintf("%s|%s", .spec_transform_name(spec$transform), spec_name(spec$skater)))
  }
  stop(sprintf("Unknown op: %s", op))
}

.spec_transform_name <- function(spec) {
  op <- spec$op
  if (op == "diff") {
    return("diff")
  }
  if (op == "frac") {
    w <- if (is.null(spec$window)) 50L else spec$window
    if (w == 50) {
      return(sprintf("frac(%s)", .spec_fmt(spec$d)))
    }
    return(sprintf("frac(%s,w=%d)", .spec_fmt(spec$d), w))
  }
  if (op == "std") {
    a <- if (is.null(spec$alpha)) 0.05 else spec$alpha
    return(sprintf("std(%s)", .spec_fmt(a)))
  }
  if (op == "ema_t") {
    return(sprintf("ema_t(%s)", .spec_fmt(spec$alpha)))
  }
  stop(sprintf("Unknown transform op: %s", op))
}

# Constructors (convenience, mirror spec.py).
#' @param k forecast horizon in steps.
#' @rdname spec
#' @export
leaf_spec <- function(k = 1L) list(op = "leaf", k = k)
#' @param alpha EMA rate.
#' @rdname spec
#' @export
ema_spec <- function(alpha = 0.05, k = 1L) list(op = "ema", alpha = alpha, k = k)
#' @param ... member specs of the ensemble.
#' @rdname spec
#' @export
ensemble_spec <- function(..., k = 1L) list(op = "ensemble", k = k, skaters = list(...))
#' @param skater_spec spec of the inner skater.
#' @param transform_spec spec of the transform wrapped around it.
#' @rdname spec
#' @export
conjugate_spec <- function(skater_spec, transform_spec) {
  list(op = "conjugate", skater = skater_spec, transform = transform_spec)
}
#' @rdname spec
#' @export
diff_spec <- function() list(op = "diff")
#' @param d fractional differencing order.
#' @param window truncation window of the fractional filter.
#' @rdname spec
#' @export
frac_spec <- function(d = 0.4, window = 50L) list(op = "frac", d = d, window = window)
#' @rdname spec
#' @export
std_spec <- function(alpha = 0.05) list(op = "std", alpha = alpha)
#' @rdname spec
#' @export
ema_t_spec <- function(alpha = 0.05) list(op = "ema_t", alpha = alpha)
