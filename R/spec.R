# Symbolic specification for skater pipelines. Port of skaters/spec.py.
#
# A spec is a plain list that fully describes how to build a skater. It
# serializes to JSON, compares, and materializes via spec_build(spec).
# Grammar (op key required):
#   leaf(k) | ema(alpha, k) | ensemble(k, skaters) | conjugate(skater, transform)
#   transforms: diff | frac(d, window) | std(alpha) | ema_t(alpha)

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
  if (op == "diff") return(difference())
  if (op == "frac") {
    w <- if (is.null(spec$window)) 50L else spec$window
    return(fractional_difference(d = spec$d, window = w))
  }
  if (op == "std") {
    a <- if (is.null(spec$alpha)) 0.05 else spec$alpha
    return(standardize(alpha = a))
  }
  if (op == "ema_t") return(ema_transform(alpha = spec$alpha))
  stop(sprintf("Unknown transform op: %s", op))
}

.spec_infer_k <- function(spec) {
  if (!is.null(spec$k)) return(spec$k)
  if (!is.null(spec$skater)) return(.spec_infer_k(spec$skater))
  if (!is.null(spec$skaters)) return(.spec_infer_k(spec$skaters[[1]]))
  stop("Cannot infer k from spec")
}

.spec_fmt <- function(x) {
  if (x == as.integer(x)) return(sprintf("%d", as.integer(x)))
  sprintf("%.6g", x)
}

spec_name <- function(spec) {
  op <- spec$op
  if (op == "leaf") return("leaf")
  if (op == "ema") return(sprintf("ema(%s)", .spec_fmt(spec$alpha)))
  if (op == "ensemble") {
    inner <- paste(vapply(spec$skaters, spec_name, ""), collapse = ",")
    return(sprintf("ensemble(%s)", inner))
  }
  if (op == "conjugate") {
    return(sprintf("%s|%s", .spec_transform_name(spec$transform),
                   spec_name(spec$skater)))
  }
  stop(sprintf("Unknown op: %s", op))
}

.spec_transform_name <- function(spec) {
  op <- spec$op
  if (op == "diff") return("diff")
  if (op == "frac") {
    w <- if (is.null(spec$window)) 50L else spec$window
    if (w == 50) return(sprintf("frac(%s)", .spec_fmt(spec$d)))
    return(sprintf("frac(%s,w=%d)", .spec_fmt(spec$d), w))
  }
  if (op == "std") {
    a <- if (is.null(spec$alpha)) 0.05 else spec$alpha
    return(sprintf("std(%s)", .spec_fmt(a)))
  }
  if (op == "ema_t") return(sprintf("ema_t(%s)", .spec_fmt(spec$alpha)))
  stop(sprintf("Unknown transform op: %s", op))
}

# Constructors (convenience, mirror spec.py).
leaf_spec <- function(k = 1L) list(op = "leaf", k = k)
ema_spec <- function(alpha = 0.05, k = 1L) list(op = "ema", alpha = alpha, k = k)
ensemble_spec <- function(..., k = 1L) list(op = "ensemble", k = k, skaters = list(...))
conjugate_spec <- function(skater_spec, transform_spec)
  list(op = "conjugate", skater = skater_spec, transform = transform_spec)
diff_spec <- function() list(op = "diff")
frac_spec <- function(d = 0.4, window = 50L) list(op = "frac", d = d, window = window)
std_spec <- function(alpha = 0.05) list(op = "std", alpha = alpha)
ema_t_spec <- function(alpha = 0.05) list(op = "ema_t", alpha = alpha)
