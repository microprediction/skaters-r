# Skater state as plain JSON data in the reference's canonical layout
# (docs/STATE.md in microprediction/skaters, skaters/state.py).
#
# state_to_json(state) returns the canonical image of an R skater state as
# plain R data in the form from_json() returns and to_json() writes: a named
# list is a JSON object, an unnamed list a JSON array, a length-one atomic
# a scalar, NULL null. state_from_json(obj) is the inverse: it takes that
# form (from from_json(), or jsonlite::fromJSON(simplifyVector = FALSE)) and
# returns an R state the skater resumes from. A checkpoint written by the
# Python, JS or R port restores in any of them.
#
# R keeps some state in a different shape from the reference (a Dist is
# list(w, m, s), sticky keeps parallel vals and wts, a period detector keeps
# lag-indexed vectors, a numeric list is a vector), and a length-one R vector
# cannot say whether it is a scalar or a list of one. Each state object is
# therefore recognised by its field names and written or read by its own
# rule; everything a rule does not mention follows the plain mapping above.

# --- plain values ---

.st_offend <- function(path, what) {
  stop(sprintf("skater state is not plain data at %s: %s", if (nzchar(path)) path else "<root>", what),
       call. = FALSE)
}

.st_is_dist <- function(x) is.list(x) && identical(names(x), c("w", "m", "s"))

.st_num <- function(x, path) {
  if (is.na(x) && !is.nan(x)) {
    return(NULL)
  }
  if (!is.finite(x)) .st_offend(path, sprintf("non-finite number %s (JSON has no literal for it)", format(x)))
  x
}

.st_scalar <- function(x, path) {
  if (is.character(x)) {
    return(if (is.na(x)) NULL else x)
  }
  if (is.logical(x)) {
    return(if (is.na(x)) NULL else x)
  }
  if (is.numeric(x)) {
    return(.st_num(as.numeric(x), path))
  }
  .st_offend(path, sprintf("a value of type %s", typeof(x)))
}

# An atomic vector as a JSON array, whatever its length.
.st_vec <- function(x, path = "") {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x)) {
    return(lapply(seq_along(x), function(i) .st_out(x[[i]], sprintf("%s[%d]", path, i - 1L))))
  }
  lapply(seq_along(x), function(i) .st_scalar(x[i], sprintf("%s[%d]", path, i - 1L)))
}

# A list of length one that is a list, not the thing it holds.
.st_list <- function(x, path = "", f = .st_out) {
  lapply(seq_along(x), function(i) f(x[[i]], sprintf("%s[%d]", path, i - 1L)))
}

.st_dist_out <- function(d, path) {
  if (isTRUE(d$spliced)) {
    return(list(
      spliced = TRUE, body = .st_dist_out(d$body, paste0(path, ".body")),
      t_lo = d$t_lo, t_up = d$t_up, zeta_lo = d$zeta_lo, zeta_up = d$zeta_up,
      g_lo = d$g_lo, s_lo = d$s_lo, g_up = d$g_up, s_up = d$s_up
    ))
  }
  if (!.st_is_dist(d)) {
    kind <- if (isTRUE(d$pushforward)) "a pushforward" else if (isTRUE(d$mixture)) "a granger mixture" else
      if (isTRUE(d$corrected)) "a corrected dist" else if (isTRUE(d$homogenized)) "a homogenized dist" else "an object"
    .st_offend(path, paste(kind, "(not plain data)"))
  }
  for (v in c(d$w, d$m, d$s)) .st_num(v, path)
  list(components = lapply(seq_along(d$w), function(i) list(d$w[i], d$m[i], d$s[i])))
}

# Generic writer: dispatch a named list on its field names, else map plainly.
.st_out <- function(x, path = "") {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.function(x) || is.environment(x)) {
    .st_offend(path, "a function (closures belong in the skater, not in state)")
  }
  if (is.list(x)) {
    if (isTRUE(x$spliced) || .st_is_dist(x) || isTRUE(x$pushforward) || isTRUE(x$mixture) ||
        isTRUE(x$corrected) || isTRUE(x$homogenized)) {
      return(.st_dist_out(x, path))
    }
    nm <- names(x)
    if (is.null(nm)) {
      return(.st_list(x, path))
    }
    rule <- .st_rule_for(nm, "r")
    if (!is.null(rule)) {
      return(rule$to(x, path))
    }
    return(.st_obj(x, path))
  }
  if (length(x) == 1) {
    return(.st_scalar(x, path))
  }
  .st_vec(x, path)
}

# A named list written field by field with the plain mapping, except the
# fields named in `as` (field -> writer function(x, path)).
.st_obj <- function(x, path, as = list(), drop = character(0)) {
  out <- list()
  for (key in names(x)) {
    if (key %in% drop) next
    v <- x[[key]]
    p <- paste0(path, ".", key)
    f <- as[[key]]
    out[key] <- list(if (is.null(f)) .st_out(v, p) else f(v, p))
  }
  if (length(out) == 0) names(out) <- character(0)
  out
}

# --- reading ---

.st_is_dist_obj <- function(o) {
  is.list(o) && !is.null(names(o)) &&
    (isTRUE(o$spliced) || (identical(names(o), "components") && is.list(o$components)))
}

.st_dist_in <- function(o) {
  if (isTRUE(o$spliced)) {
    return(spliced_new(.st_dist_in(o$body), .st_d(o$t_lo), .st_d(o$t_up), .st_d(o$zeta_lo),
                       .st_d(o$zeta_up), .st_d(o$g_lo), .st_d(o$s_lo), .st_d(o$g_up), .st_d(o$s_up)))
  }
  comps <- o$components
  w <- vapply(comps, function(c) .st_d(c[[1]]), 0.0)
  m <- vapply(comps, function(c) .st_d(c[[2]]), 0.0)
  s <- vapply(comps, function(c) .st_d(c[[3]]), 0.0)
  # Keep the weights bit for bit when they already sum to one (within 1e-9
  # of the compensated sum), as Dist.from_dict does; normalise otherwise.
  total <- .st_fsum(w)
  if (abs(total - 1.0) <= 1e-9) list(w = w, m = m, s = s) else dist_new(w, m, s)
}

# Neumaier compensated sum (math.fsum is exact; this is within an ulp).
.st_fsum <- function(v) {
  s <- 0.0
  cmp <- 0.0
  for (x in v) {
    t <- s + x
    cmp <- cmp + (if (abs(s) >= abs(x)) (s - t) + x else (x - t) + s)
    s <- t
  }
  s + cmp
}

.st_d <- function(v) if (is.null(v)) NA_real_ else as.numeric(v)

# A JSON array of numbers as a double vector (null as NA).
.st_numvec <- function(o) {
  if (is.null(o) || length(o) == 0) {
    return(numeric(0))
  }
  vapply(o, .st_d, 0.0)
}

# Generic reader: numbers become doubles, arrays of numbers vectors, arrays
# of anything else lists, objects named lists dispatched on their keys.
.st_in <- function(o) {
  if (is.null(o)) {
    return(NULL)
  }
  if (is.list(o)) {
    if (.st_is_dist_obj(o)) {
      return(.st_dist_in(o))
    }
    nm <- names(o)
    if (is.null(nm)) {
      if (length(o) > 0 && all(vapply(o, function(v) is.numeric(v) && length(v) == 1, FALSE))) {
        return(.st_numvec(o))
      }
      return(lapply(o, .st_in))
    }
    rule <- .st_rule_for(nm, "py")
    if (!is.null(rule)) {
      return(rule$from(o))
    }
    return(.st_obj_in(o))
  }
  if (is.numeric(o)) {
    return(as.numeric(o))
  }
  o
}

.st_obj_in <- function(o, as = list()) {
  out <- list()
  for (key in names(o)) {
    f <- as[[key]]
    out[key] <- list(if (is.null(f)) .st_in(o[[key]]) else f(o[[key]]))
  }
  if (length(out) == 0) names(out) <- character(0)
  out
}

# A JSON array read as a list (never collapsed to a vector), each element
# by `f`.
.st_list_in <- function(o, f = .st_in) {
  if (is.null(o)) {
    return(list())
  }
  lapply(o, f)
}

# --- rules: state objects recognised by their field names ---

.st_rules <- function() .ST_RULES

.st_rule_for <- function(nm, side) {
  key <- paste(sort(nm), collapse = ",")
  for (rule in .ST_RULES) {
    if (identical(rule[[side]], key)) {
      return(rule)
    }
  }
  NULL
}

.st_rule <- function(r, py = r, to, from) {
  list(r = paste(sort(r), collapse = ","), py = paste(sort(py), collapse = ","), to = to, from = from)
}

# Field helpers for rules: writers take (x, path), readers take the JSON value.
.w_vec <- function(x, path) .st_vec(x, path)
.w_list <- function(x, path) .st_list(x, path)
.w_listvec <- function(x, path) .st_list(x, path, .st_vec)
.w_listlist <- function(x, path) .st_list(x, path, .st_list)
.r_vec <- function(o) .st_numvec(o)
.w_mat <- function(x, path) .st_vec(c(t(x)), path)
.r_mat <- function(o) {
  v <- .st_numvec(o)
  n <- as.integer(round(sqrt(length(v))))
  matrix(v, n, n, byrow = TRUE)
}
.r_int <- function(o) if (is.null(o) || length(o) == 0) integer(0) else as.integer(vapply(o, .st_d, 0.0))
.r_chr <- function(o) if (is.null(o) || length(o) == 0) character(0) else vapply(o, as.character, "")
.r_list <- function(o) .st_list_in(o)
.r_listvec <- function(o) .st_list_in(o, .st_numvec)
.r_listlist <- function(o) .st_list_in(o, .st_list_in)
.r_listlistlist <- function(o) .st_list_in(o, .r_listlist)

# A rule that renames nothing and only fixes the shape of some fields.
.st_shape_rule <- function(fields, w = list(), r = list()) {
  .st_rule(fields,
    to = function(x, path) .st_obj(x, path, as = w),
    from = function(o) .st_obj_in(o, as = r)
  )
}

.st_vec_rule <- function(fields, vecs) {
  .st_shape_rule(fields,
    w = stats::setNames(rep(list(.w_vec), length(vecs)), vecs),
    r = stats::setNames(rep(list(.r_vec), length(vecs)), vecs)
  )
}

# multiscale keys its per-stride tables by the stride as a decimal string;
# R keeps them positionally, strides ascending. A stride is the length of
# its phase list.
.st_ms_to <- function(x, path) {
  strides <- as.character(vapply(x$phase, length, 0L))
  keyed <- function(v, f) {
    out <- list()
    for (i in seq_along(strides)) out[strides[i]] <- list(f(v[[i]], sprintf("%s[%s]", path, strides[i])))
    if (length(out) == 0) names(out) <- character(0)
    out
  }
  latest <- list()
  for (i in seq_along(strides)) {
    if (!is.null(x$latest[[i]])) latest[strides[i]] <- list(.st_list(x$latest[[i]], path))
  }
  if (length(latest) == 0) names(latest) <- character(0)
  list(
    t = x$t,
    phase = keyed(x$phase, function(v, p) .st_list(v, p)),
    pending = keyed(x$pending, function(v, p) .st_list(v, p)),
    latest = latest,
    score = keyed(as.list(x$score), function(v, p) .st_scalar(v, p))
  )
}

.st_ms_from <- function(o) {
  keys <- names(o$phase)
  keys <- keys[order(as.integer(keys))]
  list(
    t = .st_d(o$t),
    phase = lapply(keys, function(k) .r_list(o$phase[[k]])),
    pending = lapply(keys, function(k) .r_list(o$pending[[k]])),
    latest = lapply(keys, function(k) if (is.null(o$latest[[k]])) NULL else .r_list(o$latest[[k]])),
    score = vapply(keys, function(k) .st_d(o$score[[k]]), 0.0, USE.NAMES = FALSE)
  )
}

.ST_RULES <- list(
  # conjugate: an empty transform state is an empty object
  .st_rule(c("t_state", "s_state"),
    to = function(x, path) {
      ts <- if (is.list(x$t_state) && length(x$t_state) == 0) {
        structure(list(), names = character(0))
      } else {
        .st_out(x$t_state, paste0(path, ".t_state"))
      }
      list(t_state = ts, s_state = .st_out(x$s_state, paste0(path, ".s_state")))
    },
    from = function(o) list(t_state = .st_in(o$t_state), s_state = .st_in(o$s_state))
  ),
  # parade
  .st_shape_rule(c("base", "pending", "pit", "z", "skipped"),
    w = list(pending = .w_listlist, pit = .w_vec, z = .w_vec),
    r = list(pending = .r_listlist, pit = .r_vec, z = .r_vec)),
  # gpdtails, a tail entry, one tail
  .st_shape_rule(c("base", "pending", "tails"),
    w = list(pending = .w_listlist, tails = .w_list), r = list(pending = .r_listlist, tails = .r_list)),
  .st_vec_rule(c("up", "lo", "warm", "n", "skipped"), "warm"),
  .st_vec_rule(c("t", "exc", "s1", "nx", "r", "g", "s", "since", "run"), "exc"),
  # sticky: the reference keeps [value, weight] pairs in first-seen order
  .st_rule(c("base", "vals", "wts"), c("base", "counts"),
    to = function(x, path) list(
      base = .st_out(x$base, paste0(path, ".base")),
      counts = lapply(seq_along(x$vals), function(i) list(x$vals[i], x$wts[i]))
    ),
    from = function(o) list(
      base = .st_in(o$base),
      vals = vapply(o$counts, function(p) .st_d(p[[1]]), 0.0),
      wts = vapply(o$counts, function(p) .st_d(p[[2]]), 0.0)
    )
  ),
  # terminal_leaf_ensemble; adapt is null unless adaptive temperature is on.
  # An R state from before adapt existed is read as adapt = null.
  .st_rule(c("sub", "qdist", "log_w", "leaf_state", "leaf_pred", "mean_q", "adapt"),
    to = function(x, path) .st_obj(x, path, as = list(
      sub = .w_list, qdist = .w_listlist, log_w = .w_vec, leaf_state = .w_list, leaf_pred = .w_list,
      mean_q = .w_listvec, adapt = function(v, p) if (is.null(v)) NULL else .st_obj(v, p, as = list(pred_q = .w_list))))[
      c("sub", "qdist", "log_w", "leaf_state", "leaf_pred", "mean_q", "adapt")],
    from = function(o) .st_obj_in(o, as = list(
      sub = .r_list, qdist = .r_listlist, log_w = .r_vec, leaf_state = .r_list, leaf_pred = .r_list,
      mean_q = .r_listvec, adapt = function(v) if (is.null(v)) NULL else .st_obj_in(v, as = list(pred_q = .r_list))))
  ),
  .st_rule(c("sub", "qdist", "log_w", "leaf_state", "leaf_pred", "mean_q"),
           c("sub", "qdist", "log_w", "leaf_state", "leaf_pred", "mean_q", "x_unused"),
    to = function(x, path) {
      x["adapt"] <- list(NULL)
      .st_rule_for(names(x), "r")$to(x, path)
    },
    from = function(o) stop("unreachable")
  ),
  # bayesian_ensemble, precision_weighted_ensemble
  .st_shape_rule(c("sub", "queues", "log_w", "n_obs"),
    w = list(sub = .w_list, queues = function(x, p) .st_list(x, p, .w_listlist), log_w = .w_listvec),
    r = list(sub = .r_list, queues = .r_listlistlist, log_w = .r_listvec)),
  .st_shape_rule(c("sub", "queues", "stats"),
    w = list(sub = .w_list, queues = function(x, p) .st_list(x, p, .w_listvec), stats = .w_listlist),
    r = list(sub = .r_list, queues = function(o) .st_list_in(o, .r_listvec), stats = .r_listlist)),
  # multiscale
  .st_rule(c("t", "phase", "pending", "latest", "score"), to = .st_ms_to, from = .st_ms_from),
  # leaves
  .st_vec_rule(c("v", "w", "n"), "w"),
  .st_vec_rule(c("h", "s2", "n", "omega", "alpha", "beta", "buf", "w", "last_r2"), c("buf", "w")),
  # adaptive search: a pool entry is rebuilt from its recipe
  .st_rule(c("pool", "n_obs", "buffer", "pd_state", "detected_periods"),
    to = function(x, path) .st_obj(x, path, as = list(pool = .w_list, buffer = .w_vec, detected_periods = .w_vec)),
    from = function(o) .st_obj_in(o, as = list(pool = .r_list, buffer = .r_vec, detected_periods = .r_int))),
  .st_rule(c("s", "depth", "recipe", "cost", "age", "warmed", "log_w", "queues", "dists"),
    to = function(x, path) .st_obj(x, path, as = list(
      recipe = .w_vec, log_w = .w_vec, queues = .w_listlist,
      dists = function(v, p) if (is.null(v)) NULL else .st_list(v, p))),
    from = function(o) .st_obj_in(o, as = list(
      recipe = .r_chr, log_w = .r_vec, queues = .r_listlist,
      dists = function(v) if (is.null(v)) NULL else .r_list(v)))),
  # period detector: cross keyed by lag string
  .st_rule(c("buffer", "n", "mean", "var", "cross"),
    to = function(x, path) {
      cr <- list()
      for (key in names(x$cross)) cr[key] <- list(.st_scalar(x$cross[[key]], path))
      if (length(cr) == 0) names(cr) <- character(0)
      out <- .st_obj(x, path, as = list(buffer = .w_vec), drop = "cross")
      out["cross"] <- list(cr)
      out
    },
    from = function(o) {
      out <- .st_obj_in(o, as = list(buffer = .r_vec))
      out$cross <- vapply(o$cross, .st_d, 0.0)
      out
    }
  ),
  # granger (horizon 1): the reference keeps a queue of [laplace, learned]
  # pairs, cum per horizon and phase, and z per horizon
  .st_rule(c("lap", "learned", "pend", "cum", "t", "z", "repeat_frac", "abs_z_mean", "n", "last_y"),
    to = function(x, path) {
      out <- .st_obj(x, path, drop = c("pend", "cum", "z"))
      # one pair of predictive lists (laplace's, the learned expert's)
      out["pend"] <- list(if (is.null(x$pend)) list() else
        list(lapply(x$pend, function(d) list(.st_out(d, paste0(path, ".pend"))))))
      out["cum"] <- list(list(list(.st_vec(x$cum, paste0(path, ".cum")))))
      out["z"] <- list(list(.st_scalar(x$z, paste0(path, ".z"))))
      out[c("lap", "learned", "pend", "cum", "t", "z", "repeat_frac", "abs_z_mean", "n", "last_y")]
    },
    from = function(o) {
      out <- .st_obj_in(o[setdiff(names(o), c("pend", "cum", "z"))])
      np <- length(o$pend)
      out["pend"] <- list(if (np == 0) NULL else lapply(o$pend[[np]], function(ds) .st_in(ds[[1]])))
      out$cum <- .r_vec(o$cum[[1]][[1]])
      out$z <- if (length(o$z) == 0) NA_real_ else .st_d(o$z[[1]])
      out$last_y <- .st_d(o$last_y)
      out
    }
  ),
  # homogenize and its pool, residual_transform needs no rule
  .st_shape_rule(c("base", "raw", "pending", "cell", "correction"),
    w = list(pending = .w_listlist, cell = .w_list, correction = .w_list),
    r = list(pending = .r_listlist, cell = .r_list, correction = .r_list)),
  .st_shape_rule(c("dynamic", "issued", "n"), w = list(dynamic = .w_list), r = list(dynamic = .r_list)),
  .st_vec_rule(c("weights", "scales"), c("weights", "scales")),
  # transforms with vector fields
  .st_vec_rule("buffer", "buffer"),
  .st_shape_rule(c("ema", "buffer", "n"), w = list(ema = .w_list, buffer = .w_vec),
                 r = list(ema = .r_list, buffer = .r_vec)),
  .st_vec_rule(c("mu_g", "var_g", "mu", "var", "vis", "n"), c("mu", "var", "vis")),
  # ar, grouped_ar: R keeps P as a matrix, the reference as a flat row-major list
  .st_shape_rule(c("buffer", "phi", "phi_used", "P", "n"),
    w = list(buffer = .w_vec, phi = .w_vec, phi_used = .w_vec, P = .w_mat),
    r = list(buffer = .r_vec, phi = .r_vec, phi_used = .r_vec, P = .r_mat)),
  .st_shape_rule(c("buffer", "theta", "phi_used", "P", "n"),
    w = list(buffer = .w_vec, theta = .w_vec, phi_used = .w_vec, P = .w_mat),
    r = list(buffer = .r_vec, theta = .r_vec, phi_used = .r_vec, P = .r_mat)),
  # gaussianize: knots are [cs, zs, ts], ts null when unsmoothed
  .st_rule(c("buf", "n", "mu", "m2", "knots"),
    to = function(x, path) {
      out <- .st_obj(x, path, as = list(buf = .w_vec), drop = "knots")
      kn <- x$knots
      out["knots"] <- list(if (is.null(kn)) NULL else list(
        .st_vec(kn$cs), .st_vec(kn$zs), if (is.null(kn$ts)) NULL else .st_vec(kn$ts)))
      out
    },
    from = function(o) {
      out <- .st_obj_in(o, as = list(buf = .r_vec))
      kn <- o$knots
      out["knots"] <- list(if (is.null(kn)) NULL else list(
        cs = .r_vec(kn[[1]]), zs = .r_vec(kn[[2]]), ts = if (is.null(kn[[3]])) NULL else .r_vec(kn[[3]])))
      out
    }
  )
)

# --- exported API ---

# Canonical plain-data image of a skater state (or of a list of Dists, a
# skater's output). Raises an error naming the path of anything that is not
# plain data. Port of skaters.state.state_to_json.
state_to_json <- function(state) .st_out(state, "")

# Inverse of state_to_json: takes the canonical layout (as from_json() or
# jsonlite::fromJSON(simplifyVector = FALSE) parse it) and returns a state
# the R skater resumes from. Port of skaters.state.state_from_json.
state_from_json <- function(obj) .st_in(obj)

# Raise an error if state holds anything that would not survive JSON.
assert_plain_state <- function(state) {
  .st_out(state, "")
  invisible(state)
}
