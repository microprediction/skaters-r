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
  # Python: str(int(x)) for an integral value, else f"{x:.6g}" (C %.6g, so
  # 1e-05, not 1.0e-5).
  if (x == trunc(x)) {
    return(sprintf("%.0f", x))
  }
  sprintf("%.6g", x)
}

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
leaf_spec <- function(k = 1L) list(op = "leaf", k = k)
ema_spec <- function(alpha = 0.05, k = 1L) list(op = "ema", alpha = alpha, k = k)
ensemble_spec <- function(..., k = 1L) list(op = "ensemble", k = k, skaters = list(...))
conjugate_spec <- function(skater_spec, transform_spec) {
  list(op = "conjugate", skater = skater_spec, transform = transform_spec)
}
diff_spec <- function() list(op = "diff")
frac_spec <- function(d = 0.4, window = 50L) list(op = "frac", d = d, window = window)
std_spec <- function(alpha = 0.05) list(op = "std", alpha = alpha)
ema_t_spec <- function(alpha = 0.05) list(op = "ema_t", alpha = alpha)

# Python's name for spec_build.
build <- function(spec) spec_build(spec)

# --- JSON, byte-identical to Python's json.dumps(spec, separators=(",", ":")) ---
#
# Types map as in Python: an R integer is a JSON int, an R double a float
# (written as Python's repr writes it: 0.05, 1.0, 1e-05, 1e+16), a named list
# an object (in its order), an unnamed list or a vector of length other than
# one an array, a string a string (ASCII-escaped), TRUE/FALSE true/false and
# NULL or NA null. from_json reads ints back as integers and floats as
# doubles, so a round trip is exact.

# Exact decimal arithmetic for .py_float_repr: a nonnegative decimal is
# list(d = integer digit vector, e = exponent), value = digits * 10^e.
.dec_parse <- function(txt) {
  parts <- strsplit(txt, "e", fixed = TRUE)[[1]]
  mant <- parts[1]
  ex <- if (length(parts) > 1) as.integer(parts[2]) else 0L
  dot <- as.integer(regexpr(".", mant, fixed = TRUE))
  frac_len <- if (dot > 0) nchar(mant) - dot else 0L
  digits <- as.integer(strsplit(gsub(".", "", mant, fixed = TRUE), "")[[1]])
  list(d = digits, e = ex - frac_len)
}

# printf prints a double's exact binary value when asked for enough digits.
.dec_exact <- function(x) .dec_parse(sprintf("%.1100e", x))

.dec_align <- function(a, b) {
  e <- min(a$e, b$e)
  da <- c(a$d, integer(a$e - e))
  db <- c(b$d, integer(b$e - e))
  n <- max(length(da), length(db))
  list(a = c(integer(n - length(da)), da), b = c(integer(n - length(db)), db), e = e)
}

.dec_cmp <- function(a, b) {
  al <- .dec_align(a, b)
  diff <- al$a - al$b
  nz <- which(diff != 0L)
  if (length(nz) == 0) 0L else sign(diff[nz[1]])
}

.dec_add <- function(a, b) {
  al <- .dec_align(a, b)
  v <- al$a + al$b
  carry <- 0L
  for (i in rev(seq_along(v))) {
    t <- v[i] + carry
    v[i] <- t %% 10L
    carry <- t %/% 10L
  }
  list(d = c(if (carry > 0L) carry else integer(0), v), e = al$e)
}

# a - b for a >= b.
.dec_sub <- function(a, b) {
  al <- .dec_align(a, b)
  v <- al$a - al$b
  borrow <- 0L
  for (i in rev(seq_along(v))) {
    t <- v[i] - borrow
    if (t < 0L) {
      t <- t + 10L
      borrow <- 1L
    } else {
      borrow <- 0L
    }
    v[i] <- t
  }
  list(d = v, e = al$e)
}

.dec_absdiff <- function(a, b) if (.dec_cmp(a, b) >= 0) .dec_sub(a, b) else .dec_sub(b, a)

.dec_half <- function(a) {
  v <- a$d * 5L # x / 2 = 5x / 10
  carry <- 0L
  for (i in rev(seq_along(v))) {
    t <- v[i] + carry
    v[i] <- t %% 10L
    carry <- t %/% 10L
  }
  list(d = c(if (carry > 0L) carry else integer(0), v), e = a$e - 1L)
}

# The reals that round to x > 0: (lo, hi) between the midpoints to its
# neighbours, closed when x's significand is even (round half to even).
.dec_interval <- function(x) {
  ex <- .dec_exact(x)
  down <- .pf_nextafter(x, -1)
  up <- .pf_nextafter(x, 1)
  lo <- .dec_half(.dec_add(ex, .dec_exact(down)))
  # Past the largest finite double, the bound is half an ulp above it.
  hi <- if (is.finite(up)) .dec_half(.dec_add(ex, .dec_exact(up))) else .dec_add(ex, .dec_exact((x - down) / 2))
  hexs <- sprintf("%a", x)
  mant <- sub("p.*$", "", hexs)
  last <- substr(mant, nchar(mant), nchar(mant))
  even <- !grepl(".", mant, fixed = TRUE) || nchar(sub("^.*\\.", "", mant)) < 13L ||
    strtoi(last, 16L) %% 2L == 0L
  list(lo = lo, hi = hi, even = even)
}

.dec_in_interval <- function(d, iv) {
  c_lo <- .dec_cmp(d, iv$lo)
  if (c_lo < 0 || (c_lo == 0 && !iv$even)) {
    return(FALSE)
  }
  c_hi <- .dec_cmp(d, iv$hi)
  c_hi < 0 || (c_hi == 0 && iv$even)
}

.py_float_repr <- function(x) {
  if (is.nan(x)) {
    return("NaN")
  }
  if (is.infinite(x)) {
    return(if (x > 0) "Infinity" else "-Infinity")
  }
  if (x == 0) {
    return(if (1 / x < 0) "-0.0" else "0.0")
  }
  # Shortest digits that read back as x, nearest x among those (Python's
  # repr, David Gay's shortest mode). R's own string-to-double conversion is
  # not correctly rounded, so the read-back is decided exactly, against the
  # rounding interval of x. The nearest p-digit string can fall outside an
  # asymmetric interval (at a power of two) where its neighbour falls inside.
  neg <- x < 0
  ax <- abs(x)
  iv <- .dec_interval(ax)
  ex <- .dec_exact(ax)
  best <- NULL
  for (p in 1:17) {
    near <- .dec_parse(sprintf("%.*e", p - 1L, ax))
    unit <- list(d = 1L, e = near$e)
    cands <- list(near, .dec_add(near, unit))
    if (.dec_cmp(near, unit) > 0) cands[[3]] <- .dec_sub(near, unit)
    ok <- Filter(function(cd) .dec_in_interval(cd, iv), cands)
    if (length(ok)) {
      best <- ok[[1]]
      for (cd in ok[-1]) {
        if (.dec_cmp(.dec_absdiff(cd, ex), .dec_absdiff(best, ex)) < 0) best <- cd
      }
      break
    }
  }
  if (is.null(best)) best <- .dec_parse(sprintf("%.16e", ax))
  dg <- best$d
  nzs <- which(dg != 0L)
  dg <- dg[nzs[1]:nzs[length(nzs)]]
  e <- best$e + (length(best$d) - nzs[length(nzs)]) + length(dg) - 1L
  digits <- paste(dg, collapse = "")
  nd <- nchar(digits)
  out <- if (e >= -4 && e < 16) {
    if (e >= 0) {
      if (nd <= e + 1) {
        paste0(digits, strrep("0", e + 1 - nd), ".0")
      } else {
        paste0(substr(digits, 1, e + 1), ".", substring(digits, e + 2))
      }
    } else {
      paste0("0.", strrep("0", -e - 1), digits)
    }
  } else {
    mant <- if (nd > 1) paste0(substr(digits, 1, 1), ".", substring(digits, 2)) else digits
    paste0(mant, "e", if (e < 0) "-" else "+", sprintf("%02d", abs(e)))
  }
  if (neg) paste0("-", out) else out
}

.json_string <- function(x) {
  cps <- utf8ToInt(enc2utf8(x))
  out <- character(length(cps))
  for (i in seq_along(cps)) {
    cp <- cps[i]
    out[i] <- if (cp == 34L) {
      "\\\""
    } else if (cp == 92L) {
      "\\\\"
    } else if (cp == 10L) {
      "\\n"
    } else if (cp == 13L) {
      "\\r"
    } else if (cp == 9L) {
      "\\t"
    } else if (cp == 8L) {
      "\\b"
    } else if (cp == 12L) {
      "\\f"
    } else if (cp < 32L || (cp > 126L && cp <= 0xFFFF)) {
      sprintf("\\u%04x", cp)
    } else if (cp > 0xFFFF) {
      v <- cp - 0x10000
      sprintf("\\u%04x\\u%04x", 0xD800 + v %/% 0x400, 0xDC00 + v %% 0x400)
    } else {
      intToUtf8(cp)
    }
  }
  paste0("\"", paste(out, collapse = ""), "\"")
}

# A double written so any correctly rounded reader recovers it exactly:
# 17 significant digits, or the integer itself when it is integral and
# exactly representable (as JS writes it).
.json_float17 <- function(x) {
  if (!is.finite(x)) {
    return(.py_float_repr(x))
  }
  if (x == trunc(x) && abs(x) < 2^53 && !(x == 0 && 1 / x < 0)) {
    return(sprintf("%.0f", x))
  }
  sprintf("%.17g", x)
}

.json_scalar <- function(x, digits = NULL) {
  if (is.character(x)) {
    return(if (is.na(x)) "null" else .json_string(x))
  }
  if (is.na(x) && !(is.double(x) && is.nan(x))) {
    return("null")
  }
  if (is.logical(x)) {
    return(if (x) "true" else "false")
  }
  if (is.integer(x)) {
    return(as.character(x))
  }
  if (is.numeric(x)) {
    return(if (is.null(digits)) .py_float_repr(x) else .json_float17(x))
  }
  stop(sprintf("to_json cannot encode a %s", class(x)[1]))
}

# digits = NULL (the default) writes floats as Python's repr does, so a spec
# is byte-identical to json.dumps. digits = 17 writes 17 significant digits,
# exact under any correctly rounded reader and much faster, for large data
# such as a skater state.
to_json <- function(spec, digits = NULL) {
  x <- spec
  if (is.null(x)) {
    return("null")
  }
  if (is.list(x)) {
    if (!is.null(names(x))) {
      keys <- names(x)
      body <- vapply(seq_along(x), function(i) paste0(.json_string(keys[i]), ":", to_json(x[[i]], digits)), "")
      return(paste0("{", paste(body, collapse = ","), "}"))
    }
    body <- vapply(x, to_json, "", digits = digits)
    return(paste0("[", paste(body, collapse = ","), "]"))
  }
  if (length(x) == 1) {
    return(.json_scalar(x, digits))
  }
  paste0("[", paste(vapply(seq_along(x), function(i) .json_scalar(x[i], digits), ""), collapse = ","), "]")
}

# Correctly rounded decimal-to-double. R's own conversion is exact for at
# most 15 significant digits and a decimal exponent within 22 (one exact
# product or quotient); otherwise its answer is moved to the double whose
# rounding interval holds the decimal.
.dec_read <- function(txt) {
  v <- as.numeric(txt)
  if (is.na(v)) {
    return(v)
  }
  neg <- substr(txt, 1, 1) == "-"
  t2 <- tolower(sub("^[-+]", "", txt))
  parts <- strsplit(t2, "e", fixed = TRUE)[[1]]
  mant <- parts[1]
  ex <- if (length(parts) > 1) as.integer(parts[2]) else 0L
  dot <- as.integer(regexpr(".", mant, fixed = TRUE))
  frac_len <- if (dot > 0) nchar(mant) - dot else 0L
  digits <- sub("^0+", "", gsub(".", "", mant, fixed = TRUE))
  if (!nzchar(sub("0+$", "", digits))) {
    return(v) # zero, with its sign
  }
  trail <- nchar(digits) - nchar(sub("0+$", "", digits))
  sig <- nchar(digits) - trail
  e10 <- ex - frac_len + trail
  if (sig <= 15 && abs(e10) <= 22) {
    # an integer below 2^53 times an exact power of ten: one rounding
    r <- as.numeric(paste0(substr(digits, 1, sig), "e", e10))
    return(if (neg) -r else r)
  }
  d <- .dec_parse(t2)
  # A start within a few ulps: R's own conversion can be far off at extreme
  # exponents, so scale the significand by libm powers of ten in two steps.
  m <- as.numeric(paste0("0.", digits))
  e_all <- ex - frac_len + nchar(digits)
  av <- m * 10^(e_all %/% 2) * 10^(e_all - e_all %/% 2)
  if (is.infinite(av)) {
    xm <- .Machine$double.xmax
    return((if (neg) -1 else 1) * (if (.dec_in_interval(d, .dec_interval(xm))) xm else Inf))
  }
  if (av == 0) av <- 2^-1074
  for (i in 1:200) { # R's answer can be several ulps off at extreme exponents
    if (.dec_in_interval(d, .dec_interval(av))) break
    if (.dec_cmp(d, .dec_exact(av)) < 0) {
      if (av == 2^-1074) {
        av <- 0 # below half the smallest subnormal
        break
      }
      av <- .pf_nextafter(av, -1)
    } else {
      av <- .pf_nextafter(av, 1)
    }
  }
  if (neg) -av else av
}

# Parse JSON text (Python's json.loads, including NaN and Infinity) into the
# list form above. Numbers are correctly rounded, so a double written with
# 17 significant digits (or as Python's repr) reads back exactly.
from_json <- function(s) {
  chars <- strsplit(enc2utf8(s), "")[[1]]
  n <- length(chars)
  pos <- 1L
  ws <- function() {
    while (pos <= n && chars[pos] %in% c(" ", "\t", "\n", "\r")) pos <<- pos + 1L
  }
  fail <- function(msg) stop(sprintf("from_json: %s at character %d", msg, pos))
  expect <- function(lit) {
    k <- nchar(lit)
    if (pos + k - 1L > n || paste(chars[pos:(pos + k - 1L)], collapse = "") != lit) fail(paste("expected", lit))
    pos <<- pos + k
  }
  hex4 <- function() {
    if (pos + 3L > n) fail("short \\u escape")
    v <- strtoi(paste(chars[pos:(pos + 3L)], collapse = ""), 16L)
    if (is.na(v)) fail("bad \\u escape")
    pos <<- pos + 4L
    v
  }
  parse_string <- function() {
    pos <<- pos + 1L # opening quote
    cps <- integer(0)
    repeat {
      if (pos > n) fail("unterminated string")
      ch <- chars[pos]
      if (ch == "\"") {
        pos <<- pos + 1L
        break
      }
      if (ch == "\\") {
        e <- chars[pos + 1L]
        pos <<- pos + 2L
        cp <- switch(e,
          "\"" = 34L, "\\" = 92L, "/" = 47L, b = 8L, f = 12L, n = 10L, r = 13L, t = 9L,
          u = {
            v <- hex4()
            if (v >= 0xD800 && v <= 0xDBFF && pos + 1L <= n && chars[pos] == "\\" && chars[pos + 1L] == "u") {
              pos <<- pos + 2L
              lo <- hex4()
              v <- 0x10000 + (v - 0xD800) * 0x400 + (lo - 0xDC00)
            }
            v
          },
          fail("bad escape")
        )
        cps <- c(cps, cp)
      } else {
        cps <- c(cps, utf8ToInt(ch))
        pos <<- pos + 1L
      }
    }
    if (length(cps) == 0) "" else intToUtf8(cps)
  }
  parse_number <- function() {
    start <- pos
    while (pos <= n && grepl("[-+0-9.eE]", chars[pos])) pos <<- pos + 1L
    txt <- paste(chars[start:(pos - 1L)], collapse = "")
    if (grepl("[.eE]", txt)) {
      return(.dec_read(txt))
    }
    v <- .dec_read(txt)
    if (abs(v) <= .Machine$integer.max) as.integer(v) else v
  }
  parse_value <- function() {
    ws()
    if (pos > n) fail("unexpected end")
    ch <- chars[pos]
    if (ch == "{") {
      pos <<- pos + 1L
      out <- list()
      keys <- character(0)
      ws()
      if (chars[pos] == "}") {
        pos <<- pos + 1L
        names(out) <- character(0)
        return(out)
      }
      repeat {
        ws()
        if (chars[pos] != "\"") fail("expected a key")
        key <- parse_string()
        ws()
        expect(":")
        val <- parse_value()
        hit <- match(key, keys)
        if (is.na(hit)) {
          out[length(out) + 1L] <- list(val)
          keys <- c(keys, key)
        } else {
          out[hit] <- list(val) # a repeated key keeps its first position, last value
        }
        ws()
        if (chars[pos] == ",") {
          pos <<- pos + 1L
          next
        }
        expect("}")
        break
      }
      names(out) <- keys
      return(out)
    }
    if (ch == "[") {
      pos <<- pos + 1L
      out <- list()
      ws()
      if (chars[pos] == "]") {
        pos <<- pos + 1L
        return(out)
      }
      repeat {
        out[length(out) + 1L] <- list(parse_value())
        ws()
        if (chars[pos] == ",") {
          pos <<- pos + 1L
          next
        }
        expect("]")
        break
      }
      return(out)
    }
    if (ch == "\"") {
      return(parse_string())
    }
    if (ch == "t") {
      expect("true")
      return(TRUE)
    }
    if (ch == "f") {
      expect("false")
      return(FALSE)
    }
    if (ch == "n") {
      expect("null")
      return(NULL)
    }
    if (ch == "N") {
      expect("NaN")
      return(NaN)
    }
    if (ch == "I") {
      expect("Infinity")
      return(Inf)
    }
    if (ch == "-" && pos + 1L <= n && chars[pos + 1L] == "I") {
      expect("-Infinity")
      return(-Inf)
    }
    parse_number()
  }
  v <- parse_value()
  ws()
  if (pos <= n) fail("extra data")
  v
}
