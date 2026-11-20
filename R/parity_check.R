# Parity engine, called by BOTH entry points: tests/parity.R (script, repo root
# or R CMD check) and tests/testthat/test-parity.R. Scenarios, probes and the
# 1e-6 tolerance are the five-port contract; this is the original script's
# assertion engine wrapped in a function VERBATIM (a local `cat` shadow silences
# it under testthat), returning a structured result so a framework failure can
# name the component that drifted. A migration that rewrites assertions instead
# of porting them is how parity goes falsely green; this one ports them.
run_parity <- function(vec_path, verbose = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read the parity vectors")
  }
  cat <- if (verbose) base::cat else function(...) invisible(NULL)
  result <- list(scenarios = list(), values_checked = 0L,
                 structure_ok = TRUE, fails = 0L)
  v <- jsonlite::fromJSON(vec_path, simplifyVector = FALSE)
  series <- unlist(v$series)
  ATOL <- 1e-6
  RTOL <- 1e-6
  probe <- function(d, p, qlo, qhi) {
    c(
      dist_mean(d),
      dist_std(d),
      dist_logpdf(d, p),
      dist_cdf(d, p),
      dist_quantile(d, qlo),
      dist_quantile(d, qhi),
      dist_crps(d, p)
    )
  }
  scenarios <- list()
  for (k in c(1L, 3L)) {
    suf <- if (k == 1L) "" else sprintf("_k%d", k)
    add <- function(name, sk) scenarios[[paste0(name, suf)]] <<- list(k = k, sk = sk)
    add("leaf", leaf(k))
    add("diff", conjugate(leaf(k), difference(), k))
    add("ema_t", conjugate(leaf(k), ema_transform(0.1), k))
    add("standardize", conjugate(leaf(k), standardize(), k))
    add("theta", conjugate(leaf(k), theta(0.1), k))
    add("drift", conjugate(leaf(k), drift(0.05, 0.01), k))
    add("holt", conjugate(leaf(k), holt_linear(0.1, 0.05), k))
    add("garch", conjugate(leaf(k), garch(), k))
    add("seasonal", conjugate(leaf(k), seasonal_difference(7L), k))
    add("power", conjugate(leaf(k), power_transform(0.5), k))
    add("ar1", conjugate(leaf(k), ar(1L), k))
    add("ar2", conjugate(leaf(k), ar(2L, decay = 1), k))
    add("frac", conjugate(leaf(k), fractional_difference(0.4, 30L), k))
    add("grouped_ar", conjugate(leaf(k), grouped_ar(8L), k))
    add("yeojohnson_log", conjugate(leaf(k), yeo_johnson(0.0), k))
    add("yeojohnson_half", conjugate(leaf(k), yeo_johnson(0.5), k))
    add("ou", conjugate(leaf(k), ou_transform(0.1), k))
    add("ou_sqrt", conjugate(conjugate(leaf(k), ou_transform(0.1), k), yeo_johnson(0.5), k))
    add("seasonal_scale", conjugate(leaf(k), seasonal_scale(7L), k))
    add("seasonal_scale_c", conjugate(leaf(k), seasonal_scale(7L, 0.1, center = TRUE), k))
    add("power_exact", conjugate(leaf(k), power_transform(0.5, exact = TRUE), k))
    add("yeojohnson_exact", conjugate(leaf(k), yeo_johnson(0.5, exact = TRUE), k))
    add("log_lower", conjugate(leaf(k), log_transform(lower = -15.0), k))
    add("log_upper", conjugate(leaf(k), log_transform(upper = 15.0), k))
    add("gaussianize", conjugate(leaf(k), gaussianize(), k))
    add("gaussianize_smooth", conjugate(leaf(k), gaussianize(smooth = TRUE), k))
    add("gaussianize_delete", conjugate(leaf(k), gaussianize(max_knots = 15L, guard = "delete", slope_cap = 1.5), k))
    add("gaussianize_emin", conjugate(leaf(k), gaussianize(max_knots = 15L, guard = "emin", slope_cap = 1.5), k))
    add("ema_skater", ema(0.05, k))
    add("pw_ensemble", precision_weighted_ensemble(list(ema(0.05, k), ema(0.2, k)), k))
    add("multiscale", multiscale(function(kk) conjugate(leaf(kk), ema_transform(0.1), kk), k))
    add(
      "bayes_ensemble",
      bayesian_ensemble(
        list(ema(0.05, k), conjugate(leaf(k), difference(), k)),
        k = k,
        learning_rate = 0.5,
        complexity_penalty = 0.02,
        depths = c(1, 1)
      )
    )
  }
  scenarios[["scale_mixture_leaf"]] <- list(k = 1L, sk = scale_mixture_leaf(1L))
  scenarios[["crps_leaf"]] <- list(k = 1L, sk = crps_leaf(1L))
  scenarios[["garch_leaf"]] <- list(k = 1L, sk = garch_leaf(1L))
  scenarios[["scalemix_ema"]] <- list(
    k = 1L,
    sk = conjugate(scale_mixture_leaf(1L), ema_transform(0.1), 1L)
  )
  scenarios[["gpd_tails"]] <- list(
    k = 1L,
    sk = gpdtails(conjugate(leaf(1L), ema_transform(0.1), 1L), k = 1L, level = 0.9, nexc = 50L, warmup = 100L)
  )
  scenarios[["search_default"]] <- list(
    k = 1L,
    sk = adaptive_search(k = 1L, expand_interval = 50L)
  )
  scenarios[["spec_diff_ensemble"]] <- list(
    k = 1L,
    sk = spec_build(
      conjugate_spec(ensemble_spec(ema_spec(0.01, 1L), ema_spec(0.1, 1L), k = 1L), diff_spec())
    )
  )
  scenarios[["spec_ema"]] <- list(k = 1L, sk = spec_build(ema_spec(0.05, 1L)))
  scenarios[["pol_laplace"]] <- list(k = 1L, sk = laplace(k = 1L))
  scenarios[["pol_laplace_k3"]] <- list(k = 3L, sk = laplace(k = 3L))
  scenarios[["pol_flesaker"]] <- list(k = 1L, sk = flesaker(k = 1L, lower = -15.0))

  fails <- 0L
  checked <- 0L
  check_block <- function(scenarios, series, expected_block) {
    for (name in names(scenarios)) {
      sc <- scenarios[[name]]
      expected <- expected_block[[name]]$out
      st <- NULL
      row <- 0L
      # Per-scenario tally. Without it the summary line below printed "ok" for
      # EVERY scenario, failing or not, so a drifted port looked clean apart from
      # whichever scenario happened to consume the global 7-line print budget.
      sc_fails <- 0L
      for (i in seq_along(series)) {
        r <- sc$sk(series[i], st)
        st <- r$state
        if (i - 1 >= v$burn) {
          row <- row + 1L
          for (h in seq_len(sc$k)) {
            got <- probe(r$dists[[h]], v$probe, v$q_lo, v$q_hi)
            exp_ <- unlist(expected[[row]][[h]])
            exp_ <- suppressWarnings(as.numeric(exp_))
            for (j in seq_along(got)) {
              checked <<- checked + 1L
              if (is.na(exp_[j])) {
                next
              }
              if (abs(got[j] - exp_[j]) > ATOL + RTOL * abs(exp_[j])) {
                fails <<- fails + 1L
                sc_fails <- sc_fails + 1L
                # Budget per SCENARIO, not globally, so every drifted scenario
                # shows evidence instead of only the first one.
                if (sc_fails <= 3) {
                  cat(sprintf("FAIL %s row %d h %d probe %d: got %.9g want %.9g\n", name, row, h, j, got[j], exp_[j]))
                }
              }
            }
          }
        }
      }
      result$scenarios[[name]] <<- list(fails = sc_fails)
      if (sc_fails > 0) {
        cat(sprintf("FAIL %-16s %d mismatches\n", name, sc_fails))
      } else {
        cat(sprintf("ok   %-16s\n", name))
      }
    }
  }
  check_block(scenarios, series, v$scenarios)

  # Sticky/dirac on the repeat-heavy series (exact repeats + 0.25-grid jumps).
  repeat_series <- unlist(v$repeat_series)
  repeat_scenarios <- list(
    sticky_ema = list(k = 1L, sk = sticky(conjugate(leaf(1L), ema_transform(0.1), 1L), k = 1L))
  )
  check_block(repeat_scenarios, repeat_series, v$repeat_scenarios)

  # Periodicity detector: ranked (lag, acf) per step on the main series.
  {
    pd <- period_detector()
    st <- NULL
    row <- 0L
    for (i in seq_along(series)) {
      r <- pd(series[i], st)
      st <- r$state
      if (i - 1 >= v$burn) {
        row <- row + 1L
        expected <- v$periodicity[[row]]
        checked <- checked + 1L
        if (length(r$scores) != length(expected)) {
          fails <- fails + 1L
          if (fails < 8) {
            cat(sprintf("FAIL periodicity row %d: %d scores, want %d\n", row, length(r$scores), length(expected)))
          }
          next
        }
        for (j in seq_along(expected)) {
          lag_want <- as.integer(expected[[j]][[1]])
          acf_want <- suppressWarnings(as.numeric(expected[[j]][[2]]))
          got <- r$scores[[j]]
          checked <- checked + 2L
          if (as.integer(got[1]) != lag_want) {
            fails <- fails + 1L
            if (fails < 8) {
              cat(sprintf("FAIL periodicity row %d rank %d: lag %d want %d\n", row, j, as.integer(got[1]), lag_want))
            }
          }
          if (
            !is.na(acf_want) &&
              abs(got[2] - acf_want) > ATOL + RTOL * abs(acf_want)
          ) {
            fails <- fails + 1L
            if (fails < 8) cat(sprintf("FAIL periodicity row %d rank %d: acf %.9g want %.9g\n", row, j, got[2], acf_want))
          }
        }
      }
    }
    cat("ok   periodicity\n")
  }

  # Covariance estimators on the fixed multivariate series.
  vec_series <- lapply(v$vec_series, unlist)
  cov_fns <- list(running = running_cov, ema = ema_cov, ledoit = ledoit_wolf_cov)
  for (nm in names(cov_fns)) {
    fn <- cov_fns[[nm]]
    expected <- v$cov[[nm]]
    st <- NULL
    row <- 0L
    for (i in seq_along(vec_series)) {
      r <- fn(vec_series[[i]], st)
      st <- r$state
      if (i - 1 >= v$burn) {
        row <- row + 1L
        got <- c(r$mean, r$cov)
        exp_ <- suppressWarnings(as.numeric(unlist(expected[[row]])))
        for (j in seq_along(got)) {
          checked <- checked + 1L
          if (is.na(exp_[j])) {
            next
          }
          if (abs(got[j] - exp_[j]) > ATOL + RTOL * abs(exp_[j])) {
            fails <- fails + 1L
            if (fails < 8) cat(sprintf("FAIL cov %s row %d probe %d: got %.9g want %.9g\n", nm, row, j, got[j], exp_[j]))
          }
        }
      }
    }
    cat(sprintf("ok   cov_%s\n", nm))
  }

  # --- structure: the candidate population, not its numbers ---
  # Numeric probes cannot see a missing candidate. This port sat at 57 candidates
  # against the reference's 60 (no seasonal_anchor) while every individual
  # transform matched to 1e-6; only the composed laplace drifted. Count and depth
  # histogram catch that in one line, and catch a candidate added in the wrong
  # position too, since the ensemble aligns depths by index.
  if (!is.null(v$structure)) {
    for (kname in names(v$structure)) {
      kk <- as.integer(kname)
      want <- v$structure[[kname]]
      cd <- build_candidates(kk)
      got_n <- length(cd$candidates)
      want_n <- as.integer(want$n_candidates)
      if (got_n != want_n) {
        fails <- fails + 1L
        result$structure_ok <- FALSE
        cat(sprintf("FAIL structure k=%d: %d candidates, want %d\n", kk, got_n, want_n))
      }
      got_d <- as.integer(cd$depths)
      want_d <- as.integer(unlist(want$depths))
      if (length(got_d) != length(want_d) || any(got_d != want_d)) {
        fails <- fails + 1L
        result$structure_ok <- FALSE
        cat(sprintf("FAIL structure k=%d: depth vector differs (got %s want %s)\n", kk,
                    paste(table(got_d), collapse = "/"), paste(table(want_d), collapse = "/")))
      }
      if (got_n == want_n && length(got_d) == length(want_d) && all(got_d == want_d)) {
        cat(sprintf("ok   structure k=%d   %d candidates\n", kk, got_n))
      }
    }
  }
  result$values_checked <- checked
  result$fails <- fails
  cat(sprintf("%d values checked\n", checked))
  if (fails > 0L) cat(sprintf("PARITY FAILED: %d\n", fails)) else cat("PARITY OK\n")
  result
}

# Transform parity: inst/parity/transform_vectors.json, generated from the
# reference by tools/gen_transform_vectors.py. seasonal_scale, gaussianize,
# log_transform, flesaker, the exact flags, a pushforward read through the
# next transform's components view, and every component-reading transform
# over a tail-spliced predictive. Same probes and 1e-6 tolerance as run_parity.
transform_parity_scenarios <- function() {
  spliced <- function(k) {
    gpdtails(conjugate(leaf(k), ema_transform(0.1), k), k = k, level = 0.9, nexc = 50L, warmup = 100L)
  }
  s <- list()
  for (k in c(1L, 3L)) {
    suf <- if (k == 1L) "" else sprintf("_k%d", k)
    add <- function(name, series, sk) s[[paste0(name, suf)]] <<- list(series = series, k = k, sk = sk)
    add("seasonal_scale", "main", conjugate(leaf(k), seasonal_scale(7L), k))
    add("seasonal_scale_center", "main", conjugate(leaf(k), seasonal_scale(7L, center = TRUE), k))
    add("seasonal_scale_sdiff", "main", conjugate(conjugate(leaf(k), seasonal_difference(7L), k), seasonal_scale(7L), k))
    add("gaussianize", "main", conjugate(leaf(k), gaussianize(), k))
    add("gaussianize_smooth", "main", conjugate(leaf(k), gaussianize(smooth = TRUE), k))
    add("gaussianize_delete", "main", conjugate(leaf(k), gaussianize(guard = "delete", slope_cap = 2.0), k))
    add("gaussianize_emin", "main", conjugate(leaf(k), gaussianize(guard = "emin", slope_cap = 2.0), k))
    add("gaussianize_move_tight", "main", conjugate(leaf(k), gaussianize(slope_cap = 2.0, max_knots = 15L), k))
    add("power_exact", "main", conjugate(leaf(k), power_transform(0.5, exact = TRUE), k))
    add("yeojohnson_log_exact", "main", conjugate(leaf(k), yeo_johnson(0.0, exact = TRUE), k))
    add("yeojohnson_half_exact", "main", conjugate(leaf(k), yeo_johnson(0.5, exact = TRUE), k))
    add("power_over_yj_exact", "main", conjugate(conjugate(leaf(k), yeo_johnson(0.5, exact = TRUE), k), power_transform(0.5), k))
    add("yj_over_power_exact", "main", conjugate(conjugate(leaf(k), power_transform(0.5, exact = TRUE), k), yeo_johnson(0.0), k))
    add("log_lower", "positive", conjugate(leaf(k), log_transform(), k))
    add("log_lower_shifted", "main", conjugate(leaf(k), log_transform(lower = -6.0), k))
    add("log_upper", "main", conjugate(leaf(k), log_transform(upper = 10.0), k))
    add("spliced_power", "main", conjugate(spliced(k), power_transform(0.5), k))
    add("spliced_yj_log", "main", conjugate(spliced(k), yeo_johnson(0.0), k))
    add("spliced_yj_half", "main", conjugate(spliced(k), yeo_johnson(0.5), k))
    add("spliced_gaussianize", "main", conjugate(spliced(k), gaussianize(), k))
    add("spliced_ema", "main", conjugate(spliced(k), ema_transform(0.05), k))
    add("spliced_seasonal_difference", "main", conjugate(spliced(k), seasonal_difference(2L), k))
    add("spliced_seasonal_anchor", "main", conjugate(spliced(k), seasonal_anchor(2L), k))
    add("spliced_seasonal_scale", "main", conjugate(spliced(k), seasonal_scale(7L), k))
    add("spliced_log", "positive", conjugate(spliced(k), log_transform(), k))
    add("spliced_power_over_yj", "main", conjugate(conjugate(spliced(k), yeo_johnson(0.5), k), power_transform(0.5), k))
  }
  s[["flesaker"]] <- list(series = "positive", k = 1L, sk = flesaker(k = 1L))
  s
}

run_transform_parity <- function(vec_path, verbose = TRUE, only = NULL,
                                 scenarios = transform_parity_scenarios(), label = "TRANSFORM PARITY") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read the parity vectors")
  }
  cat <- if (verbose) base::cat else function(...) invisible(NULL)
  v <- jsonlite::fromJSON(vec_path, simplifyVector = FALSE)
  ATOL <- 1e-6
  RTOL <- 1e-6
  all_names <- names(scenarios)
  if (!is.null(only)) scenarios <- scenarios[only]
  result <- list(scenarios = list(), values_checked = 0L, fails = 0L, missing = character(0))
  for (name in names(scenarios)) {
    sc <- scenarios[[name]]
    want <- v$scenarios[[name]]
    if (is.null(want)) {
      result$missing <- c(result$missing, name)
      result$fails <- result$fails + 1L
      cat(sprintf("FAIL %-28s not in the vectors\n", name))
      next
    }
    series <- unlist(v$series[[want$series]])
    # A per-scenario rtol marks a float-limited reference (see the generator).
    rtol <- if (is.null(want$rtol)) RTOL else as.numeric(want$rtol)
    st <- NULL
    row <- 0L
    sc_fails <- 0L
    for (i in seq_along(series)) {
      r <- sc$sk(series[i], st)
      st <- r$state
      if (i - 1 >= v$burn) {
        row <- row + 1L
        for (h in seq_len(sc$k)) {
          d <- r$dists[[h]]
          got <- c(
            dist_mean(d), dist_std(d), dist_logpdf(d, v$probe), dist_cdf(d, v$probe),
            dist_quantile(d, v$q_lo), dist_quantile(d, v$q_hi), dist_crps(d, v$probe)
          )
          exp_ <- unlist(want$out[[row]][[h]])
          exp_ <- vapply(exp_, function(e) switch(as.character(e), inf = Inf, `-inf` = -Inf, nan = NaN, as.numeric(e)), 0.0)
          for (j in seq_along(got)) {
            result$values_checked <- result$values_checked + 1L
            ok <- if (is.nan(exp_[j])) {
              is.nan(got[j])
            } else if (is.infinite(exp_[j])) {
              identical(got[j], exp_[j])
            } else {
              is.finite(got[j]) && abs(got[j] - exp_[j]) <= ATOL + rtol * abs(exp_[j])
            }
            if (!ok) {
              sc_fails <- sc_fails + 1L
              if (sc_fails <= 3) {
                cat(sprintf("FAIL %s row %d h %d probe %d: got %.12g want %.12g\n", name, row, h, j, got[j], exp_[j]))
              }
            }
          }
        }
      }
    }
    if (row != length(want$out)) {
      sc_fails <- sc_fails + 1L
      cat(sprintf("FAIL %s: %d rows, want %d\n", name, row, length(want$out)))
    }
    result$scenarios[[name]] <- list(fails = sc_fails)
    result$fails <- result$fails + sc_fails
    cat(sprintf("%s %-28s %s\n", if (sc_fails) "FAIL" else "ok  ", name, if (sc_fails) paste(sc_fails, "mismatches") else ""))
  }
  extra <- setdiff(names(v$scenarios), all_names)
  if (length(extra) && is.null(only)) {
    result$fails <- result$fails + length(extra)
    result$missing <- c(result$missing, extra)
    cat(sprintf("FAIL vectors carry scenarios this runner lacks: %s\n", paste(extra, collapse = ", ")))
  }
  cat(sprintf("%d values checked\n", result$values_checked))
  cat(if (result$fails > 0L) sprintf("%s FAILED: %d\n", label, result$fails) else sprintf("%s OK\n", label))
  result
}

# Wave 2 parity: inst/parity/wave2_vectors.json from tools/gen_wave2_vectors.py.
# homogenize, residual_transform and its diagnostics, dantzig and search,
# pushforward affine maps, and the spec JSON and names (byte for byte).
wave2_parity_scenarios <- function() {
  spliced <- function(k) {
    gpdtails(conjugate(leaf(k), ema_transform(0.1), k), k = k, level = 0.9, nexc = 50L, warmup = 100L)
  }
  ema_leaf <- function(k) conjugate(leaf(k), ema_transform(0.1), k)
  s <- list()
  add <- function(name, series, k, sk) s[[name]] <<- list(series = series, k = k, sk = sk)
  add("homogenize_leaf", "main", 1L, homogenize(leaf(1L)))
  add("homogenize_ema", "main", 1L, homogenize(ema_leaf(1L)))
  for (couple in c("none", "weights", "state", "shared")) {
    add(paste0("homogenize_ema_k3_", couple), "main", 3L, homogenize(ema_leaf(3L), k = 3L, couple = couple, lam = 0.4))
  }
  add("homogenize_custom", "main", 1L, homogenize(ema_leaf(1L),
    candidates = homogenize_candidates(c(0.0, 0.9), 0.05, c(0.0, 0.3, 0.9), delta_max = 0.6)))
  add("homogenize_spliced", "main", 1L, homogenize(spliced(1L)))
  for (model in c("m0", "m1", "m2")) {
    add(paste0("residual_", model), "main", 1L, residual_transform(ema_leaf(1L), model = model, lr = 0.05))
  }
  add("residual_m2_spliced", "main", 1L, residual_transform(spliced(1L), model = "m2"))
  add("dantzig", "main", 1L, dantzig(k = 1L, expand_interval = 40L))
  add("search", "main", 1L, adaptive_search(k = 1L, expand_interval = 60L, max_pool = 10L))
  add("standardize_over_yj_exact", "main", 1L,
      conjugate(conjugate(leaf(1L), yeo_johnson(0.0, exact = TRUE)), standardize(0.05)))
  add("seasonal_scale_over_log", "positive", 3L,
      conjugate(conjugate(leaf(3L), log_transform(lower = -1.0), 3L), seasonal_scale(7L), 3L))
  add("gaussianize_over_yj_exact", "main", 1L,
      conjugate(conjugate(leaf(1L), yeo_johnson(0.5, exact = TRUE)), gaussianize()))
  s
}

run_wave2_parity <- function(vec_path, verbose = TRUE) {
  cat <- if (verbose) base::cat else function(...) invisible(NULL)
  res <- run_transform_parity(vec_path, verbose = verbose, scenarios = wave2_parity_scenarios(),
                              label = "WAVE2 SCENARIO PARITY")
  v <- jsonlite::fromJSON(vec_path, simplifyVector = FALSE)
  ATOL <- 1e-6
  RTOL <- 1e-6
  num <- function(e) switch(as.character(e), inf = Inf, `-inf` = -Inf, nan = NaN, as.numeric(e))
  close <- function(got, want) {
    if (is.nan(want)) return(is.nan(got))
    if (is.infinite(want)) return(identical(got, want))
    is.finite(got) && abs(got - want) <= ATOL + RTOL * abs(want)
  }
  probe <- function(d) {
    c(dist_mean(d), dist_std(d), dist_logpdf(d, v$probe), dist_cdf(d, v$probe),
      dist_quantile(d, v$q_lo), dist_quantile(d, v$q_hi), dist_crps(d, v$probe))
  }
  tally <- function(section, ok, msg) {
    res$values_checked <<- res$values_checked + 1L
    if (!ok) {
      res$fails <<- res$fails + 1L
      res$sections[[section]] <<- (if (is.null(res$sections[[section]])) 0L else res$sections[[section]]) + 1L
      cat(sprintf("FAIL %s %s\n", section, msg))
    }
  }
  res$sections <- list(affine = 0L, diagnostics = 0L, specs = 0L)

  # pushforward affine, shift and scale
  inner <- dist_new(c(0.6, 0.4), c(0.2, -0.5), c(0.7, 1.3))
  bases <- list(
    table = yeo_johnson(0.0, exact = TRUE)$inverse_k(list(inner), list())[[1]],
    analytic = log_transform()$inverse_k(list(inner), list())[[1]]
  )
  for (cs in v$affine) {
    d <- bases[[cs$kind]]
    e <- switch(cs$op,
      affine = dist_affine(d, cs$a, cs$b),
      shift = dist_shift(d, 1.5),
      scale = dist_scale(d, -2.0),
      affine2 = dist_affine(dist_affine(d, -0.5, 2.0), -3.0, 1.0)
    )
    got <- probe(e)
    want <- vapply(cs$probe, num, 0.0)
    for (j in seq_along(got)) {
      tally("affine", close(got[j], want[j]),
            sprintf("%s %s(%g, %g) probe %d: got %.12g want %.12g", cs$kind, cs$op, cs$a, cs$b, j, got[j], want[j]))
    }
  }
  cat(sprintf("%s affine (%d cases)\n", if (res$sections$affine) "FAIL" else "ok  ", length(v$affine)))

  # residual diagnostics, one summary row per update
  dg <- v$diagnostics
  rd <- residual_diagnostics(dg$K, dg$halflife)
  for (i in seq_along(dg$z)) {
    rd <- residual_diagnostics_update(rd, as.numeric(dg$z[[i]]))
    sm <- residual_diagnostics_summary(rd)
    got <- c(sm$n, sm$`E[z]`, sm$`E[z2]`, sm$`E[z3]`, sm$`E[z4]`, sm$L, sm$cross_sq)
    want <- vapply(dg$rows[[i]], num, 0.0)
    tally("diagnostics", length(got) == length(want), sprintf("row %d length", i))
    for (j in seq_along(want)) {
      tally("diagnostics", close(got[j], want[j]), sprintf("row %d col %d: got %.12g want %.12g", i, j, got[j], want[j]))
    }
  }
  cat(sprintf("%s residual_diagnostics (%d rows)\n", if (res$sections$diagnostics) "FAIL" else "ok  ", length(dg$z)))

  # specs: JSON byte for byte through a round trip, canonical names, build
  for (sp in v$specs) {
    parsed <- from_json(sp$json)
    tally("specs", identical(to_json(parsed), sp$json), sprintf("to_json(from_json(%s)) = %s", sp$json, to_json(parsed)))
    tally("specs", identical(spec_name(parsed), sp$name), sprintf("spec_name %s: got %s want %s", sp$json, spec_name(parsed), sp$name))
    f <- tryCatch(build(parsed), error = function(e) NULL)
    tally("specs", identical(!is.null(f), isTRUE(sp$builds)), sprintf("build %s: builds %s", sp$json, !is.null(f)))
    if (!is.null(f)) {
      tally("specs", identical(attr(f, "spec_name"), sp$name), sprintf("build name %s", sp$json))
    }
  }
  cat(sprintf("%s specs (%d)\n", if (res$sections$specs) "FAIL" else "ok  ", length(v$specs)))
  cat(if (res$fails > 0L) sprintf("WAVE2 PARITY FAILED: %d\n", res$fails) else "WAVE2 PARITY OK\n")
  res
}

# Cross-language checkpoints: inst/parity/state_vectors.json.gz from
# tools/gen_state_vectors.py. R restores each Python checkpoint, steps 20
# ticks and must match the reference's outputs; its own state then must
# match the reference's, field for field.
state_parity_scenarios <- function() {
  cj <- function(t, k = 1L) conjugate(leaf(k), t, k)
  s <- list()
  add <- function(name, k, sk) s[[name]] <<- list(k = k, sk = sk)
  tr <- list(
    difference = difference(), fractional_difference = fractional_difference(0.4, 30L),
    standardize = standardize(), ema_transform = ema_transform(0.1), ou_transform = ou_transform(0.1),
    theta = theta(0.1), drift = drift(0.05, 0.01), holt_linear = holt_linear(0.1, 0.05), garch = garch(),
    seasonal_difference = seasonal_difference(7L), seasonal_anchor = seasonal_anchor(7L),
    seasonal_scale = seasonal_scale(7L), seasonal_scale_center = seasonal_scale(7L, center = TRUE),
    power_transform = power_transform(0.5), power_transform_exact = power_transform(0.5, exact = TRUE),
    yeo_johnson = yeo_johnson(0.5), yeo_johnson_exact = yeo_johnson(0.5, exact = TRUE),
    ar = ar(2L), grouped_ar = grouped_ar(8L), gaussianize = gaussianize(),
    gaussianize_smooth = gaussianize(smooth = TRUE)
  )
  for (nm in names(tr)) add(paste0("conjugate_", nm), 1L, cj(tr[[nm]]))
  add("conjugate_log_transform", 1L, cj(log_transform()))
  add("laplace", 1L, laplace(1L))
  add("laplace_k3", 3L, laplace(3L))
  add("laplace_likelihood", 1L, laplace(1L, objective = "likelihood"))
  add("laplace_gaussian_tails", 1L, laplace(1L, tails = "gaussian"))
  add("laplace_repeat", 1L, laplace(1L))
  add("laplace_gaps", 1L, laplace(1L))
  add("laplace_k3_gaps", 3L, laplace(3L))
  add("granger", 1L, granger())
  add("flesaker", 1L, flesaker(1L, lower = -100.0))
  add("dantzig", 1L, dantzig(1L, expand_interval = 50L))
  add("search_k2", 2L, adaptive_search(k = 2L, expand_interval = 40L, max_pool = 12L))
  add("multiscale_k4", 4L, multiscale(function(kk) cj(ema_transform(0.1), kk), 4L))
  add("parade_spliced", 1L, parade(gpdtails(cj(ema_transform(0.1)), k = 1L, nexc = 50L, warmup = 100L), 1L))
  add("parade_laplace_k3", 3L, parade(laplace(3L), 3L))
  add("sticky_ema_repeat", 1L, sticky(cj(ema_transform(0.1)), k = 1L))
  add("homogenize", 1L, homogenize(cj(ema_transform(0.1)), k = 1L))
  add("homogenize_k2_weights", 2L, homogenize(cj(ema_transform(0.1), 2L), k = 2L, couple = "weights"))
  add("homogenize_k2_shared", 2L, homogenize(cj(ema_transform(0.1), 2L), k = 2L, couple = "shared"))
  for (model in c("m0", "m1", "m2")) {
    add(paste0("residual_transform_", model), 1L, residual_transform(cj(ema_transform(0.1)), model = model))
  }
  add("pw_ensemble_k2", 2L, precision_weighted_ensemble(list(ema(0.05, 2L), ema(0.2, 2L)), k = 2L))
  add("bayes_ensemble_k2", 2L, bayesian_ensemble(list(ema(0.05, 2L), cj(difference(), 2L)), k = 2L))
  add("terminal_adaptive", 1L, terminal_leaf_ensemble(list(ema(0.05, 1L), cj(difference())), k = 1L,
                                                      depths = c(0, 1), adaptive_temperature = TRUE))
  add("terminal_k3", 3L, terminal_leaf_ensemble(list(ema(0.05, 3L), cj(difference(), 3L)), k = 3L,
                                                depths = c(0, 1), forget = 0.99))
  add("leaf", 1L, leaf(1L))
  add("scale_mixture_leaf", 1L, scale_mixture_leaf(1L))
  add("crps_leaf", 1L, crps_leaf(1L))
  add("garch_leaf", 1L, garch_leaf(1L))
  add("ema_k2", 2L, ema(0.1, 2L))
  s
}

# Restore every reference checkpoint, step `continue` ticks, compare the
# outputs with the reference's at 1e-6, then compare R's own state with the
# reference's state at that point, field for field. If out_dir is given, R's
# state there is written as <scenario>_<t>.json for the reference to resume
# (tools/check_r_state_in_python.py). With roundtrip, R's state is also put
# through JSON and the twin stepped beside it for another `continue` ticks,
# bit for bit.
run_state_parity <- function(vec_path, verbose = TRUE, only = NULL, out_dir = NULL, roundtrip = TRUE) {
  cat <- if (verbose) base::cat else function(...) invisible(NULL)
  con <- gzfile(vec_path)
  v <- jsonlite::fromJSON(con, simplifyVector = FALSE)
  ATOL <- 1e-6
  RTOL <- 1e-6
  num <- function(e) switch(as.character(e), inf = Inf, `-inf` = -Inf, nan = NaN, as.numeric(e))
  close <- function(got, want) {
    if (is.nan(want)) return(is.nan(got))
    if (is.infinite(want)) return(identical(got, want))
    is.finite(got) && abs(got - want) <= ATOL + RTOL * abs(want)
  }
  probe <- function(d) {
    c(dist_mean(d), dist_std(d), dist_logpdf(d, v$probe), dist_cdf(d, v$probe),
      dist_quantile(d, v$q_lo), dist_quantile(d, v$q_hi), dist_crps(d, v$probe))
  }
  state_diffs <- function(a, b, path, acc) {
    if (length(acc) > 20) return(acc)
    kind <- function(x) if (is.null(x)) "null" else if (is.list(x)) (if (is.null(names(x))) "array" else "object") else "scalar"
    ta <- kind(a)
    tb <- kind(b)
    if (length(a) == 0 && length(b) == 0 && ta %in% c("array", "object") && tb %in% c("array", "object")) return(acc)
    if (ta != tb) return(c(acc, sprintf("%s: %s vs %s", path, ta, tb)))
    if (ta == "object") {
      if (!setequal(names(a), names(b))) return(c(acc, sprintf("%s: keys differ", path)))
      for (key in names(b)) acc <- state_diffs(a[[key]], b[[key]], paste0(path, ".", key), acc)
    } else if (ta == "array") {
      if (length(a) != length(b)) return(c(acc, sprintf("%s: length %d vs %d", path, length(a), length(b))))
      for (i in seq_along(a)) acc <- state_diffs(a[[i]], b[[i]], sprintf("%s[%d]", path, i - 1L), acc)
    } else if (ta == "scalar") {
      if (is.numeric(a) && is.numeric(b)) {
        if (!close(as.numeric(a), as.numeric(b))) acc <- c(acc, sprintf("%s: %.12g vs %.12g", path, a, b))
      } else if (!identical(a, b)) {
        acc <- c(acc, sprintf("%s: %s vs %s", path, format(a), format(b)))
      }
    }
    acc
  }
  scenarios <- state_parity_scenarios()
  if (!is.null(only)) scenarios <- scenarios[only]
  res <- list(scenarios = list(), values_checked = 0L, fails = 0L,
              missing = setdiff(names(v$scenarios), names(state_parity_scenarios())))
  for (name in names(scenarios)) {
    want <- v$scenarios[[name]]
    if (is.null(want)) {
      res$fails <- res$fails + 1L
      cat(sprintf("FAIL %-28s not in the vectors\n", name))
      next
    }
    series <- vapply(v$series[[want$series]], num, 0.0)
    sk <- scenarios[[name]]$sk
    sc_fails <- 0L
    deep <- 0L
    for (cp in want$checkpoints) {
      st <- state_from_json(cp$state)
      for (m in seq_len(v$`continue`)) {
        r <- sk(series[cp$t + m], st)
        st <- r$state
        for (h in seq_along(r$dists)) {
          got <- probe(r$dists[[h]])
          exp_ <- vapply(cp$out[[m]][[h]], num, 0.0)
          tail_ <- .deep_tail_probe(r$dists[[h]], v$probe)
          for (j in seq_along(got)) {
            res$values_checked <- res$values_checked + 1L
            if (tail_ && j %in% c(3L, 4L) && !close(got[j], exp_[j]) && is.finite(got[j]) &&
                abs(got[j] - exp_[j]) <= ATOL + 1e-4 * abs(exp_[j])) {
              deep <- deep + 1L
              next
            }
            if (!close(got[j], exp_[j])) {
              sc_fails <- sc_fails + 1L
              if (sc_fails <= 3) {
                cat(sprintf("FAIL %s t=%d +%d h %d probe %d: got %.12g want %.12g\n", name, cp$t, m, h, j, got[j], exp_[j]))
              }
            }
          }
        }
      }
      mine <- state_to_json(st)
      d <- state_diffs(mine, cp$state_next, "", character(0))
      if (length(d)) {
        sc_fails <- sc_fails + length(d)
        cat(sprintf("FAIL %s t=%d state after %d steps: %s\n", name, cp$t, v$`continue`, paste(d[seq_len(min(3L, length(d)))], collapse = "; ")))
      }
      text <- to_json(mine, digits = 17L)
      if (!is.null(out_dir)) {
        writeLines(text, file.path(out_dir, sprintf("%s_%d.json", name, cp$t)))
      }
      if (roundtrip) {
        # Within R a JSON checkpoint resumes bit for bit: step the live state
        # and its JSON twin side by side; outputs identical, states equal as text.
        twin <- state_from_json(jsonlite::parse_json(text, simplifyVector = FALSE))
        if (!identical(to_json(state_to_json(twin), digits = 17L), text)) {
          sc_fails <- sc_fails + 1L
          cat(sprintf("FAIL %s t=%d: re-encoded restore differs\n", name, cp$t))
        }
        for (m in seq_len(v$`continue`)) {
          y <- series[cp$t + v$`continue` + m]
          ra <- sk(y, st)
          rb <- sk(y, twin)
          st <- ra$state
          twin <- rb$state
          pa <- unlist(lapply(ra$dists, probe))
          pb <- unlist(lapply(rb$dists, probe))
          if (!identical(pa, pb)) {
            sc_fails <- sc_fails + 1L
            cat(sprintf("FAIL %s t=%d +%d: R JSON restore is not bit-exact\n", name, cp$t, v$`continue` + m))
            break
          }
        }
        if (!identical(to_json(state_to_json(st), digits = 17L), to_json(state_to_json(twin), digits = 17L))) {
          sc_fails <- sc_fails + 1L
          cat(sprintf("FAIL %s t=%d: R JSON restore state drifted\n", name, cp$t))
        }
      }
    }
    res$scenarios[[name]] <- list(fails = sc_fails, deep_tail = deep)
    res$fails <- res$fails + sc_fails
    cat(sprintf("%s %-28s %s%s\n", if (sc_fails) "FAIL" else "ok  ", name,
                if (sc_fails) paste(sc_fails, "mismatches") else "",
                if (deep) sprintf(" (%d deep-tail logpdf/cdf probes at 1e-4)", deep) else ""))
  }
  if (length(res$missing) && is.null(only)) {
    res$fails <- res$fails + length(res$missing)
    cat(sprintf("FAIL vectors carry scenarios this runner lacks: %s\n", paste(res$missing, collapse = ", ")))
  }
  cat(sprintf("%d values checked\n", res$values_checked))
  cat(if (res$fails > 0L) sprintf("STATE PARITY FAILED: %d\n", res$fails) else "STATE PARITY OK\n")
  res
}

# The probe point lies where a spliced predictive's body has less than 1e-8
# mass beyond it. There the splice's z comes from the body CDF within ~1e-12
# of 0 or 1, which both ports resolve only to the last bit of that CDF (the
# reference's 0.5 * (1 + erf), R's pnorm sum), so logpdf and cdf agree only
# to about 1e-5 relative. The deep-tail quantile was fixed in the reference
# (0ad69a5); logpdf and cdf there still read the body CDF.
.deep_tail_probe <- function(d, x) {
  if (!isTRUE(d$spliced)) {
    return(FALSE)
  }
  cm <- dist_components(d$body)
  if (is.null(cm)) {
    return(FALSE)
  }
  lower <- sum(cm$w * pnorm((x - cm$m) / cm$s))
  upper <- sum(cm$w * pnorm((x - cm$m) / cm$s, lower.tail = FALSE))
  min(lower, upper) < 1e-8
}
