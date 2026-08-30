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
