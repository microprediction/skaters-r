# Skater state as plain JSON in the reference's canonical layout
# (docs/STATE.md in microprediction/skaters). Mirrors tests/test_state_json.py
# there, across ports: R restores every Python checkpoint in
# inst/parity/state_vectors.json.gz, steps 20 ticks and must match Python's
# outputs and state to 1e-6, and a JSON round trip within R resumes bit for
# bit. With SKATERS_PY_SRC set (the reference's checkout), the R-written
# checkpoints are also resumed in Python.
.state_vectors_path <- function() {
  p <- system.file("parity", "state_vectors.json.gz", package = "skaters")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("..", "..", "inst", "parity", "state_vectors.json.gz"),
                 file.path("inst", "parity", "state_vectors.json.gz"))) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  ""
}

test_that("Python checkpoints resume in R, and R's resume in Python", {
  skip_if_not(nzchar(.state_vectors_path()), "state_vectors.json.gz not found")
  skip_if_not_installed("jsonlite")
  rsp <- if (exists("run_state_parity")) run_state_parity else skaters:::run_state_parity
  out_dir <- tempfile("rstate")
  dir.create(out_dir)
  res <- rsp(.state_vectors_path(), verbose = FALSE, out_dir = out_dir)
  for (nm in names(res$scenarios)) {
    expect_equal(res$scenarios[[nm]]$fails, 0L, info = paste("scenario", nm))
  }
  expect_equal(res$fails, 0L)
  expect_gt(res$values_checked, 20000)
  py_src <- Sys.getenv("SKATERS_PY_SRC")
  tool <- file.path(dirname(dirname(.state_vectors_path())), "..", "tools", "check_r_state_in_python.py")
  if (nzchar(py_src) && file.exists(tool) && nzchar(Sys.which("python3"))) {
    tool <- normalizePath(tool)
    root <- normalizePath(file.path(dirname(.state_vectors_path()), "..", ".."))
    owd <- setwd(root) # the tool reads inst/parity/state_vectors.json.gz
    on.exit(setwd(owd), add = TRUE)
    out <- system2("python3", c(shQuote(tool), shQuote(out_dir)), stdout = TRUE, stderr = TRUE,
                   env = sprintf("PYTHONPATH=%s:%s", file.path(py_src, "src"), file.path(py_src, "parity")))
    expect_true(any(grepl("R->PYTHON RESUME OK", out)), info = paste(tail(out, 5), collapse = "\n"))
  }
})

test_that("the walker names the offending path", {
  expect_error(state_to_json(list(a = list(b = c(1.0, Inf)))), "at \\.a\\.b\\[1\\]: non-finite")
  expect_error(state_to_json(list(f = sum)), "a function")
  pf <- yeo_johnson(0.0, exact = TRUE)$inverse_k(list(dist_gaussian(0, 1)), list())[[1]]
  expect_error(assert_plain_state(list(d = pf)), "a pushforward")
  expect_silent(assert_plain_state(list(a = 1, b = list(c = "x", d = NULL))))
})

test_that("Dist encoding is exact and unrenormalised", {
  set.seed(7)
  for (i in 1:200) {
    n <- sample(1:9, 1)
    d <- dist_new(runif(n) + 1e-3, rnorm(n, 0, 4), 0.1 + runif(n))
    back <- state_from_json(from_json(to_json(state_to_json(d), digits = 17L)))
    expect_identical(back, d)
  }
  loose <- state_from_json(list(components = list(list(1.0, 0.0, 1.0), list(3.0, 1.0, 1.0))))
  expect_identical(loose$w[1], 0.25)
  sp <- spliced_new(dist_new(c(0.4, 0.6), c(0, 1), c(1, 2)), -1.5, 1.5, 0.05, 0.04, 0.1, 0.8, 0.2, 0.9)
  back <- state_from_json(from_json(to_json(state_to_json(sp), digits = 17L)))
  expect_identical(back, sp)
})

test_that("from_json reads doubles exactly and to_json writes them exactly", {
  set.seed(3)
  x <- c(runif(300) * 10^sample(-300:300, 300, TRUE), 2^-1074, .Machine$double.xmax, 0.1, 1 / 3)
  for (xi in x) {
    expect_identical(from_json(to_json(xi, digits = 17L)), xi)
    expect_identical(from_json(to_json(xi)), xi)
  }
})

test_that("R checkpoints from before the canonical layout still resume", {
  f <- leaf(1L)
  st <- NULL
  for (y in c(1, 2, 0.5, 3)) st <- f(y, st)$state
  old <- st$var # the old R leaf state was the Welford triple itself
  expect_identical(f(1.5, old)$dists, f(1.5, st)$dists)
  pd <- period_detector()
  st <- NULL
  set.seed(1)
  ys <- rnorm(80)
  for (y in ys[1:60]) st <- pd(y, st)$state
  old <- st
  old$cross <- unname(old$cross)
  for (y in ys[61:80]) {
    a <- pd(y, st)
    b <- pd(y, old)
    expect_identical(a$scores, b$scores)
    st <- a$state
    old <- b$state
  }
})
