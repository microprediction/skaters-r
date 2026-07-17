# Contract tests: determinism and checkpoint-resume equivalence.
#
# Determinism: two independently constructed skaters fed the same series
# must emit bit-identical probe values at every step.
#
# Checkpoint-resume: the state is pure data (closures live in the skater,
# never in the state), so saving the state mid-stream with saveRDS and
# resuming with a freshly built skater must reproduce the uninterrupted
# run exactly. RDS serializes doubles losslessly, so the comparison is
# identical(), not a tolerance.
#
# Run: Rscript tests/contract.R   (also runs under R CMD check)
if (nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  library(skaters)
} else {
  for (.f in sort(list.files("R", full.names = TRUE))) {
    source(.f)
  }
}

failures <- 0L
check <- function(cond, label) {
  if (!isTRUE(cond)) {
    failures <<- failures + 1L
    cat(sprintf("FAIL %s\n", label))
  }
}

lcg <- function(seed) {
  s <- seed %% 2^32
  function() {
    s <<- (1664525 * s + 1013904223) %% 2^32
    s / 2^32
  }
}
gauss <- function(rand) {
  u <- 0.0
  v <- 0.0
  while (u == 0) {
    u <- rand()
  }
  while (v == 0) {
    v <- rand()
  }
  sqrt(-2.0 * log(u)) * cos(2.0 * pi * v)
}

make_series <- function(n, seed) {
  rand <- lcg(seed)
  ys <- numeric(n)
  lvl <- 0.0
  for (t in seq_len(n)) {
    lvl <- lvl + 0.02 + 0.3 * gauss(rand)
    ys[t] <- lvl + sin(2 * pi * t / 12) + gauss(rand)
  }
  ys
}

probe <- function(dists) {
  unlist(lapply(dists, function(d) {
    c(
      dist_mean(d),
      dist_std(d),
      dist_logpdf(d, 0.3),
      dist_cdf(d, 0.3),
      dist_quantile(d, 0.1),
      dist_quantile(d, 0.9),
      dist_crps(d, 0.3)
    )
  }))
}

run_probed <- function(f, ys, st = NULL, from = 1L) {
  out <- vector("list", length(ys) - from + 1L)
  for (i in from:length(ys)) {
    r <- f(ys[i], st)
    st <- r$state
    out[[i - from + 1L]] <- probe(r$dists)
  }
  list(probes = out, state = st)
}

# --- determinism -----------------------------------------------------------
# Two separately constructed instances, one series, exact agreement at
# every step (not just the end).
for (k in c(1L, 3L)) {
  ys <- make_series(800L, 7L)
  a <- run_probed(laplace(k = k), ys)
  b <- run_probed(laplace(k = k), ys)
  check(identical(a$probes, b$probes), sprintf("determinism k=%d: probe streams differ", k))
  check(identical(a$state, b$state), sprintf("determinism k=%d: final states differ", k))
  cat(sprintf("ok   determinism k=%d (800 ticks, %d probes/tick)\n", k, length(a$probes[[1]])))
}

# --- checkpoint-resume -----------------------------------------------------
# Run to the midpoint, saveRDS the state, readRDS it back, and continue
# with a FRESHLY BUILT skater. The resumed tail must equal the
# uninterrupted run's tail exactly. 1200 ticks puts the midpoint (600)
# beyond the GPD tail warm-up (500), so the splice and its excess buffers
# are part of the round-tripped state.
for (k in c(1L, 3L)) {
  ys <- make_series(1200L, 11L)
  full <- run_probed(laplace(k = k), ys)

  first <- run_probed(laplace(k = k), ys[1:600])
  tf <- tempfile(fileext = ".rds")
  saveRDS(first$state, tf)
  restored <- readRDS(tf)
  unlink(tf)
  check(identical(restored, first$state), sprintf("resume k=%d: state changed by RDS round-trip", k))

  resumed <- run_probed(laplace(k = k), ys, st = restored, from = 601L)
  check(
    identical(resumed$probes, full$probes[601:1200]),
    sprintf("resume k=%d: resumed probes differ from uninterrupted run", k)
  )
  check(identical(resumed$state, full$state), sprintf("resume k=%d: final states differ", k))
  cat(sprintf("ok   checkpoint-resume k=%d (save at 600/1200)\n", k))
}

# Adaptive search: its state holds only recipes and plain data (the
# Python reference keeps live callables in the state; the R port rebuilds
# them from recipes), so it must satisfy the same resume contract, across
# an expansion boundary (expand at 100 and 200 with the save at 150).
{
  ys <- make_series(300L, 13L)
  full <- run_probed(adaptive_search(k = 1L), ys)
  first <- run_probed(adaptive_search(k = 1L), ys[1:150])
  tf <- tempfile(fileext = ".rds")
  saveRDS(first$state, tf)
  restored <- readRDS(tf)
  unlink(tf)
  resumed <- run_probed(adaptive_search(k = 1L), ys, st = restored, from = 151L)
  check(identical(resumed$probes, full$probes[151:300]), "search resume: probes differ from uninterrupted run")
  check(identical(resumed$state, full$state), "search resume: final states differ")
  cat("ok   checkpoint-resume adaptive_search (save at 150/300)\n")
}

# The splice must actually be active in the resumed run at k=1, otherwise
# the checkpoint test is not exercising the tail state.
{
  ys <- make_series(1200L, 11L)
  f <- laplace(1L)
  st <- NULL
  d <- NULL
  for (y in ys) {
    r <- f(y, st)
    st <- r$state
    d <- r$dists[[1]]
  }
  check(isTRUE(d$spliced), "resume: splice not active at 1200 ticks")
  cat("ok   splice active across checkpoint\n")
}

if (failures > 0) {
  cat(sprintf("CONTRACT FAILED: %d violation(s)\n", failures))
  quit(status = 1)
}
cat("CONTRACT OK (determinism, checkpoint-resume)\n")
