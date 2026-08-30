# Locate the parity vectors and the engine in both harness contexts:
# R CMD check (installed package) and a repo-root run.
if (!nzchar(system.file("parity", "vectors.json", package = "skaters"))) {
  for (.f in sort(list.files(file.path("..", "..", "R"), full.names = TRUE))) {
    source(.f)
  }
}
if (!exists("run_parity")) run_parity <- skaters:::run_parity
skaters_vectors_path <- function() {
  p <- system.file("parity", "vectors.json", package = "skaters")
  if (nzchar(p)) return(p)
  for (cand in c(file.path("..", "..", "inst", "parity", "vectors.json"),
                 file.path("inst", "parity", "vectors.json"))) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  ""
}
