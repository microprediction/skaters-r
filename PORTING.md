# Porting map: Python skaters -> R

The discipline, inherited from the JavaScript twin: no module ships
without passing the parity vectors (`inst/parity/vectors.json`, generated
by the Python reference) at 1e-6. `Rscript tests/parity.R` is the gate.

| module | Python source | status |
|---|---|---|
| dist (mixture type, incl. prune) | src/skaters/dist.py | DONE |
| leaf | src/skaters/leaf.py::leaf | DONE |
| difference, ema_transform | src/skaters/transform.py | DONE |
| conjugate | src/skaters/conjugate.py | DONE |
| remaining transforms (ar, theta, garch, seasonal, yeo-johnson, ou, frac, ...) | transform.py | DONE |
| scale_mixture_leaf, crps_leaf, garch_leaf | leaf.py | DONE |
| ensembles (precision, bayesian) | ensemble.py, bayesian.py | DONE |
| ema skater | ema.py | DONE |
| terminal leaf ensemble | terminal.py | open |
| sticky (lattice) | sticky.py | DONE |
| multiscale | multiscale.py | DONE |
| parade (pit/z state) | parade.py | open |
| tails (GPD splice, 0.13.0 default) | tails.py | open |
| laplace (the composition) | api.py | open (last) |

Suggested order: transforms and leaves are independent and parallelize
well; ensembles need `prune` (port its ulp-tolerant pair merge exactly,
it exists to keep platforms in agreement); the composition comes last and
is then covered by the `pol_laplace` scenarios in the same vectors file.

Parity vectors provenance: inst/parity/vectors.json is copied from
microprediction/skaters at commit 9cb5d1c (parity/vectors.json there);
refresh by rerunning parity/gen_vectors.py in that repository and
copying the output here.

Release channels: R-universe from day one (automatic builds per push once
the repo is registered), CRAN manually when the surface stabilizes. CRAN
has no automated submission path and expects at most a release every month or
two; r-universe is the continuous channel meanwhile.
