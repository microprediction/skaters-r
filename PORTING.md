# Porting map: Python skaters -> R

The discipline, inherited from the JavaScript twin: no module ships
without passing the parity vectors (`inst/parity/vectors.json`, generated
by the Python reference) at 1e-6. `Rscript tests/parity.R` is the gate.
Two further gates run alongside it (standalone and under R CMD check):
`tests/robustness.R`, the adversarial streams from the reference repo's
release gate, and `tests/contract.R`, determinism plus exact
checkpoint-resume through saveRDS/readRDS.

Every row below is ported. The table stays as the map of what lives where.

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
| terminal leaf ensemble | terminal.py | DONE |
| sticky (lattice) | sticky.py | DONE |
| multiscale | multiscale.py | DONE |
| parade (pit/z state) | parade.py | DONE |
| tails (GPD splice, 0.13.0 default) | tails.py | DONE |
| laplace (the composition) | api.py | DONE |
| adaptive search (`adaptive_search`; `search` would mask base) | search.py | DONE |
| spec build path | spec.py | DONE |
| periodicity detector | periodicity.py | DONE |
| covariance estimators (running, ema, ledoit_wolf) | cov/ | DONE |

A porting note for R: factories must `force()` their arguments (and any
wrapper must `force(base)`). R promises otherwise capture loop variables
by reference, so candidates built in loops silently share the last
parameter value, and `f <- wrapper(f)` recurses into itself.

Suggested order: transforms and leaves are independent and parallelize
well; ensembles need `prune` (port its ulp-tolerant pair merge exactly,
it exists to keep platforms in agreement); the composition comes last and
is then covered by the `pol_laplace` scenarios in the same vectors file.

Parity vectors provenance: inst/parity/vectors.json is copied from
microprediction/skaters at commit 38ddfa9, the merge of
microprediction/skaters#365 (parity/vectors.json there); refresh by
rerunning parity/gen_vectors.py in that repository and copying the
output here.

Quantile contract: inst/parity/quantile_contract.json is copied unchanged
from the same commit (parity/quantile_contract.json there, described in
parity/QUANTILE_CONTRACT.md). It holds 64 inverse-CDF cases shared with
the Python, JavaScript and Rust runners, and
tests/testthat/test-quantile-contract.R runs them against `dist_quantile`
with the fixture's `atol`. `dist_quantile` and `spliced_quantile` default
to `max_iter = 2100` as the reference does. The granger mixture quantile
keeps its own bisection and budget.

Spliced maps: `dist_affine`, `dist_shift` and `dist_scale` dispatch on the
spliced predictive to `spliced_affine`, `spliced_shift` and
`spliced_scale`, ports of `SplicedDist.affine` (microprediction/skaters#358)
and its shift and scale companions. The splice lives in the body's PIT
space, so only the body moves; a negative multiplier swaps the tails and
negates the thresholds. Without them, conjugating laplace with standardize,
garch, ema_transform, ou_transform, seasonal_difference or seasonal_anchor
failed once the
GPD splice went live. tests/testthat/test-spliced-affine.R mirrors
tests/test_spliced_affine.py.

Spliced through component-reading transforms: the k>1 widening in
`ema_transform`, `seasonal_difference` and `seasonal_anchor` goes through
`.widen`, a port of `_widen` in transform.py (microprediction/skaters
commit bc67bd8). On a spliced predictive it widens the body and rebuilds
the splice with the same thresholds, tail masses and GPD parameters. The
nonlinear inverses of `power_transform` and `yeo_johnson` emit the exact
pushforward of a spliced predictive, a port of the table path of
`PushforwardDist` in pushforward.py (R/pushforward.R: `table_from_map`
over the extreme quantiles at 1e-6 padded by their span, slope from the
exact inverse, table quadrature CRPS, moments from the 257-point quantile
rule). The `dist_*` queries dispatch on `pushforward = TRUE`. Mixture
inputs take the unchanged component-wise paths. `gaussianize` and the
`exact = TRUE` options of power_transform and yeo_johnson are not yet
ported. tests/testthat/test-spliced-nonaffine.R mirrors
tests/test_spliced_nonaffine.py.

Release channels: R-universe from day one (automatic builds per push once
the repo is registered), CRAN manually when the surface stabilizes. CRAN
has no automated submission path and expects at most a release every month or
two; r-universe is the continuous channel meanwhile.
