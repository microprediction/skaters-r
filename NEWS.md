# skaters 0.16.0

The version number tracks the Python reference implementation, whose parity
vectors this port is verified against (105,798 values across 55 scenarios at
a tolerance of 1e-6).

* Caught up to the reference: `standardize` and `garch` corrected,
  AR coefficients damped into the stationary region, `seasonal_anchor` added.
* testthat harness alongside the standalone script harnesses, porting the
  original assertions rather than rewriting them.
* S3 wrapper around the closure API: `skater()`, `observe()`, `predict()`,
  `quantile()`, and `print()` methods.
* NAMESPACE lists every export explicitly instead of pattern-matching.
* roxygen2 documentation with usage, arguments, and return values for every
  exported function.

# skaters 0.0.1

* Initial port: transforms, leaves, ensembles, spec grammar, tail splice,
  parade, multiscale, sticky, periodicity, covariance, adaptive search, and
  the `laplace` forecaster, with the parity, contract, and robustness
  harnesses.
