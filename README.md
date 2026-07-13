# skaters for R

Port of [skaters](https://github.com/microprediction/skaters): online,
univariate, distributional forecasting by conjugation. Pure R, no
dependencies (jsonlite only for the parity tests).

Status: the full surface is ported, including `laplace` (the deployed
default with GPD tails and the parade wrapper), the adaptive search, the
spec builders, the periodicity detector, and the covariance estimators.
Everything is parity-verified against the Python reference at 1e-6
(105,798 probe values), and two further gates run under R CMD check:
the adversarial robustness streams and a determinism plus exact
checkpoint-resume contract. See PORTING.md for the module map.

## Installation

You can install the development version of skaters from
[GitHub](https://github.com/microprediction/skaters-r) with:

```r
# install.packages("pak")
pak::pak("microprediction/skaters-r")
```

## Example

```r
library(skaters)

sk <- conjugate(leaf(1L), ema_transform(0.1), 1L)
state <- NULL
for (y in rnorm(100)) {
  r <- sk(y, state)
  state <- r$state
}
d <- r$dists[[1]]
c(dist_mean(d), dist_quantile(d, 0.9))
```
