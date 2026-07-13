# skaters for R

Port of [skaters](https://github.com/microprediction/skaters): online,
univariate, distributional forecasting by conjugation. Pure R, no
dependencies (jsonlite only for the parity tests).

Status: core types and first transforms, parity-verified against the
Python reference at 1e-6 (11,760 probe values and counting). See
PORTING.md for the module map and how to help.

```r
source("R/dist.R"); source("R/leaf.R")
source("R/transform.R"); source("R/conjugate.R")

sk <- conjugate(leaf(1L), ema_transform(0.1), 1L)
state <- NULL
for (y in rnorm(100)) {
  r <- sk(y, state); state <- r$state
}
d <- r$dists[[1]]
c(dist_mean(d), dist_quantile(d, 0.9))
```
