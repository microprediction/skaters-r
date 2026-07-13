\name{laplace}
\alias{laplace}
\alias{build_candidates}
\title{laplace: the general forecaster}
\description{
  The deployed default: a likelihood-weighted trunk over the candidate
  population from \code{build_candidates}, a CRPS terminal leaf, the
  lattice projection, multi-scale at \code{k > 1}, the GPD tail splice,
  and the parade wrapper. Call as \code{f <- laplace(k)} then
  \code{r <- f(y, state)}; \code{r$dists} holds one mixture per
  horizon. Port of \code{skaters/api.py}.
}
