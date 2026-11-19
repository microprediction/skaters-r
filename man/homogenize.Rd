\name{homogenize}
\alias{homogenize}
\alias{homogenize_candidates}
\alias{homogenized_dist}
\title{Online homogenized scale-mixture correction}
\description{
  \code{homogenize(base, candidates, k, couple, lam)} wraps a skater
  with a Bayesian-hedged pool of two-point Gaussian scale-mixture
  hypotheses for its own PIT residual stream, one pool per horizon
  (\code{couple = "none"}) or coupled across horizons by
  \code{"weights"}, \code{"state"} or one \code{"shared"} pool.
  \code{homogenize_candidates} builds the candidate grid (Python's
  \code{make_candidates}). The forecasts are \code{homogenized_dist}
  objects, the base predictive reshaped by the pooled correction, which
  answer the \code{dist_*} probes. Port of \code{skaters/homogenize.py}.
}
