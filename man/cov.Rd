\name{cov-estimators}
\alias{running_cov}
\alias{ema_cov}
\alias{ledoit_wolf_cov}
\title{Online covariance estimation}
\description{
  Streaming estimators processing one observation vector at a time:
  \code{running_cov} (Welford), \code{ema_cov} (exponentially
  weighted), and \code{ledoit_wolf_cov} (EMA correlations shrunk
  toward identity). Each returns \code{list(mean, cov, state)} with
  \code{cov} flat row-major. Port of \code{skaters/cov/}.
}
