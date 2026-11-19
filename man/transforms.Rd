\name{transforms}
\alias{difference}
\alias{ema_transform}
\alias{standardize}
\alias{ou_transform}
\alias{theta}
\alias{drift}
\alias{holt_linear}
\alias{garch}
\alias{seasonal_difference}
\alias{seasonal_anchor}
\alias{seasonal_scale}
\alias{power_transform}
\alias{yeo_johnson}
\alias{log_transform}
\alias{gaussianize}
\alias{fractional_difference}
\alias{ar}
\alias{grouped_ar}
\title{Invertible online transforms}
\description{
  The transform grammar: each factory returns a forward/inverse pair
  that \code{\link{conjugate}} wraps around an inner skater, so the
  inner skater models the transformed residual and predictions are
  mapped back. Port of \code{skaters/transform.py}, parity-checked
  at 1e-6. Note \code{ar} masks \code{stats::ar} and \code{garch} here
  is the variance-tracking transform, not a fitter.

  \code{seasonal_difference} subtracts the value one period ago;
  \code{seasonal_anchor} subtracts a hedge between that seasonal-naive
  value and a phase-indexed EMA, which adapts as fast as the naive
  without inheriting its single-draw noise (\code{weight = 0} recovers
  \code{seasonal_difference}). \code{seasonal_scale} standardizes by
  per-phase EWMA statistics shrunk toward the global ones.

  \code{power_transform} and \code{yeo_johnson} map each component by
  the delta method; \code{exact = TRUE} emits the exact pushforward
  instead. \code{log_transform} is the log coordinate for a series with
  a one-sided bound (\code{lower}, default 0, or \code{upper}) and
  always emits the exact pushforward. \code{gaussianize} maps through
  the probit of a smoothed fence-post empirical CDF. A pushforward
  predictive answers the \code{dist_*} probes by the change of
  variables, and \code{dist_components} gives its moment-matched
  mixture view.
}
