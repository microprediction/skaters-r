\name{tails}
\alias{gpdtails}
\alias{spliced_new}
\alias{spliced_cdf}
\alias{spliced_logpdf}
\alias{spliced_quantile}
\alias{spliced_qgrid}
\alias{spliced_mean}
\alias{spliced_var}
\alias{spliced_crps}
\title{GPD tail splice}
\description{
  \code{gpdtails} wraps a skater and, after a warm-up, splices
  censored-ML generalized-Pareto tails into the body's predictive
  beyond frozen z-thresholds. The \code{spliced_*} functions probe the
  resulting spliced distribution; \code{dist_cdf} and friends dispatch
  to them when \code{d$spliced} is set. Port of \code{skaters/tails.py}.
}
