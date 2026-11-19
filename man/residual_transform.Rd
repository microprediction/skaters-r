\name{residual_transform}
\alias{residual_transform}
\alias{corrected_dist}
\alias{residual_diagnostics}
\alias{residual_diagnostics_update}
\alias{residual_diagnostics_summary}
\title{Online residual-transform correction}
\description{
  \code{residual_transform(base, model, lr, k = 1)} reshapes a k = 1
  skater's predictive by a residual model learned online from its own
  PIT residuals: \code{"m0"} identity, \code{"m1"} dynamic scale,
  \code{"m2"} scale plus leverage. Forecasts are \code{corrected_dist}
  objects (Python's \code{CorrectedDist}) answering the \code{dist_*}
  probes. \code{residual_diagnostics(K, halflife)} is Python's
  \code{ResidualDiagnostics} as a value: \code{residual_diagnostics_update}
  returns it updated with one z, and \code{residual_diagnostics_summary}
  gives the moments and the leverage statistics. Port of
  \code{skaters/residual_transform.py}.
}
