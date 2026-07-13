\name{parade}
\alias{parade}
\title{Prediction parade: PIT/z calibration state}
\description{
  Pass-through wrapper that scores each arriving observation against
  the predictive issued for it, exposing per-horizon PIT values and
  normal z-scores in \code{state$pit} and \code{state$z}. Also gates
  pathological inputs before the tree consumes them. Port of
  \code{skaters/parade.py}.
}
