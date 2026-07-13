\name{conjugate}
\alias{conjugate}
\title{Conjugate a skater by an invertible transform}
\description{
  \code{conjugate(skater, transform, k)} runs the transform forward on
  each observation, lets the inner skater forecast the transformed
  series, and pushes the predictive distributions back through the
  inverse. The composition operator of the whole package. Port of
  \code{skaters/conjugate.py}.
}
