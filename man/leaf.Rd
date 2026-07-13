\name{leaf}
\alias{leaf}
\alias{scale_mixture_leaf}
\alias{crps_leaf}
\alias{garch_leaf}
\title{Leaf skaters: residual noise models}
\description{
  The terminal noise models of the transform tree. \code{leaf} is the
  running-moments Gaussian; \code{scale_mixture_leaf} learns a mixture
  over a scale basis by likelihood; \code{crps_leaf} learns it by CRPS;
  \code{garch_leaf} refits a GARCH(1,1) window. Each returns a skater
  \code{function(y, state)} yielding \code{list(dists, state)}. Port of
  \code{skaters/leaf.py}.
}
