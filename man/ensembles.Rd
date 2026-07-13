\name{ensembles}
\alias{precision_weighted_ensemble}
\alias{bayesian_ensemble}
\alias{terminal_leaf_ensemble}
\title{Skater ensembles}
\description{
  \code{precision_weighted_ensemble} weights members by inverse
  mean-squared error; \code{bayesian_ensemble} by exponentiated
  log-likelihood with a complexity penalty; \code{terminal_leaf_ensemble}
  adds a terminal leaf fitted to the ensemble's own residuals (the
  trunk of \code{\link{laplace}}). Ports of \code{skaters/ensemble.py},
  \code{bayesian.py}, and \code{terminal.py}.
}
