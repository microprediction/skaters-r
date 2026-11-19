\name{adaptive_search}
\alias{adaptive_search}
\alias{search}
\alias{dantzig}
\title{Adaptive search over the transform tree}
\description{
  Beam search over the transform grammar: candidates are scored by
  cumulative clamped log-likelihood, top performers are expanded with
  new transforms (children replay recent history so they join warm),
  and losers are pruned. \code{search} is the reference's name for
  the same function; it masks \code{base::search} while the package is
  attached, and \code{adaptive_search} is the unmasked name.
  \code{dantzig(k, ...)} is the policy name in the reference's
  \code{api.py}, a call to \code{adaptive_search}. Port of
  \code{skaters/search.py}.
}
