\name{adaptive_search}
\alias{adaptive_search}
\alias{dantzig}
\title{Adaptive search over the transform tree}
\description{
  Beam search over the transform grammar: candidates are scored by
  cumulative clamped log-likelihood, top performers are expanded with
  new transforms (children replay recent history so they join warm),
  and losers are pruned. The reference calls it \code{search}; the R
  name is \code{adaptive_search} because \code{search} would mask
  \code{base::search}.
  \code{dantzig(k, ...)} is the policy name in the reference's
  \code{api.py}, a call to \code{adaptive_search}. Port of
  \code{skaters/search.py}.
}
