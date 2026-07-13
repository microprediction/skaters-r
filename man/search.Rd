\name{adaptive_search}
\alias{adaptive_search}
\title{Adaptive search over the transform tree}
\description{
  Beam search over the transform grammar: candidates are scored by
  cumulative clamped log-likelihood, top performers are expanded with
  new transforms (children replay recent history so they join warm),
  and losers are pruned. Named \code{adaptive_search} because
  \code{search} would mask \code{base::search}. Port of
  \code{skaters/search.py}.
}
