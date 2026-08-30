\name{quantile.skater}
\alias{quantile.skater}
\title{Quantile paths from a skater}
\usage{
\method{quantile}{skater}(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), ...)
}
\arguments{
\item{x}{a \code{"skater"} object that has seen at least one observation.}
\item{probs}{probabilities.}
\item{...}{unused.}
}
\value{a matrix with one row per horizon and one column per probability.}
\description{
Exact quantiles of each horizon's predictive mixture.
}
