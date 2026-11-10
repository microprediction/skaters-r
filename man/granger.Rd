\name{granger}
\alias{granger}
\alias{granger_params}
\title{granger: laplace's trained twin mixed with laplace}
\description{
  Two forecasters run on every observation: the reference
  \code{\link{laplace}}, and laplace's own pipeline with trained values for
  its continuous hyperparameters plus a small gate that shifts the ensemble's
  combining weights each step. The issued predictive is their
  Aggregating-Algorithm mixture under log loss with prior \code{pi_lap} on
  laplace, so its total log score is never more than \code{log(1 / pi_lap)}
  (log 4 at the shipped 0.25) behind laplace's on any series, up to
  floating-point rounding. Horizon 1 only. After Clive Granger, whose work
  with Bates (1969) founded forecast combination. Port of
  \code{skaters/granger.py}.

  \code{granger_params(frequency)} returns a shipped parameter set
  (\code{"monthly"}, \code{"weekly"} or \code{"hourly"}), stored as RDS in
  \code{inst/extdata/granger}.
}
\usage{
granger(params = granger_params("monthly"), k = 1L)
granger_params(frequency = c("monthly", "weekly", "hourly"))
}
\arguments{
  \item{params}{a parameter list from \code{granger_params()}.}
  \item{k}{forecast horizon; only 1 is trained.}
  \item{frequency}{which shipped parameter set.}
}
\value{
  \code{granger} returns a skater \code{f(y, state)} returning
  \code{list(dists, state)}; \code{dists[[1]]} is the issued mixture, which
  the \code{dist_*} functions accept. \code{granger_params} returns a list.
}
\examples{
f <- granger(granger_params("monthly"))
s <- NULL
for (y in c(0.1, -0.2, 0.05, 0.3)) {
  r <- f(y, s)
  s <- r$state
}
dist_mean(r$dists[[1]])
}
