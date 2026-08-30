\name{sandwich}
\alias{rosenblatt}
\alias{sandwich_transport}
\alias{sandwich_new}
\alias{sandwich_cdf}
\alias{sandwich_logpdf}
\alias{sandwich_quantile}
\alias{sandwich_qgrid}
\alias{sandwich_mean}
\alias{sandwich_var}
\alias{sandwich_crps}
\title{The laplace sandwich: run any model in z-coordinates}
\description{
  \code{rosenblatt} runs a skater over a series and returns the causal
  Rosenblatt stream \code{z_t = Phi^-1(F_t(y_t))}, where \code{F_t} is the
  predictive the skater issued for \code{y_t}. Under calibration that stream
  is close to i.i.d. standard normal, so an outside model can be fitted to
  \code{z} rather than \code{y} and inherits a stationarised, unit-scale
  series. Whatever it finds there is, by construction, structure the skater
  did not already capture.

  \code{sandwich_new} is the way back. Given the inner model's Gaussian in
  z-space and the transport the skater issued for that step, it assembles the
  y-space predictive with exact change-of-variables accounting,
  \code{log f_Y(y) = log f_Z(z) + log f_t(y) - log phi(z)}. The
  \code{sandwich_*} functions probe the result; \code{dist_cdf} and friends
  dispatch to them when \code{d$sandwiched} is set, ahead of the spliced
  check, so a sandwich may wrap a tail-spliced body.

  Nothing here knows what the inner model is. It takes a mean and a standard
  deviation in z-space, so fitting \code{ets}, \code{auto.arima}, or anything
  else is caller-side glue and costs this package no dependency.
}
\details{
  \code{rosenblatt} reports \code{z[1]} as \code{NA}, since no predictive had
  been issued when the first observation arrived. The Python mirror
  (\code{prophet-laplace}) substitutes 0 there because Prophet requires a
  gap-free series; R models accept \code{NA}, so this port reports the absence
  instead.

  \code{sandwich_transport(r, step)} returns the transport for a given step
  ahead. Steps past \code{k} reuse the \code{k}-step transport, which is an
  approximation rather than an exact multi-step map.

  With \code{z_mu = 0} and \code{z_sigma = 1} the sandwich is the identity on
  its body: an inner model that says nothing recovers the skater alone.
}
\examples{
set.seed(1)
y <- cumsum(rnorm(200)) * 0.1
r <- rosenblatt(y, k = 2L)

# fit anything to r$z, then hand back its z-space mean and sd
d <- sandwich_new(sandwich_transport(r, 1L), z_mu = 0.2, z_sigma = 0.9)
sandwich_quantile(d, 0.5)
dist_cdf(d, 0.0)
}
