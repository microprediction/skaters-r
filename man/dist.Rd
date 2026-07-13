\name{dist}
\alias{dist_new}
\alias{dist_gaussian}
\alias{dist_combine}
\alias{dist_mean}
\alias{dist_var}
\alias{dist_std}
\alias{dist_logpdf}
\alias{dist_cdf}
\alias{dist_crps}
\alias{dist_quantile}
\alias{dist_shift}
\alias{dist_scale}
\alias{dist_affine}
\alias{dist_prune}
\title{Gaussian mixture distribution: the distributional prediction type}
\description{
  Constructors and probes for the Gaussian mixture that every skater
  emits: \code{dist_new(w, m, s)} builds one, \code{dist_gaussian} a
  single component, \code{dist_combine} a weighted mixture of mixtures.
  Probes: \code{dist_mean}, \code{dist_var}, \code{dist_std},
  \code{dist_logpdf}, \code{dist_cdf}, \code{dist_crps},
  \code{dist_quantile}. Affine maps: \code{dist_shift},
  \code{dist_scale}, \code{dist_affine}. \code{dist_prune} reduces the
  component count by ulp-tolerant closest-pair merging. Port of the
  Python reference \code{skaters/dist.py}, parity-checked at 1e-6.
}
