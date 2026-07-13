\name{spec}
\alias{spec_build}
\alias{spec_name}
\alias{leaf_spec}
\alias{ema_spec}
\alias{ensemble_spec}
\alias{conjugate_spec}
\alias{diff_spec}
\alias{frac_spec}
\alias{std_spec}
\alias{ema_t_spec}
\title{Symbolic specification of skater pipelines}
\description{
  A spec is a plain list describing how to build a skater; it can be
  serialized, compared, and materialized with \code{spec_build}.
  \code{spec_name} derives the canonical name (transform chains read
  left to right, e.g. \code{"diff|ema_t(0.1)|leaf"}). The \code{*_spec}
  functions are constructors for the grammar. Port of
  \code{skaters/spec.py}.
}
