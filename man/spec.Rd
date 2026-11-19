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
\alias{build}
\alias{to_json}
\alias{from_json}
\title{Symbolic specification of skater pipelines}
\description{
  A spec is a plain list describing how to build a skater; it can be
  serialized, compared, and materialized with \code{spec_build}.
  \code{spec_name} derives the canonical name (transform chains read
  left to right, e.g. \code{"diff|ema_t(0.1)|leaf"}). The \code{*_spec}
  functions are constructors for the grammar. \code{build} is
  Python's name for \code{spec_build}. \code{to_json} writes a spec
  byte for byte as Python's \code{json.dumps(spec, separators=(",", ":"))}
  does: an R integer is a JSON int and a double a float in Python's
  repr form (\code{0.05}, \code{1.0}, \code{1e-05}), a named list an
  object, an unnamed list an array. \code{from_json} reads JSON back,
  ints as integers and floats as doubles, so a round trip is exact.
  Port of \code{skaters/spec.py}.
}
