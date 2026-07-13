\name{periodicity}
\alias{period_detector}
\alias{top_periods}
\title{Online periodicity detection}
\description{
  \code{period_detector} maintains exponentially weighted
  autocorrelation estimates at candidate lags, O(n_lags) per
  observation, returning ranked \code{(lag, acf)} scores;
  \code{top_periods} extracts the significant ones. Port of
  \code{skaters/periodicity.py}.
}
