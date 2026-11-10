\name{observe}
\alias{observe}
\alias{observe.skater}
\title{Feed observations to a skater}
\usage{
observe(object, y, ...)
\method{observe}{skater}(object, y, ...)
}
\arguments{
\item{object}{a \code{"skater"} object.}
\item{y}{numeric vector of observations, fed in order.}
\item{...}{unused.}
}
\value{the updated \code{"skater"} object; assign it back, as with any
functional update in R.}
\description{
Steps the wrapped closure once per element of \code{y}, in order. The latest
predictive distributions are then available through
\code{\link{predict.skater}}.
}
