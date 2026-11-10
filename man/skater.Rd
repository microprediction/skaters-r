\name{skater}
\alias{skater}
\alias{print.skater}
\title{Wrap a skater closure in an S3 object}
\usage{
skater(fn, name = NULL)
\method{print}{skater}(x, ...)
}
\arguments{
\item{fn}{a skater closure, such as the value of \code{\link{laplace}}.}
\item{name}{optional label used by \code{print()}.}
\item{x}{a \code{"skater"} object.}
\item{...}{unused.}
}
\value{an object of class \code{"skater"}.}
\description{
A skater in this package is a closure \code{f(y, state)} returning
\code{list(dists, state)}. That closure is the parity surface, verified
against the reference implementation's vectors to 1e-6, and it never changes
shape. \code{skater()} wraps any such closure in a small S3 object so the
usual R generics work: \code{\link{observe}} to feed data,
\code{\link{predict.skater}} for interval forecasts,
\code{\link{quantile.skater}} for quantile paths, and \code{print()}.

The wrapper is deliberately the only integration seam. A faster backend
producing the same closure contract plugs in here with no change to callers.
}
\examples{
sk <- skater(laplace(k = 3))
sk <- observe(sk, c(0.4, -0.2, 0.1, 0.6, -0.3))
predict(sk)
quantile(sk, probs = c(0.1, 0.5, 0.9))
}
