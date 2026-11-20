\name{state}
\alias{state_to_json}
\alias{state_from_json}
\alias{assert_plain_state}
\title{Skater state as plain JSON in the canonical layout}
\description{
  A skater is a fold, \code{(y, state) -> (predictives, state)}, and its
  state is plain data in the layout the reference documents in
  \code{docs/STATE.md} (microprediction/skaters), so a checkpoint
  written by the Python, JS or R port resumes in any of them.
  \code{state_to_json(state)} returns that layout as plain R data (a
  named list is a JSON object, an unnamed list an array, a length-one
  atomic a scalar, \code{NULL} null); \code{to_json(x, digits = 17)}
  writes it as text. \code{state_from_json(obj)} takes the layout as
  \code{from_json} or \code{jsonlite::parse_json(simplifyVector = FALSE)}
  reads it and returns a state the R skater resumes from.
  \code{assert_plain_state} raises an error, naming the path, on
  anything that would not survive JSON. Within R a JSON checkpoint
  resumes bit for bit; across ports the parity suite holds the resumed
  stream to 1e-6.
}
\examples{
f <- leaf(1L)
st <- NULL
for (y in c(1.0, 2.0, 0.5)) st <- f(y, st)$state
text <- to_json(state_to_json(st), digits = 17L)
st2 <- state_from_json(from_json(text))
identical(f(1.5, st)$dists, f(1.5, st2)$dists)
}
