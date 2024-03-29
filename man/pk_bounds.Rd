% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_bounds.R
\name{pk_bounds}
\alias{pk_bounds}
\title{Get bounds for local maxima in a vector}
\usage{
pk_bounds(v)
}
\arguments{
\item{v}{A numeric vector}
}
\value{
A list of numeric vectors, where each element corresponds to the bounds of a region containing a local maximum in \code{v}
}
\description{
This function takes a numeric vector \code{v} and returns the bounds of the regions containing the local maxima in \code{v}. Updated to include only true peaks (those with a lower value on either side; no endpoints).
}
\examples{
v <- c(1, 2, 3, 2, 1, 2, 2, 1, 2)
pk.bounds(v)
# Returns: list(c(1, 3), c(4, 5))
}
