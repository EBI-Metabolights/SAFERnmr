% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localMinima.R
\name{localMinima}
\alias{localMinima}
\title{Find local minima in a vector}
\usage{
localMinima(v)
}
\arguments{
\item{v}{a numeric vector}
}
\value{
a numeric vector containing the indices of local minima in \code{v}.
}
\description{
Given a numeric vector, returns the indices of local minima in the vector. Endpoints
can be included. Specifically, a point is considered a local minimum if there
are no lesser points adjacent.
}
