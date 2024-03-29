% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keep_inds_in_bounds.R
\name{keep_inds_in_bounds}
\alias{keep_inds_in_bounds}
\title{keep_inds_in_bounds.R}
\usage{
keep_inds_in_bounds(check, against)
}
\arguments{
\item{check}{A numeric vector of indices to check.}

\item{against}{A numeric vector of indices to check against.}
}
\value{
A numeric vector of indices within the bounds of against.
}
\description{
Given a set of indices to check and a set of indices to check against, returns
the subset of check that lies within the bounds of against.
}
\examples{
keep_inds_in_bounds(1:10, 5:15)

}
