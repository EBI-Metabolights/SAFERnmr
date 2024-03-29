% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_expand.R
\name{corr_expand}
\alias{corr_expand}
\title{Expand a peak region to the closest correlation minimums}
\usage{
corr_expand(peak, localMins, vRange)
}
\arguments{
\item{peak}{An integer representing the index of the peak in the correlation vector}

\item{localMins}{A numeric vector containing the indices of local minima in the correlation vector}

\item{vRange}{A numeric vector containing the range of valid indices for the correlation vector}
}
\value{
A list with two elements: \code{lower}, the index of the closest local minimum to the left of the peak, and \code{upper}, the index of the closest local minimum to the right of the peak

 @importFrom magrittr %>%
}
\description{
This function takes a correlation profile vector and an index of a peak, and returns the bounds of the region containing the peak, extended to the closest local minima on either side of the peak.
}
\examples{
v <- c(1, 2, 3, 2, 1)
localMins <- c(1, 3, 5)
corr_expand(peak = 3, localMins = localMins, vRange = c(1, 5))
# Returns: list(lower = 1, upper = 5)
}
