% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compression.R
\name{expand_runs}
\alias{expand_runs}
\title{Expand ranges to a vector}
\usage{
expand_runs(ranges)
}
\arguments{
\item{ranges}{position ranges (for row-wise linear inds)}
}
\value{
expanded position vector (row-wise linear inds)

Example: all.equal(expand_runs(runs2ranges(keep)), keep)
}
\description{
Expand ranges to a vector
}
