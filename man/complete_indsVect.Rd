% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complete_indsVect.R
\name{complete_indsVect}
\alias{complete_indsVect}
\title{Given a sparse (but evenly spaced) vector of indexes,
with NAs in the missing places,
expand it outwards to fill the entire vector. Basically 
linear inter/extrapolation of a perfectly linear vect.}
\usage{
complete_indsVect(indVect)
}
\arguments{
\item{indVect}{vector of indices with gaps (can be on the ends), in any order.}
}
\value{
vector with values completed
}
\description{
Given a sparse (but evenly spaced) vector of indexes,
with NAs in the missing places,
expand it outwards to fill the entire vector. Basically 
linear inter/extrapolation of a perfectly linear vect.
}
\author{
MTJ
}
