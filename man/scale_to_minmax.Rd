% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scale_to_minmax.R
\name{scale_to_minmax}
\alias{scale_to_minmax}
\title{Scale a vector to the range of a matrix, element-wise}
\usage{
scale_to_minmax(v, mat, useInds = TRUE)
}
\arguments{
\item{v}{A numeric vector to be scaled}

\item{mat}{A numeric matrix used to define the scaling range}

\item{useInds}{A logical vector or integer vector indicating the columns of `mat` or the elements of `v`
to use in the calculation. Default is to use all columns or elements.}
}
\value{
A numeric vector with the same length as `v`, scaled to the range of `mat` element-wise.
}
\description{
Scale a vector to the range of a matrix, element-wise
}
