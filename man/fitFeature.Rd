% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitFeature.R
\name{fitFeature}
\alias{fitFeature}
\title{Fit a feature to a spectrum}
\usage{
fitFeature(feature, spectrum, spectrum.position, method = "minmax scaling")
}
\arguments{
\item{feature}{A numeric vector containing the feature to fit to the spectrum}

\item{spectrum}{A numeric vector containing the spectrum to fit the feature to}

\item{spectrum.position}{A numeric value specifying the position of the spectrum}

\item{method}{A character string specifying the method to use for fitting. Can be "least squares" or "minmax scaling".}
}
\value{
A list containing the feature fit, the feature position, the ratio, the residuals, and the overfit

 @export
}
\description{
This function fits a feature to a spectrum using either least squares fitting or minmax scaling.
}
