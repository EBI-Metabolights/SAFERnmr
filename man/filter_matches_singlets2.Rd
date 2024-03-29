% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_matches_singlets2.R
\name{filter_matches_singlets2}
\alias{filter_matches_singlets2}
\title{Remove singlet matches}
\usage{
filter_matches_singlets2(match.info, feature.c, refmat.c, ncores = 1)
}
\arguments{
\item{pars}{A list of input parameters.}
}
\value{
A list of filtered and processed matched peak information, including back-fits to the original spectra.
}
\description{
Opting for a much lighter, much faster peak characterization approach. 
Characterize peak maxima locations first for features and refs. Use 
compressed matrix storage for peak models.
}
\details{
For a spectral signature (i.e. feature, or ref signature):
Determine matched points: those which were both 
- not NA in the signature
- contained within the matched region

How many of the true peak maxima in the signature exist within the matched points?

In other words, how many resonances do each of the fits contain?
}
