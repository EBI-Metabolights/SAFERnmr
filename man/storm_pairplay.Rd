% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/storm_pairplay.R
\name{storm_pairplay}
\alias{storm_pairplay}
\title{storm_pairplay: Run modified STORM on the provided spectral region and ref shape.
Built for accepting corrPocketPairs results. Notes:}
\usage{
storm_pairplay(
  xmat = NULL,
  ppm = NULL,
  b = 30,
  corrthresh = 0.8,
  q = 0.05,
  minpeak = 10,
  refSpec = NULL,
  ref.idx = NULL,
  driver = NULL
)
}
\arguments{
\item{xmat}{A matrix of spectral data (rows are spectra, columns are spectral points)}

\item{ppm}{A vector of the spectral points in ppm (optional, default is all columns of xmat)}

\item{b}{An integer giving the expansion parameter for the reference peak}

\item{corrthresh}{A numeric giving the minimum correlation value to be considered for inclusion (for both subset AND reference optimization)}

\item{q}{A numeric giving the p-value threshold for correlation significance (both subset AND reference optimization)}

\item{minpeak}{An integer giving the minimum number of points allowed in a run of significant points in the reference}

\item{refSpec}{A vector of spectral data to use as the initial reference}

\item{ref.idx}{A vector of the spectral points (columns of xmat) to use as the initial reference}
}
\value{
A list with components "reconstructed" and "status". "reconstructed" is a matrix
containing the reconstructed metabolite concentrations (rows are samples, columns are metabolites).
"status" is a character string indicating whether the method converged successfully or failed.
}
\description{
STORM: Joram Posma's STORM has been adapted and optimized to accept these
}
