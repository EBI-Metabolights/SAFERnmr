% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateMatchInfoRow.R
\name{updateMatchInfoRow}
\alias{updateMatchInfoRow}
\title{In match propagation, take a copy of a match.info row and update the relevant 
fields. Return the updated row. Accessory fcn, called by filter.matches.R.}
\usage{
updateMatchInfoRow(mi.row, fit)
}
\arguments{
\item{mi.row}{a row of match.info data frame}

\item{fit}{pre-computed fit for intracluster lag-adjusted feature (the other cluster member) to the ref feat}
}
\value{
row of match.info, updated for other cluster members.
}
\description{
Within a feature cluster, only one representative (or the mean profile) was matched
to all the refs. Using the matches to that, and the lags of each cluster member
to the key feature/mean profile, update the fit data for each cluster member. 
It's annoying to do these updates in the filter.matches function once we 
recalculate the fits, so I do it in this one and return the row to add to 
match.info.
}
