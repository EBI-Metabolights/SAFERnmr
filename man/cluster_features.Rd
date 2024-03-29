% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_features.R
\name{cluster_features}
\alias{cluster_features}
\title{The CLUSTERING part of tina ####
Ultimately, it doesn't matter much what clustering method is used, as this is
primarily a means of combining highly similar features to reduce the computational
burden of pairwise comparisons to reference spectra.
Taking this part out into its own function to allow reversion upon failure.}
\usage{
cluster_features(
  pars,
  feature,
  min.features = 1000,
  do.clustering = F,
  max.features = NULL
)
}
\arguments{
\item{pars}{A list of parameters for the TINA pipeline.}

\item{feature}{A feature obj}

\item{min.features}{minimum features required to do OPTICS}

\item{do.clustering}{actually do clustering. If FALSE, save dummy clusters objs.}
}
\value{
Nothing, saves cluster files and plots to pars$dirs$temp
}
\description{
The CLUSTERING part of tina ####
Ultimately, it doesn't matter much what clustering method is used, as this is
primarily a means of combining highly similar features to reduce the computational
burden of pairwise comparisons to reference spectra.
Taking this part out into its own function to allow reversion upon failure.
}
