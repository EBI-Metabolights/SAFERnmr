% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotly_functions.R
\name{drawHeatmap}
\alias{drawHeatmap}
\title{Makes a plotly heatmap for a matrix in the style of browse_evidence().}
\usage{
drawHeatmap(mat, dropRowNames = F, clipRowNames = NA, source.name = "heatmap")
}
\arguments{
\item{mat}{a numerical matrix of values, can have rownames and colnames.}

\item{dropRowNames}{whether or not to replace the rownames with 1:nrow(mat)}

\item{clipRowNames}{number of characters to cut off rownames to (if too long)}

\item{source.name}{tag for tracking in shiny apps}
}
\value{
interactive plotly heatmap object
}
\description{
Makes a plotly heatmap for a matrix in the style of browse_evidence().
}
