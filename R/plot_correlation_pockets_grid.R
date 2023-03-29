#' Plot correlation pockets in a grid
#'
#' Wrapper function for \code{\link{plot_correlation_pocket}} to generate a grid of plots
#' from a list of correlation pocket data. Prints the resulting grid to a PDF file
#' located in the specified directory.
#'
#' @param allvals a list of correlation pocket data (output from \code{\link{find_correlation_pockets}})
#' @param plotLoc character string specifying the directory to save the resulting PDF file
#'
#' @return None
#' @export
#'
#'
#' @import ggplot2
#' @import gridExtra
#' @importFrom stringr str_c
plot_correlation_pockets_grid <- function(allvals, plotLoc){
  
  # Wrapper for plot_correlation_pocket()
  # print to file in plotLoc
  # MTJ 2022
  myplots <- lapply(1:length(allvals), function(x) plot_correlation_pocket(allvals[[x]]) )
  
  # How big to make the page? 2 inches for each plot, and grid will be square.
  dim <- 2*round(sqrt(length(allvals)))
  pdf(file = str_c(plotLoc,"correlationPockets_grid.pdf"),   # The directory you want to save the file in
      width = dim, # The width of the plot in inches
      height = dim)
  
  gridExtra::grid.arrange(grobs = myplots)
  
  dev.off()
  
}