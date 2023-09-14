#' Plot storm Results for reference regions in grid format
#'
#' Wrapper for \code{\link{plot_storm_refRegions}} to generate plots for multiple
#' reference regions in a grid format and save the result to a pdf file.
#'
#' @param xmat a numeric matrix of NMR spectra with ppm values in the first column.
#' @param ppm a numeric vector of ppm values corresponding to the NMR spectra.
#' @param stormResults a list containing the storm analysis results for each reference region.
#' @param plotLoc a character string indicating the directory to save the output file.
#' @param filename a character string of the output filename. Defaults to "storm_results_refPlots_grid.pdf".
#' @param calcStocsy a logical indicating whether to calculate and include STOCSY scores in the plots.
#' @param n_xticks an integer indicating the number of x-axis ticks to use in the plots.
#'
#' @return a grid of plots of the storm analysis results for each reference region, saved to a pdf file.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom tictoc tic toc
#' @importFrom dplyr tibble
#' @importFrom magrittr %>%
#' @export
plot_stormRefRegions_grid <- function(xmat, ppm, stormResults, 
                                      plotLoc, filename = "storm_results_refPlots_grid.pdf", 
                                      calcStocsy = TRUE, n_xticks = 4){
  # Wrapper for plot_storm_refRegions()
  # print to file in plotLoc
  
  tictoc::tic()
  message("Generating plots")
  
  # If no xmat, assume only plotting profiles (way faster) ####
  if (is.null(xmat)){
    
    myplots <- pblapply(1:length(stormResults), function(x) {  
      plot_storm_refRegions_onlyFeat(ppm, stormResults[[x]], n_xticks = n_xticks, plot.title = x)
    

    })
    message("Printing Plots to pdf")
    gentime <- tictoc::toc()
  } else {
    
    # Otherwise, plot all the ref and spec data for each ####
      myplots <- pblapply(1:length(stormResults), function(x) {
        # print(str_c(x , " / ", length(stormResults)))
        plot_storm_refRegions(xmat, ppm, stormResults[[x]], calcStocsy = calcStocsy, n_xticks = n_xticks, plot.title = x)
        
      })
      gentime <- tictoc::toc()
      message("Printing Plots to pdf")
      message("this will take longer; maybe ~ 3-5 min / 250 plots)")
      message(str_c("Estimated: ", round(diff(c(gentime$tic,gentime$toc)) * 360/18/60, 1), " minutes"))
  }
  
    
  # How big to make the page? 2 inches for each plot, and grid will be square.
    dim <- 3*round(sqrt(length(stormResults)))
  pdf(file = str_c(plotLoc,filename),   # The directory you want to save the file in
      width = dim, # The width of the plot in inches
      height = dim)
  
  gridExtra::grid.arrange(grobs = myplots)
  
  dev.off()
  
}