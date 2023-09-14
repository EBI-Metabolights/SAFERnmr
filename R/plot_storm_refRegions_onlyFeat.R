#' Plot storm Results for reference regions in grid format - only plot feature profile (not subsets) to reduce weight.
#'
#' Wrapper for \code{\link{plot_storm_refRegions}} to generate plots for multiple
#' reference regions in a grid format and save the result to a pdf file.
#'
#' @param s storm output object (element of list containing the storm analysis results for each reference region, from fse.R)
#' @param ppm a numeric vector of ppm values corresponding to the NMR spectra.
#' @param n_xticks an integer indicating the number of x-axis ticks to use in the plots.
#'
#' @return ggplot object, e.g. to place in grid
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom tictoc tic toc
#' @importFrom dplyr tibble
#' @importFrom magrittr %>%
#' @export
plot_storm_refRegions_onlyFeat <- function(ppm,s, n_xticks = 4, plot.title = ''){
  # Plot the storm results (s) from storm_play, etc. 
  # just the profiles, no xmat data
  # 
  # MTJ 2022
   
    ppmRegion = ppm[s$finalRegion]
    
##############
    # Plot the reference in black, bold

    refInReg <- s$finalRegion %in% s$ref.idx %>% which
    ref <- s$ref.vals
      
      # Make a copy of the ref over the plot (final) region with NAs in the gap
        stretchedRef <- rep(NA, length(s$finalRegion))
        stretchedRef[refInReg] <- ref
        
      # Make the plot
        g <- simplePlot(stretchedRef, xvect = ppmRegion, linecolor = 'black', linewidth = 1)
                       
##############   

      # Plot corr boundaries and driver as vert lines  
        g <- g +
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "grey")
        #   geom_vline(xintercept = ppm[s$finalRegion[bounds]],
        #              linetype = 1, col = "grey")
    
      # Add title

      g <- g + ggplot2::ggtitle(plot.title)
      
    return(g)
  
}
