#' Plot STOCSY-like plots of corrpocketpairs results in a grid
#' 
#' 
#' @param xmat A matrix of dimensions (n, m), where each row is a pocket
#' @param ppm A vector of ppm values
#' @param pocketPairs list of corrpocket pairs results
#' @param res Resolution of the plot in ppm
#' @param plotLoc A string specifying the directory to save the plot in
#' 
#' @return A PDF file of correlation plots between the specified pocket pairs arranged in a grid
#'
#' @export
plot_corrPocketPairs_grid <- function(xmat, ppm, pocketPairs, res = 4000, plotLoc = "./"){
  # Wrapper for plot_storm_refRegions()
  # print to file in plotLoc
  # 
  # res <- 4000
    regs <- 1:ceiling(length(ppm)/res)
  specregions <- matrix (NA, nrow = res, ncol = length(regs))
    specregions[seq_along(ppm)] <- seq_along(ppm)
    
  # regs <- 1:3
  # g<- plot_corrPocketPairs(pocketPairs, 
  #                      region = specregions[,1] %>% na.omit, 
  #                      vshift = 5, show = TRUE)
  # plot(g)
  
  plots <- lapply(regs, function(x) plot_corrPocketPairs(pocketPairs,
                                                         xvect = ppm[specregions[,x] %>% na.omit],
                                                         region = specregions[,x] %>% na.omit, 
                                                         vshift = 5, show = FALSE))
  
  # How big to make the page? 4 inches for each plot, and grid will be square.
    dim <- 4*round(sqrt(length(regs)))
    
    pdf(file = str_c(plotLoc,"corrPairPlots_grid.pdf"),   # The directory you want to save the file in
        width = dim, # The width of the plot in inches
        height = dim)
  
      gridExtra::grid.arrange(grobs = plots)
  
    dev.off()
  
}
