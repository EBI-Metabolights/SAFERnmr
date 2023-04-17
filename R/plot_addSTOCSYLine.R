#' Plot STOCSY line over the current plot
#'
#'
#' @param g plot to add STOCSY line to
#' @param specRegion region of the NMR spectrum to plot
#' @param ppm vector of chemical shift values
#' @param corr matrix of correlations between reference and query spectra
#' @param covar matrix of covariances between reference and query spectra
#' @param driver chemical shift value of the driver peak, if given
#' @param bounds chemical shift values of the boundaries of the region of interest, if given
#'
#' @return plot with STOCSY line added
#'
#' @importFrom transport wasserstein1d
#' @importFrom ggplot2 geom_line geom_vline scale_colour_gradientn
#' @importFrom scales rescale
#' @importFrom colorRamps matlab.like2
#' 
#' @export
plot_addSTOCSYLine <- function(g,specRegion,ppm,
                               corr = NULL, covar = NULL, driver = NULL,
                               bounds = NULL){
  
  corr <- (corr/max(corr, na.rm = T)) %>% c %>% t
  stocsyLine <- covar %>% c %>% t
  spr <- t(specRegion)
  
  # New scaling (as of 26AUG22)
    spr <- spr - min(spr, na.rm = T)
    stocsyLine <- stocsyLine - min(stocsyLine, na.rm = T)
    stocsyLine <- stocsyLine * max(spr, na.rm = T) / max(stocsyLine, na.rm = T) + min(specRegion, na.rm = T)
  
  # Add the final (thresholded) STOCSY result in STORM to the plot
  
  # Plot the stocsy line over the current plot
    
    g <- g + ggnewscale::new_scale_color() +
      geom_line(data = data.frame(refvals = stocsyLine,
                                  refppms = ppm,
                                  corr = corr), 
                mapping = aes(x = refppms, y = refvals, colour = corr),
                size = 1) +
      scale_colour_gradientn(colours = colorRamps::matlab.like2(10),
                             limits = c(-1, 1))
  
  # Plot corr boundaries and driver as vert lines, if given
    if (!is.null(driver)){
      g <- g +
        geom_vline(xintercept = driver, linetype = 2, col = "grey")
    }
    
    if (!is.null(bounds)){
      g <- g + 
        geom_vline(xintercept = bounds,
                   linetype = 1, col = "grey")
    }
    return(g)
   
}