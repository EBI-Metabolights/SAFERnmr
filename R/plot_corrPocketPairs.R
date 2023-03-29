#' Plot Stackplot of Correlation Pocket Pairs
#'
#' This function plots a stackplot of correlation pocket peaks (same data as the heatmap, but plotted like spectral regions)
#' 
#' @param pocketPairs Result of corrPocketPairs_al(...) with elements named "peakMap", "corr", "cov", and "regions".
#' @param region A vector of integers specifying the columns of the peakMap to be included in the plot.
#' @param vshift The vertical shift to apply to the distributions for stacking.
#' @param xvect The x values at which to plot the distributions.
#' @param show Logical. Whether or not to display the plot.
#' 
#' @return The corrPocketPairs stack plot.
#' 
#' @examples
#' # Load data
#' pocketPairs <- corrPocketPairs_al(...)
#' 
#' # Plot correlation heatmaps
#' plot_corrPocketPairs(pocketPairs, region = 1000:1500)
#' 
#' 
#' @export
plot_corrPocketPairs <- function(pocketPairs, region = NULL, vshift = 5,
                                 xvect = NULL, show = TRUE){
  
  if (is.null(xvect)){xvect <- 1:ncol(pocketPairs$peakMap)}
  if (is.null(region)){region <- 1:ncol(pocketPairs$peakMap)}
  
  ##### Get data in line and filter ########

    peakMap <- pocketPairs$peakMap[,region]
    f <- peakMap %>% is.na %>% "!"(.)
    cc <- pocketPairs$corr[,region]
    cv <- pocketPairs$cov[,region]
    inds <- pocketPairs$regions[,region]
      inds <- inds - min(region) + 1
      inds[!f] <- NA
      cc[!f] <- NA
      cv[!f] <- NA
      
  ##### Plotting ########

  b <- cv
  bi <- inds
  
  peaks <- peakMap[f] %>% unique
    peakMap <- peakMap - min(peaks) + 1
    peaks <- peaks - min(peaks) + 1
  
  cmat <- matrix(NA, nrow = max(peaks), ncol = ncol(b))
  
  # Only loop through cols with peaks
    colsWithPeak <- f %>% t %>% rowSums(na.rm = TRUE) %>% ">"(.,0) %>% which
  
  for (i in colsWithPeak){
    
    # Map the indices from the compact corrmat, a, -> expanded offset matrix, cmat
      colinds_c <- bi[f[,i], i]
      rowinds_b <- f[,i] %>% which
      
      inrange <- 
        colinds_c < ncol(cmat) & colinds_c > 0 &
        rowinds_b < nrow(cmat) & rowinds_b > 0

      # Pull the corrs into their spots on the nxn matrix
      linds_cmat <- sub2indR(peakMap[f[,i], i] %>% .[inrange], # use peak number as cmat row number
                             colinds_c[inrange],
                             nrow(cmat))
      cmat[linds_cmat] <- b[rowinds_b[inrange], i]
  }
    
  # Make the plot using stackplot
    # since ggridges geom_ridgeline won't plot disconnected distributions, like geom_density_ridges does,
    # we can force separation by putting each peak region on its own row.
    
    g <- stackplot(cmat, vshift = vshift, hshift = 0, xdir = "reverse", show = show, xvect = xvect)
    
    return(g)
}