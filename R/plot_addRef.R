#' Add reference region to an existing ggplot object
#' 
#' 
#' @param g existing ggplot object (e.g. from simplePlot)
#' @param specRegion x[,ref.idx]
#' @param ref entire ref object - this function uses information from window and signature
#' @param ppm ppm[,ref$wind$inds]
#' @param scaling logical, whether to use new scaling
#' @param type character vector specifying the type of plot to make ("bold line", "line", "shaded area")
#' 
#' @return ggplot object with reference region added
#'
#' @export
plot_addRef <- function(g,specRegion,ref,ppm, scaling = TRUE, type = "bold line"){

############# New ref plot ################
 # * allows for ref gaps
 # * better scaling (0.75 x the STOCSY line scaling)
 # * uses more ref information; inputs simplified
 #  g:          existing ggplot (e.g. from simplePlot)
 #  specRegion: x[,ref.idx]
 #  ref:        entire ref object - this fn uses info from window and signature
 #  ppm:        ppm[,ref$wind$inds]
 
  
  refInReg <- ref$signature$idx.wind
  
  # Sept 2022 scaling (ripped from grid plot/plot_addstocsy)
    ref.vals <- ref$signature$vals
  
  # Make a copy of the ref over the plot (final) region with NAs in the gap
    stretchedRef <- rep(NA, length(ref$wind$inds))
    stretchedRef[refInReg] <- ref.vals
    
    # browser()
    spr <- t(specRegion[,refInReg])
  
  # New scaling (as of 26AUG22)
    if (scaling){
      spr <- spr - min(spr, na.rm = TRUE)
      stretchedRef <- stretchedRef - min(stretchedRef, na.rm = TRUE)
      stretchedRef <- stretchedRef * max(spr, na.rm = TRUE) / max(stretchedRef, na.rm = TRUE) + min(spr, na.rm = TRUE)
      stretchedRef <- stretchedRef * 0.75
    }
  
  
  # Make line plot
    if (any(type %in% "bold line")){
    g <- g + geom_path(data = data.frame(refvals = stretchedRef,
                                         refppms = ppm),
                       na.rm = TRUE,
                       mapping = aes(x = refppms, y = refvals),
                       colour = "black", linewidth = 2,
                       lineend = "round",
                       linejoin = "round",
                       linemitre = "5")
    }
    if (any(type %in% "line")){
    g <- g + geom_path(data = data.frame(refvals = stretchedRef,
                                         refppms = ppm),
                       na.rm = TRUE,
                       mapping = aes(x = refppms, y = refvals),
                       colour = "black", linewidth = g$plot_env$g$theme$line$linewidth,
                       lineend = "round",
                       linejoin = "round",
                       linemitre = "5")
    }

  # Make area plot
    if (any(type %in% "shaded area")){
    g <- g + geom_area(data = data.frame(refvals = stretchedRef,
                                         refppms = ppm),
                       na.rm = TRUE,
                       mapping = aes(x = refppms, y = refvals),
                       fill="pink",alpha=0.4)
    }  
  
  return(g)
}