#' Scale a vector to the range of a matrix, element-wise
#'
#' @param v A numeric vector to be scaled
#' @param mat A numeric matrix used to define the scaling range
#' @param useInds A logical vector or integer vector indicating the columns of `mat` or the elements of `v`
#'                to use in the calculation. Default is to use all columns or elements.
#'
#' @return A numeric vector with the same length as `v`, scaled to the range of `mat` element-wise.
#'
#' @importFrom magrittr %>%
#'
#'
#' @export
scale.to.minmax <- function(v,             # vector (e.g. feature profile)
                            mat,           # data to scale v to. MATRIX.
                            useInds = TRUE # logical vector of columns of mat/elements of v to use
                            ){ 
    # mat <- mat[,useInds]
    # mat <- mat - min(mat, na.rm = TRUE)
    # v <- vinit - min(v, na.rm = TRUE)
    # v <- v * max(mat, na.rm = TRUE) / max(v, na.rm = TRUE) + min(mat, na.rm = TRUE)

    
    # The calculation only works on points at which neither v nor mat cols are NA.
    # For mat, it's only important that not all the vals of a col are NA.
      
      useInds <- !(   is.na(v) | 
                      (mat %>% is.na %>% matrix(nrow = nrow(mat)) %>% apply(2,all))    
                      # what's with the nrow=1?? big issue..
                  ) & useInds
      
      useInds <- which(useInds) # make into integer inds
    
    vselect <- v[useInds]
    # if (is.null(vselect)){browser()}
    mmin <-  min(mat, na.rm = TRUE)
    vmin <- min(vselect, na.rm = TRUE)
    
    # Get the highest available point in the v that is also not NA in spectrum.
    # This should naturally be a local max in v (once filtered to make sure there
    # are points available for comparison in that col in mat). 
    
      vmax.pt <- vselect %>% which.max %>% useInds[.] # for use in v and mat, NOT vselect
    
    # Base the max end of the scaling on that point in the v and mat.
    
      vmax.at.pt <- v[vmax.pt, drop = FALSE]
      
      mmax.at.pt <- mat[,vmax.pt, drop = FALSE] %>% max(na.rm = TRUE)
      
      # if (any(is.infinite(c(vmax.at.pt, vmin, mmin)))) {browser()}
      # heatmap(mat, Rowv = NA, Colv = NA, scale = "row")
      
      
      # ratio <- ((v[useInds] - apply(mat[,useInds,drop=FALSE],2,mean)) ^ 2) * apply(mat[,useInds,drop=FALSE],2,mean)
      # (mean(ratio) * v) %>% rbind(.,mat) %>% simplePlot
      
    return(
             (v - vmin)/(vmax.at.pt - vmin) *
               # profile points scaled to [0,1]
             (mmax.at.pt - mmin) +
               # range of specdata
             mmin 
             # vshift to match specdata
          )
    
    # *** The max of the feature is matched to the spectrum at that point, the mins 
    # are not anchored. *** 
    
    # if (is.null(useInds)){useInds <- 1:length(v)}

    # Each point votes on a scaling 
      
      # vselect <- v[useInds]
      # ratios <- v[useInds] / mat[,useInds]
      # rbind(v/median(ratios), mat) %>% simplePlot
      # rbind(v/mean(ratios), mat) %>% simplePlot
      # rbind(v/min(ratios), mat) %>% simplePlot # this is what I tried in Matlab - really sensitive to noise
      
    # return(v/median(ratios))                                       

    
}