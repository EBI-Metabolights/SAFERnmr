#' Gives logical mask for points covered up by points rows whose index is lower. 
#' By replacing the TRUE points with NA, one can give the effect of plotted spectra
#' being layered over top of each other (such that the AUCs of lower row inds 
#' cover up points in the higher rows). In other words, make an overlay plot with
#' no criss-crossed lines without introducing unnecessary figure size and slowdown. 
#'
#' Intended that row 1 is the bottom spectrum in the plot. 
#'
#' @param mat xmat subset (row and col) that will be plotted
#' @param feats features lined up with cols of mat. Corresponding row in mat (for each feature) is given in ss.rows.
#' @param ss.rows row in mat for each feature (row in feats) - often duplicates in this vector
#' @return filter for the features matrix
#' @importFrom magrittr %>%
#' @author MTJ
#' 
#' @export 
points_covered <- function(mat){
  
  # Because a spec feature is associated with exactly one row of xmat, 
  # we just need to check to see if any spectra lower than that one cover it up in 
  # each column.
  
   # feats <- f.stack
   # mat <- xs
   # ss.rows <- lapply(ss.rows, function(ssr) which(ss.rows.unique %in% ssr)) %>% unlist
   covered <- lapply(1:nrow(mat), function(m) {
      # mat <- xmat[1:50, 1:100]
      # Additionally, we don't need to check cols for which the spec-feature is 
      # already NA:
      pt.covered <- lapply(1:ncol(mat), function(n){
                          
                            any( 
                                  mat[m, n] < mat[ 1 : (m-1), n]  # if m == 1, just gives F
                              )
    
                    }) %>% unlist 
      pt.covered
      
    }) %>% do.call(rbind,.) # the bottom spectrum is never covered

   
  return(covered)
      # # mat is spectral matrix
      # mat[points_covered(mat)] <- NA
      # simplePlot(mat, xdir = 'normal')     

}
