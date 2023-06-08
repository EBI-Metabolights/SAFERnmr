#' Indicate which points are covered up by a lower spectrum (lower row number)
#' Use to filter add NAs to a v-and-h-shifted spectral matrix, avoids line crossings 
#' in stackplots (without explicit layering of area plots). 
#' * This version also accepts a features matrix (of same size as mat), and ss.rows
#' to prevent features showing up where spectra lower in the stack are covering them. 
#'
#' Intended that row 1 is the bottom spectrum in the plot. 
#'
#' @param mat xmat subset (row and col) that will be plotted
#' @param feats features lined up with cols of mat. Corresponding row in mat (for each feature) is given in ss.rows.
#' @param ss.rows row in mat for each feature (row in feats) - often duplicates in this vector
#' @return filter for the features matrix
#' @importFrom magrittr %>%
#' 
#' @export 
points.covered.feats <- function(mat,      
                                 feats,    
                                 ss.rows
                                 ){
  
  # Because a spec feature is associated with exactly one row of xmat, 
  # we just need to check to see if any spectra lower than that one cover it up in 
  # each column.
  
   # feats <- f.stack
   # mat <- xs
   # ss.rows <- lapply(ss.rows, function(ssr) which(ss.rows.unique %in% ssr)) %>% unlist
   covered <- lapply(1:nrow(feats), function(f) {
     
      # f <- 2
      m <- ss.rows[f] # row of x that matches this feature
      v <- rep(F, ncol(feats))
      
      if (m > 1){
      # Additionally, we don't need to check cols for which the spec-feature is 
      # already NA:
      
        cols.to.assess <- which(!is.na(feats[f, ]))
        
        
        feat.covered <- lapply(cols.to.assess, function(n){
                          
                          column.mat <- mat[,n]
                          feat.val <- feats[f,n]
                          
                            any( 
                                  feat.val < column.mat[ 1 : (m-1) ]  # if m == 1, just gives F
                              )
    
                        }) %>% unlist 
        
        v[cols.to.assess] <- feat.covered
      }
      
      return(v)
        
    }) %>% do.call(rbind,.) # the bottom spectrum is never covered

   
  return(covered)
      # mat[covered] <- NA
      # simplePlot(mat, xdir = 'normal')     

}
