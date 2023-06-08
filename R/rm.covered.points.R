#' Remove points on each line which are covered by the AUC of lower rows (to give
#' the effect of spectra covering each other up when rastering instead of plotting
#' multiple geom_area curves in a loop). 
#' 
#' If only mat (spectral matrix), then calc for mat. 
#' 
#' If feats and ss.rows are supplied, assume we want to remove feat points that 
#' are covered by mat. Return feats matrix after filtering (NA addition)
#' ss.rows are the row inds in mat corresponding to each row of feats.
#' length(ss.rows) == nrow(feats), and all(unique(ss.rows) %in% 1:nrow(mat))
#'
#' @param mat xmat subset (row and col) that will be plotted
#' @param feats features lined up with cols of mat. Corresponding row in mat (for each feature) is given in ss.rows.
#' @param ss.rows row in mat for each feature (row in feats) - often duplicates in this vector
#' @return filter for the features matrix
#' @importFrom magrittr %>%
#' 
#' @export 
rm.covered.points <- function(mat, feats = NULL, ss.rows = NULL){

  mat[points.covered(mat)] <- NA
  
  if (!is.null(feats) & !is.null(ss.rows)){
  # Check if additional inputs for features are valid
    valid <- all(unique(ss.rows) %in% 1:nrow(mat)) &
    length(ss.rows) == nrow(feats)
  
  if (!valid){stop('check inputs')}

    feats[points.covered.feats(mat, feats, ss.rows)] <- NA
    mat <- feats
  }
  
  return(mat)
  # simplePlot(mat, xdir = 'normal')
}