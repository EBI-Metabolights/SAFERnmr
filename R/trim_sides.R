#' Trim the sides of a matrix to remove columns with all NA values
#'
#' @param mat a matrix or data.frame to be trimmed
#' @param out character string indicating the output type. Possible values are "mat"
#' (default), which returns the trimmed matrix, or "inds", which returns the indices
#' of the columns that were kept.
#'
#' @return a matrix with the same number of rows as \code{mat} but fewer columns,
#' with columns containing NA values taken out.
#'
#' @examples
#' mat <- matrix(c(1,2,NA,4,NA,NA,7,8,NA), nrow = 3)
#' trim_sides(mat)
#'
#' @export
trim_sides <- function(mat, out = "mat"){
  # Enforce matrix format (so columns are trimmed)
    if (!is.matrix(mat)){mat <- matrix(c(mat), nrow = 1)}
  
  # Does each column have only NAs?
    cols.withVals <- !(   colSums(is.na(mat))   == nrow(mat))
  
  # Error handling in case < 2 columns have values
    if (sum(cols.withVals)<2){
      warning('trim.sides: matrix has < 2 columns with values. Result will have 0 columns/NULL inds.')
      inds <- NULL
    } else {
      inds <- cols.withVals %>% which %>% range %>% fillbetween
    }
  
  # If all we want are the inds, give those
  if (out == "inds"){
    return( inds )
  }
  
  mat.out <- mat[,inds, drop = FALSE]
  if (ncol(mat.out) == 1){mat.out <- t(mat.out)} # take care of transposition if apply() is dumb and flips it
  return(mat.out)
}