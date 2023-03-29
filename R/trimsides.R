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
#' trim.sides(mat)
#'
#' @export
trim.sides <- function(mat, out = "mat"){
    if (!is.matrix(mat)){mat <- as.matrix(c(mat), nrow = 1)}
  if (out == "inds"){
    return( apply(mat, 2, function(x) !all(is.na(x))) %>% which %>% range %>% fillbetween)
  }
  return(
    mat[,mat %>% apply(2, function(x) !all(is.na(x))) %>% which %>% range %>% fillbetween, drop = FALSE]
    )
}