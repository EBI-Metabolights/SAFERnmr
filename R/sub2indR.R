#' Convert Subscripts to Linear Indices in Column-Major Order
#'
#' This function converts subscripts of a matrix to their corresponding linear indices in a column-major order. Intended to replicate MATLAB's sub2ind.
#'
#' @param rows The row subscripts.
#' @param cols The column subscripts.
#' @param m The number of rows in the matrix.
#'
#' @return The linear indices corresponding to the input subscripts.
#'
#' @examples
#' sub2indR(1:3, 1:3, 3)
#' # [1] 1 4 7 2 5 8 3 6 9
#'
#' @export
sub2indR <- function(rows, cols, m) {
  return((as.vector(cols) - 1) * m + as.vector(rows))
}
