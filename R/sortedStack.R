#' sortedStack
#'
#' Sorts a matrix by row sum and stacks the rows.
#'
#' @param mat The matrix to sort and stack.
#' @param ppm The x-axis values for the columns of the matrix. Defaults to NULL.
#' @param vshift The vertical shift for the stacked plots. Defaults to 10.
#' @param hshift The horizontal shift for the stacked plots. Defaults to 0.
#'
#' @return A stacked plot of the sorted matrix.
#'
#' @import ggplot2
#' @importFrom pracma repmat
#' @importFrom scales breaks_pretty
#' @importFrom reshape2 melt
#'
#' @export
#'
sortedStack <- function(mat, ppm = NULL, vshift = 10, hshift = 0){
  if(is.null(ppm)){ppm <- 1:ncol(mat)}
  mat %>% rowSums(na.rm = TRUE) %>% order %>% mat[.,] %>% stackplot(xvect = ppm, vshift = vshift, hshift = hshift)
}