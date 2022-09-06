#' ppTargetSpec
#'
#' Convert row, column subscripts to linear indices.
#' Comparable to sub2ind in MATLAB.
#'
#' MTJ mjudge\at\imperial.ac.uk
#' 9MARCH2022
#' @param rows rows object
#' @param cols columns object
#' @param m coefficient
#' @return linear indices
#' @export
sub2indR <- function(rows, cols, m) {

  return((cols-1) * m + rows)

}