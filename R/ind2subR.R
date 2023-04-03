#' Convert linear indices to row and column subscripts
#'
#' This function takes a numeric vector of linear indices \code{inds} and a positive integer \code{m}, representing the number of columns in a matrix, and returns a list with two elements: \code{rows}, a vector of row subscripts corresponding to the linear indices, and \code{cols}, a vector of column subscripts corresponding to the linear indices. Meant to replicate MATLAB's ind2sub.
#'
#' @param inds A numeric vector of linear indices from a matrix
#' @param m A positive integer representing the number of rows in the corresponding matrix
#' @return A list with two elements: \code{rows}, a vector of row subscripts corresponding to the linear indices, and \code{cols}, a vector of column subscripts corresponding to the linear indices
#' @examples
#' inds <- 1:9
#' m <- 4
#' ind2subR(inds, m)
#' # Returns: list(rows = c(1, 2, 3, 4, 1, 2, 3, 4, 1), cols = c(1, 1, 1, 1, 2, 2, 2, 2, 3))
#' @export
ind2subR <- function(inds, m) {
  inds <- as.vector(inds)
  row.inds <- ((inds - 1) %% m) + 1
  col.inds <- floor((inds - 1) / m) + 1

  return(list(rows = row.inds, cols = col.inds))
}
