#' ind2subR
#'
#' Utility function for converting linear indices (as defined in MATLAB) to row
#' and column indices. Got the solution from
#' http://cran.r-project.org/doc/contrib/Hiebeler-matlabR.pdf via
#' https://stackoverflow.com/questions/4452039/converting-between-matrix-subscripts-and-linear-indices-like-ind2sub-sub2ind-in
#'
#' MTJ 28FEB2022
#' @param inds The linear indices to be converted IE linear index (1 = [1,1], 2 = [2,1], 3 = [3,1], 4 = [1,2], 5 = [2,2])
#' @param m The number of rows IE m = nrow(matrix)
#' @return List of row and column indices.
#' @export
ind2subR <- function(inds, m) {
  row.inds <- ((inds - 1) %% m) + 1
  col.inds <- floor((inds - 1) / m) + 1

  return(list(rows = row.inds, cols = col.inds))
}
