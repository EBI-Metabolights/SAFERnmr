#' Extract values from a matrix given shifted indices
#'
#' Given a matrix and a matrix of shifted indices (with the same number of rows), extract
#' the values of the original matrix at those shifted indices. The shifted indices can be
#' created by aligning a subregion of the original matrix with a reference vector.
#'
#' @param xmat A matrix containing the values to extract
#' @param shiftedInds A matrix of shifted indices (with the same number of rows as xmat)
#'        indicating which values to extract from xmat
#' @return A matrix containing the values of xmat at the shifted indices given in shiftedInds
#' @export
#'
#' @examples
#' xmat <- matrix(1:9, nrow = 3)
#' shiftedInds <- matrix(c(1, 2, 3, 2, 3, 4, 3, 4, 5), nrow = 3)
#' shiftedInds_toVals(xmat, shiftedInds)
#'
#' # Output:
#' #      [,1] [,2] [,3]
#' # [1,]    1    2    3
#' # [2,]    2    3    4
#' # [3,]    3    4    5
shiftedInds_toVals <- function(xmat, shiftedInds) {

  # Calculate linear inds in xmat from shiftedInds matrix
  linds <- sub2indR(
    rows = (1:nrow(xmat)) %>% matrix() %>% pracma::repmat(1, ncol(shiftedInds)),
    cols = shiftedInds,
    m = nrow(xmat)
  )

  # Extract vals and reformat the vector to the original size
  specRegion <- xmat[linds] %>% pracma::Reshape(
    nrow(shiftedInds),
    ncol(shiftedInds)
  )

  return(specRegion)
}
