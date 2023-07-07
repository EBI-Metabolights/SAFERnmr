#' keep_inds_in_bounds.R
#'
#' Given a set of indices to check and a set of indices to check against, returns
#' the subset of check that lies within the bounds of against.
#'
#' @param check A numeric vector of indices to check.
#' @param against A numeric vector of indices to check against.
#'
#' @return A numeric vector of indices within the bounds of against.
#'
#' @examples
#' keep_inds_in_bounds(1:10, 5:15)
#'
#' @export
keep_inds_in_bounds <- function(check,against){
  return(max(  (min(check)),min(against)  ) : min(  (max(check)),max(against)  ))
  
}