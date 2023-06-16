#' fillbetween
#'
#' Calculate a sequence of values between two numeric indices (e.g., like ':' ), with optional reversal.
#'
#' @param inds a numeric vector of length 2 indicating the starting and ending indices
#' @return a numeric vector representing the sequence of indices between the starting and ending indices
#' @examples
#' fillbetween(c(1, 5)) # returns 1 2 3 4 5
#' fillbetween(c(5, 1)) # returns 5 4 3 2 1
#' fillbetween(c(1, NA)) # throws an error
#' @export
fillbetween <- function(inds){
  # Check inds (surp)
  if (length(inds) <= 1 | any(is.na(inds)) | any(is.infinite(inds)) | any(is.null(inds))){stop('fillbetween: inds contain non-numeric!')} #return(inds)
  flip <- inds[1] > inds[2]
  # if (is.null(flip)){stop('fillbetween: inds are not numeric.')}
  if (flip){
    return(rev(inds[2]:inds[1]))
  } else {
    return(inds[1]:inds[2])
  }
}