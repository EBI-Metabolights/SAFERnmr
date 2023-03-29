#' Return the index/indices in a vector closest to given value(s).
#'
#' Given a vector of values, find the index/indices in another vector which is closest to each value. 
#'
#' @param vals a numeric vector of values to find the closest index/indices to.
#' @param vect a numeric vector of values to search for closest index/indices.
#' 
#' @return a numeric vector of indices in vect which are closest to each value in vals.
#'
#' @examples
#' vectInds(1:5, c(3.5, 4.6, 2.1))
#' 
#' @export
vectInds <- function(vals,vect){
  
  vals <- as.vector(vals)
  inds <- vals
  
# Find closest val(s):
  for (i in 1:length(vals)){
    inds[i] <- which.min(abs(vect - vals[i]))
  }
  
# Return indices
  return(inds)
  
}