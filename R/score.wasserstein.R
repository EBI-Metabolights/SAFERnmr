#' Calculate the Wasserstein distance between two vectors using transport::wasserstein1d.
#'
#' @param v1 numeric vector of values.
#' @param v2 numeric vector of values.
#' 
#' @return The Wasserstein distance between v1 and v2.
#'
#' @importFrom transport wasserstein1d
#'
#' @examples
#' v1 <- c(1,2,3,4,5)
#' v2 <- c(3,5,7,9,11)
#' score.wasserstein(v1,v2)
#'
#' @export
score.wasserstein <- function(v1,v2){
    require(transport)
  # Zero-fill NA vals
    na.vals <- is.na(v1 + v2)
      v1[na.vals] <- 0
      v2[na.vals] <- 0
      
    return(transport::wasserstein1d(v1,v2))
}