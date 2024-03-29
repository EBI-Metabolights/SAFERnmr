#' Calculate the intersection or union of two numeric ranges
#'
#' @param a a numeric vector representing a range
#' @param b a numeric vector representing a range
#' @param operation a character string indicating whether to calculate the "intersection" or "union" of the two ranges
#'
#' @return a numeric vector representing the intersection or union of the two ranges, or a vector of NA if the ranges do not overlap
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' range_intersect(c(1, 5), c(3, 6), operation = "intersection")
#' range_intersect(c(1, 5), c(3, 6), operation = "union")
#'
#' @export
range_intersect <- function(a,b,operation = "intersection"){
  # a <- c(1, 2)
  # b <- c(  3,    5)
  if (min(a) <= max(b) & 
        min(b) <= max(a))
    {
    # from the 4 bounds, take the middle two (must be the overlap)
      if(operation == "intersection"){
        overlap.bounds <- c(a,b) %>% sort %>% .[c(2,3)]
      }
      if(operation == "union"){
        overlap.bounds <- c(a,b) %>% range
      }
     
     return(overlap.bounds)
    }
  return(c(NA,NA))

}

#' Calculate pairwise intersection or union between ranges
#'
#' Given a matrix of ranges (rows: features, columns: ranges),
#' calculates the pairwise intersection or union between each range.
#' Returns a matrix with dimensions (ncol(ranges), ncol(ranges)).
#'
#' @param ranges a matrix with dimensions (n x m), where n is the number of features
#' and m is the number of ranges
#' @param operation a string indicating the type of operation to perform. Can be either
#' "intersection" or "union". Defaults to "intersection".
#'
#' @importFrom magrittr %>%
#'
#' @return a matrix with dimensions (ncol(ranges), ncol(ranges)) containing the pairwise
#' intersection or union between each range.
#' @export
range_intersect_all <- function(ranges,operation = "intersection"){
      if(ncol(ranges) == 1){return(diff(ranges))}
      intersections <- matrix(NA, ncol(ranges), ncol(ranges))
      combns <- combn(1:ncol(ranges), 2)
      overlap <- apply(combns, 2, function(x){
                      range_intersect(ranges[,x[1]], ranges[,x[2]], operation = operation) # overlap (size)
                    }) %>% diff

      intersections[sub2indR(rows = combns[1,], 
                             cols = combns[2,], 
                             m = nrow(intersections))] <- overlap
      
      intersections[lower.tri(intersections)]<-t(intersections)[lower.tri(intersections)]
      
      diag(intersections) <- diff(ranges)
      
      return(intersections)
}