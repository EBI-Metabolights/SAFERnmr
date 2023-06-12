#' Calculate distances between two sets of ranges (each given as 2 x n columns)
#'
#' Given a pair of ranges, return distance between them ( for integer vectors, 
#' this should be the number of elements between them). Overlapping
#' ranges return negative distance. 
#'
#' @param a range 1 (given as column, or multiple columns)
#' @param b range 2 (given as column, or mulitple columns) - must be same length as a!
#' @return distance between them
#' @export
#' @importFrom magrittr %>%
#'
#'
range_dist <- function(a,b){
    # a <- c(1, 3)
    # b <- c(4, 5)
    # b <- c(1, 2)
    # a <- spec.rngs %>% t
    # b <- ref.rngs %>% t
    
    a <- sortPairs(a) # sort each column 
    b <- sortPairs(b) # sort each column 
    
    # from the 4 bounds, if overlapping, take the middle two (must be the overlap)
      
    # overlaps <- min(a) <= max(b) & min(b) <= max(a)
      
    # length(union) - (length(a) + length(b))
      # union 
        
        outer.bnds <- Rfast::colMinsMaxs(rbind(a, b))
        
      # Diff between ranges (negative = overlaps)
        if (is.integer(a) & is.integer(b)){
          # need different behavior for integers?
          # d <- span(outer.bnds) - (span(a) + span(b))
          warning('Integer range diff under development - for indices, span instead of diff may be needed')
          d <- (outer.bnds[2, ]-outer.bnds[1, ]) - ( (a[2, ]-a[1, ]) + (b[2, ]-b[1, ]) )
        } else {
          d <- (outer.bnds[2, ]-outer.bnds[1, ]) - ( (a[2, ]-a[1, ]) + (b[2, ]-b[1, ]) )
        }

  return(d)
}