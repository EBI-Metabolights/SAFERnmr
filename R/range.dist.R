#' Calculate range distances
#'
#' Given a pair of ranges, return distance between them (number of elements). Overlapping
#' ranges return negative distance. 
#'
#' @param a range 1
#' @param b range 2
#' @return distance between them
#' @export
#' @importFrom magrittr %>%
#'
#'
range.dist <- function(a,b){
    # a <- c(1, 3)
    # b <- c(4, 5)
    # b <- c(1, 2)
    
    a <- sort(a)
    b <- sort(b)
    # from the 4 bounds, if overlapping, take the middle two (must be the overlap)
      
    # overlaps <- min(a) <= max(b) & min(b) <= max(a)
      
    # length(union) - (length(a) + length(b))
      # union 
        outer.bnds <- c(a,b) %>% range
      # Diff between ranges (negative = overlaps)
        if (is.integer(a) & is.integer(b)){
          # need different behavior for integers?
          # d <- span(outer.bnds) - (span(a) + span(b))
          d <- diff(outer.bnds) - (diff(a) + diff(b))
        } else {
          d <- diff(outer.bnds) - (diff(a) + diff(b))
        }

  return(d)
}