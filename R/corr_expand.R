#' Expand a peak region to the closest correlation minimums
#'
#' This function takes a correlation vector and an index of a peak and returns the bounds of the region containing the peak, extended to the closest local minima on either side of the peak.
#'
#' @param peak An integer representing the index of the peak in the correlation vector
#' @param localMins A numeric vector containing the indices of local minima in the correlation vector
#' @param vRange A numeric vector containing the range of valid indices for the correlation vector
#' @return A list with two elements: \code{lower}, the index of the closest local minimum to the left of the peak, and \code{upper}, the index of the closest local minimum to the right of the peak
#' @examples
#' v <- c(1, 2, 3, 2, 1)
#' localMins <- c(1, 3, 5)
#' corr_expand(peak = 3, localMins = localMins, vRange = c(1, 5))
#' # Returns: list(lower = 1, upper = 5)
#' @export
corr_expand <- function(peak, localMins, vRange){
    
    upperbound <- (localMins >= peak) %>%  # Take local mins >= peak (allow peak = border)
      which() %>% localMins[.] %>% c(.,vRange[2]) %>% min()  # get the ind of the leftmost one --- inds in wind (if empty, use NA)
                                                             # take lesser between that and right bound of wind. ----  
    lowerbound <- (localMins <= peak) %>%  # Take the local mins <= peak (allow for peak is on border)
      which() %>% localMins[.] %>% c(.,vRange[1])  %>% max()      # get the ind of the rightmost one to left of peak (if empty, use NA)
                                                          # make sure it's not less than the left wind bound 
  
  # Store data within list
  return(list(lower = lowerbound,
              upper = upperbound))
}