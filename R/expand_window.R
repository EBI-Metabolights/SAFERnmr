#' Expand a window by a certain distance while keeping indices within a given range.
#'
#' @param window A numeric vector representing the original window to be expanded.
#' @param within A numeric vector representing the range within which the expanded window should be contained.
#' @param by A numeric scalar representing the amount by which to expand the window. If not specified, the default is to expand the window by the span of its original values.
#' 
#' @return A numeric vector representing the expanded window.
#' 
#' @importFrom magrittr %>%
#' @import keep_inds_in_bounds
#' 
#' @examples
#' # Expand a window with default arguments
#' expand_window(c(1, 3), 1:5)
#'
#' # Expand a window with custom expansion distance
#' expand_window(c(1, 3), 1:5, by = 2)
#' 
#' @export
expand_window <- function(window, within, by = NULL){
  if (is.null(by)){by <- (max(window) - min(window))+1} # default is span
  wexp <- (min(window)-by):(max(window)+by) %>% 
    keep_inds_in_bounds(check = ., against = within)
  return(wexp)
}