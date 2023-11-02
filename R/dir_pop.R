#' dir_pop
#'
#' Return the parent directory of a filepath
#' 
#' @param d directory path
#' @return d/..
#' 
#' @importFrom magrittr %>%
#'
#' @export
dir_pop <- function(d, times=1){d %>% strsplit('/') %>% .[[1]] %>% rev %>% .[-(1:times)] %>% rev %>% paste(collapse = '/')}