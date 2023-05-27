#' Give logical indices to evenly subset "select" elements of "from", an array of something.
#'
#' @param select target length of the evenly spaced subset of "from"
#' @param from length of an array (assumption if length == 1), or the array itself
#'
#' @return logical indices of "from" or the thing that "from" is the length of
#'
#' @export
every.nth <- function(select, from){
  # select = pars$storm$number.of.plots
  # from = sum(filt)
    if (length(from) > 1){
      # Assume you want to select elements of from
      from <- length(from)
    }
      
      log.inds <-
                  c(T,rep(F, 
                          max(
                            c(floor(from / select ) - 1,  0)
                          ) # at least 0 times (no -1 args to rep...)
                        )
                    )
      return(log.inds)
  
  
  
}