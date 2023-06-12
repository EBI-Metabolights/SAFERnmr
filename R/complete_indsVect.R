#' Given a sparse (but evenly spaced) vector of indexes,
#' with NAs in the missing places,
#' expand it outwards to fill the entire vector. Basically 
#' linear inter/extrapolation of a perfectly linear vect.
#'
#' @param indVect vector of indices with gaps (can be on the ends), in any order. 
#' 
#' @return vector with values completed
#'
#'
#' @importFrom magrittr %>%
#'
#' @author MTJ
#' @export
complete_indsVect <- function(indVect){
  iv.min <- which.min(indVect)
  startval <- (1 - iv.min) + indVect[iv.min]
  endval <- (startval + length(indVect)) -1
  return(startval:endval)

  # indVect <- new.feat.pos[1,]
  # # Fill in the ranges so the entire ref shape is covered
  #   freg <- which( !is.na( indVect))
  #   return((c(1, length(indVect)) - 
  #               range(freg) + 
  #               range(indVect[freg])) %>% fillbetween)
}

