#' Identify True Peaks
#'
#'
#' @param v a vector of intensities
#' @return a logical vector indicating the true peaks and tail regions
#' @importFrom fillbetween fillbetween
#' @importFrom magrittr %>%
#' @importFrom runs runs.labelBy.lengths
#' @importFrom tidyr fill
#' @importFrom utils extractPeaks_corr
#' @examples
#' v <- rnorm(100)
#' is.true.peak(v)
#'
#' @export
is.true.peak <- function(v){
  # v <- res$profile$profile
  pks <- extractPeaks_corr(v,plots = F)
  not.tail <- which(pks$truePeak) %>% 
    lapply(function(x) pks$bounds[[x]] %>% unlist %>% fillbetween) %>% 
    unlist %>% unique # %>% sort # should already be in order
  filt <- rep(F, length(v))
  filt[not.tail] <- T
  
  # run.lens <- x %>% is.na %>% "!"(.) %>% runs.labelBy.lengths
  # rl.pass <- any(run.lens >= min.runlength)
  return(filt)
}