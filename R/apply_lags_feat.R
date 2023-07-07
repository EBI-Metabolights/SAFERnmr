#' Given lags and feature object for single feature (including positions and subset),
#' make a lagged spectral matrix that is aligned to the feature in its region. 
#'
#' @param feat feature object for single feature (including positions and subset)
#' @param xmat spectral matrix
#' @param lags vector of lag values (feat$position + lag gives spec pos)
#' 
#' @return list with lagged spectral region indices, values, subset, and range of columns in xmat
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @author MTJ
apply_lags_feat <- function(feat, xmat, lags){
  
  
  indsmat <- lapply(1:length(lags), function(x) {
    # x <- 1
    sub2indR(rows = feat$ss[x], 
             cols = feat$position + lags[x], 
             m = nrow(xmat))
    
  }) %>% do.call(rbind,.)
  
  valsmat <- xmat[indsmat] %>% Reshape(length(feat$ss), length(feat$profile))
  
  return(list(inds = indsmat,
              vals = valsmat,
              ss = feat$ss,
              cols.range = range(indsmat, na.rm = T) %>% ind2subR(nrow(xmat)) %>% .['cols'] %>% unlist %>% as.numeric))
}