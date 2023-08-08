#' Helper function to take fit coeffs in match.info, with feat and ref stacks,
#' to expand to the fit features and spec region. 
#' 
#' apply_fit2 does not assume compression. apply_fit does. 
#' 
#'
#' @param mi.row match.info row (mi.row)
#' @param feat.stack feature matrix, e.g. feature$stack (feats are rows)
#' @param ref.stack ref.mat (refs are cols) 
#' 
#' @return fit feature and spec (ref) regions
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @author MTJ
apply_fit2 <- function(mi.row, v1, v2){
  # mi.row <- match.info[1,]
  
  v1 <- v1[mi.row$feat.start:mi.row$feat.end] %>% scale_between
  v2 <- v2[mi.row$ref.start:mi.row$ref.end] %>% scale_between
  
  # NA-gap
    na.vals <- is.na(v1 + v2)
      v1[na.vals] <- NA
      v2[na.vals] <- NA
      
  # if (scaled){
  #  v1 <- v1 %>% scale_between
  #  v2 <- v2 %>% scale_between
  # }
  return(
         list(feat.fit = v1*mi.row$fit.scale + mi.row$fit.intercept,
              spec.fit = v2)
         )
}
