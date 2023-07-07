#' Helper function to take fit coeffs in match.info, with feat and ref stacks,
#' to expand to the fit features and spec region. 
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
apply_fit <- function(mi.row, feat.stack, ref.stack){
  # mi.row <- match.info[1,]
  
  v1 <- feat.stack[mi.row$feat, mi.row$feat.start:mi.row$feat.end] %>% scale_between
  v2 <- ref.stack[mi.row$ref.start:mi.row$ref.end, mi.row$ref] %>% scale_between
  # if (scaled){
  #  v1 <- v1 %>% scale_between
  #  v2 <- v2 %>% scale_between
  # }
  return(
         list(feat.fit = v1*mi.row$fit.scale + mi.row$fit.intercept,
              spec.fit = v2)
         )
}
