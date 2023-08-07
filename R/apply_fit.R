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
apply_fit <- function(mi.row, feat.cstack, ref.cstack){
  # mi.row <- match.info[1,]
  feat <- feat.cstack %>%
    cstack_selectRows(mi.row$feat) %>% 
    cstack_expandRows
  
  ref <- ref.cstack %>%
    cstack_selectRows(mi.row$ref) %>% 
    cstack_expandRows
  
  v1 <- feat[, mi.row$feat.start:mi.row$feat.end] %>% scale_between
  v2 <- ref[mi.row$ref, mi.row$ref.start:mi.row$ref.end] %>% scale_between
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
