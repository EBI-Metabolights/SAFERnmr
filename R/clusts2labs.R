#' Convert clusters to feature-sorted cluster labels
#'
#' Given a list of clustered feature inds in list elements, this function returns feature-sorted 
#' cluster labels. 
#'
#' @param clusters A list of cluster-feature pairs.
#' 
#'  @importFrom magrittr %>%
#'
#' @return A numeric vector of feature-sorted cluster labels.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @author MTJ
clusts2labs <- function(clusters){
  # Take cluster-feature pairs, return feature-sorted cluster labels:
  cluster.labs <- lapply(1:length(clusters), 
                           function(x) cbind(rep(x, length(clusters[[x]])),
                                             clusters[[x]] %>% as.numeric)
                           ) %>% do.call(rbind,.)
  return(cluster.labs[order(cluster.labs[,2]),] %>% .[,1]) # sort points
  
}