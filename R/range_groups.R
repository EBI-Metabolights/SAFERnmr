#' Divide Ranges into Groups with Overlapping Sections
#'
#' Given a set of ranges, return a list of groups where each group contains
#' ranges with overlapping sections.
#'
#' @param ranges A matrix with two rows, representing ranges in the form c(min,max).
#' @param operation The type of intersection operation to use when determining overlap.
#' @return A list of groups, where each group is a vector of indices
#'   corresponding to overlapping ranges.
#' @export
#' @importFrom magrittr %>%
#'
#'
range_groups <- function(ranges, operation = "intersection", clust.by = 'comps'){
  
  if (clust.by == 'comps'){
    x <- range_intersect_all(ranges) %>% is.na %>% "!"(.) %>% 
           igraph::graph.adjacency(.,mode="undirected", weighted=NULL) %>%
           igraph::components(., mode='weak') 
    group.labs <- lapply(1:length(unique(x$membership)), function(i) which(x$membership == i))
    
  }
  
  if (clust.by == 'cliques'){
    
    group.labs <- range_intersect_all(ranges) %>% is.na %>% "!"(.) %>% 
                    igraph::graph.adjacency(.,mode="undirected", weighted=NULL) %>%
                    igraph::max_cliques(., min = 0, max = NULL) %>%
                    lapply(function(x) (as.numeric(x)))
  }      
  return(group.labs)
}