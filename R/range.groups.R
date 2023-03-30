#' Divide Ranges into Groups with Overlapping Sections
#'
#' Given a set of ranges, return a list of groups where each group contains
#' ranges with overlapping sections.
#'
#' @param ranges A matrix with two rows, representing ranges in the form c(min,max).
#' @param operation The type of intersection operation to use when determining overlap. "intersection" or "union"
#' @return A list of groups, where each group is a vector of indices
#'   corresponding to overlapping ranges.
#' @export
#'
#'
range.groups <- function(ranges, operation = "intersection"){
  return(
          range.intersect.all(ranges, operation) %>% is.na %>% "!"(.) %>% 
          igraph::graph.adjacency(.,mode="undirected", weighted=NULL) %>%
            igraph::max_cliques(., min = 0, max = NULL) %>%
            lapply(function(x) (as.numeric(x)))
    )
}