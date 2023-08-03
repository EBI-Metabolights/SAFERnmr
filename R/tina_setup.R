#' TINA Setup
#'
#' Pre-compute a bunch of useful stuff for TINA and store in an organized way ("feature" object)
#' Given a list of peak regions and a matrix of spectra, creates a matrix to house
#' all the features and returns it in a list format along with additional 
#' information about the subsets and regions of interest.
#'
#' @param a A list of peak regions.
#' @param xmat A matrix of spectra.
#'
#' @return A list containing the feature stack, position stack, subsets, and 
#' region information.
#'
#' @examples
#' # Create a mock input list of peak regions and a matrix of spectra
#' peaks <- list(list(finalRegion = c(1, 2, 3), ref.vals = c(1, 2, 3), 
#'                    ref.idx = c(1, 2, 3), subset = c(1, 2)),
#'               list(finalRegion = c(4, 5, 6), ref.vals = c(4, 5, 6), 
#'                    ref.idx = c(4, 5, 6), subset = c(3, 4)))
#' spectra <- matrix(1:12, ncol = 4)
#' 
#' # Run the function
#' tina_setup(a = peaks, xmat = spectra)
#' @importFrom magrittr %>%
#'
#' @export
tina_setup <- function(a, xmat){
      
      message("TINA Setup...")

      # Build matrices to house all the features and positions
        # Extract region inds for each feature
          relinds <- lapply(1:length(a), function(x){
            tryCatch(
              expr = {
                       a[[x]]$finalRegion %in% a[[x]]$ref.idx %>% which
                
              }, error = function(cond) NULL)
          })

        # Build the matrices
          v <- rep(NA, relinds %>% unlist %>% max)
          
          featureStack <- lapply(1:length(a), function(i){
            v[relinds[[i]]] <- a[[i]]$ref.vals
            return(v)
          }) %>% do.call(rbind,.)

          positionStack <- lapply(1:length(a), function(i){
            v[relinds[[i]]] <- a[[i]]$ref.idx
            return(v)
          }) %>% do.call(rbind,.)
               
          v <- rep(FALSE, nrow(xmat))
          subsets <- lapply(1:length(a), function(i){
            v[a[[i]]$subset] <- TRUE
            return(v)
          }) %>% do.call(rbind,.)
          
          subset.sizes <- subsets %>% rowSums

        # Get relative driver within feat
        
            drivers.relative <- lapply(1:length(a), function(x) {
              
              if (is.null(a[[x]]$peak)){
                return(NA)
              } else {
                driver.ind <- which((a[[x]]$ref.idx %>% range %>% fillbetween) %in% a[[x]]$peak)
              }
              if (length(driver.ind) == 0)
              {
                warning('feature ', x , ' driver not found')
                driver.ind <- NA
              } 
              
              return(driver.ind)
              
            }) %>% unlist
            
    return(list(stack = featureStack,
                position = positionStack,
                subset = list(ss.all = subsets,
                              sizes = subset.sizes),
                driver.relative = drivers.relative
                )
           )
}