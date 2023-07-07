#' Min-max scale the feature profiles (rows) to highest true peak in each.
#' Uses peaks which were already extracted. 
#'
#'
#' @param featureStack feature profile matrix (each row is a feature)
#' @param peaks peaks for each, from feature_extractPeaks
#'
#' @return scaled feature matrix
#' 
#' @importFrom magrittr %>%
#'
#'
#' @export
feature.scale <- function(featureStack, peaks)
{
  scaledFeatures <- lapply(1:nrow(featureStack), function(x) {
           
        # Extract the peaks from the feature:
            v <- featureStack[x, ,drop = FALSE]
            
          # Focus on the highest (soon-to-be center) peak
            pks <- peaks[[x]]$pks
            maxpk <- peaks[[x]]$pkind.max
            bounds <- pks$truePeak %>% which %>% .[maxpk] %>% pks$bounds[[.]]
          
          # Use the intensity of the bound with the lower valley 
            smallerBound <- bounds %>% unlist %>% v[.] %>% min
            
          # Use the peak max intensity
            pkmax <- peaks[[x]]$pks$peaks[maxpk] %>% v[.]
            
          # Scale so peak is between 0 and 1
            v <- (v - smallerBound)/(pkmax-smallerBound) # 
                  
        # Return scaled feature (or 0 if no peak found)
            return(v)
    }) %>% do.call(rbind, .)
  return(scaledFeatures)
}
