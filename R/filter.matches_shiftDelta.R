#' Filter matches based on ppm shift difference
#'
#' This function filters the matches in match.info based on the difference in
#' ppm shift between the feature and the match. Only matches with a ppm shift
#' difference less than or equal to ppm.tol are retained.
#'
#' @param match.info A data frame containing the match information; produced in filter.matches()
#' @param feature Object containing features' intensity matrix, position (column inds of spectral matrix) matrix, and other information about the features
#' @param ppm A numeric vector of ppm values; feature$position values index this
#' @param fits.feature List of feature fits to reference spectral regions
#' @param ppm.tol The maximum allowed ppm shift difference between the feature and the match
#'
#' @return A list with two elements: the filtered `match.info` data frame, and the filtered `fits.feature` list which have been filtered by ppm shift difference.
#'
#' @export
filter.matches_shiftDelta <- function(match.info,
                                      feature,
                                      ppm,
                                      fits.feature,
                                      ppm.tol = 0.5){
  
  
  
      ######################### Calculate deltappm distance (specppm - featureppm)  #############################    
    
        
          # Each feature has a delta range/distribution ####
          # Each spec match with a feature has a delta range. 
          # Each feature has a delta range/distribution. Compare them:  ####
            
            match.info[,'ppm.difference'] <- pblapply(1:nrow(match.info), function(m){
              
              fnum <- match.info$feat[m]
              
              f.inds.trim <- match.info[m, c("feat.start","feat.end")] %>% as.numeric
              feat.range <- f.inds.trim %>% feature$position[fnum,.] %>% range(na.rm = T) %>% ppm[.]
              ref.range <- match.info[m, c("ref.start", "ref.end")] %>% as.numeric %>% ppm[.]
              return(mean(feat.range) - mean(ref.range))
              # match.info$feat.end %>% sort %>% plot
              
            }) %>% unlist
            
            # scattermore::scattermoreplot(sort(match.info[,'ppm.difference']), 
            #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
            # 
            # scattermore::scattermoreplot(sort(abs(match.info[,'ppm.difference'])), 
            #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
            # later, potentially: as fraction of known shift range?
            
        
        keep <- abs(match.info$ppm.difference) <= ppm.tol
        match.info <- match.info[keep,]
        fits.feature <- fits.feature[keep]

       #####   
          
          
            return(list(match.info = match.info,
                  fits.feature = fits.feature))
    
}