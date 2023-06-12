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
#' @importFrom magrittr %>%
#' @importFrom Rfast rowmeans
#'
#' @export
filter_matches_shiftDelta <- function(match.info,
                                      feature.positions,
                                      ppm,
                                      ppm.tol = 0.5){
  
  
  
      ######################### Calculate deltappm distance (specppm - featureppm)  #############################    
    
          # Get the xmat column of each feature's first element (even if NA)
          
          # This allows much much slimmer data. 
            fpos.el1.xmat <- lapply(1:nrow(feature.positions), function(r) {
              
                            fp <- feature.positions[r, ]
                            first.nonNA <- which.min(fp)
                            
                            return(   fp[first.nonNA] - first.nonNA +1  )
                            
                          }) %>% unlist
              
          # Calc the ppms that correspond to each feature range:
            mean.ppm.feat <- Rfast::rowmeans(
                                              cbind(ppm[fpos.el1.xmat[match.info$feat] + match.info$feat.start],
                                                    ppm[fpos.el1.xmat[match.info$feat] + match.info$feat.end])
                                              )
            
          # Simply convert the ref ranges to ppms:
            mean.ppm.ref <- Rfast::rowmeans(
                                              cbind(ppm[match.info$ref + match.info$ref.start],
                                                    ppm[match.info$ref + match.info$ref.end])
                                              )
          # Take the difference
            match.info$ppm.difference <- mean.ppm.feat-mean.ppm.ref
            
            
            # scattermore::scattermoreplot(sort(match.info[,'ppm.difference']), 
            #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
            # 
            # scattermore::scattermoreplot(sort(abs(match.info[,'ppm.difference'])), 
            #                              1:nrow(match.info), xlab = "ppm difference", ylab = "match")
            # later, potentially: as fraction of known shift range?
            
        
        keep <- abs(match.info$ppm.difference) <= ppm.tol
        match.info <- match.info[keep,]

       #####   
          
          
            return(match.info)
    
}