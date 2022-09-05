#' matchPeaksToRef
#'
#' Function to match target and reference spectral peaks by chemical shift
#' Goncalo Graca, 16 February 2021, g.gomes-da-graca-at-imperial.ac.uk
#' modified on 13 May 2021
#' modifications added:
#' matches only count if driver peak is found in the reference
#' score is the Jaccard index for the target and reference
#' matched ppms counted only once
#' result output modified to match with other functions
#' @param target target peak.
#' @param driver_ppm parts per million of driver peak.
#' @param reference reference spectral peaks as db file.
#' @param tol tolerance in ppm.
#' @param Itol
#' @param matchMethod matching methodology as a string, defaults to 'basic'
#' @param intensity flag to indicate whether to factor in intensity to matching.
#' @return list of resulting matches.
#' @export
matchPeaksToRef <- function(target, driver_ppm, reference, tol = 0.02, Itol = 20,
                            intensity = FALSE){
  # extract relevant metadata from reference
  hmdb_name <- reference$metadata$hmdb_name
  hmdb_id <- reference$metadata$hmdb_id
  spectrum_id <- reference$metadata$spectrum_id
  reference <- reference$spectrum_peaks
  # transform target into matrix if target contains only one peak
  if(is.null(dim(target))) {
    t <- target[1]
    target <- as.matrix(t(target))
  } else t <- target[,1]
  # normalise target for intensity comparison 
  target[,2] <- 100*target[,2]/max(target[,2])
  # transform reference into matrix if target contains only one peak
  if(is.null(dim(reference))) {
    r <- reference[1]
    reference <- as.matrix(t(reference))
  } else r <- reference[,1]
  # initialize matches and matched_ppms reference and target peaks vectors
  matches <- NULL
  matches <- 0
  matched_ppms <- NULL
  reference_peaks <- length(r)
  target_peaks <- length(t)
  # match using chemical shift and intensity
  if(intensity){
    if(any(abs(driver_ppm - r) <= tol)) {
      for(i in 1:length(t)) {
        if(any(abs(t[i] - r) <= tol)){
          j <- which(abs(t[i] - r) <= tol)
          if(length(j) > 1) {
            k <- which.min(abs(t[i] - r[j]))
            j <- j[k]
          } else NULL
          # matched ppms are removed from the reference vector
          # to prevent multiple matches
          if(abs(target[i,2] - reference[j,2]) < Itol){
            r <- r[-j]   
            matched_ppms <- c(matched_ppms, t[i]) 
          }
        } else NULL
      }
    } else NULL 
    if(!is.null(matched_ppms)) {
      matches <- length(matched_ppms)
    } else NULL
  } else {
  # match using chemical shift  
  if(any(abs(driver_ppm - r) <= tol)) {
    for(i in 1:length(t)) {
      if(any(abs(t[i] - r) <= tol)){
        j <- which(abs(t[i] - r) <= tol)
        if(length(j) > 1) {
          k <- which.min(abs(t[i] - r[j]))
          j <- j[k]
        } else NULL
        # matched ppms are removed from the reference vector
        # to prevent multiple matches
        r <- r[-j]   
        matched_ppms <- c(matched_ppms, t[i])
      } else NULL
    }
 } else NULL 
  if(!is.null(matched_ppms)) {
    matches <- length(matched_ppms)
  } else NULL
  }  
 
  score <- matches / (reference_peaks + target_peaks - matches)
  #browser()
  result <- list(driver_peak_ppm = driver_ppm,
                 hmdb_name = hmdb_name,
                 hmdb_id = hmdb_id,
                 spectrum_id = spectrum_id,
                 matches = matches,
                 reference_peaks = reference_peaks,
                 score = score)
  return(result)
}
