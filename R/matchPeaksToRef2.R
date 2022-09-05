#' matchPeaksToRef2
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
#' @param Itol tbc
#' @param matchMethod matching methodology as a string, defaults to 'basic'
#' @param intensity flag to indicate whether to factor in intensity to matching.
#' @return list of resulting matches.
#' @export
matchPeaksToRef2 <- function(target, driver_ppm, reference, tol = 0.02,
                            Itol = 20, matchMethod = "basic", intensity = FALSE) {

  # Clean up params
  
    # To keep same argument form as before, if intensity arg is used, switch matchMethod
      if (intensity == TRUE){ matchMethod <- "intensity"}
  
  # extract relevant metadata from reference
    hmdb_name <- reference$metadata$hmdb_name
    hmdb_id <- reference$metadata$hmdb_id
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # print(hmdb_name)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    spectrum_id <- reference$metadata$spectrum_id
    reference <- reference$spectrum_peaks # NOTE: this overwrites reference
    
  # transform target into matrix if target contains only one peak
    if (is.null(dim(target))) {
    t <- target[1]
    target <- as.matrix(t(target))
  } else {
    t <- target[, 1]
  }

  # normalise target for intensity comparison
    target[, 2] <- 100 * target[, 2] / max(target[, 2])

  # transform reference into matrix if target contains only one peak
    if (is.null(dim(reference))) {
    r <- reference[1]
    reference <- as.matrix(t(reference))
  } else {
    r <- reference[, 1]
  }
    
  #########################################################

  # initialize matches and matched_targets reference and target peaks vectors
    matches <- NULL
    matches <- 0
    matched_targets <- NULL
      target_inds <- NULL
    matched_refs <- NULL
      ref_inds <- NULL
    reference_peaks <- length(r)
    target_peaks <- length(t)
    dists <- NULL
    
    matchList = list()
   
  #########################################################     
## There are several matching methods that can be used:
    
    # Match using chemical shift and intensity  
    if (matchMethod == "intensity") {
      if (any(abs(driver_ppm - r) <= tol)) {
        for (i in 1:target_peaks)
        {
          # If any of the reference peaks are within tolerance tol of this target peak
          distances <- abs(t[i] - r) 
          if (any(distances <= tol)) {
              bestMatch <- which.min(distances)   # determine which ref peaks match (best = closest)
    
              # Check if the matched peak falls within the intensity tolerance Itol
                if (abs(target[i, 2] - reference[bestMatch, 2]) < Itol) {
                  r <- r[-bestMatch] # best-matched ppms are removed from the reference vector to prevent multiple matches
                  matched_targets <- c(matched_targets, t[i])
                }
              
          # If not, do nothing
            } else {
              NULL
            }
        }
  
      # If any peaks in target were matched to reference, report number
        if (!is.null(matched_targets)) {
          matches <- length(matched_targets)
        }
      
      } 
    }
    
    #########################################################
    # Match using chemical shift only (2021 algorithm)
    if (matchMethod == "basic"){ 
      # browser()
      # Since the actual ref list will be modified, keep track of the indices in a parallel vector
        indslist <- 1:length(r)
      
      # Do matching if the driver peak is within tolerance of any ref peaks (I suppose the driver is most important)
      if (any(abs(driver_ppm - r) <= tol)) # otherwise, no real point in matching the others?
      {
        
        # Loop through all the target peaks for this ref
          for (i in 1:target_peaks) # go through each target peak
          {
            # Calculate all distances from this target to each ref
              distances <- abs(t[i] - r)
            
            # Check to see if this target is close to any refs (within tolerance)
              if (any(distances <= tol)) # 
              {
                bestMatch <- which.min(distances)   # determine which ref peaks match (best = closest)
                
                # Record the info about the match pair
                  target_inds <- c(target_inds,i)
                    # matched_targets <- c(matched_targets, t[i])
                  
                  ref_inds <- c(ref_inds, indslist[bestMatch])
                    # matched_refs <- c(matched_refs,r[bestMatch])
                  
                  dists <- c(dists,distances[bestMatch])
                  
                  r <- r[-bestMatch] # best-matched ppm is removed from the reference vector to prevent matching to the same ref peak multiple times
                  # NOTE: this means that lower ppm values (first in the target list) will be prioritized.
                  # *** Important***: also need to update the indslist so we can keep track of ref inds:
                  indslist <- indslist[-bestMatch]
              }
          }
        
        # Record number of matches
          matches <- length(target_inds)
        
      } # If no peaks were within tol, then don't modify matches, target_inds, or ref_inds
      
    } 
    
    #########################################################
    # Match using iterative minimum distance or global cluster optimization (MTJ 2022)
    if (matchMethod == "itmin" | matchMethod == "hungarian" | matchMethod == "hungarian_scaled"){  
      
      # Using matchCluster
          
        if (any(abs(driver_ppm - r) <= tol)) # if there are any refs close to driver peak (within tolerance) - otherwise, no real point in matching the others?
          {
            matchList <- matchCluster(t, # target peaks
                                       r, # reference peaks
                                       tol = 0.02, # tolerance (ppm)
                                       method = matchMethod)

# debug when match is obtained
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              # if (!is.null(matchList)){browser()}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            matched_targets <- matchList$matchedTargets
            target_inds <- matchList$matchedTargetInds
            matched_refs <- matchList$matchedRefs
            ref_inds <- matchList$matchedRefInds
            
            dists <- matchList$distances
            
            # If any peaks in target were matched to reference, report number
              if (!is.null(matched_targets)) {
                matches <- length(matched_targets)
            }
          }
    }
  
  ########################################################
  # Gather up scores
    
    score <- matches / (reference_peaks + target_peaks - matches)
    
  ########################################################
  # Produce result struct(s)
    
      # Troubleshooting 
        
        # if (is_null(matched_targets)){browser()}  
    
      # include initial peaklists (result will be duplicated (refs in particular))
    
        peaklists <- list(reference = reference,
                          target = target)
      
      # Format the match pair data
        
        matchlist <- list(ref    = list(ppms = reference[ref_inds,"chemical-shift"],
                                        intensities = reference[ref_inds,"intensity"],
                                        inds = ref_inds),
                          target = list(ppms = target[target_inds,"chemical-shift"],
                                        intensities = target[target_inds,"intensity"],
                                        inds = target_inds),
                          distances = dists
                          )
        
      # Format results struct
        
        result <- list(driver_peak_ppm = driver_ppm,
                        hmdb_name = hmdb_name,
                        hmdb_id = hmdb_id,
                        spectrum_id = spectrum_id,
                        score = score,
                        matches = matches,
                        reference_peaks = reference_peaks,
                        target_peaks = target_peaks,
                        peaklists = peaklists,
                        matchlist = matchlist)
        

        
      # Old format:
                  # result <- list(
                  #   driver_peak_ppm = driver_ppm,
                  #   hmdb_name = hmdb_name,
                  #   hmdb_id = hmdb_id,
                  #   spectrum_id = spectrum_id,
                  #   matches = matches,
                  #   target_peaks = target,
                  #   reference_peaks = reference,
                  #   list(matched_refs = matched_refs,
                  #        matched_targets = matched_targets,
                  #        score = score
                  #   )
                  #   
  return(result)

}