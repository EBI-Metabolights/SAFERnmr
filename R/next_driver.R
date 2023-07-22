#' Select a new driver based on the ref profile
#' (accessory function for storm_pairplay() ) 
#' 
#' Option 1 - 'simpleMax': 
#'   - standard option, traditionally used in STORM
#'   - select max intensity in the ref.profile (usually covariance w/current driver)
#'   - no bounds
#'   
#'   However, this can lead to STORM walking around a lot, to the point that the ref 
#'   no longer contains the initial driver at all (therefore not extracting the relic
#'   for the peak of interest). We can remedy this situation by doing peak-based 
#'   driver selection. Using feat.extractPeaks():
#'   
#' Option 2 - 'pkMax': 
#'   - select max intensity in the ref.profile (usually covariance w/current driver), 
#'     but must belong to a true peak (lower values on both sides)
#'   - returns bounds (for good measure)
#'   
#' Option 3 - 'samePk': 
#'   - pick driver within bounds of true ref peak that contains previous (current) driver.
#'   - return bounds as well
#' 
#' All options express driver and bounds indices according to ref.idx, as this is
#' how STORM usually wants them and that's the most basic indexing level. 
#' 
#' For peak-based options, if no peak is found which satisfies the criteria, a NULL
#' result is returned for both list values.
#'
#' MTJ 18MAY2023
#' 
#' 
#' @param current.driver current driver position, should match an element in ref.idx
#' @param ref.profile covar profile for ref. can have NAs or not (usually omitted)
#' @param ref.idx indices for each value in ref.profile, usually ppm indices
#' @param behavior see options. Default is samePk.
#' 
#' @return list with driver index and bounds (using ref.idx vals, typically cols of xmat) 
#'
#' @importFrom magrittr %>%
#'
#' @export
next_driver <- function(current.driver, 
                       ref.profile,
                       ref.idx,
                       behavior = 'samePk')
{
  # Option 1: No bounds, just return max in ref signature
  
    if (behavior == 'simpleMax'){
        new.driver <- which.max(ref) %>% ref.idx[.]    
        return(list(idx = new.driver,
                    bounds = NULL))
    }

  # Otherwise, do peak-based (only true peaks) driver selection: 
  
    # Extract the peaks from the NA-expanded ref profile:
      ref.profile.expanded <- ref.idx.expanded <- rep(NA, span(ref.idx))
        # if (any(is.na((ref.idx-min(ref.idx, na.rm = T)+1)))){browser()}
        if (any(is.na(ref.idx))){ref.idx <- na.omit(ref.idx)}
        ref.profile.expanded[ref.idx-min(ref.idx, na.rm = T)+1] <- ref.profile
        ref.profile <- ref.profile.expanded
        ref.idx.expanded[ref.idx-min(ref.idx, na.rm = T)+1] <- ref.idx
        ref.idx <- ref.idx.expanded
        # simplePlot(ref.profile, xvect = ref.idx)
        
      ref.pks <- feature_extractPeaks(list(stack = matrix(ref.profile, nrow = 1))) %>% .[[1]]

    # Option 2: Just pick the max and bounds of the highest peak. 
    
      if (behavior == 'maxPk'){
        
            new.driver <- ref.pks$maxpk %>% ref.idx[.]
            driver.bounds <- ref.pks$pks$bounds %>% .[[ref.pks$pkind.max]] %>% as.numeric %>% ref.idx[.]
        
      }

    # Option 3: Keep driver within ref profile peak bounds of previous driver (to keep from walking)
    
      if (behavior == 'samePk'){
          
        
          # Index from true peak(s) containing driver 
            peak.locs <- ref.pks$pks$peaks[ref.pks$pks$truePeak]
            # store as matrix 
            peak.bnds <- ref.pks$pks$bounds[ref.pks$pks$truePeak] %>% unlist %>% ref.idx[.] %>% matrix(nrow = 2) 

          driver.pk <- (current.driver >= peak.bnds[1, ] & 
                        current.driver <= peak.bnds[2, ] ) %>% which
          
          if (length(driver.pk) > 1){
            
            # If old driver fell on a pk boundary, choose the higher peak
              
              driver.pk <- (peak.locs[driver.pk] %>% ref.profile[.] %>% c(-Inf,.) %>% which.max - 1) %>% driver.pk[.]
          }
          
          driver.bounds <- peak.bnds[, driver.pk] 
          
          new.driver <- peak.locs[driver.pk] %>% ref.idx[.]
       
      }
    
      
    # Check results, and return list:
      if (length(new.driver) == 0 | length(driver.bounds) == 0){
        new.driver <- NULL
        driver.bounds <- NULL
      }
      
      return(list(idx = new.driver,
                  bounds = driver.bounds))
      
}


# Other approaches I've tried
        # Allow ref center to travel to the nearest uphill local maximum
          # ref.max <- ref.idx[climb(pt = ref.max, pks = pks)]
          
        # ref.max = localMaxima(ref)
