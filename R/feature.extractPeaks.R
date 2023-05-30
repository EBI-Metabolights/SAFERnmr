#' Extract true peaks from feature using extractPeaks_corr. Formats with peak max
#' inds easily accessed. 
#'
#' @param feature feauture profile
#' 
#'
#' @return Peak information for true peaks, including
#'   - pks (all peak info)
#'   - maxpk highest intensity true peak
#'   - pkind.max index of maxpk
#' @importFrom magrittr %>%
#'
#' @export
feature.extractPeaks <- function(feature){
      peaks <- lapply(1:nrow(feature$stack), function(x) {
      
      # Extract the peaks from the feature:
        thisrow <- feature$stack[x,]
        pks <- extractPeaks_corr(thisrow)
        
      # Take max true peak index (not detected on boundaries):
        passpks <- pks$peaks[pks$truePeak]
        maxpk <- thisrow[passpks] %>% c(-Inf,.) %>% which.max - 1 # this is the index in passpks
          # c(-Inf,.) => allows empties in which.max with default = 1. 
          # -1 => if no true peaks, returns 0 (i.e. no lag; desirable behavior)
        
      # Return peak ind (or 0 if no peak found)
        return(list(maxpk = max(c(0,passpks[maxpk])),
                    pks = pks,
                    pkind.max = max(c(0,maxpk))))
          # likewise, default is 0
          
    })
  return(peaks)
}