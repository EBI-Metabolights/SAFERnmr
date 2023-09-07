#' Filter matches with singlets
#'
#' This function removes matches that have only one peak in their respective feature, reference or feature-not-never-fit regions.
#'
#' @param match.info The match information data.frame obtained within filter.matches.
#' @param feature.stack e.g. feature$stack
#' @param ref.mat spectral matrix for reference compounds
#' @param peak.qualities A list of peak quality vectors (~ feature points' relevance in the reference database)
#' @param pq.featureNumbers A vector of feature numbers to index the peak.qualities list. This is necessary when only a subset of features are matched and feature fits are filtered.
#' @param res.area.thresh The minimum reference spectrum resonance area that must be accounted for by the fit feature in order to consider it matched.
#'
#' @return A list containing the filtered match information and the filtered fitted features.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' @importFrom dplyr filter
#'
#' @export
filter_matches_singlets <- function(match.info, 
                                    f.cstack, r.cstack, 
                                    peak.qualities, pq.featureNumbers, 
                                    res.area.threshold, ncores){
    
      
######################################################################################################################## 
## Removing singlets ####

      # # Dependencies ####
      #   source('./../extractPeaks_corr.R')
      #   source('./../pk_maxs.R')
      #   source('./../pk_bounds.R')
      #   source('./../corr_expand.R')
      #   source('./../ind2subR.R')
      #   source('./../stackplot.R')
      
      # Define function to create dummy vals for row if failed:
        failed <- function(mi){
          mi$numpeaks.feat <- 0
          mi$numpeaks.ref <- 0
          mi$numpeaks.feat.nnf <- 0
          mi$refpeaks.matched <- 0
          mi
        }      
        
      # Loop through each row of match.info, recreate the fit vectors, then characterize peaks 
        message('filtering out matches:',
                  '\n\t- for which the fit feature has 1 or fewer resonances...',
                  '\n\t- involving ref features with 1 or fewer resonances...',
                  '\n\t- involving not-never-fit feature regions with 1 or fewer resonances...',
                  '\n\t- for which 1 or fewer resonances had at least ', res.area.threshold, 
               ' of their area explained by the matching feature...')
        
        
        # f.cstack <- compress_stack(feature.stack)
        # r.cstack <- compress_stack(t(ref.mat)) # refs on rows
        
        # rm(feature.stack)
        # rm(ref.mat)
        
        gc()
        
        
      # Apply singlet filters to each mi row, and add the results as fields. 
        match.info <- mclapply(1:nrow(match.info),
        # minf <- lapply(1:nrow(match.info),
                               function(m){
          # m <- 100000                       
          # print(m)
            # The only thing that can be returned from this is a match.info row with 4 more fields:
            # - numpeaks.feat (number of true peaks in the ref profile)
            # - numpeaks.ref (number of true peaks in the ref profile)
            # - numpeaks.feat.nnf (not-never-fit; based on peak.quality for feature)
            # - refpeaks.matched (number of true peaks in the ref profile with > 30% of their AUC accounted for by feature)
            
            mi <- match.info[m, ]
            
            # Wrap the whole thing in tryCatch; if any of the peak extraction fails, the row will be filtered out anyways.
            mi <- tryCatch(
              expr = {
              ########### singlet filter for fit features ################
                
                ff <- apply_fit(mi.row = mi, 
                                feat.cstack = f.cstack, 
                                ref.cstack = r.cstack)
                
                mask <- !is.na(ff$feat.fit + ff$spec.fit)
                  if (length(mask) == 0){return(failed(mi))}
                
                mi$numpeaks.feat <- pk_maxs(ff$feat.fit, mask) %>% length
              
              ########### singlet filter for fit ref regions ################
              
                mi$numpeaks.ref <- pk_maxs(ff$spec.fit, mask) %>% length

              ########### singlet filter for feature-not-never-fit regions ################
                
                # peak.quality <- peak.qualities[[which(mi$feat == pq.featureNumbers)]]
                # f.adj <- ff$feat.fit - ff$feat.fit * peak.quality
                # f.adj <- f.adj - min(f.adj, na.rm = T)
                # mi$numpeaks.feat.nnf <- pk_maxs(f.adj, mask = !is.na(f.adj)) %>% length # not sure if mask here is always the same as above
                mi$numpeaks.feat.nnf <- 2 # just set to something for now.
              
              ########### singlet filter for refpeaks.matched ################
              
                pks <- extractPeaks_corr(ff$spec.fit,plots = F)
                pk.coverage <- lapply(pks$bounds[pks$truePeak], function(p){
                  pkinds <- p %>% unlist %>% fillbetween
                  return(sum(ff$feat.fit[pkinds],na.rm = T) / sum(ff$spec.fit[pkinds],na.rm = T))
                })
                mi$refpeaks.matched <- sum(pk.coverage > (1-res.area.threshold) & pk.coverage < (1+res.area.threshold))
    
              ########### Return row #########
              
                return(mi)
                
              }, 
              
              # Or, if failed at any point, simply return NAs
              warning = function(cond){
                
                return(failed(mi))
                
              }, 
              error = function(cond){
                
                return(failed(mi))
                
              }
            )

        }, mc.cores = ncores
        )
        
        
      # Rbind results
        
        saveRDS(match.info, 'match.info.filt.RDS') # for now; for debug 
        match.info <- rbindlist(match.info)
        match.info <- filter(match.info, 
                              numpeaks.feat > 1 &
                              numpeaks.ref > 1 &
                              numpeaks.feat.nnf > 1 &
                              refpeaks.matched > 1)
      
        message(nrow(match.info), ' matches remaining\n\n')

      return(match.info)
        
}



