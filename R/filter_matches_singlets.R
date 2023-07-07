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
                                    feature.stack, ref.mat, 
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
      
        message('filtering out matches:',
                  '\n\t- for which the fit feature has 1 or fewer resonances...',
                  '\n\t- involving ref features with 1 or fewer resonances...',
                  '\n\t- involving not-never-fit feature regions with 1 or fewer resonances...',
                  '\n\t- for which 1 or fewer resonances had at least ', res.area.threshold, 
               ' of their area explained by the matching feature...')
        
        
      # Apply singlet filters to each mi row, and add the results as fields. 
        match.info <- mclapply(1:nrow(match.info), 
                               function(m){
          # print(m)
          
          ########### singlet filter for fit features ################
            mi <- match.info[m, ]
            
            ff <- apply_fit(mi.row = mi, feat.stack = feature.stack, ref.stack = ref.mat)
            
            mask <- !is.na(ff$residuals)
            mi$numpeaks.feat <- pk_maxs(ff$feat.fit, mask) %>% length
          
          ########### singlet filter for fit ref regions ################
          
            mi$numpeaks.ref <- pk_maxs(ff$spec.fit, mask) %>% length
          
          ########### singlet filter for feature-not-never-fit regions ################
            
            peak.quality <- peak.qualities[[which(match.info$feat[m] == pq.featureNumbers)]]
            f.adj <- ff$feat.fit - ff$feat.fit * peak.quality
            f.adj <- f.adj - min(f.adj, na.rm = T)
            mi$numpeaks.feat.nnf <- pk_maxs(f.adj, mask = !is.na(f.adj)) %>% length # not sure if mask here is always the same as above
          
          ########### singlet filter for refpeaks.matched ################
          
            pks <- extractPeaks_corr(ff$spec.fit,plots = F)
            pk.coverage <- lapply(pks$bounds[pks$truePeak], function(p){
              pkinds <- p %>% unlist %>% fillbetween
              return(sum(ff$feat.fit[pkinds],na.rm = T) / sum(ff$spec.fit[pkinds],na.rm = T))
            })
            mi$refpeaks.matched <- sum(pk.coverage > (1-res.area.threshold) & pk.coverage < (1+res.area.threshold))

          ########### Return row #########
          return(mi)
            
        }, mc.cores = ncores)
        
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