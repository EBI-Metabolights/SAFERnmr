#' Calculate fft-based alignment of spectrum to feature.
#' Report lag table.
#'
#'
#' @param feat feature object for single feature (with profile, position, and subset vects, and driver position)
#' @param xmat spectral matrix
#' @param r.thresh correlation threshold
#' 
#' @return lag table for each spectrum in xmat
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @author MTJ
align.spec2feat <- function(feat, xmat, r.thresh = 0.8){
  
  # feat <- list(profile = profile,
  #                  position = f.pos,
  #                  ss = ss,
  #                  driver.relative = driver)
                  
    
    cols.of.x <- c(1,ncol(xmat))
    # simplePlot(feat$profile) + geom_vline(xintercept = feat$driver.relative)
    # simplePlot(xmat[, feat$position]) + geom_vline(xintercept = feat$driver.relative)
    # stackplot(xmat[, feat$position], vshift = 3) + geom_vline(xintercept = feat$driver.relative)
    
    # Expand the position window so the xcorr includes external spectral data
    # - profile padded by NAs so we look for profile only
    # - spec data window is simply expanded so we can actually look around ####
      if(length(feat$driver.relative) == 0) {stop('Driver must be defined in feature object.')}
      feat.expanded <- feat
        
        feat.expanded$position <- feat$position %>% expand_window(within = cols.of.x, 
                                                                  # by = length(feat$profile)-1)
                                                                  by = (length(feat$profile)*0.5) %>% round)
        
        f.pos.in.exp <- (feat.expanded$position %in% na.omit(feat$position)) %>% which %>% range %>% fillbetween
        feat.expanded$profile <- rep(NA, length(feat.expanded$position))
        feat.expanded$profile[f.pos.in.exp] <- feat$profile
        feat.expanded$driver.relative <- feat.expanded$driver.relative + f.pos.in.exp[1]
          # simplePlot(feat.expanded$profile) + geom_vline(xintercept = feat.expanded$driver.relative)
        
      profile.exp <- feat.expanded$profile
      
      nspec <- nrow(xmat)
    
    # Correlate each subset spectrum to the feature profile ####
      
      # Extract the full xmat.ss using profile shape ####
        
        xmat.ss <- xmat[, feat.expanded$position, drop = F]
          # simplePlot(xmat.ss) + geom_vline(xintercept = feat.expanded$driver.relative)
          # simplePlot(feat$profile) + geom_vline(xintercept = feat$driver.relative)
          
      # Align spectra for this feature: ####
      
        # Before: simplePlot(xmat.ss[ss, f.pos.in.exp])
          
          use.inds <- !is.na(feat$profile) # for internal NAs, we'll need to exclude sometimes

          # Get the current subset selection score for each spectrum ####
          
            r.current <- cor(feat$profile[use.inds] %>% c, 
                             xmat[, feat$position[use.inds], 
                                  drop = F] %>% t
                             )
          
          # Calculate best alternative lags ####
            lags <- feat.align.to(align = xmat.ss, to = profile.exp, max.hits = 1)
            feat$ss <- 1:nspec # set this to all for now, in case we pick up any spectra :)
          
          # Try out each lag to see if it improves its spectrum's subset selection score ####
            
            lags.tested <- lapply(feat$ss, function(x) {
              
              # Get the subset selection score for this lag
                r.proposed <- cor(feat$profile[use.inds], xmat[x, feat$position[use.inds] + lags$lag.in.f2[x]])
              
              # only apply lag if improves similarity to initial profile AND passes cutoff
                if (r.proposed > r.current[x] &  
                    r.proposed >= r.thresh){     
                  
                  return(lags$lag.in.f2[x])
                
              # otherwise, don't make any changes. We only make changes if improvement is clear.
                } else {
                  return(0)
                }
  
            }) %>% unlist
          
          # Make the new matrices using our lag choices from above ####
            al.fts <- apply.lags.feat(feat = feat,
                                      xmat = xmat, 
                                      lags = lags.tested)
          
            indsmat <- al.fts$inds
            valsmat <- al.fts$vals
        
        
        # After: simplePlot(valsmat)
        
      # Recalc profile based on feature-aligned spectra: ####

        rvals.updated <- cor(c(feat$profile[use.inds]), t(valsmat[, use.inds]))
        rvals.passed <- rvals.updated >= r.thresh & !is.na(rvals.updated)# subset will probably expand here
                # browser()
        valsmat <- valsmat[rvals.passed, ,drop =F]
        indsmat <- indsmat[rvals.passed, ,drop =F]

      # Update the feature ####
        if (sum(rvals.passed) > 2){
          feat$profile <- cov(valsmat[, feat$driver.relative], valsmat)
          feat$corr <- cor(valsmat[, feat$driver.relative], valsmat)
        } else {
          feat$corr <- rep(NA, length(feat$profile))
        }
        
        
        feat$ss <- which(rvals.passed)
        
    # Format results ####
    
      return(
             list(feat = feat,
                  lags = lags.tested[rvals.passed],
                  rvals = rvals.updated[rvals.passed],
                  valsmat = valsmat,
                  indsmat = indsmat)
             )
}
        
