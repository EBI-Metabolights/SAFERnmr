#' Filters for relevant ref, match.info, and ref-feature info depending on 
#' selected from the resultsViewer panes. Returns a list that can be used for 
#' project_features.stackplot.
#' Filters evidence based on:
#' ref (row in heatmap)
#' columns (samples in scores scatterplot)
#' ref region (ppm of selection box in ref spectrum plot) 
#' ppm.tolerance between ref feats and spec, and ref feats and selected ref region
#' cutoff.residuals.feat 
#' cutoff.residuals.spec
#'
#' @param ref.ind index of the ref in the ss-ref scores mat
#' @param spec.ind spectra indices, from scatterplot of scores for selected ref 
#' @param match.info match.info table used for subsetting
#' @param backfits ref feature fit information
#' @param rfs.used lists of which ref features were used in each score, allows only best evidence
#' @param lib.data.processed dataset-interpolated ref library. contains compound names + data + more
#' @param xmat spectral matrix
#' @param ppm ppm for xmat
#' @param ppm.tolerance maximum allowed distance (ppm) between ref feature range in subset spectra and the selected ref region.
#' @param cutoff.residuals.feat currently not used
#' @param cutoff.residuals.spec currently not used
#'
#' @return list of evidence for given compound match to ref signature, limited to specific samples and ref spectrum region;
#' - ref.ind: index of ref spectrum
#' - ref.info: library info for ref spectrum
#' - rf.fits: list of fit objects for relevant ref features
#'
#'
#' @examples
#' v1 <- c(1, 2, 3, 4, 5)
#' v2 <- c(3, 5, 7, 9, 11)
#' score.wasserstein(v1, v2)
#'
#' @export
select.evidence_refmet <- function(ref.ind = 1, 
                                   spec.ind = NULL,
                                     # Big objects to subset using ref.ind:
                                       match.info,          
                                       backfits,
                                       rfs.used,
                                       lib.data.processed,  
                                     # Spectral data:
                                       xmat, ppm,           
                                     # Filtering thresholds:
                                       ppm.tolerance = 0.2, 
                                       cutoff.residuals.feat = .5,
                                       cutoff.residuals.spec = .5){
    
############ Accessory functions ################          

    # Provide a way out if no evidence is selected:
          no.evidence <- function(ref.ind, ld){
            list(ref.ind = ref.ind,
                 ref.info = ld,
                 rf.fits = NULL)
          }
          
          
    if (is.null(spec.ind)){spec.ind <- 1:nrow(xmat)}
  
    
          
############ Function body ######################     

    # Subset the evidence for annotation of the given ref in all samples: ####

      # Pull data for our selected ref ####
          
          ld <- lib.data.processed[[ref.ind]]
          
      # Find all backfits which involve this ref ####
          
          # inds <- match.info$ref == ref.ind
          # rfs.used is the best evidence rfs for each score coordinate
          
          rfs.selection <- rfs.used$score.mat.coords$ref == ref.ind & 
                              rfs.used$score.mat.coords$ref %in% spec.ind
          if(!any(rfs.selection)){return(no.evidence(ref.ind, ld))}
          
          inds <- match.info$id %in% (rfs.used$tot[rfs.selection] %>% unlist) # inds not always = ids 
          
          rf.specFits <- backfits[inds] %>% lapply(function(rf) rf$fits) %>% unlist(recursive = F)
          
          # Find all those which involve the spectra we're into
            specs <- lapply(rf.specFits, function(x) x$ss.spec) %>% unlist
            rf.specFits <- rf.specFits[specs %in% spec.ind]
            if(length(rf.specFits) == 0){return(no.evidence(ref.ind, ld))}
            
      # Build a slimmed-down object that project_features.stackplot() will accept ####
        # Filter out features whose origins in the data are too far from their ref spec match: ####
        
          # Get the ranges of the ref features in the spectra and in the ref spec: ####
            ref.rngs <- lapply(rf.specFits, function(x) range(x$ref.region, na.rm = T) %>% ld$mapped$ppm[.])
            spec.rngs <- lapply(rf.specFits, function(x) range(x$spec.region, na.rm = T) %>% ld$mapped$ppm[.])
            
          # Filter for fits whose dataset spectrum range is within tolerance of ref spec range: ####
            spec.dist.from.rf <- lapply(1:length(spec.rngs), function(x) range.dist(spec.rngs[[x]], ref.rngs[[x]]) ) %>% unlist
            close.enough <- spec.dist.from.rf <= ppm.tolerance

            if(!any(close.enough)){return(no.evidence(ref.ind, ld))}
            
            rf.specFits <- rf.specFits[close.enough]
              ref.rngs <- ref.rngs[close.enough]
              spec.rngs <- spec.rngs[close.enough]
            
        # Unlist all the ref feats into spectrum-fit ref feats, and expand their spectrum positions: ####  
          fit.feats <- lapply(rf.specFits, function(rff) {
            # Compute on the fly
              rf <- rff$ref.region %>% fillbetween %>% ld$mapped$data[.]
              fit.ref <- as.numeric(rff$fit.ref.coef[1]) + (rf * as.numeric(rff$fit.ref.coef[2]))
              return(fit.ref)
          })
            
          
          fit.positions <- lapply(1:length(rf.specFits), function(rff.ind) 
            {
                rff <- rf.specFits[[rff.ind]]
              # Get the positions in the xrow
                pos <- rff$spec.region %>% fillbetween
                pos[is.na(fit.feats[[rff.ind]])] <- NA
                return(pos)
            }) 
        browser()
        # Get only the best-scoring evidence, to help limit the amount plotted
          # fit.positions <- lapply(rf.specFits, function(rff) 
          #   {
          #     rff$
          #   })
            
        # Make matrices from those lists, add to bestfits list object ####
          
          maxlen.ff <- lapply(fit.feats, length) %>% unlist %>% max
        
          ffint <- lapply(1:length(fit.feats), function(i){
                      c(fit.feats[[i]], rep(NA, maxlen.ff-length(fit.feats[[i]])))
                    }) %>% do.call(rbind,.)
          
          ffpos <- lapply(1:length(fit.positions), function(i){
            
                      c(fit.positions[[i]], rep(NA, maxlen.ff-length(fit.positions[[i]])))
                    }) %>% do.call(rbind,.)
          
          fit.xrow <- lapply(1:length(rf.specFits), function(x) rf.specFits[[x]]$ss.spec) %>% unlist
        
          
        # Build the bestfits object. This will get subsetted within the plotting function during runtime ####
          
            # First, filter based on ...
             
            # Report: 
              bestfits <- list(fit.feats = ffint,        
                               fit.positions = ffpos, # in indices
                               fit.xrow = fit.xrow,
                               fit.rngs = ref.rngs) # in ppm 
              
          # Clean up our data mess
            # rm('fit.feats', 'fit.positions', 'maxlen.ff', 'ffint', 'ffpos', 'fit.xrow', 'fit.rngs', 'spec.rngs')
       

      
    # Make a data pack for spectral viewer: ####
      metab_evidence = list(ref.ind = ref.ind,
                            ref.info = ld,
                            rf.fits = bestfits)
    return(metab_evidence)
      
  }
  
