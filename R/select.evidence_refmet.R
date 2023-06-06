# Filters for relevant ref and 
# - ref list
# - backfit list
# - match.info table
# 
# The idea here is to figure out which rfs, of those used for scoring, apply
# to the ref and subset spectra selected by ref.ind and spec.ind, respectively.
# These rfs can then be passed onto the plotting function.

select.evidence_refmet <- function(ref = NULL, 
                                   sample = NULL,
                                     # Big objects to subset using ref.ind:
                                       match.info,          
                                       backfits,
                                       rfs.used,
                                       lib.data.processed,  
                                     # # Spectral data:
                                     #   xmat, ppm,           
                                     # Filtering thresholds:
                                       ppm.tolerance = pars$matching$filtering$ppm.tol, 
                                       cutoff.residuals.feat = .5,
                                       cutoff.residuals.spec = .5){
    
############ Accessory functions ################          

    # Provide a way out if no evidence is selected:
          no.evidence <- function(ref, ld){
            list(ref.ind = ref$id,
                 ref.info = ld,
                 rf.fits = NULL)
          }
          
          
    # if (is.null(spec.ind)){spec.ind <- 1:nrow(xmat)}
  
    
          
############ Function body ######################     

    # Subset the evidence for annotation of the given ref in all samples: ####

      # Pull data for our selected ref ####
          
          ld <- lib.data.processed[[ref$number]]
          
      # Find all backfits which involve this ref ####
          
          # inds <- match.info$ref == ref.ind
          # rfs.used is the best evidence rfs for each score coordinate
          
          # One of these lists for each sample selected
          rfs.selection <- rfs.used$score.mat.coords$ref == ref$id & 
                              rfs.used$score.mat.coords$ss.spec %in% sample$number
          
          if(!any(rfs.selection)){return(no.evidence(ref$number, ld))}
          
          inds <- match.info$id %in% (rfs.used$tot[rfs.selection] %>% unlist) # inds not always = ids 
          match.info <- match.info[inds, ]
          backfits <- backfits[inds ]
          
        # Ensure the rfs are all specific to this ref..
      
          correct.refnum <- match.info$ref == ref$id
          match.info <- match.info[correct.refnum, ]
          
          backfits <- backfits[correct.refnum ]
        
          # Add ref start and end as fields for each spec feature (depending on corresponding match info row)
          
            rf.specFits <- lapply(1:length(backfits), function(x) 
            {
              rf <- backfits[[x]]
              rf$ref.start <- match.info$ref.start[x]
              rf$ref.end <- match.info$ref.end[x]
              rf$id <- match.info$id[x]
              rf
            }) %>% do.call(rbind, .)
            
          # Find all those which involve the spectra we're into
          
            rf.specFits <- rf.specFits[rf.specFits$ss.spec %in% sample$number, ]
            if(nrow(rf.specFits) == 0){return(no.evidence(ref$number, ld))}
            
      # Build a slimmed-down object that project_features.stackplot() will accept ####
        # Filter out features whose origins in the data are too far from their ref spec match: ####
        
          # Get the ranges of the ref features in the spectra and in the ref spec: ####
            ref.rngs <- rbind(ld$mapped$ppm[rf.specFits$ref.start], ld$mapped$ppm[rf.specFits$ref.end])
            spec.rngs <- rbind(ld$mapped$ppm[rf.specFits$spec.start], ld$mapped$ppm[rf.specFits$spec.end])
            
          # Filter for fits whose dataset spectrum range is within tolerance of ref spec range: ####
            spec.dist.from.rf <- range.dist(spec.rngs, ref.rngs)
            close.enough <- spec.dist.from.rf <= ppm.tolerance

            if(!any(close.enough)){return(no.evidence(ref$number, ld))}
            
            rf.specFits <- rf.specFits[close.enough, ]
              ref.rngs <- ref.rngs[, close.enough]
              spec.rngs <- spec.rngs[, close.enough]
              
   ########## Expand the fits? ##########
   
        # Unlist all the ref feats into spectrum-fit ref feats, and expand their spectrum positions: ####  
          fit.feats <- lapply(1:nrow(rf.specFits), function(x) {
            # Compute on the fly
              rff <- rf.specFits[x, ]
              rf <- rff$ref.start:rff$ref.end %>% ld$mapped$data[.]
              fit.ref <- as.numeric(rff$fit.intercept) + (rf * as.numeric(rff$fit.scale))
              return(fit.ref)
          })
            
          
          fit.positions <- lapply(1:nrow(rf.specFits), function(x) 
            {
                rff <- rf.specFits[x,]
              # Get the positions in the xrow
                pos <- rff$spec.start:rff$spec.end
                pos[is.na(fit.feats[[x]])] <- NA
                return(pos)
            }) 
        
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
          
          fit.xrow <- rf.specFits$ss.spec
        
          
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
      metab_evidence = list(ref.ind = ref$number,
                            ref.info = ld,
                            rf.fits = bestfits)
    return(metab_evidence)
      
  }
  
