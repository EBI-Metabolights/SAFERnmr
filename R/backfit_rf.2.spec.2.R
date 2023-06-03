#' Back-fit reference features to a subset of spectra
#' 
#' This has gone through several iterations.
#' A list of backfits (full fit objects) to each subset spectrum was kept
#' for each ref feature. The backfit object was really all that was needed to 
#' estimate a score.
#' 
#' Then, I moved to slimming them down by cutting out data - just storing the fit
#' coefficients instead of fit rfs explicitly. This was a great reduction. However, 
#' there was still massive data duplication (e.g. ref, feature, match inds apply 
#' to all subset spectra - nspec times!).
#' 
#' So, currently, the match.info contains the data that applies at the rf level - 
#' feat ind, ref ind, rmse (feature), feat region, refspec region - this would all
#' have been duplicated. The backfits are a list of nrow(match.info) data frames,
#' each storing the ss-level information (fit coeffs, ss ind, bff scores) for that
#' rf. The combination of these two tables contains a truly ridiculous amount of
#' information, as they describe all the shape relationships between the PRCSs and 
#' dataset spectra. 
#' 
#' The match.info table can be used to filter for rfs which apply to a given
#' PRCS and/or spectral region, or to a given feature. Once the row (rf) inds are 
#' obtained, one can easily access the ss-level information for each of the rfs.
#' This can be used, in combination with xmat and refmat to extract the shapes,
#' and fits can be applied for plotting, scoring, etc.. Example:
#'   # 1: select ref
#'    - this narrows down matches without ref
#'    - apply to backfits (indexed by match)
#'    2: select ref region
#'      - this narrows down matches further
#'      - still don't need anything from the backfits
#'      3: select samples
#'        - now backfits needed
#'        - plotting can take place. Access:
#'         - spec region (really just the start)
#'          - ref region
#'          - ref spec
#'          - fit data
#' 
#' Alternatively, one can loop through the rfs relevant to a given PRCS, then pull
#' the ss-level table for each, and use that for scoring.
#'
#' @param m.inds A vector of indices corresponding to the matches between reference features and spectra. Necessary to allow for subsetting of matches (i.e. because of filtering)
#' @param fits.feature A list of fitted features (to reference spectra)
#' @param match.info A data frame containing information about the matches between features and reference spectra
#' @param feature List containing feature intensities and positions
#' @param xmat Full spectral matrix from which the features were derived
#' @param ppm A vector containing the ppm axis of the spectra
#' @param plots A logical value indicating whether to generate plots of the fits and back-fit feasibility scores - only use for a couple of backfits at a time as plots will double the weight of the result
#'
#' @return A list containing the back-fits and back-fit feasibility scores for each spectrum in the subset,
#'         as well as a grid plot of the fits and back-fit feasibility scores.
#'
#' @importFrom ggnewscale new_scale_color
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#'
#'
#' @export
backfit.rfs <- function(match.info, 
                        feature, 
                        xmat, ppm){

  # Chunk the data by ref (more or less) ####
  
    # Sort by ref, so we can group chunks by ref
      mi.order <- order(match.info$ref)
      match.info <- match.info[mi.order, ]
    
    # Assign each row of match.ind to a chunk based on number of cores
      
      ncores <- 5 # pars$par$ncores
      chunk.size <- max(1, nrow(match.info) / ncores)
      m.grp <- ceiling((1:nrow(match.info)) / chunk.size)
      
    # Assign feature, ref.mat, and match.info subsets to each chunk 
      chunks <- lapply(unique(m.grp), function(x) 
      {
        mi <- match.info[m.grp == x, ]
        ref.numbers <- unique(mi$ref)
        feat.numbers <- unique(mi$feat)
        
        # Index these to this chunk's refs and features
          # Record initial numbers
            mi$feat.init <- mi$feat
            mi$ref.init <- mi$ref
          
          # Index for this chunk
            mi$feat.init <- lapply(mi$feat, function(x) which(feat.numbers == x)) %>% unlist
            mi$ref.init <- lapply(mi$ref, function(x) which(ref.numbers == x)) %>% unlist
        
        list(match.info = mi,
             ref.numbers = ref.numbers,
             refs = ref.mat[,ref.numbers, drop = F],
             feat.numbers = feat.numbers,
             feature = list(stack = feature$stack[feat.numbers, , drop=F],
                            position = feature$position[feat.numbers, , drop=F],
                            sfe = feature$sfe[feat.numbers]
                            )
        )
      })
      
    # Clean up big objects
      rm(match.info)
      rm(feature)
      rm(ref.mat)
      
      gc()
      
  # Run through each chunk ####
    t1 <- Sys.time()
        
      backfits.by.chunk <- mclapply(chunks, function(chunk) {
      # backfits <- lapply(chunks, function(chunk) {
        # chunk <- chunks[[1]]
  
        backfits.chunk <- lapply(1:nrow(chunk$match.info), 
                           function(m) {
          # Each core gets:
          # - 1 chunk (~1Mb / 50 rows @ 150 spectra and 130K points)
          # - 1 xmat
          # - 1 ppm vector
          ############# For each match ###############
          
          # Get data, expand fit ####
            
            mi <- chunk$match.info[m, ]
            
            # Rebuild the fit locally
            
              fit <- apply.fit(mi, feat.stack = chunk$feature$stack, ref.stack = chunk$refs)
              ref.region <- c(mi$ref.start, mi$ref.end)
              
            # If sfe has not been done:
            
              feat.cols <- c(mi$feat.start, mi$feat.end) %>% as.numeric %>% chunk$feature$position[mi$feat,.] %>% fillbetween
  
            # If sfe HAS been done, wait to do this at the spectrum level using lags
            
              sfs <- data.frame(lag = chunk$feature$sfe[[mi$feat]]$lags,
                                spec.number = chunk$feature$sfe[[mi$feat]]$feat$ss)
  
          # Calculate rf fit for each spec-feature: ####
          
            fits <- lapply(1:nrow(sfs), function(s){
              
              sf <- sfs[s,]
              
              ss.spec <- sf$spec.number
              spec.cols <- feat.cols + sf$lag
              
              
              # Back-fit the ref region to the filled spec data ####
              
                spec.region <- xmat[ss.spec, spec.cols]
              
                fit.feat2spec <- fit.batman(fit$feat.fit, spec.region, 
                                            exclude.lowest = .5, 
                                            ppm = ppm[spec.cols])
                
                  # plot.fit(fit.feat2spec, type = "simple", ppm = ppm[spec.cols]) %>% plot
                
                fit.ref <- fit.feat2spec$ratio * fit$spec.fit + fit.feat2spec$intercept # here spec.fit is the ref.feature!
                  # plot(fit.ref)
                  # plot(spec.region)
                  # rbind(fit.ref, spec.region) %>% simplePlot(linecolor = "black")
                
                residuals <- fit.ref - spec.region
                  
  
              # Get pct. overshoot vector ####
                posres <- residuals > 0
                  posres[is.na(posres)] <- 0
                  
              # Assume each run of positive overshoot is a resonance. ####
              #    what % of its prominence is overshoot?             ####
              
                runs <- run.labels(posres)
                
                if (!any(runs > 0))
                  {
                    ovs.res <- 0
                    ovs.tot <- 0
                  }
                else{
                  peaks <- extractPeaks_corr(fit.ref, plots = F)
                  valleys <- peaks$bounds %>% unlist %>% unique
                  
                  # res.proms <- prominences(peaks, v2 = fit.ref)
                  max.prom <- range(fit.ref, na.rm = T) %>% diff
                  
                  # For each run, expand outwards to nearest local minima in ref spec
                  
                    worst.res <- lapply(1:max(runs), function(r) {
                        # r <- 1
                        # print(r)
                        inds <- which(runs == r)
                          if(length(inds) < 3 | sum(peaks$truePeak) < 1)
                          {
                            return(list(ratio.res = 0,
                                        ratio.total = 0))
                          }
                          
                        # Which ref peak are we looking at? What's its local prominence?
                        # Expand out from the boundaries of the peak in the positive residuals
                        # until you hit the next ref valley.
                          valleys.below <- valleys <= min(inds)
                          
                          if (any(valleys.below))
                            {lower.bound <- max(valleys[valleys.below])}
                          else{lower.bound <- min(inds)}
  
                          valleys.above <- valleys >= max(inds)
                          
                          if (any(valleys.above))
                            {upper.bound <- min(valleys[valleys.above])}
                          else{upper.bound <- max(inds)}
                          
                          ref.prom <- range(fit.ref[lower.bound:upper.bound], na.rm = T) %>% diff
                          
                        # What's the height of the positive residual peak here?
                          pos.res.prom <- range(residuals[inds], na.rm = T) %>% diff 
                          # (these will only be the positive residuals)
                        
                        # What's the ratio?
                          ratio.res <- pos.res.prom/ref.prom
                          
                        # or, weight using the ref peak's relative importance in the ref feature
                        # * note: this could increase the score if pos.res straddles >1 resonance
                          ratio.total <- pos.res.prom / max.prom
                          
                          return(list(ratio.res = ratio.res,
                                      ratio.total = ratio.total))
                          
                    })
                    
                    ovs.res <- lapply(1:length(worst.res), function(x) worst.res[[x]]$ratio.res) %>% unlist
                    ovs.tot <- lapply(1:length(worst.res), function(x) worst.res[[x]]$ratio.total) %>% unlist
                }
                
              # Return results ####
                                  # match = m,
                                  # ref = mi$ref,
                                  # feat = mi$feat,
                return(data.frame(ss.spec = ss.spec,
                                  fit.intercept = fit.feat2spec$intercept, 
                                  fit.scale = fit.feat2spec$ratio,
                                  spec.start = spec.cols[1],
                                  spec.end = tail(spec.cols,1),
                                  # These aren't necessary - held in the match.info
                                  # ref.start = ref.region[1],
                                  # ref.end = ref.region[2],
                                  bffs.res = max(ovs.res), # tmp "bff"
                                  bffs.tot = max(ovs.tot))) # tmp "bff"
            }) %>% do.call(rbind,.)
            
          # The bff score (backfit feasibility score) is defined as:
          #   1 - overshoot viability score
          #   ovs:
          #   Each protrusion in the positive residuals (ref feat - spec) is scored:
          #     -find the resonance(s) that the protrusion overlaps
          #     -% intensity of the protrusion compared to the overall ref feature height?
          #     -take the max of this for the ref feature fit
               
          # Calculate bff from ovs score ####
            
            fits$bffs.res <- 1-fits$bffs.res
              fits$bffs.res[fits$bffs.res<0] <- 0
            fits$bffs.tot <- 1-fits$bffs.tot
              fits$bffs.tot[fits$bffs.tot<0] <- 0
            
        return(fits)
        })
        
        return(backfits.chunk)
        
      }, mc.cores = ncores)
      
    print(Sys.time()-t1)  
    
  # Unchunk match.info ####
  
    match.info <- lapply(chunks, function(chunk) {
      # chunk <- chunks[[1]]
      # Undo the feature and ref number changes for each chunk - this is stored 
      # in the match table, and the corrections are different for each chunk.
      # (make sure to put the feature and ref indices back the way they were)
      chunk$match.info$feat <- chunk$match.info$feat %>% chunk$feat.numbers[.]
      chunk$match.info$ref <- chunk$match.info$ref %>% chunk$ref.numbers[.]
      return(chunk$match.info)
    }) %>% do.call(rbind,.)
  
  # Unlist the backfit tables so each one matches a row in match.info ####
  # (they are currently split across chunks): ####
  
    backfits <- backfits %>% unlist(recursive = F)
    
    # Temporarily remove match + level - just keep the fits
    # backfits <- lapply(backfits.by.chunk, function(bfc) {
    #   # bfc <- backfits.by.chunk[[1]]
    #     bfc <- lapply(bfc, function(bfc.row) {
    #       bfc.row$fits
    #     })
    #     
    #   return(bfc)
    # }) %>% unlist(recursive = F)
    
    backfits <- lapply(backfits, function(bf.df) {
      # bf.df <- backfits[[1]]
      bf.df$ref.start <- NULL
      bf.df$ref.end <- NULL
      return(bf.df)
    })
     
  # Undo mi.order ####
    match.info <- match.info[]
      # saveRDS(backfits, paste0(tmpdir, '/backfits.init.RDS'))
      # backfits.by.chunk <- readRDS(paste0(tmpdir, '/backfits.init.RDS'))
      
      
      browser()
      
  return(backfits)
}

  # Data format criteria ####
  # We want to be able to quickly filter each of these for: 
  # - ref spec region - comes from match.info
  # - ref number - comes from match.info 
  # - ss sample number - need to access backfits.chunk$ss.spec
  
  # Will need the ability to make this table (or just make it now):
  # Fast way (should be fine) 
          # fit <- bf$fits[[1]]
   
   # fastest if the ref regions are lumped, so only calc. pct.ref once
   # even though it's duplicated a bunch, it only adds one column, and it also
   # allows essentially vectorized ss x ref access
   # - what really matters is organization at the following levels:
   #    - by ss x ref
   #      - rf region in ref
   #      - region in spec
   #    - by rf:
   #    
          #   pct.ref <- sum(fit$ref.region %>% as.numeric %>% fillbetween %>% refmat[fit$ref, .], na.rm = T)
          #     # no need to sum the whole spectrum again; already normed to 1.
          #   
          # # return df (expanded this score to all ss.spec x rf combinations) 
          #   data.frame(match = fit$match, # match # = backfit #
          #              ref.start = fit$ref.region[1],
          #              ref.end = fit$ref.region[2],
          #              ref = fit$ref,
          #              # ref.start = fit$ref.region$ref.start,
          #              # ref.end = fit$ref.region$ref.end,
          #              feat = fit$feat,
          #              ss.spec = lapply(bf$fits, function(f) f$ss.spec) %>% unlist,
          #              bff.res = bf$bffs.res,
          #              bff.tot = bf$bffs.tot,
          #              pct.ref = pct.ref)
  
  # 1: select ref
  #   - this narrows down matches without ref
  #   - apply to backfits (indexed by match)
  #   2: select ref region
  #     - this narrows down matches further
  #     - still don't need anything from the backfits
  #     3: select samples
  #       - now backfits needed
  #       - plotting can take place. Access:
  #         - spec region (really just the start)
  #         - ref region
  #         - ref spec
  #         - fit data
   
   
    
     

