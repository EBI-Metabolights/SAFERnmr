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
#' Note: as of 26JUL2023, this function will handle the pct.ref as well (instead of during score_matches().)
#'
#' @param match.info A data frame containing information about the matches between features and reference spectra
#' @param feature List containing feature intensities and positions
#' @param xmat Full spectral matrix from which the features were derived
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
backfit_rfs2 <- function(match.info, 
                        feature, 
                        xmat,
                        ref.mat,
                        ncores){
  
  success = TRUE
  emptyRow <- function(){data.frame(ss.spec = NA,
                                        fit.intercept = NA, 
                                        fit.scale = NA,
                                        spec.start = NA,
                                        spec.end = NA,
                                        bffs.res = NA, # tmp "bff"
                                        bffs.tot = NA,
                                        rmse = NA,
                                        rmse.biased = NA)}
  
  # Chunk the data by ref (more or less) ####
    message('\tchunking match.info table, features objects, and ref spectra for distribution to cores...')
    # Sort by ref, so we can group chunks by ref
      mi.order <- order(match.info$ref)
      match.info <- match.info[mi.order, ]
    
    # Assign each row of match.ind to a chunk based on number of cores
      
      chunk.size <- max(1, nrow(match.info) / ncores)
      m.grp <- ceiling((1:nrow(match.info)) / chunk.size)
      refsums <- Rfast::colsums(ref.mat, na.rm = T)
        
    # Compress the features and ref stack
      
      feature.c <- compress_features(feature)
      refmat.c <- compress_stack(ref.mat %>% t)
      
    # Assign feature, ref.mat, and match.info subsets to each chunk 
      chunks <- lapply(unique(m.grp), function(x) 
      {
        mi <- match.info[m.grp == x, ]
        ref.numbers <- unique(mi$ref)
        feat.numbers <- unique(mi$feat)
        
        # Index these to this chunk's refs and features
          # the mi data will be used to index in the ref.mat and feat.stack matrices
          mi$feat <- lapply(mi$feat, function(x) which(feat.numbers == x)) %>% unlist
          mi$ref <- lapply(mi$ref, function(x) which(ref.numbers == x)) %>% unlist
        
        list(match.info = mi,
             ref.numbers = ref.numbers,
             refs = refmat.c %>% cstack_selectRows(ref.numbers),
             ref.sums = refsums[ref.numbers],
             feat.numbers = feat.numbers,
             feature = feature.c %>% select_features(feat.numbers),
             sfe = feature$sfe[feat.numbers])
             # feature = list(stack = feature$stack[feat.numbers, , drop=F],
             #                position = feature$position[feat.numbers, , drop=F],
             #                sfe = feature$sfe[feat.numbers]
                            # )
        # )
      })
      
    # Clean up big objects
      rm(match.info)
      rm(feature)
      rm(ref.mat)
      rm(feature.c)
      rm(refmat.c)
      
      gc()
      
  # Compute over each chunk in parallel ####
    
    t1 <- Sys.time()
      
      # Each core gets:
      # - 1 chunk (~1Mb / 50 rows @ 150 spectra and 130K points)
      # - 1 xmat
      # - 1 ppm vector        
      message('\tcomputing backfits over ', ncores, ' cores...')
      message('\tlarge numbers of matches or large datasets will take some time.')
      # message('\tgo eat or get a coffee...\n\n')
  
      # backfits.by.chunk <- lapply(chunks, function(chunk) {
      backfits.by.chunk <- mclapply(chunks, function(chunk) {
      ############# For each chunk (in parallel): ###############
        # chunk <- chunks[[2]]
        
        backfits.chunk <- lapply(1:nrow(chunk$match.info),
                           function(m) {
          tryCatch({
            
          #############        For each match:         ###############
              # m <- 1
              # print(m)
            # Get data, expand fit ####
              
              mi <- chunk$match.info[m, ]
              
              # Rebuild the feature fit to ref region
                
                # Get the positions out
                  # use apply_fit2 as version with stacks passed instead of cstacks
                  # compressed.feature <- chunk$feature
                  # row.nums <- mi$feat
                lil.feat <- chunk$feature %>% expand_features(row.nums = chunk$feat.numbers[mi$feat])
                lil.ref <- chunk$refs %>% cstack_selectRows(chunk$ref.numbers[mi$ref]) %>% cstack_expandRows
              
                feat.cols <- lil.feat$position[mi$feat.start:mi$feat.end]

                # Apply feature fit to ref region (just using this to extract profiles, really)
                # the ref is the thing getting fit to, but it is scaled
                
                  fit <- apply_fit2(mi, v1 = lil.feat$profile, v2 = lil.ref) #%>% plot_fit()
                    
                  # Calculate % of ref signature covered by this ref feat:
                    ref.reg <- lil.ref[mi$ref.start:mi$ref.end]
                      ref.reg[is.na(fit$feat.fit)] <- NA
                    pct.ref = sum(ref.reg, na.rm=T) / chunk$ref.sums[mi$ref]
                  
              # get info for the individual spec-features
              
                sfs <- data.frame(lag = chunk$sfe[[mi$feat]]$lags,
                                  spec.number = chunk$sfe[[mi$feat]]$feat$ss)
    
            # Calculate rf fit for each spec-feature: ####
            
              fits <- lapply(1:nrow(sfs), function(s){
                tryCatch(
                  {
                    # Get the ref region and spec data: ####
                      # s <- 1
                      sf <- sfs[s,]
                      ss.spec <- sf$spec.number
                      spec.cols <- feat.cols + sf$lag
      
                    # Back-fit the ref region to the filled spec data using the feature ####
                    
                      spec.region <- xmat[ss.spec, spec.cols %>% range(na.rm = T) %>% fillbetween]
                    
                          # Propagate fit ####
                          
                            
                            fit.ref2spec <- fit_batman(fit$spec.fit, spec.region, 
                                                  exclude.lowest = .5)
                            
                            # fit.ref <- fit.feat2spec$ratio * fit$spec.fit + fit.feat2spec$intercept # here spec.fit is the ref.feature!
                              fit.ref <- fit.ref2spec$feat.fit
                              # plot_fit(fit.ref2spec, type = "simple", ppm = ppm[spec.cols]) %>% plot
                              # plot_fit(fit.ref2spec, type = "simple") %>% plot
                              
                            # rvalue is not useful here
                            
                          # Calculate scores ####
                            
                            srf <- rbind(fit.ref, spec.region) %>% scale_between()
                            
                            residuals <- srf[1,] - srf[2,]
                            use <- !is.na(fit.ref - spec.region)
                            if (any(use)){
                              
                            
                              # rbind(srf[1,],srf[2,], resid.biased, rep(0, length(resid.biased))) %>% simplePlot
                              # plot(resid.biased)
                              rmse <- Metrics::rmse(srf[1,use], srf[2,use]) # this would penalize underfits equally
                              
                              not.neg <- residuals >= 0  &  use
                              if(sum(not.neg) == 0){
                                # This will catch NULL and empty vectors. NA, NaNs could get through. 
                                rmse.biased <- 0
                              }else{
                                # sqrt(sum(residuals[not.neg]^2)/length(residuals))
                                # rmse.pos <- Metrics::rmse(srf[1,not.neg], srf[2,not.neg]) # this is too generous
                                # rmse.pos <- Metrics::rmse(srf[1,not.neg], srf[2,not.neg]) # this is too generous
                                # rmse.biased <- sqrt(sum(resid.biased[not.neg])/length(residuals)) # too generous
                                resid.sq.biased <- residuals[not.neg]*srf[1,not.neg]
                                rmse.biased <- sqrt(sum(resid.sq.biased)/sum(not.neg))
                              }
                            }
                            
                          # Get pct. overshoot vector ####
                            posres <- not.neg
                              residuals <- fit.ref - spec.region # NOT scaled
                              # posres[is.na(posres)] <- 0 # no longer any NAs
                              
                          # Assume each run of positive overshoot is a resonance. ####
                          #    what % of its prominence is overshoot?             ####
                          
                            runs <- run.labels(posres)
                            
                            if (!any(runs > 0))
                              {
                                ovs.res <- 0
                                ovs.tot <- 0
                            }else{
                              # Calculate the backfit overlaps (res and tot)
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
                            # Failed peak extraction does NOT error here, results in 0s for ovs scores
                        
                        # Return results ####

                            return(data.frame(ss.spec = ss.spec,
                                              fit.intercept = fit.ref2spec$intercept, 
                                              fit.scale = fit.ref2spec$ratio,
                                              spec.start = spec.cols[1],
                                              spec.end = tail(spec.cols,1),
                                              bffs.res = max(ovs.res), # tmp "bff"
                                              bffs.tot = max(ovs.tot),
                                              rmse = rmse,
                                              rmse.biased = rmse.biased)) # tmp "bff"
                  }, 
                  error = function(cond){
                    emptyRow()
                  }
                )
                
              }) %>% do.call(rbind,.)
              
              fits <- fits[!is.na(rowSums(fits)),]
              
            # Calculate bff from ovs score ####
            # 
            # The bff score (backfit feasibility score) is defined as:
            #   1 - overshoot viability score
            #   ovs:
            #   Each protrusion in the positive residuals (ref feat - spec) is scored:
            #     -find the resonance(s) that the protrusion overlaps
            #     -% intensity of the protrusion compared to the overall ref feature height?
            #     -take the max of this for the ref feature fit
              
              fits$bffs.res <- 1-fits$bffs.res
                fits$bffs.res[fits$bffs.res<0] <- 0
              fits$bffs.tot <- 1-fits$bffs.tot
                fits$bffs.tot[fits$bffs.tot<0] <- 0

              fits$pct.ref <- pct.ref
              
              
            # Return the fits dataframe rows (minimal data)  ####
              return(fits)
              
          }, error = function(cond){
            # If the whole match fails, return an NA fits df row
              return(emptyRow())
          })
        })
        
        return(backfits.chunk)
        
      }, mc.cores = ncores
      )
      
     
    
  message('\tbackfit calculations finished. \n')
  print(round(Sys.time()-t1))
  message('\tun-chunking and formatting results...\n')
  
  # Unchunk match.info ####
  
    match.info <- lapply(chunks, function(chunk) {
      # Undo the feature and ref number changes for each chunk - this is stored 
      # in the match table, and the corrections are different for each chunk.
      # (make sure to put the feature and ref indices back the way they were)
      chunk$match.info$feat <- chunk$match.info$feat %>% chunk$feat.numbers[.]
      chunk$match.info$ref <- chunk$match.info$ref %>% chunk$ref.numbers[.]
      
      return(chunk$match.info)
      
    }) %>% do.call(rbind,.)
  
  # Unlist the backfit tables so each one matches a row in match.info ####
  # (they are currently split across chunks): ####
    # saveRDS(backfits.by.chunk, paste0(pars$dirs$temp, '/backfits.init.RDS'))
      # backfits.by.chunk <- readRDS(paste0(tmpdir, '/backfits.init.RDS'))
    backfits <- backfits.by.chunk %>% unlist(recursive = F)

  # Undo mi.order for both objects ####
    message('\tunsorting the match.info and backfits...\n')
    match.info <- match.info[order(mi.order),]
    backfits <- backfits[order(mi.order)]

  # Remove NA rows
    # backfits
      backfits <- lapply(backfits, function(x) x %>% is.na %>% rowSums %>% "=="(0) %>% x[.,])
      bad <- lapply(backfits, function(x) x %>% is.null) %>% unlist
      backfits <- backfits[!bad]
    # matches
      match.info <- match.info[!bad]
      
    message('\tbackfitting completed on ', length(backfits), ' ref-features.')
    # message('\thopefully your lunch was nice and you are now properly caffeinated...\n')
    if(length(backfits) != nrow(match.info)){
      success = FALSE
      warning('not all matches had backfits')
    }

  # Return list ####
    return(list(match.info = match.info,
                backfits = backfits,
                all.succeeded = success))
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
          #     # no need to sum the whole spectrum again; already sums to 1.
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
   
   
    
     

