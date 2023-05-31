#' Match Filtering
#'
#' Filter and process matched peaks based on user-specified criteria, such as singlet removal and ppm distance filtering.
#' Also calculate backfits of ref-feats (reference subsignatures fished out by feature matches to reference spectra) to each individual dataset spectrum.
#' Note: ref-feat fit is obtained using the original feature, then backfit feasibility scores are calculated. BFF scores indicate the extent to which ref-feat resonances do not exceed actual spectral signal. Specifically, each single ref-feat resonance positive residual is evaluated as a fraction of the total ref-feat height. Because the absence of a feasible fit for any single reference resonance invalidates a match, the worst-violating resonance gives the score for the whole ref-feat in a particular spectrum. Note that a ref-feat receives a BFF score for each spectrum it is fit to.
#'
#' @param pars A list of input parameters.
#' @return A list of filtered and processed matched peak information, including back-fits to the original spectra.
#' @import yaml
#' @importFrom magrittr %>%
#' 
#' @import pbapply
#' 
#' @export filter.matches
filter.matches <- function(pars){

  message('--------------------------------------------------------------')
  message('-------------------     Match Filtering    -------------------')
  message('--------------------------------------------------------------')
  message('\n\n\n')
  
################ Read parameters file ##################
  
  
  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)

##################################################################################################################
 # Read data and set up ####

    message("Loading data from files...\n\n\n")
    
    fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
    feature <- readRDS(paste0(this.run, "/feature.final.RDS"))
    ref.mat <- readRDS(paste0(this.run, "/temp_data_matching/ref.mat.RDS"))
    matches <- readRDS(paste0(this.run, "/matches.RDS"))
    cluster <- readRDS(paste0(this.run, "/cluster.final.RDS"))

    # Format matches ####
      
      matches.split <- split(matches, names(matches))
        rm(matches)
      match.info <- do.call(rbind, matches.split$matches)
        rownames(match.info) <- NULL
        
      # peak.qualities aligns with fit and match.info now, but feat number does not index it (there are missing features)!
        pq.featureNumbers <- unique(match.info[,'feat']) # this does not sort (just for good measure)
        
      peak.qualities <- matches.split$peak.quality
      fits.feature <-  matches.split$fits %>% unlist(recursive = F)
        rm(matches.split)
      
        
######################### Remove singlets ############################################

      # Do the filtering (functionalized)
        res <- filter.matches_singlets(match.info, fits.feature, 
                                       peak.qualities, pq.featureNumbers, 
                                       pars$matching$filtering$res.area.threshold)
        match.info <- res$match.info
        fits.feature <- res$fits.feature
        
######################### Propagate matches to feature clusters  ############

        # For each match, produce other matches based on clusters
        
        # Option 1: 
        # - loop through clusters
        # - copy all matches from key to other cluster members
        # - what gets copied?
        #   - match.info row = ref.feature
        #   - all that's needed is to test the ref feature in the local regions 
        #     x spectra from which the cluster member was extracted, as gained 
        #     from the sfe
        #   - fstack.row is the row of the match matrix (subset of feature matrix)
        
        message('\n Propagating matches to cluster members...')
        matched.feats <- match.info$feat %>% unique
        new.data <- mclapply(matched.feats, function(fstack.row) {
            # fstack.row <- matched.feats[1]
            # print(fstack.row)
          
          # Which cluster did it belong to?
            # Behind each row of fstack is 1 or more feature indices
              
              # For now, just assume these are the row numbers in the original featureStack (with all cluster members)

                # fstack.row <- matched.feats[[1]]
                clust.number <- which(cluster$keys %in% fstack.row) # just an index, not the key feature
          
            
          # Pull cluster info
          
            clust.key <- cluster$keys[clust.number] # actual name of the key
            cluster.members <- cluster$groups[[clust.number]]
            clust.info <- cluster$info[[clust.number]]
            
            # Put lag table in terms of matching to key
            
              clust.lags <- clust.info$lag.table %>% pw.lags.relative.to(clust.key)
            
              # lag.features(feature$stack, clust.lags, to = clust.key) %>% trim.sides %>% simplePlot(linecolor = 'black')

          # Match info for key feature
            clust.matches <- match.info[match.info$feat == clust.key, ]
            
 ############ # Specific to each cluster member: ######################################################################
 
            nonkey.members <- cluster.members[cluster.members != clust.key]
          
            if (length(nonkey.members) > 0){
              
              # Copy match info for each cluster member (recopy first one), and recalculate:
              # - match info (adjust position)
              # - fit to ref 
              
                member.data <- lapply(nonkey.members, function(cluster.member){
                  # cluster.member <- nonkey.members[1]
                  
                    # Copy the match info to other member ####
                      new.matches <- clust.matches
                      new.matches$feat <- cluster.member
                      
                    # Update the lag based on cluster.info ####
    
                      # Are the lags additive? Yes, subtract lag (f1 to f2) from lag (feat to ref)
                      # Does the feat.start, feat.end need to be adjusted?
                      lag <- clust.lags$lag.in.f2[ clust.lags$f1 == cluster.member ]
                      new.matches$lag <- new.matches$lag + lag
                      
                      # Default is the lag which puts this feature along the ref at the same spot as the key feature.
                      # Update all:
                      
                        new.matches[, c("ref.start", "ref.end")] <- new.matches[, c("ref.start", "ref.end")] - lag
                    
     # ####                   
     ############ # For each ref feat for that cluster member: ####################################################
     
                    rfs.new <- lapply(1:nrow(new.matches), function(nm.row){
                      rf <- new.matches[nm.row,]
                      # rf <- 2
                      ref.num <- rf$ref
                      ref.reg <- rf[c("ref.start", "ref.end")] %>% unlist %>% fillbetween
                      feat.reg <- rf[c("feat.start", "feat.end")] %>% unlist %>% fillbetween
                      
                    # Get the lags and fits for the new matches to rfs
                    
                      # Extract out the feature
                      
                        feat <- feature$stack[cluster.member, ]
                        
                      # Update the feature bounds to match THIS feature, not the key feature ####
                        # First, expand to use the whole feature (so we know how to ss the ref region):
                          ref.feat <- ref.inds <- rep(NA, length(feat))
                            ref.inds[feat.reg] <- ref.reg
                            ref.inds <- complete.indsVect(ref.inds)
                            
                        # Then, re-subset the feature and ref inds to match new feature:
                          f.inds <- feat %>% trim.sides(out = "inds") %>% range
                          rf[c("feat.start", "feat.end")] <- f.inds
                          rf[c("ref.start", "ref.end")] <- ref.inds[f.inds]
    
                        # Finally, update our temp vars for these two:
                          feat.reg <- rf[c("feat.start", "feat.end")] %>% unlist %>% fillbetween
                          ref.reg <- rf[c("ref.start", "ref.end")] %>% unlist %>% fillbetween
                          
                      # Make ref feat same size as full feature ####
                        # fill in the values
                          ref.feat[feat.reg] <- ref.mat[ref.reg, ref.num] # remember that ref mat is transposed
    
                      # Calculate the fit between the initial rf and the cluster.member ####
                        
                        fit1 <- fit.leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
                          # fit1$plot %>% plot
                        
                      # # Try to optimize the alignment a bit more ####
                      # 
                      #   small.adj <- rbind(feat, ref.feat) %>% feat.align(max.hits = 1, counter = F)
                      #     feat.rf.new.al <- rbind(feat, ref.feat) %>% lag.features(small.adj, to = 2) # "to" doesn't matter for just 2; it's bidirectional 
                      #                                                                                 # and lag.features will sort it out
                      #     # fit2 <- fit.leastSquares(feat.rf.new.al[1,], feat.rf.new.al[2,], plots = T, scale.v2 = T)
                      #     #   fit2$plot %>% plot
                      #   
                      #   # Remake the ref from the new region using the adjustment (don't update rf yet) ####
                      #     ref.reg <- (rf[c("ref.start", "ref.end")] %>% 
                      #                   unlist %>% fillbetween) + 
                      #                   small.adj$lag.in.f2[1]
                      #     
                      #     ref.feat <- ref.inds <- rep(NA, length(feat))
                      #       ref.inds[feat.reg] <- ref.reg
                      #       ref.feat <- ref.mat[ref.inds, ref.num] # remember that ref mat is transposed
                      # 
                      #     fit2<- fit.leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
                      #     # fit2$plot %>% plot
                      #       
                      #   # If the adjusted aligment is a better fit, keep the adjusted lag ####
                      #   
                      #     if (fit2$rmse < fit1$rmse){
                      #       
                      #       # Make the adjustments to rf
                      #         rf$lag <- rf$lag + small.adj$lag.in.f2[1] # just like above, use first of two rows (bidirectional)
                      #         rf[c("ref.start", "ref.end")] <- rf[c("ref.start", "ref.end")] + small.adj$lag.in.f2[1]
                      # 
                      #       # Use fit2
                      #         fit <- fit2
                      #       
                      #     } else {
                      #       
                      #       # Otherwise, don't update the postions, and use fit1
                      #         fit <- fit1
                      #       
                      #     }
                          
                        # Skip alignment
                          
                          fit <- fit1
                          
                        # Either way, we now have a fit. Propagate that information:
                        
                          rf <- updateMatchInfoRow(rf, fit)
                          
                        # Return the data in a list
                        
                          return(list(rf = rf,
                                      fit = fit))
                    })
                        
                    # Build updated new.matches
    
                      new.matches <- rfs.new %>% lapply(function(rf) rf$rf) %>% do.call(rbind,.)
    
                    # Build fits list
    
                      new.fits <- rfs.new %>% lapply(function(rf) rf$fit)
    
                      rm(rfs.new)
    
                  return(list(match.info = new.matches,
                              fits = new.fits))                # # Build updated new.matches
    
                  # return(rfs.new)
                }) 
                
              # Extract updated match.info
              
                member.matches <- member.data %>% lapply(function(member) member$match.info) %>% do.call(rbind,.)
                    
              # Extract updated fits 
              
                member.fits <- member.data %>% lapply(function(member) member$fits)
              
            } else {
              
              member.matches <- NULL
              member.fits <- NULL
              
            }

              
  ###########################################################################################################                    
            
            return(list(match.info = member.matches,
                        fits = member.fits))
            
        }, mc.cores = pars$par$ncores) %>% unlist(recursive = F)
        new.data <- split(new.data, names(new.data))
          match.info <- rbind(match.info, 
                              new.data$match.info %>% do.call(rbind,.))
          row.names(match.info) <- NULL
          fits.feature <- c(fits.feature, new.data$fits %>% unlist(recursive = F, use.names = F) %>% unlist(recursive = F))
          
        # Re-filter for corr, pval
        
          keep <- match.info$rval >= pars$matching$r.thresh & 
            match.info$pval <= pars$matching$p.thresh
          
          match.info <- match.info[keep, ]
          fits.feature <- fits.feature[keep]

 ######################### Calculate deltappm distance (specppm - featureppm)  #############################

        # source('./../span.R')
        # source('./../filter.matches_shiftDelta.R')

        message('\nFiltering out matches > ', pars$matching$filtering$ppm.tol, ' ppm away...')
        res <- filter.matches_shiftDelta(match.info, feature, ppm = fse.result$ppm, fits.feature,
                                         ppm.tol = pars$matching$filtering$ppm.tol)
        match.info <- res$match.info
        fits.feature <- res$fits.feature

 ######################### Back-fit reference to spectra  #############################    
    
      message('Back-fitting ref-feats to each spectrum in the relevant subset...\n\n')
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      
    # Back-fit each matched reference region to the subset spectra
      # adjusted for sfe
        m.inds <- 1:nrow(match.info)
        match.info$id <- m.inds
        t1 <- Sys.time()
        backfits <- backfit_ref.feats.2.subset.specs(m.inds, fits.feature, match.info, 
                                                      feature, # now carries sfe data
                                                      xmat, ppm, plots = F) # plots are heavy and time-expensive!
        
        Sys.time() - t1
        message('Saving backfits...\n\n\n')
        saveRDS(backfits, paste0(this.run,"/backfits.RDS"))

 # ########## save filtered data ########################################################################
         
        message('Saving split and filtered match data...\n\n')

        saveRDS(match.info, paste0(this.run, "/match.info.RDS"))

        saveRDS(fits.feature, paste0(this.run, "/fits.RDS"))

        saveRDS(peak.qualities, paste0(this.run, "/peak.qualities.RDS"))
   
  message('-----------------------------------------------------------------')
  message('-----------------  Matching Filtering Complete ------------------')
  message('-----------------------------------------------------------------')
  
  
  # Plot a subset of the backfits to check them out ####
                     
                    # # Generate plot result for a given match
                    #       m <- 100
                    #       ss.fit <- backfit_ref.feats.2.subset.specs(m.inds[m], fits.feature, match.info,
                    #                                   feature,
                    #                                   xmat, ppm, plots = T) %>% unlist(recursive = F)
                    #       ss.fit$gridplot %>% plot
                    #       m <- m + 1
                    # 

       ######################### Plotting #############################
        # m.inds <- order(match.info[,"rmse.weighted"]/match.info[,"refpeaks.matched"], decreasing = F)
        # # m.inds <- order(match.info[,"rmse"], decreasing = F)
        # 
        # 
        # x <- 0 # x <- x - 1
        # x <- x + 1
        # m <- m.inds[[x]]
        # # m <- nrow(match.info)
        # match <- fits.feature[[m]]
        #   
        # peak.quality <- peak.qualities[[  which(feature.inds.skipEmpties == match.info[m, c("feat")] )  ]]
        # 
        # thisref.ppm <- match.info[m, c("ref.start","ref.end")] %>%
        #                   as.numeric %>% fillbetween %>% fse.result$ppm[.]
        # # match %>% plot.fit(., type = "simple", ppm = thisref.ppm) %>% plot
        # # match %>% plot.fit(., type = "auc", ppm = thisref.ppm) %>% plot
        # match %>% plot.fit(., type = "color.line", ppm = thisref.ppm, 
        #                    color = peak.quality) %>% plot
        # # plot(sort(match.info[,"rmse.weighted"])) 
        # #   cutoff <- median(match.info[,"rmse.weighted"])
        # #   abline(h = cutoff, col = 'red')
        # # peak.quality %>% plot
}
  