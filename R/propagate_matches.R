#' Propagate matches from cluster keys to the rest of each cluster.
#'
#' For each match, produce other matches based on clusters
#' - loop through clusters
#' - copy all matches from key to other cluster members
#' - what gets copied?
#'   - match.info row = ref.feature
#'   - all that's needed is to test the ref feature in the local regions 
#'     x spectra from which the cluster member was extracted, as gained 
#'     from the sfe
#'   - fstack.row is the row of the match matrix (subset of feature matrix)
#' - what gets added?
#'   - new feature, lagged and fit to ref
#'   - just store the coefficients and positions
#'   
#' This boils down to a match.info row for every feature. Fits can be rebuilt
#' on the fly from this and the feature/spectral matrix/ref matrix.
#'
#' @param match.info match info table (each row describes an rf)
#' @param cluster object (list) with cluster information, as provided by tina
#' 
#' @return updated and expanded match.info object. * Note: weighted.rmse and peak-specific fields may be duplicated/missing, as they cannot be computed here (or don't matter going forward)
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' 
#' @export
propagate_matches <- function(match.info, cluster, feature.stack, ref.mat, ncores, r.thresh, p.thresh){
        
        message('\n Propagating matches to cluster members...')
        matched.feats <- match.info$feat %>% unique
        n.matches.before <- nrow(match.info)
        
        # Prep for basic load-balancing based on cluster size ####
          cluster.size <- lapply(matched.feats, function(x) {
            clust.number <- which(cluster$keys %in% x)
            cluster$groups[[clust.number]] %>% length
          }) %>% unlist
        
          feats.by.size <- matched.feats[order(cluster.size)]
        
        # Compute new match.info for cluster members ####
          t1 <- Sys.time()
          new.data <- mclapply(feats.by.size, function(fstack.row) {
            # print(fstack.row)
            # fstack.row <- feats.by.size[2]
            # Which cluster does this feature it belong to? ####
              # Behind each row of fstack is 1 or more feature indices
                
                # For now, just assume these are the row numbers in the original featureStack (with all cluster members)
  
                  # fstack.row <- matched.feats[[1]]
                  clust.number <- which(cluster$keys %in% fstack.row) # just an index, not the key feature
            
              
            # Pull cluster info ####
            
              clust.key <- cluster$keys[clust.number] # actual name of the key - should be fstack.row
              cluster.members <- cluster$groups[[clust.number]]
              clust.info <- cluster$info[[clust.number]]
              
              # Put lag table in terms of matching to key ####
              
                clust.lags <- clust.info$lag.table %>% pw_lags_relative_to(clust.key)
              
                # lag_features(feature$stack, clust.lags, to = clust.key) %>% trim_sides %>% stackplot #%>% simplePlot(linecolor = 'black')
  
            # Match info for key feature
              clust.matches <- match.info[match.info$feat == clust.key, ]
              
  
              ############ 
              ############             
              ############             
              ############ # For each non-key cluster member: ############
              ############ 
              ############ 
              
              nonkey.members <- cluster.members[cluster.members != clust.key]
              
              if (length(nonkey.members) > 0){
                
                # Copy match info for each cluster member (recopy first one), and recalculate:
                # - match info (adjust position)
                # - fit to ref 
                
                member.matches <- lapply(nonkey.members, function(cluster.member){
                  
                  # print(cluster.member)
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
                        new.matches$ref.start <- new.matches$ref.start - lag
                        new.matches$ref.end <- new.matches$ref.end - lag
                        # new.matches[, c("ref.start", "ref.end")] <- new.matches[, c("ref.start", "ref.end")] - lag
                    
     # ####                   
     ############ # For each ref feat for that cluster member: ####################################################
     
                    rfs.new <- lapply(1:nrow(new.matches), function(nm.row){
                    # Extract data for each new match ####
                      rf <- new.matches[nm.row,]
                      # rf <- new.matches[1,]
                      ref.num <- rf$ref
                      ref.reg <- rf$ref.start:rf$ref.end
                      feat.reg <- rf$feat.start:rf$feat.end
                      
                    # Get the lags and fits for the new matches to rfs ####
                    
                      # Extract out the feature
                      
                        feat <- feature.stack[cluster.member, ] %>% scale_between
                        
                      # Update the feature bounds to match THIS feature, not the key feature ####
                        # First, expand to use the whole feature (so we know how to ss the ref region):
                          ref.feat <- ref.inds <- rep(NA, length(feat))
                            ref.inds[feat.reg] <- ref.reg
                            ref.inds <- complete_indsVect(ref.inds)
                            
                        # Then, re-subset the feature and ref inds to match new feature:
                          f.inds <- feat %>% trim_sides(out = "inds") %>% range
                          # rf[c("feat.start", "feat.end")] <- f.inds
                          # rf[c("ref.start", "ref.end")] <- ref.inds[f.inds]
                          rf$feat.start <- f.inds[1]
                          rf$feat.end <- f.inds[2]
                          rf$ref.start <- ref.inds[f.inds[1]]
                          rf$ref.end <- ref.inds[f.inds[2]]
    
                        # Finally, update our temp vars for these two:
                          feat.reg <- rf$feat.start:rf$feat.end
                          ref.reg <- rf$ref.start:rf$ref.end
                          
                      # Make ref feat same size as full feature ####
                        # fill in the values
                          ref.feat[feat.reg] <- ref.mat[ref.reg, ref.num] # remember that ref mat is transposed
    
                      # Calculate the fit between the initial rf and the cluster.member ####
                        
                        fit <- fit_leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
                          # fit1$plot %>% plot
                        
                      # # Try to optimize the alignment a bit more (adds time!) ####
                      # 
                      #   small.adj <- rbind(feat, ref.feat) %>% feat_align(max.hits = 1, counter = F)
                      #     feat.rf.new.al <- rbind(feat, ref.feat) %>% lag_features(small.adj, to = 2) # "to" doesn't matter for just 2; it's bidirectional 
                      #                                                                                 # and lag.features will sort it out
                      #     # fit2 <- fit_leastSquares(feat.rf.new.al[1,], feat.rf.new.al[2,], plots = T, scale.v2 = T)
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
                      #     fit2<- fit_leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
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
                          
                        # Skip alignment, just update the row with fit data ####
                          
                          # fit <- fit1
  
                        # Either way, we now have a fit. Propagate that information:
                        
                          rf <- updateMatchInfoRow(rf, fit)
                            # rbind(fit1$feat.fit, fit1$spec.fit) %>% simplePlot(linecolor = 'black')
                            
                            # clust.matches[1,] # match this is copied from: initial feature 351 to ref 5
                            # ff <- apply_fit(mi.row = clust.matches[1,], feat.stack = feature$stack, ref.stack = ref.mat)
                            # ff <- apply_fit(mi.row = rf, feat.stack = feature$stack, ref.stack = ref.mat)
                            # rbind(ff$feat.fit, ff$spec.fit) %>% simplePlot(linecolor = 'black')
                          
                        # Return the data in a list
                        
                          return(rf)
                          
                    })
                        
                  return(rfs.new)                
    
              })
                
                return(member.matches)
                  # distribution of fit info can be interesting
                  # plot(member.matches$fit.intercept, member.matches$fit.scale)
  
              } else {
                # If single feature, don't expand.
                return(NULL)
                
              }
              
          }, mc.cores = ncores)
          
          a <- new.data %>% unlist(recursive = F) %>% unlist(recursive = F) %>% rbindlist
            rm(new.data)
            
          match.info <- rbind(match.info, a)

            row.names(match.info) <- NULL
            
          added.feats <- match.info$rmse.weighted %>% is.na %>% match.info$feat[.] %>% unique %>% length
          message('\n\t', n.matches.before, ' matches were propagated to ', nrow(match.info), ' matches:')
                  print(Sys.time()-t1)
                  message('\n\tfeatures before: ', length(matched.feats))
                  message('\n\tfeatures added : ', added.feats)
                  message('\n\tfeatures after : ', length(unique(match.info$feat)))
          
        # Re-filter for corr, pval
          message('\n\tfiltering new matches for rval > ', r.thresh, ' and pval < ', p.thresh, ' ...')
          keep <- match.info$rval >= r.thresh & 
            match.info$pval <= p.thresh
          
          match.info <- match.info[keep, ]
          # scattermore::scattermoreplot(x = 1:nrow(match.info), y = match.info$rval %>% sort)
          message('\n\t', sum(!keep), ' matches excluded by rval/pval filter (',  round(sum(!keep)/length(keep)*100), ' %)')
  return(match.info)
}