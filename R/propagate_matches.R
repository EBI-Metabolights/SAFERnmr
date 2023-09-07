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
#' @param feature.stack feature profiles on rows
#' @param ref.mat refs
#' @param ncores pars
#' @param r.thresh pars
#' @param p.thresh pars
#' @param pad.size from matching
#' @param this.run pars
#' 
#' @return updated and expanded match.info object. * Note: weighted.rmse and peak-specific fields may be duplicated/missing, as they cannot be computed here (or don't matter going forward)
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' 
#' @export
propagate_matches <- function(match.info, cluster, feature.stack.c, refmat.c, ncores, r.thresh, p.thresh, pad.size, this.run){
    
  
        # Set up default empty row:
          emptyRow <- match.info[1,]
          emptyRow[,1:ncol(emptyRow)] <- NA

        message('\n Propagating matches to cluster members...')
        matched.feats <- match.info$feat %>% unique
        n.matches.before <- nrow(match.info)
        
        # Prep for basic load-balancing based on cluster size ####
          cluster.size <- lapply(matched.feats, function(x) {
            clust.number <- which(cluster$keys %in% x)
            cluster$groups[[clust.number]] %>% length
          }) %>% unlist
          
          plan <- distribute(workloads = cluster.size, across = ncores)
          feats.by.size <- matched.feats[order(cluster.size)]
        
        # Partition data ####
          # browser()
          chunks <- lapply(unique(plan$core.ids), function(chunk.number) {
            matched.feats[plan$core.ids == chunk.number]
          })
          browser()
          
        # Compute new match.info for cluster members ####
          t1 <- Sys.time()
          new.data <- mclapply(chunks, function(chunk) {
          # new.data <- pblapply(chunks, function(chunk) {
              # chunk <- chunks[[1]]
              chunk.matches <- pblapply(chunk[1:10], function(fstack.row){
            
              
              # Which cluster does this feature it belong to? ####
                # Behind each row of fstack is 1 or more feature indices
                  
                  # For now, just assume these are the row numbers in the original featureStack (with all cluster members)
    
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
                ############ # For each non-key cluster member: ############
  
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
                        
       # ####                   
       ############ # For each ref feat for that cluster member: ####################################################
       
                      rfs.new <- lapply(1:nrow(new.matches), function(nm.row){
                      # Extract data for each new match ####
                        
                        rf <- new.matches[nm.row,]
  
                      # Get the lags and fits for the new matches to rfs ####
                      
                        # Extract out the feature
                          feat <- feature.stack.c %>% cstack_selectRows(cluster.member) %>% cstack_expandRows %>% scale_between
                          # feat <- feature.stack[cluster.member, ] %>% scale_between
                          # feat <- feature.stack[508, ] %>% scale_between
                            # simplePlot(feat)
                            
                        # Update the feature bounds to match THIS feature, not the key feature ####
                          # First, expand to use the whole feature (so we know how to ss the ref region):
                          # Make the ref feat (full length)
                            ref.inds <- rf$lag - pad.size + (0:(length(feat)-1))
                            ref.feat <- refmat.c %>% cstack_selectRows(rf$ref) %>% cstack_expandRows %>% .[ref.inds]
                            # ref.feat <- ref.mat[ref.inds, rf$ref]
                              # simplePlot(ref.feat)
                              
                          # Then, re-subset the feature and ref inds to match new feature:
                            f.inds <- feat %>% trim_sides(out = "inds") %>% range
                            rf$feat.start <- f.inds[1]
                            rf$feat.end <- f.inds[2]
                            rf$ref.start <- ref.inds[f.inds[1]]
                            rf$ref.end <- ref.inds[f.inds[2]]
      
                          # Finally, update our temp vars for these two:
                            feat.reg <- rf$feat.start:rf$feat.end
                            
                            
                          # Calculate the fit between the initial rf and the cluster.member ####
                            # there are a lot of ways this can fail, so for now, putting in a tryCatch:
                            # if (sum(is.na(feat[feat.reg] + ref.feat[feat.reg])) < 4){return(emptyRow)}
                            
                            rf <- tryCatch(
                                    {
                                      fit <- fit_leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
                                                      # fit$plot %>% plot
                                      # Recalculate match fields for new fit information: ####
                                      updateMatchInfoRow(rf, fit)
                           
                                    },
                                    error=function(cond) {
                                        return(emptyRow)
                                    },
                                    warning=function(cond) {
                                        return(emptyRow)
                                    }
                            )
                        
                            return(rf)
                            
                      })
                      
                    return(rfs.new)                
                    # this gives a list where each element holds an individual df row
                })
                  
                  if(is.null(member.matches)){
                    
                    return(emptyRow)
                    
                  }
                  
                  member.matches <- member.matches %>% unlist(recursive = F) %>% rbindlist
                  
                  return(member.matches)
                    # this gives a list of those lists
                    # distribution of fit info can be interesting
                    # plot(member.matches$fit.intercept, member.matches$fit.scale)
                    
                } else {
                  # If single feature, don't expand.
                  return(emptyRow)
                }
              }) %>% rbindlist
            }, mc.cores = ncores
          )
          
          saveRDS(new.data, paste0(this.run,'/new.data.RDS'))
          # new.data <- readRDS('/Users/mjudge/Documents/new.data.RDS')
          
          a <- new.data %>% rbindlist
            rm(new.data)
            
          # Rbind the new matches to the end of match.info
            match.info <- rbind(match.info, a)
            # remove any null rows (indicated by NA in the feat column)
            match.info <- match.info[!is.na(match.info$feat),]

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

