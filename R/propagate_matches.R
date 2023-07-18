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
propagate_matches <- function(match.info, cluster, feature.stack, ref.mat, ncores, r.thresh, p.thresh, pad.size, this.run){
    
  
        # Set up default empty row:
          emptyRow <- function(){
            data.frame( 
              feat = NA,
              ref = NA,
              lag = NA, 
              rval = NA,
              pval = NA,
              pts.matched = NA,
              pts.feat = NA,
              feat.start = NA,
              feat.end = NA,
              ref.start = NA,
              ref.end = NA,
              fit.intercept = NA,
              fit.scale = NA,
              wasserstein.score = NA,
              sum.residuals = NA, 
              rmse = NA, 
              rmse.weighted = NA,
              numpeaks.feat = NA, 
              numpeaks.ref = NA, 
              numpeaks.feat.nnf = NA, 
              refpeaks.matched = NA
            )
          }
  
        message('\n Propagating matches to cluster members...')
        matched.feats <- match.info$feat %>% unique
        n.matches.before <- nrow(match.info)
        
        # Prep for basic load-balancing based on cluster size ####
          cluster.size <- lapply(matched.feats, function(x) {
            clust.number <- which(cluster$keys %in% x)
            cluster$groups[[clust.number]] %>% length
          }) %>% unlist
        
        # Partition data ####
          
          matched.feats <- matched.feats[cluster.size > 1]
          
          chunks <- lapply(matched.feats, function(f) {
              # f <- feats.by.size[1]
              # Pull cluster info ####
                clust.number <- which(cluster$keys %in% f)
                clust.key <- cluster$keys[clust.number] # actual name of the key - should be fstack.row
                cluster.members <- cluster$groups[[clust.number]]
                clust.info <- cluster$info[[clust.number]]
    
                # Put lag table in terms of matching to key ####
    
                  clust.lags <- clust.info$lag.table %>% pw_lags_relative_to(clust.key)
    
              # Match info for key feature and cluster members
    
                clust.matches <- match.info[match.info$feat == clust.key, ]
                nonkey.members <- cluster.members[cluster.members != clust.key]
    
              # Feature profile for key feature and cluster members
                member.feats <- feature.stack[nonkey.members, ,drop=F]
    
                # balance by members x matches
                chunk <- list(features = member.feats,
                              matches = clust.matches,
                              clust.number = clust.number,
                              clust.key = clust.key,
                              clust.lags = clust.lags,
                              clust.members = cluster.members,
                              clust.info = clust.info,
                              nonkey.members = nonkey.members,
                              amount.of.work = length(nonkey.members) * nrow(clust.matches))
                return(chunk)
          })
    
          workloads <- lapply(chunks, function(x) x$amount.of.work) %>% unlist

          plan <- distribute(workloads, across = ncores)
          
          chonks <- lapply(1:ncores, function(core.number){
            chonk <- chunks[plan$core.ids == core.number]
            return(chonk)
          })
          rm(chunks)
          
          browser()
        # Compute new match.info for cluster members ####
          t1 <- Sys.time()
          new.data <- mclapply(chonks, function(chonk) {
            # chonk <- chonks[[1]]
            
            # For each chunk (cluster)
            chonk.matches <- lapply(chonk, function(chunk){
              
              # chunk <- chonk[[1]]

            # Which cluster does this feature it belong to? ####
              # Behind each row of fstack is 1 or more feature indices
                
                # For now, just assume these are the row numbers in the original featureStack (with all cluster members)
  
                clust.number <- chunk$clust.number
            
              
            # Pull cluster info ####
            
              clust.key <- chunk$clust.key
              cluster.members <- chunk$clust.members
              clust.info <- chunk$clust.info
              
              # Put lag table in terms of matching to key ####
              
                clust.lags <- chunk$clust.lags
              
                # lag_features(feature$stack, clust.lags, to = clust.key) %>% trim_sides %>% stackplot #%>% simplePlot(linecolor = 'black')
  
            # Match info for key feature
              clust.matches <- chunk$matches
              
              ############ # For each non-key cluster member: ############

              nonkey.members <- chunk$nonkey.members
              
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
                      
                    # Extract out the feature
                      feat <- chunk$features[nonkey.members %in% cluster.member, ] %>% scale_between  
                      f.inds <- feat %>% trim_sides(out = "inds") %>% range
                      # simplePlot(feat)                      
     # ####                   
     ############ # For each ref feat for that cluster member: ####################################################
     
                    rfs.new <- lapply(1:nrow(new.matches), function(nm.row){
                      # nm.row <- 1
                    # Extract data for each new match ####
                      
                      rf <- new.matches[nm.row,]

                    # Get the lags and fits for the new matches to rfs ####
                          
                      # Update the feature bounds to match THIS feature, not the key feature ####
                        # First, expand to use the whole feature (so we know how to ss the ref region):
                        # Make the ref feat (full length)
                          ref.inds <- rf$lag - pad.size + (0:(ncol(feature.stack)-1))
                          ref.feat <- ref.mat[ref.inds, rf$ref]
                            # simplePlot(ref.feat)
                            
                        # Then, re-subset the feature and ref inds to match new feature:
                          
                          rf$feat.start <- f.inds[1]
                          rf$feat.end <- f.inds[2]
                          rf$ref.start <- ref.inds[f.inds[1]]
                          rf$ref.end <- ref.inds[f.inds[2]]
    
                        # Finally, update our temp vars for these two:
                          feat.reg <- rf$feat.start:rf$feat.end
                          
                          
                      # Calculate the fit between the initial rf and the cluster.member ####
                          # there are a lot of ways this can fail, so for now, putting in a tryCatch:
                          # if (sum(is.na(feat[feat.reg] + ref.feat[feat.reg])) < 4){return(emptyRow())}
                          
                          rf <- tryCatch(
                                  {
                                    fit <- fit_leastSquares(feat[feat.reg], ref.feat[feat.reg], plots = F, scale.v2 = T)
                                                    # fit$plot %>% plot
                                    # Recalculate match fields for new fit information: ####
                                    updateMatchInfoRow(rf, fit)
                         
                                  },
                                  error=function(cond) {
                                      return(emptyRow())
                                  },
                                  warning=function(cond) {
                                      return(emptyRow())
                                  }
                          )
                      
                          return(rf)
                          
                    })
                    
                  return(rfs.new)                
                  # this gives a list where each element holds an individual df row
              })
                
                if(is.null(member.matches)){
                  
                  return(emptyRow())
                  
                }
                
                member.matches <- member.matches %>% unlist(recursive = F) %>% rbindlist
                
                return(member.matches)
                  # this gives a list of those lists
                  # distribution of fit info can be interesting
                  # plot(member.matches$fit.intercept, member.matches$fit.scale)
                  
              } else {
                # If single feature, don't expand.
                return(emptyRow())
              }
              
            }) %>% unlist(recursive = F)
            
            return(chonk.matches)
            
          }, mc.cores = ncores)
          
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

