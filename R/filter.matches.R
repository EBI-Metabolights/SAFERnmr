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
#' @importFrom data.table rbindlist
#' 
#' @import pbapply
#' 
#' @export
filter.matches <- function(pars){

  message('--------------------------------------------------------------')
  message('-------------------     Match Filtering    -------------------')
  message('--------------------------------------------------------------')
  printTime()
  message('\n\n\n')
  
################ Read parameters file ##################
  
  
  tmpdir <- pars$dirs$temp
  this.run <- paste0(tmpdir)

##################################################################################################################
 # Read data and set up ####

    message("Loading data from files...\n\n\n")
    
    fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      rm(fse.result)

    feature <- readRDS(paste0(this.run, "/feature.final.RDS"))
    ref.mat <- readRDS(paste0(this.run, "/temp_data_matching/ref.mat.RDS"))
    # matches <- readRDS(paste0(this.run, "/matches.initial.RDS"))
    matches <- readRDS(paste0(this.run, "/matches.RDS"))
      errors <- matches[names(matches) %in% c('call', 'message')]
      matches <- matches[names(matches) %in% c('matches', 'peak.quality')]
      
###########################################################################################  
     # No full feature fits beyond this point. Just storing coefficients and positions. Try
     # to move this line up further. 
###########################################################################################        
      
    cluster <- readRDS(paste0(this.run, "/cluster.final.RDS"))

    # Format matches ####
      
      
      matches.split <- split(matches, names(matches))
        rm(matches)
        
      match.info <- rbindlist(matches.split$matches)
        rownames(match.info) <- NULL
        # saveRDS(match.info, paste0(this.run, "/match.info.RDS"))
        if (!is.data.frame(match.info)){
          stop('match info is not a dataframe.')
        } else {
          message('\n\n\tmatch.info is a dataframe with ', nrow(match.info), ' rows...\n\n')
        }
        
        
      # peak.qualities aligns with match.info now, but feat number does not index it (there are missing features)!
        pq.featureNumbers <- unique(match.info$feat) # this does not sort (just for good measure)
        
      peak.qualities <- matches.split$peak.quality
        rm(matches.split)
          
######################### Remove singlets ############################################
    printTime()
      # Do the filtering (functionalized)
        match.info <- filter.matches_singlets(match.info, feature$stack, ref.mat,
                                               peak.qualities, pq.featureNumbers, 
                                               pars$matching$filtering$res.area.threshold,
                                               pars$par$ncores)
        
        saveRDS(match.info, paste0(tmpdir, "/match.info.filtered.RDS"))
    
    
######################### Propagate matches to feature clusters  ##########################

        # For each match, produce other matches based on clusters
        # - loop through clusters
        # - copy all matches from key to other cluster members
        # - what gets copied?
        #   - match.info row = ref.feature
        #   - all that's needed is to test the ref feature in the local regions 
        #     x spectra from which the cluster member was extracted, as gained 
        #     from the sfe
        #   - fstack.row is the row of the match matrix (subset of feature matrix)
        # - what gets added?
        #   - new feature, lagged and fit to ref
        #   - just store the coefficients and positions
        #   
        # This boils down to a match.info row for every feature. Fits can be rebuilt
        # on the fly from this and the feature/spectral matrix/ref matrix.
        
        # Only run if there are actually clusters!
          
          if (cluster$method == 'none'){
            # do nothing
          } else {
            match.info <- propagate.matches(match.info, cluster, feature$stack, 
                                            ref.mat, pars$par$ncores, pars$matching$r.thresh, pars$matching$p.thresh)
            saveRDS(match.info, paste0(tmpdir, "/match.info.propagated.RDS"))
          }
        
        
        # match.info <- readRDS(paste0(tmpdir, "/match.info.propagated.RDS"))
        
######################### Calculate deltappm distance (specppm - featureppm)  #############################

        # source('./../span.R')
        # source('./../filter.matches_shiftDelta.R')
  printTime()
        message('\nFiltering out matches > ', pars$matching$filtering$ppm.tol, ' ppm away...')
        res <- filter.matches_shiftDelta(match.info, feature$position, ppm = ppm,
                                         ppm.tol = pars$matching$filtering$ppm.tol)
        # scattermore::scattermoreplot(x = 1:nrow(res), y = res$ppm.difference %>% sort)

        message('\n\t', nrow(res),' / ', nrow(match.info), ' matches survived (', round(nrow(res)/nrow(match.info)*100), ' %).')

        match.info <- res
        
        # Savepoint
          saveRDS(match.info, paste0(this.run, "/match.info.propagated.filtered.RDS"))
          # match.info <- readRDS(paste0(this.run, "/match.info.propagated.filtered.RDS"))

#########################################################################################################
    # At this point, match.info is set. Assign IDs
                          
                            match.info$id <- 1:nrow(match.info)

######################### Back-fit reference to spectra  #############################    
    
      message('\n\nBack-fitting ref-feats to each spectrum in the relevant subset...\n')
    printTime()
    
    # Back-fit each matched reference region to the subset spectra
      # adjusted to account for sfe

        backfit.results <- backfit.rfs(match.info, 
                                       feature, # has sfe data 
                                       xmat,
                                       ref.mat, 
                                       pars$par$ncores)
        
        message('Saving backfits...\n\n\n')
        saveRDS(backfit.results, paste0(this.run,"/backfit.results.RDS"))
        # backfit.results <- readRDS(paste0(this.run, "/backfit.results.RDS"))

      printTime()
      
  message('-----------------------------------------------------------------')
  message('-----------------  Matching Filtering Complete ------------------')
  message('-----------------------------------------------------------------')
  message('\nWe just generated ', length(backfit.results$backfits %>% unlist(recursive = F)), ' pieces of annotation evidence.')
  message("\nNow, let's turn these into PCRS x sample scores.")
}
  