#' Match Filtering
#'
#' Filter and process matched peaks based on user-specified criteria, such as singlet removal and ppm distance filtering.
#' Also calculate backfits of ref-feats (reference subsignatures fished out by feature matches to reference spectra) to each individual dataset spectrum. 
#' Note: backfits are now calculated directly from ref features and position in 
#' ss.spec is optimized using updated batman-style fitting. This strategy has proven
#' to be much more accurate and produces meaningful scores.
#' Note: BFF (backfit feasibility) scores are no longer used. Instead, updated batman fits
#' yield a fraction.spec.accounted score which quantifies the fraction of relevant spectral
#' signal accounted for by the ref feature. positive overshoot is penalized by 
#' subtracting from the score. This is much more simple and intuitive. 
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
filter_matches <- function(pars){

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

    feature.c <- readRDS(paste0(this.run, "/feature.final.RDS"))
    refmat.c <- readRDS(paste0(this.run, "/temp_data_matching/ref.mat.RDS"))
    pad.size <- readRDS(paste0(this.run, "/temp_data_matching/pad.size.RDS"))
    # matches <- readRDS(paste0(this.run, "/matches.initial.RDS"))
    
    # Weed out empty matches or those which failed
    matches <- readRDS(paste0(this.run, "/matches.RDS"))
      # if any invalid matches for the feature
        matches <- matches[!is_nullish(matches)]
      # per feature, was NA returned, or was ('matches', 'peak.quality')?
        nomatch <- (lapply(matches, length) %>% unlist) == 1 
        matches <- matches[!nomatch]
      
      # For each feature: 
        # Check if there was an error message
          errors <- matches[names(matches) %in% c('call', 'message')]
        # Check if both of the expected fields are present
          matches <- matches[names(matches) %in% c('matches', 'peak.quality')]

###########################################################################################        
      
    cluster <- readRDS(paste0(this.run, "/cluster.final.RDS"))

    # Format matches ####
      # At this point, matches have been thoroughly validated and can be rbinded.
        matches.split <- split(matches, names(matches))
        rm(matches)
        
      match.info <- rbindlist(matches.split$matches)
        rownames(match.info) <- NULL
        
        match.info %>% debug_write("match.info.initial.RDS", pars)
        # match.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/match.info.initial.RDS"))

      # peak.qualities aligns with match.info now, but feat number does not index it (there are missing features)!
        pq.featureNumbers <- unique(match.info$feat) # this does not sort (just for good measure)
        
      peak.qualities <- matches.split$peak.quality
        if (length(pq.featureNumbers) != length(peak.qualities)) { stop(' peak qualities not of same length as pq.featureNumbers!')}
        rm(matches.split)
        
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
          # don't use pars$tina$do.clustering here, because even if that's on, if there are < 1K features it won't cluster.
          # After this step, clusters are no longer relevant. 
          
          if (cluster$method == 'none'){
            # do nothing
          } else {
            match.info <- propagate_matches(match.info, cluster, feature.c$stack, 
                                            refmat.c, pars$par$ncores, pars$matching$r.thresh, 
                                            pars$matching$p.thresh, pad.size, this.run)
            
            match.info %>% debug_write("match.info.propagated.RDS", pars)
            # match.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/match.info.propagated.RDS"))

          }
        
######################### Calculate deltappm distance (specppm - featureppm) #############################

        # source('./../span.R')
        # source('./../filter.matches_shiftDelta.R')
  printTime()
        message('\nFiltering out matches > ', pars$matching$filtering$ppm.tol, ' ppm away...')
        res <- filter_matches_shiftDelta(match.info, feature.c %>% expand_features %>% .[['position']], 
                                         ppm = ppm, ppm.tol = pars$matching$filtering$ppm.tol)
        # scattermore::scattermoreplot(x = 1:nrow(res), y = res$ppm.difference %>% sort)

        message('\n\t', nrow(res),' / ', nrow(match.info), ' matches survived (', round(nrow(res)/nrow(match.info)*100), ' %).')

        match.info <- res
        
        match.info %>% debug_write("match.info.shiftFiltered.RDS", pars)
        # match.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/match.info.shiftFiltered.RDS"))

######################### Remove singlets ############################################
    printTime()
      # Do the filtering (functionalized)
          
        # For each match, how many peaks are contained in the ref and feature match regions
        # (not including NA gaps)
        message('Filtering out singlet matches (either ref or feature)...')
          match.info <- filter_matches_singlets2(match.info, feature.c, refmat.c, pars$par$ncores)
          
        # match.info <- filter_matches_singlets(match.info, feature.c$stack, refmat.c,
        #                                        pars$matching$filtering$res.area.threshold,
        #                                        pars$par$ncores)
        
        match.info %>% debug_write("match.info.filtered.RDS", pars)
        # match.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/match.info.filtered.RDS"))

#########################################################################################################
    # At this point, match.info is set. Assign IDs
                          
       match.info$id <- 1:nrow(match.info)

######################### Back-fit reference to spectra  #############################    
    
      message('\n\nBack-fitting ref-feats to each spectrum in the relevant subset...\n')
    printTime()
    
    # Back-fit each matched reference region to the subset spectra
      # adjusted to account for sfe
        # match.info <- match.info[1:10000,]
        backfit.results <- backfit_rfs3(match.info = match.info,
                                        feature.c = feature.c, # has sfe data 
                                        xmat = xmat,
                                        refmat.c = refmat.c, 
                                        ncores = pars$par$ncores)
        
        message('Saving backfits...\n\n\n')
        saveRDS(backfit.results, paste0(this.run,"/smrf.RDS"))
        # backfit.results <- readRDS(paste0(this.run, "/smrf.RDS"))
        unlink(paste0(tmpdir, "/matches.RDS")) # 

      printTime()
      
  message('-----------------------------------------------------------------')
  message('-----------------  Matching Filtering Complete ------------------')
  message('-----------------------------------------------------------------')
  message('\nWe just generated ', lapply(backfit.results$backfits, nrow) %>% unlist %>% sum, ' pieces of annotation evidence.')
  message("\nNow, let's turn these into PCRS x sample scores.")
}
  