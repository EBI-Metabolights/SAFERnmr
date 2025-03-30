#' log_storm_core: 
#' Locally Optimized Global STORM
#' 
#' Run modified STORM on the provided spectral region and ref shape.
#' Built for accepting corrPocketPairs results. Notes:
#'
#' STORM: Joram Posma's STORM has been adapted and optimized to accept these
#   protofeatures (corrPocketPairs) in the following ways:
# - first, since many of the protofeatures are noise, we provide failure modes
#   and reporting for the following cases:
#   "empty subset",          # empty subset (no spectrum contains signature)
#   "subset degenerated",    # 1-3 spectra in the subset (not enough spectra to
#                              get a reliable correlation)
#   "reference degenerated", # signature degenerates to include < 3 points (not
#                              meaningful to correlate shapes)
#   "did not converge"       # subset continues to change after 24 iterations
#
# - additionally, the correlation r and p-value cutoff q are both used during
#   both the subset selection and reference update steps.
# - we also remove any regions of the reference for which there are fewer than
#   minpeak values after r and p value thresholding. This helps avoid noise.
#
# STORM extracts meaningful features using protofeatures to define the region of
# interest and a rough sketch of the feature shape highly correlated with each
# spectral point. In the future, HCA could be used to cluster potential starting
# feature shapes correlated with each driver, or the nonoptimal subset for each
# point could be re-STORMed to detect any other feature shapes present. It is
# also perfectly reasonable to combine feature shapes from different STORM runs
# for a given dataset, as these comprise a list of somewhat independently tested
# feature shapes, and duplication is not an issue.

#'
#'
#' @param xmat A matrix of spectral data (rows are spectra, columns are spectral points)
#' @param ppm A vector of the spectral points in ppm (optional, default is all columns of xmat)
#' @param b An integer giving the expansion parameter for the reference peak
#' @param corrthresh A numeric giving the minimum correlation value to be considered for inclusion (for both subset AND reference optimization)
#' @param q A numeric giving the p-value threshold for correlation significance (both subset AND reference optimization)
#' @param minpeak An integer giving the minimum number of points allowed in a run of significant points in the reference
#' @param refSpec A vector of spectral data to use as the initial reference
#' @param ref.idx A vector of the spectral points (columns of xmat) to use as the initial reference
#' @param range.limit A vector of the spectral points (columns of xmat) to use as the initial reference
#'
#' @return A list with components "reconstructed" and "status". "reconstructed" is a matrix
#' containing the reconstructed metabolite concentrations (rows are samples, columns are metabolites).
#' "status" is a character string indicating whether the method converged successfully or failed.
#'
#' @export log_storm_core
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_path geom_line geom_vline geom_hline ggtitle xlab ylab scale_y_continuous scale_x_continuous
#' @importFrom stringr str_pad
log_storm=function(xmat=NULL, ppm=NULL, b=30, corrthresh = .8,
                        q=0.05, minpeak = 10, refSpec=NULL, ref.idx=NULL,
                        driver = NULL, range.limit=400){

############ Setup ##################################################  
    
    only.region.between <- pars$corrpockets$only.region.between
    if (is.null(only.region.between))                       # which ppms to run fse between
      {only.region.between <- range(ppm)}                   #   (default is all)
    correlation.r.cutoff <- pars$storm$correlation.r.cutoff   # rvalue cutoff for both subset selection (ref shape) and ref update (STOCSY)
    q <- pars$storm$q                                         # q param from storm (pval cutoff after mhtc)
    b <- pars$storm$b                                         # number of peak widths to expand ref by on each side

############ Initialize for the loop ########################################

    bounds <- vectInds(only.region.between, ppm)
    testregion <- bounds[2]:bounds[1]
    
    cc.col <- cc.peaks[[1]]
    
    # Derive the inds for the storm seeds
      protofeatures <- unlist_ccpeaks(cc.peaks)
   
      unlist_ccpeaks <- function(cc.peaks){
        peak.list.all <- lapply(cc.peaks, function(cc.col){
          
          peak.list.col <- lapply(cc.col$secondary, function(secondary.peak){
          # These are currently window indices, but will be converted to "reg" indices
          # (indices of ppm vector/xmat provided to correlation_pocket_pairs())
  
            list(driver.index = cc.col$index,
                 primary.peak = cc.col$index - (cc.col$res.center-cc.col$primary%>%unlist),
                 secondary.peak = cc.col$index - (cc.col$res.center-secondary.peak%>%unlist))
          }) 
        }) %>% unlist(recursive = FALSE)
        
      }
             
      
    # Batch the mclapply so each core gets only one copy of xmat. 
    
    # Feat data is split between nodes:
        message('\nSplitting feature matrices for distribution across nodes...')
        pars$par$ncores <- 10
        chunk.size <- max(1, length(protofeatures) / pars$par$ncores)
        f.grp <- ceiling((1:nfeats) / chunk.size)
        split.scheme <- lapply(unique(f.grp), function(g) {
            list(
                f.inds = which(f.grp == g),
                f.subset = which(f.grp == g) %>% f.subset[.] 
            )
        })
          split.scheme %>% test_nullish
        f.stack.split <- lapply(unique(f.grp), function(x) f.stack[, f.grp == x, drop = F])
          f.stack.split %>% test_nullish
          rm(f.stack)
      
    
    storm_rnd1 = list()
    message("Running storm on ",length(regions_subset), " provided protofeatures between ",ppm[bounds[1]]," and ",ppm[bounds[2]]," ppm.")
   
#################################################
        storm_rnd1 <- 
              mclapply(protofeatures,
              # pblapply(regions_subset,
                  function (x) {
                    # TryCatch will 
                    tryCatch(
                      expr = {
                              # Set up the region
                              
                              # We'll limit the STORM operation to a small region of the xmat. No need to allow gratuitous expansion.
                              
                              # We'll also derive the covariance signal for the protofeature individually. A function to spin up the 
                              # protofeature profile given xmat would be helpful. 
                              
                                protofeature_profile <- function(protofeature, xmat){
                                  
                                  # return: NA-filled covariance intensity values for the protofeature
                                }
                                
                        
                                
                                # x <- regions_subset[1793]
                                
                                driver <- colwithPair[x]
                                peakPos <- pocketPairs$peakMap[,driver] %>% is.na %>% "!"(.) %>% which
                                pair.region <- pocketPairs$regions[peakPos,driver]
          
                              # Do storm_pairplay
                    
                                # Set params
                                  pw <- span(peakPos)/2 %>% ceiling
                                  wind <- pair.region
                                  shift <- range(ppm[pair.region])
                                  
                                # Use original covariance signal within corr bounds as shape seed
                                # (could also use best spectrum index)
                                  shape <- pocketPairs$cov[peakPos,driver]
                                  bestSpec = cor( xmat[ ,wind] %>% t, shape ) %>% which.max
                                  
                                # Do the storm
                                  
                                  res <- storm_pairplay(xmat, ppm,
                                                              b = (pw * b) %>% ceiling, corrthresh = correlation.r.cutoff, q = q,
                                                              minpeak = noisewidth, refSpec = shape, ref.idx = pair.region,
                                                              driver = driver)
                                  res$cpp.driver <- driver
                                  
                                  return(res)
                        
                      },
                      error = function(cond){
                        return('setup error')
                        }
                      )
                  }, mc.cores = pars$par$ncores
              )
        # Note: errors in the loop are captured and passed out as strings.
        # NULL elements are not possible, although parts of an element could be.
        # Those are checked below.  
          
################ Report run stats  ######
        
          fmodes <- lapply(storm_rnd1, 
                           function(x) {
                             if (is.character(x)){return(x)} # this will get any errors from setup or storm
                             if (x$status == 'succeeded'){
                               if (any(is_nullish(x))){
                                 # Even if STORM succeeded, it may contain NULLs in some return value elements. 
                                 return(  paste0(  is_nullish(x) %>% which %>% names, " contains NULL")) 
                               }
                             }
                             return(x$status)
                           })
          
          succeeded <- lapply(fmodes, function(x) x == 'succeeded') %>% unlist
          failed <- !succeeded
          
          message(str_c("Failed iterations (count): ", sum(failed), " (",
                        (sum(failed)/length(regions_subset) * 100) %>% round, " %)"))
          
          message(str_c("Succeeded iterations (count): ", sum(succeeded), " (",
                        (sum(succeeded)/length(regions_subset) * 100) %>% round, " %)"))
          
          # Print out breakdown of statuses
            fm <- fmodes %>% plyr::ldply(rbind)
            colnames(fm) <- 'status'
            fm <- fm %>% group_by(status) %>% count %>% as.data.frame
            
            print(fm)

################################################
}
