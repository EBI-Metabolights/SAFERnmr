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
#' Note: v3 of this uses batman_local_opt() to locally optimize the fit position and coefficients. Now doing away
#' with the bff score.
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
backfit_rfs3 <- function(match.info, 
                         feature.c, 
                         xmat,
                         refmat.c,
                         ncores){
  
  success = TRUE
  emptyRow <- function(){
    data.frame(ss.spec = NA,
              fit.intercept = NA, 
              fit.scale = NA,
              spec.start = NA,
              spec.end = NA,
              fit.fsa = NA,
              fit.rval = NA
    )
  }
  
  log_error <- function(message, chunk.num, m, error_cond) {
  
  
  cat(sprintf("%s - chunk$number: %d', chunk$match.info row: %d, Error: %s\n",
              Sys.time(), chunk.num, m, conditionMessage(error_cond)),
      file = "error_log.txt", append = TRUE)
  }
  
  # Chunk the data by ref (more or less) ####
    message('\tchunking match.info table, features objects, and ref spectra for distribution to cores...')
    # Sort by ref, so we can group chunks by ref
      mi.order <- order(match.info$ref)
      match.info <- match.info[mi.order, ]
    
    # Assign each row of match.ind to a chunk based on number of cores ####
      
      chunk.size <- max(1, nrow(match.info) / ncores)
      m.grp <- ceiling((1:nrow(match.info)) / chunk.size)
      refsums <- Rfast::rowsums(refmat.c %>% cstack_expandRows, na.rm = TRUE) # needed for normalizing ref features by total spectral signal
    
    # Assign feature, ref.mat, and match.info subsets to each chunk ####
      chunks <- lapply(unique(m.grp), function(x) 
      {
        # message('building chunk ', x)
        mi <- match.info[m.grp == x, ]
        ref.numbers <- unique(mi$ref)
        feat.numbers <- unique(mi$feat)
        
        # Index these to this chunk's refs and features
          # the mi data will be used to index in the ref.mat and feat.stack matrices
          mi$feat <- lapply(mi$feat, function(x) which(feat.numbers == x)) %>% unlist
          mi$ref <- lapply(mi$ref, function(x) which(ref.numbers == x)) %>% unlist
          
          fcss <- feature.c %>% select_features(feat.numbers)
          
          
        list(match.info = mi,
             ref.numbers = ref.numbers,
             refs = refmat.c %>% cstack_selectRows(ref.numbers),
             ref.sums = refsums[ref.numbers],
             feat.numbers = feat.numbers,
             feature = list(stack = fcss$stack, position = fcss$position),
             sfe = feature.c$sfe[feat.numbers],
             number = x)
             # Could keep feature object together and just subset, but for now just split sfe and stack info
             # feature = list(stack = feature$stack[feat.numbers, , drop=F],
             #                position = feature$position[feat.numbers, , drop=F],
             #                sfe = feature$sfe[feat.numbers]
                            # )
        # )
      })
      
    # Clean up big objects ####
      # rm(match.info)
      # rm(feature.c)
      # rm(refmat.c)
      # 
      # gc()
      
  # Compute over each chunk in parallel ####
    
    t1 <- Sys.time()
      
      # Each core gets:
      # - 1 chunk (~1Mb / 50 rows @ 150 spectra and 130K points)
      # - 1 xmat
      # - 1 ppm vector        
      message('\tcomputing backfits over ', ncores, ' cores...')
      message('\tlarge numbers of matches or large datasets will take some time.')
      # message('\tgo eat or get a coffee...\n\n')
  
      backfits.by.chunk <- mclapply(chunks, function(chunk) {
      ############# For each chunk (in parallel): ###############
        # chunk <- chunks[[1]]
        
        backfits.chunk <- lapply(1:nrow(chunk$match.info),
                           function(m) {
          tryCatch({
            
                      if (m==1){
                        # print('triggered on s =',s,' in row ',m)
                        stop("Simulated error")
                      }
          #############        For each match:         ###############
              # m <- 1
              
            # Get data, expand fit ####
              
              mi <- chunk$match.info[m, ]
              
              # Rebuild the feature fit to ref region
                
                # Get the feature and ref out
                  
                  lil.feat <- chunk$feature %>% expand_features(row.nums = chunk$feat.numbers[mi$feat])
                  lil.ref <- chunk$refs %>% cstack_selectRows(chunk$ref.numbers[mi$ref]) %>% cstack_expandRows
              
                  feat.cols <- lil.feat$position[mi$feat.start:mi$feat.end]

                # Apply feature fit to ref region (just using this to extract profiles, really)
                # the ref is the thing getting fit to, but it is scaled
                
                  fit <- apply_fit2(mi, v1 = lil.feat$stack, v2 = lil.ref) #%>% plot_fit(fit)
                    
                  # Calculate % of ref signature covered by this ref feat:
                    ref.reg <- lil.ref[mi$ref.start:mi$ref.end]
                      ref.reg[is.na(fit$feat.fit)] <- NA
                    pct.ref = sum(ref.reg, na.rm=T) / chunk$ref.sums[mi$ref]
                  
              # get info for the individual spec-features
              
                sfs <- data.frame(lag = chunk$sfe[[mi$feat]]$lags,
                                  ss.spec = chunk$sfe[[mi$feat]]$feat$ss)
                
                sfs$spec.start <- lil.feat$position[mi$feat.start] + sfs$lag
                sfs$spec.end <- lil.feat$position[mi$feat.end] + sfs$lag
                sfs$feat.start <- mi$feat.start
                sfs$feat.end <- mi$feat.end
                sfs$ref.start <- mi$ref.start
                sfs$ref.end <- mi$ref.end
                      
            # Calculate rf fit for each spec-feature: ####
              
              fits <- lapply(1:nrow(sfs), function(s){
                tryCatch(
                  {
                    # Get the ref region and spec data: ####
                      # s <- 21

                      sf <- sfs[s,]
                          
                    # Back-fit the ref region to the filled spec data using the feature ####
                    
                          # Propagate fit ####
                          
                            # optimize for single spectrum, single feature, single refspec
                              
                              # sf <- match.info[x, ] %>% [propagated to ss.specs]
                              # alternatively, can just use the sf$lag
                              # but opt_... perfects the lag
                              
                              res <- opt_specFit_slim(sf, 
                                                      feat.model = lil.feat$position, 
                                                      refspec = lil.ref, 
                                                      spec = xmat[ sf$ss.spec, ])
                              
                              fit <- res$fit
                              
                              
                              # ff <- lil.ref[sf$ref.start:sf$ref.end] * fit$ratio + fit$intercept
                              # fs <- xmat[ sf$ss.spec, res$sf$spec.start:res$sf$spec.end]
                              # 
                              # plot_fit(list(feat.fit = ff,
                              #               spec.fit = fs), type = 'auc')
                              # plot_fit(list(feat.fit = res$feat,
                              #               spec.fit = res$spec), type = 'auc')
                              
                        # Return results ####

                            return(data.frame(ss.spec = sf$ss.spec,
                                              fit.intercept = fit$intercept, 
                                              fit.scale = fit$ratio,
                                              spec.start = res$sf$spec.start, # was changed by opt_specFit_slim()
                                              spec.end = res$sf$spec.end, # but ONLY if using the actual result!
                                              fit.fsa = fit$fraction.spec.accounted,
                                              fit.rval = fit$rval
                                              )
                                   ) 
                  }, 
                  warning = function(cond){
                    ## Fit-level warning
                    # stop('backfit_rfs3: warning in
                    #         chunk$match.info row ', m,', ',
                    #         'sfs row ', s)
                    emptyRow()
                  }, 
                  error = function(cond){
                    
                    # stop('backfit_rfs3: error in
                    #         chunk$match.info row ', m,', ',
                    #         'sfs row ', s
                    #         )
                    emptyRow()
                  }
                )
                
              }) %>% do.call(rbind,.)
              
              # fits <- fits[!is.na(rowSums(fits)),] # get rid of NA rows # potential failure point; all NA rows removed means empty df
              
              fits$pct.ref <- pct.ref
              
            # Return the fits dataframe rows (minimal data)  ####
              return(fits)
              
          }, warning = function(w){
            ## Match-level warning
            # 
            # stop('backfit_rfs3: error in second layer loop, iteration: chunk$match.info row ', m)
            # If the whole match fails, return an NA fits df row
            fits <- emptyRow()
            fits$pct.ref <- NA
            log_error("backfit_rfs3 warning: ", chunk$number, m, w)
      
            return(fits)
            
          }, error = function(e){
            
            # stop('backfit_rfs3: error in second layer loop, iteration: chunk$match.info row ', m)
            # If the whole match fails, return an NA fits df row
            fits <- emptyRow()
            fits$pct.ref <- NA
            log_error("backfit_rfs3 error: ", chunk$number, m, e)
            
            return(fits)
            
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
      # chunk <- chunks[[1]]
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

  # Remove NA rows and validate data frames
    # backfits
    browser()
    message('\tcleaning backfits data...\n')
      backfits_valid <- function(backfits) {
        # backfits <- list(data.frame(a = c(1, NA, 3), b = c(NA, 2, 3)),
        #          data.frame(a = c(4, 5, 6), b = c(7, 8, 9)),
        #          data.frame(a = c(NA,NA,NA), b = c(NA,NA,NA)),
        #          c(10, 20, 30))

        # Filter to ensure all elements are data frames
        valid_dfs <- purrr::keep(backfits, is.data.frame)
        
        # Remove rows with NAs from each data frame
        cleaned_dfs <- purrr::map(valid_dfs, ~ dplyr::filter_all(.x, dplyr::all_vars(!is.na(.))))
     
        # Check each data frame to make sure it still has data
        has.rows <- sapply(cleaned_dfs, function(x) nrow(x) > 0) %>% which
        
        return(list(cleaned = cleaned_dfs,
                    inds.good = has.rows))
      }
      
      bfs <- backfits_valid(backfits)
      backfits <- bfs$cleaned[bfs$inds.good]
    if (length(backfits)==0){
      
      warning('backfit_rfs3: All backfits failed validation.')
    }
      
    # matches
      match.info <- match.info[bfs$inds.good,]      
      
    
    message('\tbackfitting completed on ', length(backfits), ' ref-features.')
    # message('\thopefully your lunch was nice and you are now properly caffeinated...\n')
    if(length(backfits) != nrow(match.info)){
      success = FALSE
      no.bfs <- match.info$id[-bfs$inds.good]
      warning('not all matches had backfits, including IDs:\n', paste(no.bfs,collapse = ' '))
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
          #              feat = fit$feat,
          #              ss.spec = lapply(bf$fits, function(f) f$ss.spec) %>% unlist,
          #              pct.ref = pct.ref,
          #              + fit scores)
  
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
   
   
    
     

