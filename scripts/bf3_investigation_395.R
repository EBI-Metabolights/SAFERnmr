# Examine the cause of backfit failure for 395 (1710856857)
# 

devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1710856857'

# Make sure the tmpdir is actually local
  pars <- yaml::yaml.load_file(paste0(tmpdir,'/params.yaml'))
    pars$dirs$temp <- tmpdir
    
# Simulate backfit function call

  filter_matches(pars)
  
  feature.c <- readRDS(paste0(tmpdir, "/feature.final.RDS"))
  fse.result <- readRDS(paste0(tmpdir, "/fse.result.RDS"))
    xmat <- fse.result$xmat
    ppm <- fse.result$ppm

  # *** these should be in temp data matching..
    refmat.c <- readRDS(paste0(tmpdir, "/ref.mat.RDS"))
    pad.size <- readRDS(paste0(tmpdir, "/pad.size.RDS"))
  
# Matches ####
###########################################################################################        
     
  # Weed out empty matches or those which failed
    matches <- readRDS(paste0(tmpdir, "/matches.RDS"))
      
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
          match.info <- filter_matches_singlets2(match.info, feature.c, refmat.c, 2)
          
        match.info %>% debug_write("match.info.filtered.RDS", pars)
        filtered.matches <- nrow(match.info)
        # match.info <- readRDS(paste0(pars$dirs$temp, "/debug_extra.outputs", "/match.info.filtered.RDS"))

  # match.info %>% debug_write("match.info.beforeBackfitting.RDS", pars) # this should go into filter_matches
#########################################################################################################
    # Estimate the number of backfits, and jettison matches to keep below computational limit
    
      # This requires that the sfe field contains ss information for each feature, and that
      # the number of features in the features object actually matches the indices listed
      # in match.info. 
    
      ss.lengths <- feature.c$sfe %>% lapply(function(x){
        x$feat$ss %>% length
      }) %>% unlist
      contribution <- ss.lengths[match.info$feat]
      est.bfs <- contribution %>% sum
      message('Estimated backfits: ', est.bfs, ' at match rval cutoff of ', pars$matching$r.thresh, '.')

      if (est.bfs > pars$matching$filtering$max.backfits){
        # limit the number of matches - only take the top n such that they don't exceed max number
        
        message('\n\tBackfit limit is set to ', pars$matching$filtering$max.backfits, '.' )
          pars$matching$filtering$select <- 'random'
          
          if (pars$matching$filtering$select == 'rval'){
            # Sort matches by rval, then take the top n until # estimated backfits < limit
            
              sort.order <- order(match.info$rval, decreasing = TRUE)

          } else {
            # Randomly select matches to satisfy
              
              sort.order <- runif(nrow(match.info)) %>% order

          }

        # Figure out what the default bf.limit would be for this set of features (min subset) ####
          smallest.subset <- which.min(contribution[sort.order])
          cutoff.msg <- paste0('\n\tmatching$filtering$max.backfits (',
                               pars$matching$filtering$max.backfits,
                               ') may be set too low.',
                               '\n\tAt this setting, no matches would be kept. ',
                               '\n\tDefaulting to smallest subset (',
                               contribution[sort.order[smallest.subset]],
                               ') to preserve at least 1 match.')
          
        # Figure out how many matches we can keep, drawn in order from sort.order so we don't exceed max.backfits ####
          cutoff <- tryCatch({max(which(cumsum(contribution[sort.order]) < pars$matching$filtering$max.backfits))}, 
                              error = function(cond){message(cutoff.msg); 1}, 
                              warning = function(cond){message(cutoff.msg); 1})
              
              keep <- sort.order[1:cutoff]
              match.info <- match.info[keep, ]
              match.info <- match.info[order(keep),] #resort so mi is in same order as before
              lost <- length(sort.order)-length(keep)
            
        message('\n\t', lost, ' matches were jettisoned (', round(lost/length(sort.order)*100),'%)')
        message('\n\tThe new effective match rval cutoff is ', min(match.info$rval), '.')
        
      }

            # At this point, match.info is set. Assign IDs ####
                          
       match.info$id <- 1:nrow(match.info)
       # saveRDS(match.info,paste0(tmpdir,'/match.info.filtered.RDS'))
       match.info <- readRDS(paste0(tmpdir,'/match.info.filtered.RDS'))
      
######################### Back-fit reference to spectra  #############################    
    
      message('\n\nBack-fitting ref-feats to each spectrum in the relevant subset...\n')
      
      
    printTime()
    
    # Back-fit each matched reference region to the subset spectra
      # adjusted to account for sfe
        
        # backfit.results <- backfit_rfs3(match.info = match.info,
        #                                 feature.c = feature.c, # has sfe data 
        #                                 xmat = xmat,
        #                                 refmat.c = refmat.c, 
        #                                 ncores = pars$par$ncores)
      
# plot_grid_scattermore ####
# plot_grid_scattermore = function(dataframe) {
#   plot_list <- list()
#   n <- nrow(dataframe)
# 
#   # Loop through each column in the dataframe
#   for (var in names(dataframe)) {
#     print(var)
#     # Check if the column is numeric
#     if (any(!sapply(dataframe[[var]], is.numeric))) {
#       # Find the first non-numeric value
#       first_non_numeric <- which(!sapply(dataframe[[var]], is.numeric))[1]
#       # Prepare a message plot
#       message_plot <- ggplot() +
#         labs(title = paste("Non-numeric value found in column", var, "at row", first_non_numeric)) +
#         theme_void()
#       plot_list[[length(plot_list) + 1]] <- message_plot
#     } else {
#       # Ensure the column is numeric and prepare the plot
#       sorted_values <- sort(dataframe[[var]], na.last = TRUE)
#       plot <- scattermoreplot(x = 1:length(sorted_values), y = sorted_values, 
#                               main = paste("Scatter Plot of", var), xlab = "Index", ylab = var)
#       plot_list[[length(plot_list) + 1]] <- plot
#     }
#   }
# 
#   # Arrange the plots in a grid
#   # do.call(grid.arrange, c(plot_list, ncol = 2))
#   return(plot_list)
# }
# 
# Example usage:
# Assuming `your_dataframe` is your dataframe
# plot_grid_scattermore(your_dataframe)

# plotlist <- plot_grid_scattermore(match.info)
# ####

ncores = 2

# #####
  success = TRUE
  emptyRow <- function(){
    data.frame(ss.spec = NA,
              fit.intercept = NA, 
              fit.scale = NA,
              spec.start = NA,
              spec.end = NA,
              fit.fsa = NA,
              fit.rval = NA,
    )
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
             sfe = feature.c$sfe[feat.numbers])
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
        
        backfits.chunk <- lapply(seq(1,1000,100),#nrow(chunk$match.info),
                           function(m) {
          
          tryCatch({
            
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
                # end result of this loop will always be in the form: emptyRow()
                tryCatch(
                  {
                    # Get the ref region and spec data: ####
                      # s <- 21
                      if (s>10 & s < 15){
                        # print('triggered on s =',s,' in row ',m)
                        warning()
                      }
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
              
              fits <- fits[!is.na(rowSums(fits)),] # get rid of NA rows
              
              fits$pct.ref <- pct.ref
              
            # Return the fits dataframe rows (minimal data)  ####
              return(fits)
              
          }, warning = function(cond){
            
            # stop('backfit_rfs3: error in second layer loop, iteration: chunk$match.info row ', m)
            # If the whole match fails, return an NA fits df row
            return(emptyRow())
          }, error = function(cond){
            
            # stop('backfit_rfs3: error in second layer loop, iteration: chunk$match.info row ', m)
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
    
    # matches
      match.info <- match.info[bfs$inds.good,]
    
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

    
    
    
    
    
     
#############################        
        message('Saving backfits...\n\n\n')
        saveRDS(backfit.results, paste0(tmpdir,"/smrf.RDS"))
 
  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

  backfit.results <- backfit_rfs3(match.info = match.info[seq(1,nrow(match.info), 1000),],
                                        feature.c = feature.c, # has sfe data 
                                        xmat = xmat,
                                        refmat.c = refmat.c, 
                                        ncores = 2)
  
  