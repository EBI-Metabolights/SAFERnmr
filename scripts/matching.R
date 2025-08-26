## On slim SAFER, trigger matching on a specific SAT

# SAT expansion
  # s
  pexp <- expand_protofeature(p = data.frame(driver = s$peak), 
                              xmat, ppm, half.window)
  
  # pexp$cv
  # pexp$specRegion.inds
  # pexp$specRegion
  
# Ref dataset 

    printTime()
    message("\n\n\n")
    
    ##################################################################################################################
    message("")
    message("Loading data from files...\n\n\n")

      ppm <- fse.result$ppm
      rm(fse.result)
      
    # Put features in a matrix ####
      
      message('Building feature matrix...')
      
          featureStack <- pexp$cv %>% length %>% matrix(NA, nrow = 1, ncol = .)
          # plot(pexp$specRegion.inds, pexp$cv)
          mask <- !(pexp$specRegion.inds %in% s$ref.idx)
          featureStack[, mask] <- pexp$cv[mask]
          
      nfeats <- nrow(featureStack)
      f.stack <- featureStack
      
    ##################################################################################################################
    ## Ref data import ####
    message("Loading and processing reference spectrum data...\n")
    
    lib.data <- lib.data.processed <- NULL
  
    # Import and process the spectra for this dataset ####
          
          lib.data <- tryCatch(
            {
              message('Looking in files$lib.data (',pars$files$lib.data,') for lib.data...')
          
              # Try to get lib.data from files$lib.data:
                
                readRDS(pars$files$lib.data)
                
            }, 
            error = function(cond){NULL}
          )
          
    
      # If lib.data was read, process it. Default is always reprocess: ####
        if (!is.null(lib.data)){
             
            # Process the data for the dataset: ####
              reg <- pars$corrpockets$only.region.between
              ppm.reg <- c(min(reg) - pars$matching$filtering$ppm.tol, max(reg) + pars$matching$filtering$ppm.tol)
              message(" - interpolating ref data to study ppm axis...\n\n")
              lib.data.processed <- prepRefs_for_dataset(lib.data,
                                                         ppm.dataset = ppm,
                                                         ref.sig.SD.cutoff = pars$matching$ref.sig.SD.cutoff,
                                                         ppm.range = pars$corrpockets$only.region.between,
                                                         n.cores = pars$par$ncores
              )
            
            message('\nsaving processed ref library to file...')
            saveRDS(lib.data.processed, paste0(tmpdir, "/lib.data.processed.RDS"))
        
            rm(lib.data)
        
        } else {
          
      # If lib.data.RDS was NOT read, check for/read in lib.data.processed directly (e.g. if re-running from results file): ####
      
         lib.data.processed <- tryCatch(
            {
               message('\nReading processed ref library data...')
               readRDS(paste0(tmpdir, "/lib.data.processed.RDS"))
            },
            error = function(cond){
               stop('No "lib.data.RDS" or "lib.data.processed.RDS" file found. ',
                    '\nCheck the following settings in params.yaml: ',
                    '\n\t- files$lib.data   or ',
                    '\n\t- galaxy$gissmo_location (for Galaxy runs) ',
                    '\n"lib.data.RDS" should exist in galaxy$gissmo_location (for Galaxy runs) or files$lib.data (for local/HPC runs)')
            }
          )
  
        }
        nrefs <- length(lib.data.processed)
      
    ##################################################################################################################

      
    # Scale the feature matrix rows ####
      message('\tscaling feature matrix...\n') 
      f.stack <- f.stack %>% apply(1, scale_between) # this will also transpose it, so no need to do later

    # Put ref spectra in a matrix ####
      message('Building reference matrix...')
      
      ref.mat <- lapply(lib.data.processed, function(x)
        {
          x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]]
        }
      ) %>% do.call(rbind, .)
  
        rm(lib.data.processed)

    # Pre-compute fts for refs, since they apply many times to each node ####
    # - create r.mat (padded, ft'd ref matrix)

        message("Pre-computing fts for refs (will take a minute)...\n")
        # Pad size for ref needs to be max.length(features)
    
        message('\tpadding ref matrix using feature size - 1...')
        pad.size <- nrow(f.stack) - 1
    
        # Make the padded ref mat matrix
          r.mat <- ref.mat %>% padmat(use = 0, col.by = pad.size)
          
          # Transpose original matrix so columns are spectra, save and remove it to clear memory ####
            # message("\nTransposing reference matrix (takes a few seconds)...\n\n")
            # ref.mat <- t(ref.mat)
              ref.mat %>% test_nullish
          
            ref.mat <- ref.mat %>% compress_stack(sparse.val = NA)
            
            message('\nWriting compressed ref matrix to file...')
            saveRDS(ref.mat, paste0(tmpdir, "/temp_data_matching/ref.mat.RDS"))
          
          
          r.mat[is.na(r.mat)] <- 0
    
        # List-format the matrices to facilitate parallel
          message('\tSplitting ref matrix to lists...')
          r.mat <- lapply(1:nrow(r.mat), function(r) r.mat[r,])
        
        # Loop through spec matrix, compute fftw::fft()
          message('\tReference matrix fft...')
          gc() # garbage collect before parallel operation
          r.mat <- mclapply(r.mat, function(ref) fftw::FFT(ref), mc.cores = pars$par$ncores) %>% do.call(cbind,.)

    # Save ref data:
    
      
      r.mat %>% test_nullish('r.mat')
      message('\nWriting transformed ref data to file...')
      saveRDS(r.mat, paste0(tmpdir, "/temp_data_matching/rmat.RDS"))
      
        
        rm(r.mat)
        rm(ref.mat)

##################################################################################################################
# Split the feature matrices for distribution across nodes ####
    
  # Feat data is split between nodes:
      message('\nSplitting feature matrices for distribution across nodes...')
      chunk.size <- max(1, nfeats / pars$par$ncores)
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

  # Save feature data
    message('\nWriting split feature data to file...')
    saveRDS(f.stack.split, paste0(tmpdir, "/temp_data_matching/f.stack.split.RDS"))

    rm(f.stack.split)
    
    saveRDS(pad.size, paste0(tmpdir, "/temp_data_matching/pad.size.RDS"))
    saveRDS(split.scheme, paste0(tmpdir, "/temp_data_matching/split.scheme.RDS"))

  message('\n--------------------------------------------------------------------------')
  message('-------------------  Parallel Matching Setup complete. -------------------')
  message('--------------------------------------------------------------------------')
  
  
  
  
  
  
  
feat_align_to(align = xmat.ss, to = profile.exp, max.hits = 1)
feat_align_to(align = , to = , max.hits = 5, max.lag = 100)