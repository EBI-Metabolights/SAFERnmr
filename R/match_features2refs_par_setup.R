#' match.features2refs.par.explicit.setup function
#'
#' Chunks feature data and precomputes ref ffts for match.features2refs.par.explicit. 
#'
#' @param pars a list of input parameters
#'
#' @return writes necessary files for matching script
#' @importFrom magrittr %>%
#' @import fftw
#' @importFrom parallel mclapply
#'
#' @export
match_features2refs_par_setup <- function(pars) {
    message("-------------------------------------------------------------")
    message("-------------------     Matching Setup    -------------------")
    message("-------------------------------------------------------------")
    printTime()
    message("\n\n\n")
    
    ##################################################################################################################
    # Parse file paths ####

    tmpdir <- pars$dirs$temp
    dir.create(paste0(tmpdir, "/temp_data_matching"), showWarnings = F)

    ##################################################################################################################
    # Read data and set up ####
    message("### Setup for Parallel Matching")
    message("### MTJ 2023")
    message("")
    message("Loading data from files...\n\n\n")

    fse.result <- readRDS(paste0(tmpdir, "/fse.result.RDS"))
      fse.result %>% test_nullish('fse.result did not pass nullish check. Cannot proceed with matching setup.')
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      rm(fse.result)
       
      
    # This runs even if there was no clustering (clusters of 1 feature each):
    cluster.final <- readRDS(paste0(tmpdir, "/cluster.final.RDS"))
        cluster.final %>% test_nullish('cluster.final did not pass nullish check. Cannot proceed with matching setup.')
    
      c.labs <- cluster.final$labels

    # Put features in a matrix ####
      
      
      message('Building feature matrix...')
      
          feature.final <- readRDS(paste0(tmpdir, "/feature.final.RDS")) %>% expand_features
            feature.final %>% test_nullish('feature.final')
          
          featureStack <- feature.final$stack
          
        # Check that cluster labels actually correspond to rows of feature matrix. If not, redo clustering as trivial clustering
          if(!all(c.labs %in% 1:nrow(featureStack))){
            message('Cluster labels do not correspond with feature rows. Reverting to comprehensive feature matching...')
            cluster_features(pars, feature.final, min.features = 1000, do.clustering = FALSE)
          }
          
          rm(feature.final)

    # If throttling comparisons, keep just a random subset of features
    
        max.features <- length(cluster.final$keys)
        
        if ( 
             pars$debug$enabled == TRUE &
             nrow(featureStack) > pars$debug$throttle_features
           ) 
        {

            max.features <- pars$debug$throttle_features
          
        }
      
      # Allow use of representative.feature profiles, or cluster weighted.means:
    
      if (pars$matching$cluster.profile == 'representative.feature'){
        
        # If using a representative feature, just get the stack. 
         
          f.subset <- cluster.final$keys %>% .[1:max.features]
          
          rm(cluster.final)
          
      } else {
  
        # If using the recalculated cluster profiles, get those:
        # - be aware that matches MUST be propagated to cluster members later on.
        stop('using the recalculated cluster profiles is still under development. Sorry!')
        p.width <- lapply(cluster.final$info, function(x) x$profile %>% length) %>% unlist %>% max
        
        featureStack <- lapply(cluster.final$info, function(x) c(x$profile, rep(NA, p.width-length(x$profile)))) %>% do.call(rbind,.)
        f.subset <- 1:nrow(featureStack) %>% .[1:max.features]
        
      }
    
      nfeats <- length(f.subset)
      f.stack <- featureStack[f.subset, , drop = F]
      
    ##################################################################################################################
    ## Ref data import ####
    message("Loading and processing reference spectrum data...\n")
    
    lib.data <- lib.data.processed <- NULL
  
    # Import and process the spectra for this dataset ####
        if (pars$galaxy$enabled) {
          
          message('Detected Galaxy usage. Looking in galaxy$gissmo_location (',pars$galaxy$gissmo_location,') for lib.data...')
          
          # Try to get lib.data from galaxy$gissmo_location. If that fails, try files$lib.data. If that fails, return NULL
          lib.data <- tryCatch(
            {          
                message('\tReading lib.data from ',galaxy$gissmo_location,'...')
                
                readRDS(pars$galaxy$gissmo_location)
            },
            error = function(cond){
              
              warning('In match...par_setup: galaxy$gissmo_location ',pars$galaxy$gissmo_location,' does not exist...')
              
              tryCatch(
                {
                  message('Looking in files$lib.data (',pars$files$lib.data,') for lib.data...')
              
                  # Try to get lib.data from files$lib.data:
                    
                    lib.data <- readRDS(pars$files$lib.data)
                    lib.data %>% test_nullish
                    
                    return(lib.data)
                }, 
                error = function(cond){NULL}
              )
  
            }
          )
            
        } else {
          
          # If not using Galaxy, just 
          lib.data <- tryCatch(
            {
              message('Looking in files$lib.data (',pars$files$lib.data,') for lib.data...')
          
              # Try to get lib.data from files$lib.data:
                
                readRDS(pars$files$lib.data)
                
            }, 
            error = function(cond){NULL}
          )
          
        }
    
      # If lib.data was read, process it. Default is always reprocess: ####
        if (!is.null(lib.data)){
             
            # Process the data for the dataset: ####
              message(" - interpolating ref data to study ppm axis...\n\n")
              lib.data.processed <- prepRefs_for_dataset(lib.data,
                                                         ppm.dataset = ppm,
                                                         ref.sig.SD.cutoff = pars$matching$ref.sig.SD.cutoff,
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
        
      
    ##################################################################################################################

      
    # Scale the feature matrix rows ####
      message('\tscaling feature matrix...\n') 
      f.stack <- f.stack %>% apply(1, scale_between) # this will also transpose it, so no need to do later

    # Put ref spectra in a matrix ####
      message('Building and scaling reference matrix...')
      
      ref.mat <- lapply(lib.data.processed, function(x)
        {
          x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]] %>% scale_between(0,1)
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
}
