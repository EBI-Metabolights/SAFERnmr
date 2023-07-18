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
    this.run <- paste0(tmpdir)
    dir.create(paste0(this.run, "/temp_data_matching"), showWarnings = F)

    ##################################################################################################################
    # Read data and set up ####
    message("### Setup for Parallel Matching")
    message("### MTJ 2023")
    message("")
    message("Loading data from files...\n\n\n")

    fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
      xmat <- fse.result$xmat
      ppm <- fse.result$ppm
      rm(fse.result)
      
    # This runs even if there was no clustering (clusters of 1 feature each):
    cluster.final <- readRDS(paste0(this.run, "/cluster.final.RDS"))
    
      c.labs <- cluster.final$labels
      

    # Put features in a matrix ####
    
      message('Building feature matrix...')
      
          feature.final <- readRDS(paste0(this.run, "/feature.final.RDS")) 
          featureStack <- feature.final$stack
          
          rm(feature.final)

    # If throttling comparisons, keep just a random subset of features
    
        max.features <- length(cluster.final$keys)
        
        if ( pars$debug$enabled == TRUE &
             nrow(featureStack) > pars$debug$throttle_features) 
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

    # Import and process the spectra for this dataset ####
      if (pars$galaxy$enabled) {
        lib.data <- readRDS(pars$galaxy$gissmo_location)
      } else {
        lib.data <- readRDS(pars$files$lib.data)
      }
      message(" - interpolating ref data to study ppm axis...\n\n")
      lib.data.processed <- prepRefs_for_dataset(lib.data,
          ppm.dataset = ppm,
          ref.sig.SD.cutoff = pars$matching$ref.sig.SD.cutoff,
          n.cores = pars$par$ncores
      )
      
        rm(lib.data)
        
      message('\nsaving processed ref library to file...')
      saveRDS(lib.data.processed, paste0(this.run, "/lib.data.processed.RDS"))

    ##################################################################################################################

      
    # Scale the feature matrix rows
      message('\tscaling feature matrix...\n') 
      f.stack <- f.stack %>% apply(1, scale_between) # this will also transpose it, so no need to do later

    # Put ref spectra in a matrix ####
      message('Building and scaling reference matrix...')
      ref.mat <- lapply(
          lib.data.processed,
          function(x) x$mapped$data %>% scale_between(0,1)
      ) %>% do.call(rbind, .)
  
        rm(lib.data.processed)

    # Pre-compute fts for refs, since they apply many times to each node ####
    # - create r.mat (padded, ft'd ref matrix)

        message("Pre-computing fts for refs (will take a minute)...\n")
        # Pad size for ref needs to be max.length(features)
    
        message('\tpadding ref matrix...')
        pad.size <- ncol(f.stack) - 1
    
        # Make the padded ref mat matrix
          r.mat <- ref.mat %>% padmat(use = 0, col.by = pad.size)
          
          # Transpose original matrix so columns are spectra, save and remove it to clear memory ####
            message("\nTransposing reference matrix (takes a few seconds)...\n\n")
            ref.mat <- t(ref.mat)
            message('\nWriting ref matrix to file...')
            saveRDS(ref.mat, paste0(this.run, "/temp_data_matching/ref.mat.RDS"))
          
          
          r.mat[is.na(r.mat)] <- 0
    
        # List-format the matrices to facilitate parallel
          message('\tSplitting ref matrix to lists...')
          r.mat <- lapply(1:nrow(r.mat), function(r) r.mat[r,])
        
        # Loop through spec matrix, compute fftw::fft()
          message('\tReference matrix fft...')
          gc() # garbage collect before parallel operation
          r.mat <- mclapply(r.mat, function(ref) fftw::FFT(ref), mc.cores = pars$par$ncores) %>% do.call(cbind,.)

    # Save ref data:
    
      
      
      message('\nWriting transformed ref data to file...')
      saveRDS(r.mat, paste0(this.run, "/temp_data_matching/rmat.RDS"))
      
        
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

      f.stack.split <- lapply(unique(f.grp), function(x) f.stack[, f.grp == x, drop = F])
        rm(f.stack)

  # Save feature data
    message('\nWriting feature data to file...')
    saveRDS(f.stack.split, paste0(this.run, "/temp_data_matching/f.stack.split.RDS"))

    rm(f.stack.split)
    
    saveRDS(pad.size, paste0(this.run, "/temp_data_matching/pad.size.RDS"))
    saveRDS(split.scheme, paste0(this.run, "/temp_data_matching/split.scheme.RDS"))

  message('\n--------------------------------------------------------------------------')
  message('-------------------  Parallel Matching Setup complete. -------------------')
  message('--------------------------------------------------------------------------')
}
