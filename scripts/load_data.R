# Load Data
# - import from RDS in spectral matrix format
# - adjust resolution if necessary
# 
# Follows:
# - setup

################ Read parameters file ##################
  
  tmpdir <- pars$dirs$temp

################ Get Data from RDS ##################

    X_raw <- readRDS(pars$files$spectral.matrix)
    # X_raw <- readRDS(paste0(tmpdir, "/spectral.matrix.RDS"))
      xmat <- X_raw[-1,]          # spectral matrix (each row is a spectrum; 
                                  # doesn't require alignment. normalization ok
                                  # but not necessary. scaling, no.)
      ppm <- X_raw[1,]            # ppm vector (corresponding to cols of xmat)
      digital.res <- ppm %>% diff %>% mean %>% abs # ppm per element
      
      message('digital resolution of xmat is ', digital.res, ' ppm.')
      if (!is.null(pars$opts$npoints)){
        if (pars$opts$npoints < length(ppm)){
          
          # Re-interpolate dataset spectra to a lower number of points to save compute
            rs <- resample_spectra(xmat, ppm, npoints = pars$opts$npoints, cores = pars$par$ncores)
            ppm <- rs$ppm
            xmat <- rs$spectra
            if (nrow(rs$spectra) < pars$tina$min.subset){stop('The minimum number of spectra is: ', pars$tina$min.subset,'. The submitted dataset only has ', nrow(rs$spectra),'. SAFER terminated.')}
            
            rm(rs)
          
          digital.res <- ppm %>% diff %>% mean %>% abs # ppm per element
          message('Using opts:npoints @ ',pars$opts$npoints, ' points. 
                  \nNew xmat digital resolution: ', digital.res)

        }
      }
      # if not set, do nothing (warning is printed)

