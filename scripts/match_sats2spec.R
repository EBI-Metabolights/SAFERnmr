# Submit: an SAT
# Get: list of fit plots

fit_matches <- function(allmatches.feat, feat, ref.mat){
  fits <- lapply(1:nrow(allmatches.feat), function(m)
  {
      
    # Get f and r indices for this row
      f <- allmatches.feat[m, 'feat']
      r <- allmatches.feat[m, 'ref']
      feat.pos <- allmatches.feat[m, c('feat.start','feat.end')] %>% as.numeric %>% fillbetween
      ref.pos <- allmatches.feat[m, c('ref.start','ref.end')] %>% as.numeric %>% fillbetween

    # Get spectral signatures which matched
      ref <- ref.mat[,r,drop = F]
      
    # Fit
      
      fit <- fit_leastSquares(feat[feat.pos], ref[ref.pos], plots = F, scale.v2 = T)

      return(fit)
  })
  return(fits)
}


get_full_lib <- "/Documents/GitHub/SAFERnmr/scripts/slim.lib.R"

# Prep feature(s)

i <- i + 500
s <- sat.list[[i]]
plot_protofeature(p = data.frame(driver = s$peak),
          half.window = half.window, ppm = data$ppm,
          xmat = xmat[s$subset,],
          # bgplot = 'stack', line.shape = 'covar', line.color = "corr",
          bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
          showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)

pexp <- expand_protofeature(p = data.frame(driver = s$peak), 
                               xmat[s$subset,], ppm, half.window)

# simplePlot(ymat = pexp$cv, xvect = ppm[pexp$specRegion.inds])

featureStack <- pexp$cv %>% length %>% matrix(NA, nrow = 1, ncol = .)
mask <- (pexp$specRegion.inds %in% s$ref.idx)
featureStack[, mask] <- pexp$cv[mask]
# plot(featureStack%>%t)

load_lib_data <- function(pars){
  # lib.data <- readRDS(pars$files$lib.data)
  lib.data <- lib.data.700
  
  ppm.range <- range(pexp$ppmRegion)
  ppm.tol <- 1
  ref.range <- c(ppm.range[1]-ppm.tol, ppm.range[2]+ppm.tol)
  
  lib.data.processed <- prepRefs_for_dataset(lib.data,
                                             ppm.dataset = ppm,
                                             ref.sig.SD.cutoff = 0,#pars$matching$ref.sig.SD.cutoff,
                                             ppm.range = ref.range,
                                             n.cores = pars$par$ncores
  )
  return(lib.data.processed)
}

lib.data.processed <- load_lib_data(pars)

# Scale the feature matrix rows ####
  f.stack <- featureStack %>% apply(1, scale_between) # this will also transpose it, so no need to do later

# Put ref spectra in a matrix ####
  ref.mat <- lapply(lib.data.processed, function(x)
    {
      x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]]
    }
  ) %>% do.call(rbind, .)

  ref.mat[1,] %>% simplePlot
  
  
  # Pad the ref spectra to feature size
  pad.size <- nrow(f.stack) - 1
  r.mat <- ref.mat %>% padmat(use = 0, col.by = pad.size)
  r.mat[is.na(r.mat)] <- 0
  
  # Transpose original matrix so columns are spectra, save and remove it to clear memory ####
    ref.mat <- ref.mat %>% compress_stack(sparse.val = NA)
  
  # List-format the matrices to facilitate parallel
    message('\tSplitting ref matrix to lists...')
    r.mat <- lapply(1:nrow(r.mat), function(r) r.mat[r,])
  
  # Loop through spec matrix, compute fftw::fft()
    message('\tReference matrix fft...')
    r.mat <- mclapply(r.mat, function(ref) fftw::FFT(ref), mc.cores = pars$par$ncores) %>% do.call(cbind,.)
    
    ref.mat <- ref.mat %>% cstack_expandRows %>% t

# At this point, we have:
# - features:         f.stack (vertical)
# - references:       ref.mat (vertical; compressed)
# - padded refs:      r.mat
# - fft'd references: r.mat   (vertical)
  
    f.ind = 1
    feat = f.stack[,f.ind, drop = F]
      plot(feat)
    f.num = f.ind
    
    padded.feat <- feat %>% c(rep(0, nrow(r.mat) - length(feat)))
    padded.feat[is.na(padded.feat)] <- 0

    feat.ft <- Conj(fftw::FFT(padded.feat))

    message("    - cross-correlating to refs...")
    
    # Locate best positions in all available refs  ####
    
      # r.mat # ftd refs
      # ref.mat # refs
      # feat
      # feat.ft
      # pad.size
      
      allmatches.feat <-  foreach(r.num = 1:ncol(r.mat),
                                  ref = ref.mat,
                                  ref.ft = r.mat,
                                  .combine = 'rbind',
                                  .errorhandling="pass") %do%
      {
        # r.num = 1
        # ref = ref.mat[,r.num, drop = F]
        # ref.ft = r.mat[,r.num, drop = F]
        # 
        # simplePlot(feat %>% t %>% trim_sides(out = "inds") %>% feat[.])
        # simplePlot(ref %>% t %>% trim_sides(out = "inds") %>% ref[.])
        
        # Cross-correlate to find locations and scores:
          matches <- feature_match2ref_slim(f.num, r.num,
                                            feat, ref,
                                            pad.size = pad.size,
                                            feat.ft, ref.ft,
                                            max.hits = 5,#pars$matching$max.hits,
                                            r.thresh = .7,#pars$matching$r.thresh,
                                            p.thresh = .01)#pars$matching$p.thresh)
          
          return(matches)
      }
    
    # # There is a distribution of matches
    # allmatches.feat$rval %>% sort %>% plot
    # 
    # # The whole feature width is matched
    # # (allmatches.feat$ref.end-allmatches.feat$ref.start) %>% sort %>% plot
    # 
    # # However, some points are not matched
    # (allmatches.feat$pts.matched) %>% sort %>% plot
    # 
    # # Does the Number of points > 0.8 * correlation indicate much?
    # (allmatches.feat$pts.matched * allmatches.feat$rval) %>% sort %>% plot
     
    # Escape and return nothing if there are no matches to evaluate: ####
        # if (is.null(allmatches.feat)){return(NA)}
    
        ######################################################################

      # Evaluate fits for top-scoring positions (regardless of ref) ####
        message("    - calculating ", nrow(allmatches.feat), " fits...")
    
        allmatches.fits <- fit_matches(allmatches.feat, feat, ref.mat)
     
        # i <- 0
        # i <- i+1
        # plot_fit(allmatches.fits[[i]],type = "auc")
        
        
        # Extract out minimal fit data ####
          fit.data <- lapply(allmatches.fits, function(f) f$fit %>% as.numeric) %>% do.call(rbind,.)
          allmatches.feat[,"fit.intercept"] <- fit.data[,1]
          allmatches.feat[,"fit.scale"] <- fit.data[,2]
  
        # Add some different scores from the fits ####
          message("    - calculating additional scores...")
          allmatches.feat[,'sum.residuals'] <-
            lapply(allmatches.fits, function(x) x$sum.residuals) %>% unlist
          allmatches.feat[,'rmse'] <-
            lapply(allmatches.fits, function(x) x$rmse) %>% unlist
        
              # Format results ####
      
  # return(allmatches.feat)

  
  n.passing <- sum(allmatches.feat$rval >= pars$matching$r.thresh)
  n.refs <- ncol(r.mat)
  specificity.score <- 1-(n.passing/n.refs)
  

# Functionalize the matching, fitting, and match-fitting. 

# Feature -> Refs

  # Inputs: 
  #   feature(s) [align to] 
  #   ref(s). 
  #   Both are JUST MATRICES, NAs okay. 
  
  # Process both. 
  # Things to keep in mind:
  #   features have the same sizes, and refs have the same sizes
  #   should refs be region-limited?
  #   should refs be compressed, then unpacked for a set of features?
  #   ref compression is pointless if no repeated vals
  #   fft operations should NOT be repeated for feature-pair comparisons unless absolutely necessary (for memory purposes)
  #   so, the initial strategy should be to compute ffts and handle the vectors in full.
  #   
  #   Alternatively, trim the refmat to the region, compute the xcorr, and time it. 
  #     thought: if the ffts lived on backend nodes and were available to respond to requests, it would be way faster..
  #     
  # Compute pairwise cross-correlations
  #   get xcorr: positions, pseudo.rvals for top N hits for pseudo.rval > threshold
  #   format the compare vectors for corr, fit calcs from the match table
  #   *recompute real rvals for top N hits - this is done in feature_match2ref_slim()
  #   
  #   get fits (for later): leastSquares
  #   * match table to plot, match table to vectors (NA-filled)
  #     i.e. expand_match(match_table, features, refs, feature.ppms = NULL, ref.ppms = NULL)
  #     i.e. expand_fit_match(match, fit.type)
  #     
  #   * NOTE: record the total number of hits: pseudo.rval > threshold

  # Return a table of match results
  
# For each feature, align to dataset spectra
# Align to refs

selected <- seq(1,length(allmatches.fits), by=1)
lib.data.processed <- load_lib_data(pars)


plots <- pbapply::pblapply(selected, function(f){
  fit <- allmatches.fits[[f]]
  match <- allmatches.feat[f, ]
  # Set ppms
  # range(not.na(spec.fit))
  ref.ppm <- lib.data.processed[[match$ref]]$mapped$data.compressed %>% expand_stacklist(which.stacks = 'ppm') %>% .[[1]]
  ppm.match <- ref.ppm[match$ref.start:match$ref.end]
  
  # Set the width to be equal to original feature
  
  plot_fit(fit,type = "auc", ppm = ppm.match)
}) 

plots %>% grid_pdf(plotLoc=tmpdir, filename="/sat_ref_matches.pdf")

grid_pdf <- function(plots=NULL, plotLoc="./", filename="grid_plot.pdf"){
  # How big to make the page? 2 inches for each plot, and grid will be square.
  dim <- 3*round(sqrt(length(plots)))
  pdf(file = str_c(plotLoc,filename),   # The directory you want to save the file in
      width = dim, # The width of the plot in inches
      height = dim)
  
  gridExtra::grid.arrange(grobs = plots)
  
  dev.off()  
}





          