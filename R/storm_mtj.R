#' My attempt at a pure STORM as described in the manuscript's SI.
#'
#' @param X A matrix of spectral data (rows are spectra, columns are spectral points)
#' @param ppm A vector of the spectral points in ppm (optional, default is all columns of X)
#' @param b An integer giving the expansion parameter for the reference peak
#' @param q A numeric giving the p-value threshold for correlation significance (both subset AND reference optimization)
#' @param ref integer (if providing best spectrum) or vector (if providing feature shape)
#' @param ref.index A vector of the spectral points (columns of X) to use as the initial reference
#'
#' @return A list 
#'
#' @importFrom magrittr %>%
#' 
#' @export
storm_mtj <- function( X=NULL, 
                       ppm=NULL, 
                       b=30, 
                       q=0.05, 
                       ref=NULL,
                       ref.index=NULL, 
                       driver = NULL ){

############ Setup ##################################################  
    
  # Provide checks for input parameters
    # X:
    # - must be a matrix with ppm on columns, spectra on rows
    # - must not contain NAs
    # - ideally not scaled
    # 
    # ppm: 
    # - vector of ppm values corresponding to columns of X
    # 
    # b: 
    # - integer
    # - number of points to expand the window by on either side of the driver
    # 
    # How is the region selected?
    # - ref index is supplied as xcols (can be discontinuous)
    # 
    # driver:
    # - x column used to drive the STOCSY; default is max
    
  
############ Initialize for the loop ########################################

    n <- nrow(X)
    
    ns <- 5
    
    b.half <- 15
 
    itlimit = 25
    
    if (length(ref) == 1){
      # Assume spectrum number
      ref <- X[ref, ref.index]
    }
    
############ Run storm loop ###################################################################
    
    i=1       
    
    subset = list() # cannot be same as .previous or loop won't run
      subset[[1]] <- 1:n # first subset includes all spectra

    while( (subset %>% not_seen_yet) & (i < itlimit)){ 
      
  ## Update the subset ########################################################################
      print(i)
      subset.i <- subset[[1]]
      
      Xr=X[subset.i, ref.index]
      
      r=cor(t(Xr), ref)
      pval <- two_t_cdf(r)
      
      sspass <- (pval<q & r>0) %>% which
      r <- r[sspass]
      pval <- pval[sspass]

      # Pick the spectra with the lowest ns pvals #####

        index <- order(pval, decreasing = FALSE)
        
        subset.i <- subset.i[index]
        subset.i <- subset.i[1:ns]
        subset.i <- subset.i[!is.na(subset.i)]
        
        
  ## Update the reference ########################################################################
       
    # update driver
      
        driver <- which.max(ref) %>% 
          ref.index[.] # put back in terms of xcols

    # STOCSY the new driver within subsets.previous and the widened window to get new ref ##############
      
        ref.index.expanded <- expand_window(ref.index, 
                                            within = c(1:ncol(X)), 
                                            by = b.half + 1, 
                                            keep.nas = FALSE)
        
        Xprime <- X[subset.i, ref.index.expanded]
        
        covar <- cov(Xprime, X[subset.i,driver])
        r <- cor(Xprime, X[subset.i,driver])
        pval <- two_t_cdf(r)

        ref.pass <- r > 0 & pval < q
        
        ref <- covar[ref.pass]
        ref.index <- ref.index.expanded[ref.pass]
        
      # Finish the loop by updating the counter
        subset[[i]] <- subset.i
        i <- i+1
        
        # Plot 
          # ppm.region <- ref.index %>% range %>% fillbetween
          # stackplot(X[ref, ppm.region])
          simplePlot(ref, xvect = ref.index) + geom_vline(xintercept = driver)

        Sys.sleep(1)
    }
  
############ Finish up and return results ########################################################

  return(list(subsets = subset,
              ref.index = ref.index,
              ref = ref,
              corr = r,
              pval = pval,
              covar = covar,
              driver = driver,   # index in ref.expanded$wind
              pass = ref.pass,  # indices in ref.expanded$wind
              driver.initial = driver.init,
              iterations = i-1) # (completed iterations only) 
         )
 
}

not_seen_yet  <- function(subset){
  
  if (length(subset) < 2){return(TRUE)}
  i <- length(subset)
  subset.i <- subset[[i]]
  
  lapply(subset[1:(i-1)], function(ss) 
    
    !identical(subset[[ss]], subset.i)
    
  ) %>% any
}

two_t_cdf <- function(corr){
  a=-abs(corr * sqrt((length(corr)-2)/(1-corr^2)))
  pval=2*pt(a,(length(corr)-2))
  return(pval)
}
