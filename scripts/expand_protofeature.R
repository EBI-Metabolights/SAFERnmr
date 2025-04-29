# Protofeatures are just a sketch of a correlation signature fragment
# The actual shapes, xmat segment,nad 

expand_protofeature <- function(p, xmat, ppm, half.window){
  
    driver <- p$driver

  # Driver locates the index, everything else can be built around it
    
    p.abs <- driver + p
    
    fullView <- (driver - half.window):(driver + half.window)
    
    in.bounds <- !(fullView < 1 | fullView > length(ppm))
    
    specreg.inds <- fullView[in.bounds]
    
    specRegion = xmat[,
                      specreg.inds]
    
    ppmRegion = ppm[specreg.inds]
    
      # simplePlot(specRegion, xvect = ppmRegion)
      
  # Recalculate cov and corr
  
    cv <- cr <- rep(NA, length(fullView))
    cv[in.bounds] <- cov(xmat[,driver], specRegion)
    cr[in.bounds] <- cor(xmat[,driver], specRegion)

  # Apply peaks
  
    # If this isn't a protofeature, but just a driver: let "peak" bounds = spec Region
    if (is.null(p.abs$primary.lower)){
      p.abs$primary.lower <- min(specreg.inds)
      p.abs$primary.upper <- max(specreg.inds)
      p.abs$secondary.lower <- p.abs$primary.lower
      p.abs$secondary.upper <- p.abs$primary.upper
    }
    
    primary.bounds <- c(p.abs$primary.lower, p.abs$primary.upper) %>% sort
    secondary.bounds <- c(p.abs$secondary.lower, p.abs$secondary.upper) %>% sort
    
    
    pk.mask <- rep(0, length(fullView))
    pk.mask[which(fullView %in% primary.bounds) %>% fillbetween] <- 1
    pk.mask[which(fullView %in% secondary.bounds) %>% fillbetween] <- 2
    peak.inds <- specreg.inds[pexp$peak.mask > 0]
    
    return(list(driver = driver,
                cv = cv, # covariance between driver and specRegion.inds - mask with [peak.mask>0] for ref.idx in log_storm
                cr = cr, # correlation between driver and specRegion.inds
                specRegion.inds = specreg.inds, # ppm indices for region of interest - mask with [peak.mask>0] for ref.idx in log_storm
                specRegion = specRegion, # xmat[, specRegion.inds] 
                ppmRegion = ppmRegion, # ppm vals in specRegion.inds
                peak.mask = pk.mask, # labels peaks 1, 2 >0 in specRegion.inds
                primary.bounds = primary.bounds, # ppm indices, bounds of peak
                secondary.bounds = secondary.bounds)) # ppm indices, bounds of peak
}
