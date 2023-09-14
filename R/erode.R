erode <- function(v, vis = F){
  # v <- feat

  
  previous <- chop_peaks(v, chop = FALSE) # just characterize v
  current <- previous
  
  # simplePlot(rbind(current$v, v))
  
  rval <- cor(current$v, v)
  peaks.lost <- previous$n.pks - current$n.pks
  n.pks <- current$n.pks
  iteration <- 0
  real.pts <- current$valley.locs %>% length
  v.list <- list(current$v)
  n <- 1

  while (previous$n.pks > 0)
  {
      print(n)
      current <- chop_peaks(previous$v, pks = previous$pks, plots = FALSE)
      
    # Calculate stuff
      
      rval <- c(rval, cor(current$v, v))
      n.pks <- c(n.pks, current$n.pks)
      iteration <- c(iteration, n)
      real.pts <- c(real.pts, current$valley.locs %>% length)
      peaks.lost <- c(peaks.lost, previous$n.pks - current$n.pks)
      v.list[[n+1]] <- current$v
      
    # Update iterator
      
      n <- n + 1
      previous <- current
  }
    
  if (vis){
    stackplot(do.call(rbind, v.list[length(v.list):1]))
  }

  
  return(data.frame(iteration = iteration,
             n.pks = n.pks,
             r2 = rval^2,
             pts = real.pts,
             peaks.lost = peaks.lost)
         )
  # If we get rid of 
}

chop_peaks <- function(v, chop = TRUE, pks = NULL, plots = FALSE){
  
  if (is.null(pks)){
    pks <- extractPeaks_corr(v, plots = plots)
  }
  
  # Try eliminating the bottom fraction of prominences
    
    proms <- prominences(pks, feat) %>% .[pks$truePeak]
    
  # Calculate peaks info
  
    p.i <- pks_info(v, pks, keep.proms = which(proms > median(proms)))

  # Chop peaks away by interpolating between current peak bounds
    
    if(chop){
      v <- pracma::interp1(x = p.i$valley.locs, 
                           y = p.i$valley.heights, 
                           xi = 1:length(v))
    }
    
  # Recalculate peaks info
  
    pks <- extractPeaks_corr(v, plots = plots)
    
    p.i <- pks_info(v, pks, keep.proms = TRUE)
    
  
  current <- list(n.pks = p.i$truepks %>% length,
                  v = v,
                  pks = pks,
                  peak.locs = p.i$peak.locs,
                  peak.heights = p.i$peak.heights,
                  valley.locs = p.i$valley.locs,
                  valley.heights = p.i$valley.heights)
  return(current)
}

pks_info <- function(v, pks, keep.proms = TRUE){
  
  truepks <- which(pks$truePeak)
  
  peak.locs <- pks$peaks[-keep.proms]
    peak.heights <- feat[peak.locs]
  
  valley.locs <- pks$bounds[truepks] %>% unlist %>% as.numeric %>% 
    c(., 1, length(v), pks$peaks[keep.proms]) %>% unique %>% sort
    valley.heights <- feat[valley.locs]

  p.i <- list(
    truepks = truepks,
    peak.locs = peak.locs,
    peak.heights = peak.heights,
    valley.locs = valley.locs,
    valley.heights = valley.heights
  )
  return(p.i)
}
