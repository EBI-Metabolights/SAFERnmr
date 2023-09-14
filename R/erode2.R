erode2 <- function(v, vis = F){
  # v <- feat
#####################################################################################

rm.smallest.pk <- function(v, pks = NULL){
  
  if (is.null(pks)){
    pks <- extractPeaks_corr(v) %>% only_true
  }
  
  # Eliminate the smallest prominence
    
    proms <- prominences(pks, v)
    
    rep.idx <- which.min(proms) %>% pks$bounds[.] %>% unlist
    fill.pts <- rep.idx %>% fillbetween
    
  # Chop peak away by interpolating between current peak bounds
    
      v[fill.pts] <- pracma::interp1(x = rep.idx, 
                           y = v[rep.idx], 
                           xi = fill.pts)
    
  # Recalculate peaks info
  
    pks <- extractPeaks_corr(v) %>% only_true
  
  current <- list(v = v,
                  pks = pks,
                  # npoints = c(pks$peaks, 
                  #             pks$bounds %>% unlist, 
                  #             1, length(v)) %>% unique %>% sort %>% length,
                  n.pks = length(pks$peaks),
                  rm.pk = which.min(proms))
  return(current)
}

#####################################################################################

only_true <- function(pks){
  pks$bounds <- pks$bounds[pks$truePeak]
  pks$peaks <- pks$peaks[pks$truePeak]
  pks
}

#####################################################################################
  pks <- extractPeaks_corr(v) %>% only_true
  
  previous <- list(
    v = v,
    pks = pks,
    n.pks = length(pks$peaks)
  )
  
  current <- previous
  
  rval <- cor(current$v, v)
  peaks.lost <- previous$n.pks - current$n.pks
  n.pks <- current$n.pks
  iteration <- 0
  v.list <- list(current$v)
  rm.pk <- list(NULL)
  
  n <- 0
  
  while (previous$n.pks > 0)
  {
      # print(n)
    # Remove the smallest peak
      
      current <- rm.smallest.pk(v = previous$v, 
                                pks = previous$pks)
    
    # Calculate stuff
      
      rval <- c(rval, cor(current$v, v))
      n.pks <- c(n.pks, current$n.pks)
      iteration <- c(iteration, n)
      peaks.lost <- c(peaks.lost, previous$n.pks - current$n.pks)
      v.list[[n+1]] <- current$v
      rm.pkbounds <- current$rm.pk %>% previous$pks$bounds[.] %>% unlist
      rm.pk[[n+1]] <- rm.pkbounds

    # Update iterator and obj
      
      n <- n + 1
      previous <- current    

  }
  


  if (vis){
    stackplot(do.call(rbind, v.list[length(v.list):1]))
  }


  return(
      list( 
            it.info = data.frame(iteration = iteration,
                                     n.pks = n.pks,
                                     r2 = rval^2,
                                     pts = real.pts,
                                     peaks.lost = peaks.lost),
             rm.pk = rm.pk
          )
  
         )

}

