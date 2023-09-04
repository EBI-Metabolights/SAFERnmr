climb <- function(start.pt, v){
  
  tryCatch(
    {
      pks <- extractPeaks_corr(v)
      bnds <- pks$bounds[pks$truePeak] %>% unlist %>% matrix(nrow = 2)
      start.pt.peak <- which((bnds[1, ] <= start.pt) & (start.pt <= bnds[2, ]))
      
      if (length(start.pt.peak) == 1){
        peak.max <- pks$peaks[pks$truePeak] %>% .[start.pt.peak]
        
        return(peak.max)
      }
      
      if (length(start.pt.peak) == 0){
        
        return(start.pt)
      }
  
      if (length(start.pt.peak) > 1){
        # go with the closer one?
        # go with the higher one?
        # do nothing:
        
        return(start.pt)
      }
          
    },
    error = function(cond){
      # warning('in climb: peak extraction failed for start.pt', start.pt, '. length(v) = ', length(v), '...')
      
      return(start.pt)
    }
  )
  
}
