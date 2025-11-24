feature_match2ref_pcc <- function(f.num, r.num,
                                  feat, ref,
                                  r.thresh = 0.8,
                                  max.hits = 100) {
  
  
  # 1. compute NA-masked sliding Pearson correlations
  xc.res <- compute_pearson_sliding_overlap(feat, ref)
  xc   <- xc.res$xcorr
  lags <- xc.res$lags
  pts  <- xc.res$pts
  
  # 2. find local maxima
  peaks <- localMaxima(xc)
  
  # restrict to complete-overlap region
  valid <- peaks >= xc.res$range.complete.overlap[1] &
           peaks <= xc.res$range.complete.overlap[2]
  peaks <- peaks[valid]
  
  if (length(peaks) == 0) return(NULL)
  
  # restrict by r threshold
  peaks <- peaks[xc[peaks] >= r.thresh]
  if (length(peaks) == 0) return(NULL)
  
  # sort by peak height
  peaks <- peaks[order(xc[peaks], decreasing = TRUE)]
  peaks <- peaks[seq_len(min(max.hits, length(peaks)))]
  
  # Extract vectors
  lag_vec <- lags[peaks]
  rval    <- xc[peaks]
  npts    <- pts[peaks]
  
  # Vectorized p-values
  tstat <- abs(rval * sqrt((npts - 2) / pmax(1e-12, (1 - rval^2))))
  pvals <- 2 * pt(-tstat, df = npts - 2)
  
  # reference/feature lengths
  f.len <- xc.res$feat_len
  r.len <- xc.res$ref_len
  
  # padded reference location
  ref_start_pad <- f.len
  ref_end_pad   <- f.len + r.len - 1
  
  # feature padded coordinates
  feat_start_pad <- ref_start_pad + lag_vec
  feat_end_pad   <- feat_start_pad + f.len - 1
  
  # convert padded → reference coordinates
  ref_start <- feat_start_pad - (ref_start_pad - 1)
  ref_end   <- ref_start + f.len - 1
  
  # clip to reference
  ref_start <- pmax(ref_start, 1)
  ref_end   <- pmin(ref_end, r.len)
  
  # assemble output
  data.frame(
    feat = f.num,
    ref  = r.num,
    lag  = as.numeric(lag_vec),
    rval = as.numeric(rval),
    pval = as.numeric(pvals),
    pts.matched = as.numeric(npts),
    pts.feat = f.len,
    
    # corrected semantics
    feat.start = 1,
    feat.end   = f.len,
    ref.start  = ref_start,
    ref.end    = ref_end,
    
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

compute_pearson_sliding_overlap <- function(feat, ref) {
  
  feat <- as.double(feat)
  ref  <- as.double(ref)
  
  f.len <- length(feat)
  r.len <- length(ref)
  N <- f.len + r.len - 1
  
  # lag indexing identical to FFT Pearson
  lags <- -(f.len - 1):(r.len - 1)   # length N
  
  pearson <- rep(NA_real_, N)
  pts_vec <- integer(N)
  
  # Reference occupies padded indices:
  ref_start_pad <- f.len
  ref_end_pad   <- f.len + r.len - 1
  
  # Feature padded positions (vectorized)
  feat_start_pad <- ref_start_pad + lags
  feat_end_pad   <- feat_start_pad + f.len - 1
  
  # MAIN LOOP (just slicing + mask)
  for (i in seq_len(N)) {
    
    fs_pad <- feat_start_pad[i]
    fe_pad <- feat_end_pad[i]
    
    # enforce bounds
    if (fs_pad < 1 || fe_pad > N)
      next
    
    # local indices for feature
    local_feat_idx <- 1:f.len
    
    # global padded indices for this alignment
    global_f <- fs_pad:fe_pad
    
    # overlap with reference padded region
    valid_ref <- (global_f >= ref_start_pad) & (global_f <= ref_end_pad)
    
    if (!any(valid_ref))
      next
    
    # convert to reference indices
    ref_idx <- global_f[valid_ref] - (ref_start_pad - 1)
    feat_idx <- local_feat_idx[valid_ref]
    
    x <- feat[feat_idx]
    y <- ref[ref_idx]
    
    use_mask <- !(is.na(x) | is.na(y))
    n <- sum(use_mask)
    
    pts_vec[i] <- n
    
    if (n < 3)
      next
    
    x <- x[use_mask]
    y <- y[use_mask]
    
    # Fast Pearson
    xm <- mean(x)
    ym <- mean(y)
    
    cov_xy <- sum((x - xm) * (y - ym))
    var_x  <- sum((x - xm)^2)
    var_y  <- sum((y - ym)^2)
    
    pearson[i] <- cov_xy / sqrt(var_x * var_y)
  }
  
  list(
    feat = feat,
    ref = ref,
    xcorr = pearson,
    lags = lags,
    pts = pts_vec,
    feat_len = f.len,
    ref_len = r.len,
    range.complete.overlap = c(f.len, r.len)
  )
}

convert_for_plot <- function(feat, ref, xc, match_row) {
  
  f.len <- length(feat)
  r.len <- length(ref)
  N <- f.len + r.len - 1
  
  # padded matrix
  vals <- matrix(NA_real_, 3, N)
  rownames(vals) <- c("corr", "ref", "feat")
  
  # place reference in padded space
  ref_start_pad <- f.len
  ref_end_pad   <- f.len + r.len - 1
  vals["ref", ref_start_pad:ref_end_pad] <- ref
  
  # correlation spans entire padded domain
  vals["corr", ] <- xc
  
  # convert ref.start to padded feature placement
  fs <- match_row$ref.start + (ref_start_pad - 1)
  fe <- fs + f.len - 1
  
  vals["feat", fs:fe] <- feat
  
  use_mask <- !(is.na(vals["feat", ]) | is.na(vals["ref", ]))
  
  list(
    f.len = f.len,
    r.len = r.len,
    N = N,
    inds = list(
      feat_start = 1,
      feat_end   = f.len,
      ref_start  = match_row$ref.start,
      ref_end    = match_row$ref.end,
      use        = use_mask,
      lag        = match_row$lag,
      peak_loc   = match_row$lag + f.len
    ),
    vals = vals
  )
}

plot_conv_match <- function(match){

    fit.feat.ref <- fit_leastSquares(match$vals["feat",], 
                                     match$vals["ref",], plots = FALSE, scale.v2 = FALSE)#; fit.feat.ref$plot

    match$vals["feat",] <- match$vals["feat",] * fit.feat.ref$fit[2] + fit.feat.ref$fit[1]

    range.vals <- match$vals[c("feat", "ref"),] %>% range(na.rm = TRUE)
    match$vals["corr",] <- match$vals["corr",] %>% scale_between(range.vals[1], range.vals[2]) + range.vals[2]
    
    match$vals <- match$vals[c("corr", "ref", "feat"),]
    
    simplePlot(match$vals, 
               linecolor = c('black', 'gray','blue')) + geom_vline(xintercept=match$inds$peak_loc, color='blue', alpha=0.4)      
              
}

  