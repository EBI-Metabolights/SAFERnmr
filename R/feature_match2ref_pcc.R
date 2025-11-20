feature_match2ref_pcc <- function(f.num, r.num, feat, ref,
                                  r.thresh = 0.8,
                                  max.hits = 100) {
  
  # 1. Compute all masked Pearson correlations
  xc.res <- compute_pearson_sliding_overlap(feat, ref)
  xc <- xc.res$xcorr
  lags <- xc.res$lags
  pts  <- xc.res$pts
  
  # 2. Find local maxima
  peaks <- localMaxima(xc)
  
  # restrict to complete overlap region
  valid <- peaks >= xc.res$range.complete.overlap[1] &
           peaks <= xc.res$range.complete.overlap[2]
  peaks <- peaks[valid]
  
  if (length(peaks) == 0) return(NULL)
  
  # restrict by r threshold
  peaks <- peaks[xc[peaks] >= r.thresh]
  if (length(peaks) == 0) return(NULL)
  
  # sort by correlation height
  peaks <- peaks[order(xc[peaks], decreasing = TRUE)]
  peaks <- peaks[seq_len(min(max.hits, length(peaks)))]
  
  # Vectorized extraction
  lag_vec <- lags[peaks]
  rval    <- xc[peaks]
  npts    <- pts[peaks]
  
  # Vectorized p-value
  tstat <- abs(rval * sqrt((npts - 2) / pmax(1e-12, (1 - rval^2))))
  pvals <- 2 * pt(-tstat, df = npts - 2)
  
  f.len <- xc.res$feat_len
  r.len <- xc.res$ref_len
  
  ref_start <- f.len
  ref_end   <- f.len + r.len - 1
  
  feat_start <- ref_start + lag_vec
  feat_end   <- feat_start + f.len - 1
  
  data.frame(
    feat = f.num,
    ref  = r.num,
    lag  = as.numeric(lag_vec),
    rval = as.numeric(rval),
    pval = as.numeric(pvals),
    pts.matched = as.numeric(npts),
    pts.feat = f.len,
    feat.start = as.numeric(feat_start),
    feat.end   = as.numeric(feat_end),
    ref.start = ref_start,
    ref.end   = ref_end,
    row.names = NULL,             
    stringsAsFactors = FALSE)
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
  pts_matched <- integer(N)
  
  # Reference occupies padded indices:
  ref_start <- f.len
  ref_end   <- f.len + r.len - 1
  
  # ---------------------------------------------------------
  # PURE VECTOR LOGIC
  # feat_start = ref_start + lag
  # feat_end   = feat_start + f.len - 1
  # ---------------------------------------------------------
  feat_start <- ref_start + lags
  feat_end   <- feat_start + f.len - 1
  
  # Outer loop is already minimal (just slicing + mask)
  for (i in seq_len(N)) {
    
    fs <- feat_start[i]
    fe <- feat_end[i]
    
    # enforce bounds
    if (fs < 1 || fe > N)
      next
    
    # Compute overlap
    # local feat indices
    local_f <- 1:f.len
    
    # global indices of feat in padded frame
    global_f <- fs:fe
    
    # valid ref region
    valid_ref <- (global_f >= ref_start) & (global_f <= ref_end)
    
    if (!any(valid_ref))
      next
    
    # reference local indices
    ref_idx <- global_f[valid_ref] - ref_start + 1
    feat_idx <- local_f[valid_ref]
    
    x <- feat[feat_idx]
    y <- ref[ref_idx]
    
    use_mask <- !(is.na(x) | is.na(y))
    n <- sum(use_mask)
    
    pts_matched[i] <- n
    
    if (n < 3)
      next
    
    x <- x[use_mask]
    y <- y[use_mask]
    
    # Fast correlation formula
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
    pts = pts_matched,
    feat_len = f.len,
    ref_len = r.len,
    range.complete.overlap = c(f.len, r.len)
  )
}

convert_for_plot <- function(feat, ref, xc, match_row) {
  
  f.len <- length(feat)
  r.len <- length(ref)
  N <- f.len + r.len - 1
  
  vals <- matrix(NA_real_, 3, N)
  rownames(vals) <- c("corr", "ref", "feat")
  
  ref_start <- f.len
  ref_end   <- f.len + r.len - 1
  
  vals["ref", ref_start:ref_end] <- ref
  vals["corr", ] <- xc
  
  fs <- match_row$feat.start
  fe <- match_row$feat.end
  vals["feat", fs:fe] <- feat
  
  use_mask <- !(is.na(vals["feat", ]) | is.na(vals["ref", ]))
  
  list(
    f.len = f.len,
    r.len = r.len,
    N = N,
    inds = list(
      feat_start = fs,
      feat_end   = fe,
      ref_start  = ref_start,
      ref_end    = ref_end,
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

  