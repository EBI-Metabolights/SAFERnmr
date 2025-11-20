#' Match a feature against a reference using FFT-based convolution and correlation
#' picks the top max.hits peaks in the convolution, then calculates correlation and
#' pvalue of the PCC at those points in the ref.
#' Returns a dataframe with the xcorr information.
#'
#' @param f.num numeric: Feature number
#' @param r.num numeric: Reference number
#' @param feat numeric: Vector of feature intensities
#' @param ref numeric: Vector of reference intensities
#' @param feat.ft.c numeric: FFT of feature intensities, conjugated
#' @param ref.ft numeric: FFT of reference intensities
#' @param pad.size numeric: Size of the padding applied during FFT-based convolution
#' @param max.hits numeric: Maximum number of candidate lags to consider
#' @param r.thresh numeric: Correlation threshold for considering a match
#' @param p.thresh numeric: P-value threshold for considering a match
#'
#' @return A data frame with the following columns:
#'   - feat: Feature number or identifier
#'   - ref: Reference number or identifier
#'   - lag: The lag that produces the highest correlation
#'   - rval: The correlation coefficient
#'   - pval: The p-value of the correlation
#'   - pts.matched: Number of data points used in the correlation
#'   - pts.feat: Number of data points in the feature vector
#'   - feat.start: The starting index of the feature vector used in the correlation
#'   - feat.end: The ending index of the feature vector used in the correlation
#'   - ref.start: The starting index of the reference vector used in the correlation
#'   - ref.end: The ending index of the reference vector used in the correlation
#'
#' @importFrom fftw FFT
#' @importFrom magrittr %>%
#'
#' @export
feature_match2ref_slim_pcc <- function(f.num, r.num, feat, ref,
                                   feat.ft.c, ref.ft, 
                                   max.hits = 100, 
                                   r.thresh = 0.8, p.thresh = 0.01,
                                   trim = FALSE){
  
                                  ## To debug, run these:
                                  # pad.size <- length(feat)-1
                                  # feat.ft.c <- feat.padded.ft.c

    # Just do Pearson-based xcorr ####
      
      xc.res <- compute_pearson_sliding_overlap(feat, ref)
      xc <- xc.res$xcorr
      
    # Get maxima (candidate lags) ####
      
      lmxs <- localMaxima(xc)
      
      # plot(xc, type='l')
      # points(lmxs,xc[lmxs], col='blue')
      
      not.partial <- lmxs >= xc.res$range.complete.overlap[1] & lmxs <= xc.res$range.complete.overlap[2]
      lmxs <- lmxs[not.partial]
      
      # plot(fft.res$ncc, type='l')
      # points(lmxs,fft.res$ncc[lmxs], col='blue')
      
    # Sort maxima
      
      lmxs <- lmxs[xc[lmxs] >= r.thresh]
      if (length(lmxs) == 0){
        return(NULL)
      }
      peaks <- lmxs[order(xc[lmxs], decreasing = TRUE)] # sort by xcorr peak height
        # xc[lmxs] %>% sort(decreasing = T)
        
      # plotly::plot_ly(data.frame(x=1:length(ncc), y=ncc), x = ~x, y= ~y)
      
    # Restrict lags to those not inside padding (padding is really just for end 
    # effects in the FT, not for actual comparison. Perhaps it's necessary to 
    # pad the matrix twice?). Also only take the top n hits (to keep calculations
    # reasonable).
    
      # peaks are the max posns in ncc[valid]
      # these lags are actually peaks + pad.size+1
      
      peaks <- peaks[1:(min(max.hits, length(peaks)))]
      lags <- xc.res$lags[peaks]
      
      # plot(xc, type='l')
      # points(peaks,xc[peaks], col='red')
      #   max(xc, na.rm=TRUE)
      
    # Loop though candidate lags and evaluate fit at each one ####
      
      feat <- as.double(feat)
      ref <- as.double(ref)
      
      
      mapped <- map_ref_xcorr(xc.res, ref)
        # a<-mapped$vals
        
      
      fits <- lapply(lags, function(lag){
        # Calculate corr at each shift
        # Note: to look at these, uncomment [PLOT] sections:
        
            # # [PLOT] # # # # # # # #
            # i <- 0
            # i <- i + 1
            # lag<- lags[i]
            # # [PLOT] # # # # # # # #
            
            match <- map_feat_xcorr(mapped, feat, lag)

            # # [PLOT] # # # # # # # #
            #   g <- plot_conv_match(match)
            #   g + ggtitle(xc[match$inds$peak_loc] %>% round(4))
            # # [PLOT] # # # # # # #
            
            return(data.frame(feat.start = match$inds$feat_start,
                              feat.end = match$inds$feat_end,
                              ref.start = match$inds$ref_start,
                              ref.end = match$inds$ref_end,
                              pts.matched = sum(match$inds$use),
                              rval = xc[match$inds$peak_loc]))
      }) %>% do.call(rbind,.)

        # In case of infinite rvals, set to zero:
          fits$rval[is.infinite(fits$rval)] <- 0
          r.pass <- fits$rval > r.thresh
        
        # Calculate pvals using t-distribution
          a <- -abs(fits$rval * sqrt( (fits$pts.matched-2) /(1-fits$rval^2)))
          pvals <- 2*pt(a,(fits$pts.matched-2))
          p.pass <- TRUE #pvals < p.thresh

      # average fit intensity as a fraction of the feature signal (want to fit parts that are dominant)
      
      r.p.pass <- which(r.pass & p.pass)
      
      if (!any(r.p.pass)){return(NULL)}
      
      
     matches.ranked <- order(pvals[r.p.pass]) %>% r.p.pass[.]
     if (length(matches.ranked) == 0) {
        return(NULL)
     }
     matches <- data.frame( feat = f.num,
                            ref = r.num,
                            lag = lags[matches.ranked],
                            rval = fits$rval[matches.ranked],
                            pval = pvals[matches.ranked],
                            pts.matched = fits$pts.matched[matches.ranked],
                            pts.feat = length(feat),
                            feat.start = fits$feat.start[matches.ranked],
                            feat.end = fits$feat.end[matches.ranked],
                            ref.start = fits$ref.start[matches.ranked],
                            ref.end = fits$ref.end[matches.ranked],
                            row.names = NULL)
           
  # Record results
    return(matches)

     
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
  compute_pearson_sliding_overlap <- function(feat, ref) {
    
    # 3. Convert everything to plain double early
    feat <- as.double(feat)
    ref  <- as.double(ref)
    
    f.len <- length(feat)
    r.len <- length(ref)
    N <- f.len + r.len - 1
    
    # lag indexing identical to FFT Pearson / NCC
    lags <- -(f.len - 1):(r.len - 1)
    
    pearson <- rep(NA_real_, N)
    
    # place reference at canonical FFT alignment:
    ref_start <- f.len
    ref_end   <- f.len + r.len - 1
    
    # 1. Precompute base vals matrix ONCE per reference
    base_vals <- matrix(NA_real_, 3, N)
    rownames(base_vals) <- c("feat", "ref", "corr")
    
    base_vals["ref", ref_start:ref_end] <- ref
    
    # 2. Precompute ref mask once (never changes)
    ref_mask <- !is.na(base_vals["ref", ])
    
    # Reuse a single vals matrix; only overwrite the feat row region
    vals <- base_vals
    last_feat_start <- NULL
    last_feat_end   <- NULL
    
    # 4. Avoid S3 dispatch by taking a local handle to stats::cor
    .cor <- stats::cor
    
    # ---------------------------------------------------------------
    # Slide the feature across the reference, reproduce your overlay logic
    # ---------------------------------------------------------------
    for (i in seq_along(lags)) {
      
      lag <- lags[i]
      
      feat_start <- ref_start + lag
      feat_end   <- feat_start + f.len - 1
      
      # Only compute if feature is within the padded frame
      if (feat_start < 1 || feat_end > N)
        next
      
      # clear previous feature region only (no full copy)
      if (!is.null(last_feat_start)) {
        vals["feat", last_feat_start:last_feat_end] <- NA_real_
      }
      
      # write new feature region
      vals["feat", feat_start:feat_end] <- feat
      
      last_feat_start <- feat_start
      last_feat_end   <- feat_end
      
      # 2. Vectorized overlap mask using precomputed ref_mask
      feat_mask <- !is.na(vals["feat", ])
      use_mask  <- feat_mask & ref_mask
      
      # Need at least 3 valid points
      if (sum(use_mask) < 3)
        next
      
      pearson[i] <- suppressWarnings(
        .cor(
          vals["feat", use_mask],
          vals["ref",  use_mask],
          use    = "pairwise.complete.obs",
          method = "pearson"
        )
      )
    }
    
    range.complete.overlap <- c(f.len, r.len)
    
    list(
      feat = feat,
      ref = ref,
      xcorr = pearson,   # matches FFT version field name
      lags = lags,
      feat_len = f.len,
      ref_len = r.len,
      range.complete.overlap = range.complete.overlap
    )
  }

  # compute_pearson_sliding_overlap <- function(feat, ref) {
  #   
  #   # keep NA exactly as-is
  #   feat <- as.numeric(feat)
  #   ref  <- as.numeric(ref)
  #   
  #   f.len <- length(feat)
  #   r.len <- length(ref)
  #   N <- f.len + r.len - 1
  #   
  #   # lag indexing identical to FFT Pearson / NCC
  #   lags <- -(f.len - 1):(r.len - 1)
  #   
  #   pearson <- rep(NA_real_, N)
  #   
  #   # place reference at canonical FFT alignment:
  #   ref_start <- f.len
  #   ref_end   <- f.len + r.len - 1
  #   
  #   # Base template for vals matrix (same structure as your mapping code)
  #   base_vals <- matrix(NA_real_, 3, N)
  #   rownames(base_vals) <- c("feat", "ref", "corr")
  #   
  #   base_vals["ref", ref_start:ref_end] <- ref
  #   
  #   # ---------------------------------------------------------------
  #   # Slide the feature across the reference, reproduce your overlay logic
  #   # ---------------------------------------------------------------
  #   for (i in seq_along(lags)) {
  #     
  #     lag <- lags[i]
  #     
  #     feat_start <- ref_start + lag
  #     feat_end   <- feat_start + f.len - 1
  #     
  #     # Only compute if feature is within the padded frame
  #     if (feat_start < 1 || feat_end > N)
  #       next
  #     
  #     vals <- base_vals
  #     vals["feat", ] <- NA
  #     vals["feat", feat_start:feat_end] <- feat
  #     
  #     # EXACT same mask as your plotting code
  #     use_mask <- !(is.na(vals["feat", ]) | is.na(vals["ref", ]))
  #     
  #     # Need at least 3 valid points
  #     if (sum(use_mask) < 3)
  #       next
  #     
  #     # EXACTLY your desired Pearson logic
  #     pearson[i] <- suppressWarnings(
  #       cor(
  #         vals["feat", use_mask],
  #         vals["ref",  use_mask],
  #         use    = "pairwise.complete.obs",
  #         method = "pearson"
  #       )
  #     )
  #   }
  #   
  #   range.complete.overlap <- c(f.len, r.len)
  #   
  #   list(
  #     feat = feat,
  #     ref = ref,
  #     xcorr = pearson,   # matches FFT version field name
  #     lags = lags,
  #     feat_len = f.len,
  #     ref_len = r.len,
  #     range.complete.overlap = range.complete.overlap
  #   )
  # }
