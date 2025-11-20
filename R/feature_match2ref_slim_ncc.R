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
feature_match2ref_slim_ncc <- function(f.num, r.num, feat, ref,
                                   feat.ft.c, ref.ft, 
                                   pad.size,
                                   max.hits = 5, 
                                   r.thresh = 0.8, p.thresh = 0.01,
                                   trim = FALSE){
  
                                  ## To debug, run these:
                                  # pad.size <- length(feat)-1
                                  # feat.ft.c <- feat.padded.ft.c

    # Do the FFT-based conv/xcorr ####
      
      # r.conv <- (feat.ft.c*ref.ft) %>% fftw::FFT(.,inverse = TRUE) %>% Re %>% c
      
      xc.res <- compute_ncc_fft(feat, ref)
      
    # Get maxima (candidate lags) ####
      
      lmxs <- localMaxima(xc.res$xcorr)
      xc <- xc.res$xcorr
      
      # plot(xc.res$xcorr, type='l')
      # points(lmxs,fft.res$ncc[lmxs], col='blue')
      
      not.partial <- lmxs >= xc.res$range.complete.overlap[1] & lmxs <= xc.res$range.complete.overlap[2]
      lmxs <- lmxs[not.partial]
      # plot(fft.res$ncc, type='l')
      # points(lmxs,fft.res$ncc[lmxs], col='blue')
      
    # Sort maxima
    
      peaks <- lmxs[order(xc.res$xcorr[lmxs], decreasing = TRUE)] # sort by xcorr peak height
      
      # plotly::plot_ly(data.frame(x=1:length(ncc), y=ncc), x = ~x, y= ~y)
      
    # Restrict lags to those not inside padding (padding is really just for end 
    # effects in the FT, not for actual comparison. Perhaps it's necessary to 
    # pad the matrix twice?). Also only take the top n hits (to keep calculations
    # reasonable).
    
      # peaks are the max posns in ncc[valid]
      # these lags are actually peaks + pad.size+1
      
      # peaks <- peaks[1:max.hits]
      lags <- xc.res$lags[peaks]
      plot(xc.res$xcorr, type='l')
      points(peaks,xc.res$xcorr[peaks], col='red')
       
    # Loop though candidate lags and evaluate fit at each one ####
      
      feat.inds <- 1:length(feat)
      feat <- t(c(feat))
      ref <- t(c(ref))
      
      mapped <- map_ref_xcorr(xc.res, ref)
      
      fits <- lapply(lags, function(lag){
        # Calculate corr at each shift
        # Note: to look at these, uncomment [PLOT] sections:
        
            # # [PLOT] # # # # # # # #
            i <- 0
            i <- i + 1
            lag<- lags[i]
            # # [PLOT] # # # # # # # #
            
            match <- map_feat_xcorr(mapped, feat, lag)
            
            # [PLOT] # # # # # # # #
              g <- plot_conv_match(match)
            # [PLOT] # # # # # # # #
            
  
          # Get the overlapping, non-NA values of ref and feat
          
            # Make sure there are enough points to do a correlation:
            if (sum(match$inds$use) < 3){return(NULL)}
            r <- suppressWarnings( 
                                   cor(match$vals["feat",match$inds$use], 
                                       match$vals["ref",match$inds$use],
                                       use = "pairwise.complete.obs",
                                       method = "pearson")            
                                   )       
            
            # [PLOT] # # # # # # # #
              g + ggtitle(r %>% round(4))
            # [PLOT] # # # # # # # #
            
            return(data.frame(feat.start = match$inds$feat.start,
                              feat.end = match$inds$feat.end,
                              ref.start = match$inds$ref.start,
                              ref.end = match$inds$ref.end,
                              pts.matched = sum(match$inds$use),
                              rval = r))
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
      
      lag <- match$inds$lag
      peak_loc <- lag + match$f.len

      fit.feat.ref <- fit_leastSquares(match$vals["feat",], 
                                       match$vals["ref",], plots = FALSE, scale.v2 = FALSE)#; fit.feat.ref$plot

      match$vals["feat",] <- match$vals["feat",] * fit.feat.ref$fit[2] + fit.feat.ref$fit[1]

      range.vals <- match$vals[c("feat", "ref"),] %>% range(na.rm = TRUE)
      match$vals["corr",] <- match$vals["corr",] %>% scale_between(range.vals[1], range.vals[2]) + range.vals[2]
      
      match$vals <- match$vals[c("corr", "ref", "feat"),]
      
      simplePlot(match$vals, 
                 linecolor = c('black', 'gray','blue')) + geom_vline(xintercept=peak_loc, color='blue')      
                
  }

