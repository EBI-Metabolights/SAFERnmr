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
feature_match2ref_slim <- function(f.num, r.num, feat, ref,
                                   feat.ft.c, ref.ft, 
                                   pad.size,
                                   max.hits = 5, 
                                   r.thresh = 0.8, p.thresh = 0.01,
                                   trim = FALSE){
  
                                  ## To debug, run these:
                                  # pad.size <- length(feat)-1
                                  # feat.ft.c <- feat.padded.ft.c
                                  
    # Do the FFT-based conv/xcorr ####
      
      r.conv <- (feat.ft.c*ref.ft) %>% fftw::FFT(.,inverse = TRUE) %>% Re %>% c

    # Get maxima (candidate lags) ####
      
      lmxs <- localMaxima(r.conv)

    # Sort maxima
      lags <- lmxs[order(r.conv[lmxs], decreasing = TRUE)] # sort by xcorr peak height
      # plot(r.conv, type='l')
      # points(lmxs,r.conv[lmxs], col='blue')
      # plotly::plot_ly(data.frame(x=1:length(r.conv), y=r.conv), x = ~x, y= ~y)
      
    # Restrict lags to those not inside padding (padding is really just for end 
    # effects in the FT, not for actual comparison. Perhaps it's necessary to 
    # pad the matrix twice?). Also only take the top n hits (to keep calculations
    # reasonable).
    
      lags <- lags[lags>=pad.size] %>% .[1:max.hits] 
      # plot(r.conv, type='l')
      # points(lags,r.conv[lags], col='red')
      
    # Loop though candidate lags and evaluate fit at each one ####
      
      # inds.trim.feat <- trim_sides(feat, out = "inds")
      feat.inds <- 1:length(feat)
      feat <- t(c(feat)) %>% rev # *** must reverse
      ref <- t(c(ref))
      
      fits <- lapply(lags, function(lag){
        # Calculate corr at each shift
        # Note: to look at these, uncomment [PLOT] sections:
        
            # # [PLOT] # # # # # # # #
            # i <- 0
            # i <- i + 1
            # lag<- lags[i]
            # # [PLOT] # # # # # # # #
            
            ref.pos <- lag - pad.size - feat.inds
            
            # # [PLOT] # # # # # # # #
            #   g <- plot_conv_match(feat, ref, r.conv, ref.pos, pad.size, lag)
            # # [PLOT] # # # # # # # #
            
  
          # Get the overlapping, non-NA values of ref and feat
          
            valid <- which(ref.pos >= 1 & ref.pos <= length(ref))
            
            if (length(valid) < 3) return(NULL)
            
            feat.pos <- feat.inds[valid]
            ref.pos  <- ref.pos[valid]
            use <- !is.na(feat[feat.pos] + ref[ref.pos])

            # use <- !is.na(feat + ref[ref.pos])
            # rbind(feat[feat.pos]%>% scale_between(),ref[ref.pos]%>% scale_between())  %>% simplePlot

            # Make sure there are enough points to do a correlation:
            if (sum(use) < 3){return(NULL)}
            r <- suppressWarnings( 
                                   cor(feat[feat.pos[use]], 
                                       ref[ref.pos[use]],
                                       use = "pairwise.complete.obs",
                                       method = "pearson")            
                                   )       
            
            # # [PLOT] # # # # # # # #
            #   g + ggtitle(r %>% round(4))
            # # [PLOT] # # # # # # # #
            
            return(data.frame(ref.start = min(ref.pos),
                              ref.end = max(ref.pos),
                              pts.matched = sum(use),
                              rval = r))
      }) %>% do.call(rbind,.)

        # In case of infinite rvals, set to zero:
          fits$rval[is.infinite(fits$rval)] <- 0
          r.pass <- fits$rval > r.thresh
        
        # Calculate pvals using t-distribution
          a <- -abs(fits$rval * sqrt( (fits$pts.matched-2) /(1-fits$rval^2)))
          pvals <- 2*pt(a,(fits$pts.matched-2))
          p.pass <- pvals < p.thresh

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
                            pts.feat = length(use),
                            feat.start = 1,
                            feat.end = length(feat),
                            ref.start = fits$ref.start[matches.ranked],
                            ref.end = fits$ref.end[matches.ranked],
                            row.names = NULL)
           
  # Record results
    return(matches)

}

plot_conv_match <- function(feat, ref, r, ref.pos, pad.size, lag){
  
            r.inds <- 1:length(r)
            feat.inds.in.r <- ref.pos + pad.size
            ref.inds.in.r <- pad.size + 1:length(ref)
  
            unified.inds <- c(r.inds, feat.inds.in.r, ref.inds.in.r) %>% range %>% fillbetween
            
            fit.feat.ref <- fit_leastSquares(feat, ref[ref.pos], plots = T, scale.v2 = FALSE); fit.feat.ref$plot
            # fit.feat.ref <- fit_batman(feat, ref[ref.pos], plots = T); fit.feat.ref$plot
            
            
            feat.filled <- ref.filled <- r.filled <- matrix(NA, 1, length(unified.inds))

            f <- fit.feat.ref$fit
            fr <- fit.r.ref$fit
            
            feat.filled[feat.inds.in.r]<- (feat) * f[2] + f[1]
            ref.filled[ref.inds.in.r]<- ref
            
            allvals <- c(feat.filled, ref.filled)
            range.vals <- range(allvals, na.rm = TRUE)
            r <- r %>% scale_between(range.vals[1], range.vals[2])
            r <- r + range.vals[2]
            
              simplePlot(rbind(ref.filled,
                               r,
                               feat.filled), 
                         linecolor = c('black','gray', 'blue')) + geom_vline(xintercept=lag, color='blue')      
              
}
# r <- convolve(feat.long,rev(ref.long), conj = T, type = c("circular", "open", "filter"))
# lag <- which.max(r)
  
# simplePlot(res)
    # plot(t(r))
    #   points(lag, r[lag], col = 'red', cex = 1)
    #   ref.hit <- ref.long[(1:length(feat)) + lag - r.pad.size]
    #   fit <- fit_leastSquares(feat, ref.hit, plots = T)
    #     fit$plot %>% plot
      