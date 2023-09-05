#' Fit a two-component model to a feature and a reference spectrum
#'
#' Batman fit seeks to stretch the feature signal up from the spectral minimum
#' (or zero) until peaks start exceeding spectral signal. Ebbels' BATMAN uses a 
#' more nuanced optimization strategy that penalizes overshoots; I'm just seeking 
#' a rough optimum. The algorithm is:
#' 
#' - assume the positions of the vectors are matched
#' - divide the spectral signal by the spectral signal by the feature signal. 
#'   - perfectly matched signals will have a ratio of 1 at each point
#'   - if feature is too low in a spot (could be raised before exceeding spectral
#'     signal), that point will have a high (>1) ratio
#'   - if feature is too high (the extreme case to be prevented), the ratio there 
#'     will be low (<1).
#'   - *** if the feature is zero, no information is retained. currently excluded.
#' - the point with the lowest ratio is the key point, i.e. the point which needs 
#'   will cross the spectral signal first if scaling up from zero; the point which
#'   should determine the scaling of the feature profile. Thus, a heuristic for 
#'   good fit is to make sure this point is perfectly matched. 
#'   - multiplying feature by a point's ratio will match the signals at that point
#' - multiply feature by that ratio
#' 
#' This algorithm makes some assumptions:
#' 
#' - no negative values in either vector (these will always be chosen as the min)
#' - where is zero? since scaling is linear, it shouldn't actually matter how far
#'   away the feature/spectra are from zero. However, unlike least squares fitting,
#'   the intercept is not optimized here. 
#' - for our purposes (fitting spectral signal), negative values make no sense. 
#'   We must assume that values are positive. Large negative spectral values will
#'   disrupt the shape of the fit if zero-bottoming (the features will be stretched
#'   and appear long and thin), but this would be 
#' - what to do with zero vals? exclude?
#' 
#' The idea behind excluding the lowest n% of points is to minimize the effect 
#' of noise on the fit. It doesn't matter if noise gets exceeded from time to time,
#' but we're most concerned about the highest feature points exceeding the spectral
#' signal being fit to. Thus we ignore a fraction of the lower feature points for
#' fitting purposes. However, once we have the ratio, we still apply it to all 
#' points in v1, and our fit quality is assessed using all non-NA points. 
#'
#' Notes:
#' - exclude.lowest around 50% seems to work well, and after a certain point (~0.2)
#'   raising this param doesn't appear to have much practical effect
#' - zero-bottoming the vectors helps to avoid having negative fits (e.g. from negative
#'   baseline effects in one of the vectors)
#' - v1 * ratio = v1 fit to v2
#'
#' @param feat a numeric vector of the feature spectrum
#' @param spec a numeric vector of the reference spectrum
#' @param exclude.lowest the fraction of the lowest values to exclude from the ratio calculation (default is 0.5)
#' @param ppm a numeric vector of the parts per million (ppm) values for the spectra (default is NULL)
#' @param plots logical, whether to generate plots (default is FALSE)
#'
#' @return a list with the following elements:
#' \itemize{
#' \item \code{feat.fit}: a numeric vector of the fitted feature spectrum
#' \item \code{spec.fit}: a numeric vector of the fitted reference spectrum
#' \item \code{keypoint}: the index of the keypoint where the ratio of the two spectra is minimum
#' \item \code{ratio}: the ratio of the feature and reference spectra at the keypoint
#' \item \code{intercept}: the intercept of the reference spectrum
#' \item \code{residuals}: a numeric vector of the residuals between the fitted feature and reference spectra
#' \item \code{plot}: a plot of the fitted feature and reference spectra (if plots = TRUE)
#' }
#'
#' @examples
#' # Example usage:
#' feat <- c(1, 2, 3, 4, 5)
#' spec <- c(1, 2, 3, 4, 5)
#' fit.batman(feat, spec)
#'
#' @importFrom ggplot2 geom_line scale_colour_gradientn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggnewscale new_scale_color
#' @importFrom colorRamps matlab.like2
#' @importFrom magrittr %>%
#'
#'
#' @export
fit_batman <- function(feat, spec, 
                       exclude.lowest = .5, # in practice I see convergence about here
                       ppm = NULL,
                       plots = FALSE){

  # Get the overlap between the two vects ####
    # can't compute where one is NAs
    
    use <- !is.na(feat + spec)

  # Zero-bottom the vects ####
    # negative values do bad things. Start both @ zero for scaling purposes
    # since we're fitting to v2, we'll want to add sub.v2 back to the fit after scaling
    
    sub.v1 <- min(feat[use])
    v1 <- feat[use] - sub.v1
    
    sub.v2 <- min(spec, na.rm = TRUE)
    v2 <- spec[use] - sub.v2
    
  # Exclude low intensity points ####
    v1.use.sortorder <- order(v1)
    use.exclude <- use # for keeping track of non-excluded, used points
    
    if (!is.null(exclude.lowest)){
      
      # Further subset the points 
      
      too.low <- v1.use.sortorder[1:floor(exclude.lowest * sum(use))]
      pairedRatios <- (v2[-too.low]) / (v1[-too.low]);
      
      use.exclude[too.low] <- FALSE # make sure we note which ones in 'use' were too low
      
    } else {
      
      too.low <- FALSE
      pairedRatios <- v2 / v1

    }
      # Plot ####
                  # plot(v2)
                  # points(too.low, v2[too.low], col = 'red')
                  # plot(v1)
                  # points(too.low, v1[too.low], col = 'red')
                  # plot(x = c(1:length(v1), 1:length(v2)), y = c(v1, v2))
                  # points(1:length(v1), v1, col = 'red')
                  # 
                  # simplePlot(rbind(v1,v2) %>% trim_sides)
  
  # Identify key point and ratio ####
  
    # Watch out for calculations involving values <= 0 - these can't offer much help
    
      pairedRatios[pairedRatios <= 0 | is.infinite(pairedRatios)] <- Inf 
      # never choose these points, but preserve indices in use
      
    keypoint <- which.min(pairedRatios)
      # use the index as-is to pull the ratio; currently indexes across TRUE vals in use.exclude
    ratio <- pairedRatios[keypoint]
    
    # now that the ratio is pulled, correct the keypoint index for the full vector
      keypoint <- which(use.exclude) %>% .[keypoint]
    
  # Calculate fit ####
    # v1.fit = (feat - sub.v1)*ratio + sub.v2 # zero-bottomed feature, scaled, then slid up to match v2.
      # = (feat * ratio) - (sub.v1 * ratio) + sub.v2
      # = (feat * ratio) + sub.v2 - (sub.v1*ratio)
      # = (feat * ratio) + sub.v2 - (sub.v1*ratio)
      # = (  m  *   x  ) + (           b         )
      
    v1.fit = (feat * ratio) + (sub.v2 - sub.v1*ratio) # zero-bottomed feature, scaled, then slid up to match v2.
    v2.fit = spec # subtracting both to zero causes high-scoring zero fits! 
    
  
  # Score fit ####
    
    # Want to count the overaccounted as negative..
    
    v1.f.s <- v1.fit-sub.v2
    v2.f.s <- v2.fit-sub.v2
    
    overshoot <- v1.f.s - v2.f.s
     pos <- overshoot>0
     pos.overshoot <- overshoot[pos]
     
    fraction.v2.accounted <- (sum(v1.f.s, na.rm = TRUE) - sum(pos.overshoot, na.rm = TRUE)) / sum(v2.f.s, na.rm = TRUE)

  # Plot ####
    g <- NULL
  
  if (plots){
    
    if (is.null(ppm)){ppm <- 1:length(v1.fit)}
    # Plot ####
      plot.new()
        plot(v2.fit, type = 'l',
             ylim = c(0, max(v2.fit)*1.2),
             col = 'black')
        points(which(use), v2.fit[use], col = 'black', cex = .2)
        lines(v1.fit, type = 'l', col = 'blue')
        points(which(use), v1.fit[use], col = 'blue', cex = .2)
          inds <- which(use) %>% .[too.low]
        points(c(inds, inds), c(v1.fit[inds], v2.fit[inds]), col = 'red')
        points(keypoint, v1.fit[keypoint], col = 'purple', pch = 10, cex = 5)
        g <- recordPlot()
      
      
      # g <- simplePlot(v2.fit,
      #                 xvect = ppm,
      #                 linecolor = "gray",
      #                 opacity = .9, 
      #                 linewidth = 1)
      # 
      # g <- g + new_scale_color() +
      #         geom_line(data = data.frame(vals = v1.fit,
      #                                     ppm = ppm,
      #                                     corr = rep(-1, length(v1.fit))),
      #                   na.rm = TRUE,
      #                   mapping = aes(x = ppm, y = vals, colour = corr),
      #                   linewidth = .5) +
      #         scale_colour_gradientn(colours = matlab.like2(10),
      #                                limits = c(-1, 1))
  }
    
  # return fit obj ####
  
    return(list(feat.fit = v1.fit,
                spec.fit = v2.fit,
                keypoint = keypoint,
                ratio = ratio,
                intercept = sub.v2 - sub.v1*ratio, # see above
                residuals = v2.fit - v1.fit,
                fraction.spec.accounted = fraction.v2.accounted,
                plot = g))
  
    
}

#' Fit least squares model between two signals and calculate statistics
#'
#'
#' @param v1 a numeric vector representing one of the signals
#' @param v2 a numeric vector representing the other signal
#' @param ppm a numeric vector representing the ppm values for the signals
#' @param plots a logical value indicating whether to produce a plot (default is \code{FALSE})
#'
#' @return A list with the following elements:
#' @import Metrics
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom magrittr %>%
#'
#' @export
fit_leastSquares <- function(v1 = NA, v2 = NA, ppm = NULL, plots = FALSE, scale.v2 = TRUE){
  
  fit <- 
        tryCatch(expr = {
          
          use <- !is.na(v1+v2)
                # simplePlot(rbind(v1,v2) %>% trim_sides)
                # simplePlot(rbind(v1[use],v2[use]) %>% trim_sides)
        
          # Check if there's even a need to fit. 
          # If not, skip it, scale if needed, and make a dummy fit.
            
            if (!any(use)){
              return(fit_obj())
            } else {
              if (all(v1[use] == v2[use])){
                  # If scaling v2, also do v1
                  if (scale.v2){
                    v2 <- v2 %>% scale_between
                    v1 <- v1 %>% scale_between #
                  }
                
                  res <- dummy_fit(v1[use], v2[use])
                
              } else {
                
                # Do the least squares fit
                
                  if (scale.v2){v2 <- v2 %>% scale_between}
                  
                  res <- lm(v2[use]~v1[use], na.action = na.exclude)
                  
              }
            }
            
          # Either way:
            v1.fit <- residuals <- rep(NA, length(v1))
          
          # rather than subset [use] for every calc in this section, don't calc the whole v1.fit yet
            v1.fit[use] <- res$fitted.values
        
            fit.neg <- v1.fit < 0
            # v1.fit[fit.neg] <- 0
            
            residuals[use] <- -res$residuals # (lm does spec-ref; we want ref-spec)
            resid.frac.fit <- residuals / v1.fit
          
            sum.residuals <- sum(abs(residuals), na.rm = T)
            mean.resid <- sum.residuals/sum(use)
            overshoot <- sum(residuals[residuals > 0], na.rm = T)
            undershoot <- sum(residuals[residuals < 0], na.rm = T)
            rmse <- Metrics::rmse(v1.fit[use], v2[use])
            # fraction.v2.unaccounted <- sum(v2[use] - res$fitted.values[use]) / sum(v2[use], na.rm = T)
            
            # Residual portions as fraction of each signal
              res.pos <- residuals
                res.pos[res.pos*v1.fit <= 0] <- NA
                # simplePlot(rbind(res.pos,res.neg))
              
                pos.res.pct.feat <- 1 -
                                     (v1.fit - res.pos) /
                                      max(v1.fit, na.rm = T)
                pos.res.pct.feat[pos.res.pct.feat > 1] <- 1
                # plot(pos.res.pct.feat)
                
              res.neg <- residuals
                res.neg[res.neg*v2 >= 0] <- NA
                res.neg <- abs(res.neg) # should this be?
                # simplePlot(rbind(res.pos,res.neg))
              
                neg.res.pct.spec <- 1 -
                                     (v2 - res.neg) /
                                      max(v2, na.rm = T)
                neg.res.pct.spec[neg.res.pct.spec > 1] <- 1
                # neg.res.pct.spec <- neg.res.pct.spec ^ 2
                # # plot(neg.res.pct.spec)
                # spec.missing <- v2 - neg.res.pct.spec * v2
                # Express as a pct of relevant spec signal...
                # simplePlot(rbind(v1.fit, newspec, newspec, newspec))
                # simplePlot(rbind(v2, newspec, newspec, newspec))
        
          # Now that the metrics have been calculated, replace v1.fit with the whole v1.fit
            v1.fit <- v1 * res$coefficients[2] + res$coefficients[1]
                # simplePlot(rbind(v1.fit,v2))
          
          
          g <- NULL
          # # Old plots
          #   if (plots == "simple"){
          #     g <- simplePlot(rbind(v1.fit, v2,v2,v2) %>% trim_sides) # ,-res$residuals
          #     # g %>% plot
          #   }
          
          if (plots){
            mat <-  rbind(v1.fit, v2)
            mat.cols <- mat %>% trim_sides(out = "inds")
            mat <- mat[,mat.cols]
            
            if (is.null(ppm)){ppm <- 1:ncol(mat)} else {ppm <- ppm[mat.cols]}
            
              g <- simplePlot(mat[2,],
                              linecolor = "gray",
                              xvect = ppm,
                              opacity = .9, 
                              linewidth = 1) # ,-res$residuals
        
              
              g <- g + new_scale_color() +
                      geom_line(data = data.frame(vals = mat[1,],
                                                  ppm = ppm,
                                                  corr = rep(-1, length(ppm))), 
                                mapping = aes(x = ppm, y = vals, colour = corr),
                                linewidth = .5,
                                na.rm = T) +
                      scale_colour_gradientn(colours = matlab.like2(10),
                                             limits = c(-1, 1))
          }
        
         list(feat.fit = v1.fit, # this is now the updated, full v1.fit
              spec.fit = v2,
              ratio = NA,
              overshoot = overshoot,
              undershoot = undershoot,
              fit = res$coefficients,
              residuals = residuals,
              resid.frac.fit = resid.frac.fit,
              sum.residuals = sum.residuals,
              mean.residual = mean.resid,
              pts.matched = sum(use),
              rmse = rmse,
              overshoot.pct.feat = pos.res.pct.feat,
              spec.missing.pct = sum(neg.res.pct.spec, na.rm = T),
              fit.neg = fit.neg,
              # fraction.spec.unaccounted = (v2 - res$fitted.values)/sum(v2, na.rm = T),
              plot = g)
        }, 
        error = function(cond){
          return(fit_obj())
        })
  
}

dummy_fit <- function(v1,v2){
  list(coefficients = c(0,1),
       residuals = rep(0, length(v1)),
       fitted.values = v1)
}

fit_obj <- function(){
  list(
        feat.fit = NA, # this is now the updated, full v1.fit
        spec.fit = NA,
        ratio = NA,
        overshoot = NA,
        undershoot = NA,
        fit = NA,
        residuals = NA,
        resid.frac.fit = NA,
        sum.residuals = NA,
        mean.residual = NA,
        pts.matched = NA,
        rmse = NA,
        overshoot.pct.feat = NA,
        spec.missing.pct = NA,
        fit.neg = NA,
        plot = NULL
  )
}
  
