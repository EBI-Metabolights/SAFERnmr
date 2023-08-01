#' Fit a two-component model to a feature and a reference spectrum
#'
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
                       exclude.lowest = .5,
                       ppm = NULL,# in practice I see convergence about here
                       plots = FALSE){
   # v1 * ratio = v1 fit to v2
   # 
   # feat <- featureStack[2002,]
   # spec <- featureStack[2000,] / 2
   # feat <- featureStack[comb[1],]
   # spec <- featureStack[comb[2],]
  
  use <- which(!is.na(feat + spec))
  
  sub.v1 <- min(feat[use])
  sub.v2 <- min(spec[use])          
                
  v1 <- feat[use] - sub.v1
  v2 <- spec[use] - sub.v2
  
  v2.use.sortorder <- order(v2)
  
  if (!is.null(exclude.lowest)){
    exclude.pts <- v2.use.sortorder[1:floor(exclude.lowest * length(use))]
    pairedRatios <- v2[-exclude.pts] / v1[-exclude.pts];
  } else {
    pairedRatios <- v2 / v1;
  }
  
  # plot(v2)
  # points(exclude.pts, v2[exclude.pts], col = 'red')

  # simplePlot(rbind(v1,v2) %>% trim_sides)
  
  
  keypoint <- which.min(pairedRatios)
  ratio <- pairedRatios[keypoint]
  
  v1.fit = (feat - sub.v1)*ratio + sub.v2
  v2.fit = spec #- sub.v2
  
  # fraction.v2.unaccounted <- sum(v2.fit[use]  -  v1.fit[use], na.rm = TRUE) / sum(v2.fit[use], na.rm = TRUE)
  g <- NULL
  # simplePlot(rbind(v1 * ratio,v2) %>% trim_sides)
  if (plots){
    
    if (is.null(ppm)){ppm <- 1:length(v1.fit %>% t %>% trim_sides)}
    
      v2.inds <- v2 %>% t %>% trim_sides(out = "inds")
      g <- simplePlot(v2[v2.inds],
                      xvect = ppm[v2.inds],
                      linecolor = "gray",
                      opacity = .9, 
                      linewidth = 1)

      v1.inds <- v1 %>% t %>% trim_sides(out = "inds")
      g <- g + new_scale_color() +
              geom_line(data = data.frame(vals = v1.fit[v1.inds],
                                          ppm = ppm[v1.inds],
                                          corr = rep(-1, length(v1.inds))), 
                        mapping = aes(x = ppm, y = vals, colour = corr),
                        linewidth = .5) +
              scale_colour_gradientn(colours = matlab.like2(10),
                                     limits = c(-1, 1))
  }
  # 
  
  return(list(feat.fit = v1.fit,
              spec.fit = v2.fit,
              keypoint = use[keypoint],
              ratio = ratio,
              intercept = sub.v2,
              residuals = v2.fit - v1.fit,
              # mean.feat.fit = mean(v1.fit, na.rm = TRUE),
              # mean.spec.fit = mean(v2.fit, na.rm = TRUE),
              # fraction.spec.unaccounted = fraction.v2.unaccounted,
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
fit_leastSquares <- function(v1, v2, ppm = NULL, plots = FALSE, scale.v2 = TRUE){
  
  require(Metrics)
  fit <- 
        tryCatch(expr = {
          use <- !is.na(v1+v2)
                # simplePlot(rbind(v1,v2) %>% trim_sides)
                # simplePlot(rbind(v1[use],v2[use]) %>% trim_sides)
        
          # Check if there's even a need to fit. 
          # If not, skip it, scale if needed, and make a dummy fit.
            
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
  
