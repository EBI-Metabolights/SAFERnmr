#' Back-fit reference features to a subset of spectra
#'
#' @param m.inds A vector of indices corresponding to the matches between reference features and spectra. Necessary to allow for subsetting of matches (i.e. because of filtering)
#' @param fits.feature A list of fitted features (to reference spectra)
#' @param match.info A data frame containing information about the matches between features and reference spectra
#' @param feature List containing feature intensities and positions
#' @param xmat Full spectral matrix from which the features were derived
#' @param ppm A vector containing the ppm axis of the spectra
#' @param plots A logical value indicating whether to generate plots of the fits and back-fit feasibility scores - only use for a couple of backfits at a time as plots will double the weight of the result
#'
#' @return A list containing the back-fits and back-fit feasibility scores for each spectrum in the subset,
#'         as well as a grid plot of the fits and back-fit feasibility scores.
#'
#' @importFrom ggnewscale new_scale_color
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#'
#'
#' @export
backfit_ref.feats.2.subset.specs <- function(m.inds, fits.feature, match.info, 
                                    feature, 
                                    xmat, ppm, 
                                    plots = F){

      backfits <- lapply(m.inds, function(m)
      {
        # print(m)
        ############# For each match ###############
        # Get data ####
          # m <- 4
          # m <- m.inds[3]
          
          fit <- fits.feature[[m]]
          ref.region <- match.info[m, c("ref.start","ref.end")] %>% as.numeric
          feat <- match.info$feat[m]
          
          # If sfe has not been done:
            spec.cols <- match.info[m, c("feat.start","feat.end")] %>% as.numeric %>% feature$position[feat,.] %>% fillbetween
            thisref.ppm <- ref.region %>% fillbetween %>% ppm[.]
          
          # If sfe HAS been done, wait to do this at the spectrum level using lags
            
          # fit %>% plot.fit(., type = "simple", ppm = thisref.ppm) %>% plot
          # feat <- f.subset[match.info$feat[m]]
          
          lags <- feature$sfe[[feat]]$lags
          ss <- feature$sfe[[feat]]$feat$ss
           
           
           # lapply(1:length(ss), function(s){
           #   pos <- spec.cols + lags[s]
           #   vals <- xmat[ss[s], pos]
           #   return(vals)
           # }) %>% do.call(rbind,.) %>% simplePlot
           
          # spec.cols <- feature$position[feat,] %>% t %>% trim.sides(out = "inds") %>% complete.indsVect
                        # simplePlot(xmat[ss, spec.cols], xdir = "n")
                        # simplePlot(fit$spec.fit, xdir = "n", linecolor = "black")
                        # simplePlot(fit$feat.fit, xdir = "n", linecolor = "blue")
            
        # Calculate fits to subset spectra: ####
          fits <- lapply(1:length(ss), function(s.ind){
            
            ss.spec <- ss[s.ind]
            spec.cols <- spec.cols + lags[s.ind]
            
            # ss.spec <- 1 31  39  97 132
            # print(ss.spec)
            # Back-fit the ref region to the filled spec data ####
              spec.region <- xmat[ss.spec, spec.cols]
              fit.feat2spec <- fit.batman(fit$feat.fit, spec.region, 
                                          exclude.lowest = .5, 
                                          ppm = ppm[spec.cols])
                # plot.fit(fit.feat2spec, type = "simple", ppm = ppm[spec.cols]) %>% plot
              
              fit.ref <- fit.feat2spec$ratio * fit$spec.fit + fit.feat2spec$intercept # here spec.fit is the ref.feature!
                # plot(fit.ref)
                # plot(spec.region)
                # rbind(fit.ref, spec.region) %>% simplePlot(linecolor = "black")
              residuals <- fit.ref - spec.region
                

            # Get pct. overshoot vector ####
              posres <- residuals > 0
                posres[is.na(posres)] <- 0
              # pct.overshoot <- sum(fit.ref[posres], na.rm=T) / sum(fit.ref, na.rm=T)
              # pct.overshoot <- residuals / fit.ref
              #   pct.overshoot[is.infinite(pct.overshoot)] <- 0
              #   pct.overshoot[!posres] <- 0
                
            # Assume each run of positive overshoot is a resonance. ####
            #    what % of its prominence is overshoot?             ####
            
              runs <- run.labels(posres)
              
              if (!any(runs > 0))
                {
                  ovs.res <- 0
                  ovs.tot <- 0
                }
              else{
                peaks <- extractPeaks_corr(fit.ref, plots = F)
                valleys <- peaks$bounds %>% unlist %>% unique
                
                # res.proms <- prominences(peaks, v2 = fit.ref)
                max.prom <- range(fit.ref, na.rm = T) %>% diff
                
                # For each run, expand outwards to nearest local minima in ref spec
                
                  worst.res <- lapply(1:max(runs), function(r) {
                      # r <- 1
                      # print(r)
                      inds <- which(runs == r)
                        if(length(inds) < 3 | sum(peaks$truePeak) < 1)
                        {
                          return(list(ratio.res = 0,
                                      ratio.total = 0))
                        }
                        
                      # Which ref peak are we looking at? What's its local prominence?
                      # Expand out from the boundaries of the peak in the positive residuals
                      # until you hit the next ref valley.
                        valleys.below <- valleys <= min(inds)
                        
                        if (any(valleys.below))
                          {lower.bound <- max(valleys[valleys.below])}
                        else{lower.bound <- min(inds)}

                        valleys.above <- valleys >= max(inds)
                        
                        if (any(valleys.above))
                          {upper.bound <- min(valleys[valleys.above])}
                        else{upper.bound <- max(inds)}
                        
                        ref.prom <- range(fit.ref[lower.bound:upper.bound], na.rm = T) %>% diff
                        
                      # What's the height of the positive residual peak here?
                        pos.res.prom <- range(residuals[inds], na.rm = T) %>% diff 
                        # (these will only be the positive residuals)
                      
                      # What's the ratio?
                        ratio.res <- pos.res.prom/ref.prom
                        
                      # or, weight using the ref peak's relative importance in the ref feature
                      # * note: this could increase the score if pos.res straddles >1 resonance
                        ratio.total <- pos.res.prom / max.prom
                        
                        return(list(ratio.res = ratio.res,
                                    ratio.total = ratio.total))
                        
                  })
                  
                  ovs.res <- lapply(1:length(worst.res), function(x) worst.res[[x]]$ratio.res) %>% unlist
                  ovs.tot <- lapply(1:length(worst.res), function(x) worst.res[[x]]$ratio.total) %>% unlist
              }
              
            # Plot result ####
                  # AUC ####
                  g <- NULL
                
                  if(plots){
                    
                        # For area plot to work correctly, no negative values allowed
                        pl.adj <- min(spec.region, na.rm = T)
                        ff <- fit.ref - pl.adj
                        sr <- spec.region - pl.adj
                        
                          # ff[fit.ref < 0] <- 0
                        g <- simplePlot(sr,
                                        xvect = ppm[spec.cols],
                                        linecolor = "gray",
                                        opacity = .9,
                                        linewidth = 1) # + geom_hline(yintercept = pl.adj)
                  
                        g <- g + new_scale_color() +
                              geom_area(data = data.frame(vals = ff,
                                                          ppms = ppm[spec.cols]),
                                        na.rm = TRUE,
                                        mapping = aes(x = ppms, y = vals),
                                        fill="blue", alpha=0.1,
                                        # fill="pink", alpha=0.6,
                                        color = "black",    # line color
                                        lwd = 0.25,    # line width
                                        linetype = 1)
    
                        
                  # Simple ####
                  # g <- simplePlot(spec.region,
                  #                 xvect = ppm[spec.cols],
                  #                 linecolor = "gray",
                  #                 opacity = .9, 
                  #                 linewidth = 1)
                  # 
                  # g <- g + new_scale_color() +
                  #         geom_line(data = data.frame(vals = fit.ref,
                  #                                     ppm = ppm[spec.cols],
                  #                                     corr = rep(-1, length(spec.cols))), 
                  #                   mapping = aes(x = ppm, y = vals, colour = corr),
                  #                   linewidth = .5) +
                  #         scale_colour_gradientn(colours = matlab.like2(10),
                  #                                limits = c(-1, 1))
                  }
            # Return results ####
              
              return(list(match = m,
                          ref = match.info[m, "ref"],
                          feat = match.info[m, "feat"],
                          ss.spec = ss.spec,
                          fit.ref = NA,
                          fit.ref.coef = data.frame(intercept = fit.feat2spec$intercept, 
                                                    ratio = fit.feat2spec$ratio),
                          spec.region = range(spec.cols),
                          ref.region = ref.region,
                          # residuals = residuals, # save on storage
                          # pct.overshoot = pct.overshoot, # save on storage
                          ovs.res = max(ovs.res),
                          ovs.tot = max(ovs.tot),
                          plot = g))
          })
          
        # The bff score (backfit feasibility score) is defined as:
        #   1 - overshoot viability score
        #   ovs:
        #   Each protrusion in the positive residuals (ref feat - spec) is scored:
        #     -find the resonance(s) that the protrusion overlaps
        #     -% intensity of the protrusion compared to the overall ref feature height?
        #     -take the max of this for the ref feature fit
        #     
        # Calculate bff from ovs score ####
          ovs.res <- lapply(fits, function(f) f$ovs.res) %>% unlist
          ovs.tot <- lapply(fits, function(f) f$ovs.tot) %>% unlist
          bffs.res <- 1-ovs.res
            bffs.res[bffs.res<0] <- 0
          bffs.tot <- 1-ovs.tot
            bffs.tot[bffs.tot<0] <- 0
          
        # Plotting (individual ref feat matches to spectra) ####
          gridplot <- NULL
          
          if (plots){
            # Generate a plot for with all the backfits, sorted by bff
              # Sort fits and bffs by bff score ####
                plot.order <- order(bffs.tot)
                fits <- fits[plot.order]
                bffs.tot <- bffs.tot[plot.order]
                bffs.res <- bffs.res[plot.order]

              # Add titles and bffs to the individual plots, do formatting ####
                plots <- lapply(1:length(fits), function(x){
                # plots <- lapply(which(abs(lags) > 3), function(x){
                  fits[[x]]$plot +
                    ggtitle(paste0("Spectrum ", fits[[x]]$ss.spec, 
                                   ": bff.res = ",round(bffs.res[x],2),
                                   ": bff.tot = ",round(bffs.tot[x],2))) +
                    theme(# axis.title.x=element_blank(),
                          # axis.text.x=element_blank(),
                          # axis.ticks.x=element_blank(),
                          plot.title = element_text(size=10))
                })
              
              # Make lattice plot ####
                gridplot <- gridExtra::grid.arrange(grobs = plots)
          }
            
      return(list(match = m,
                  fits = fits,
                  bffs.res = bffs.res,
                  bffs.tot = bffs.tot,
                  gridplot = gridplot))
      })#, mc.cores = pars$par$ncores)
      
      
  return(backfits)
}