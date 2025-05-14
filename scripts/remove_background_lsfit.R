

            i <- i + 10
            i
            
            p <- pfs[[i]]
            plot_protofeature(p,
                      half.window = half.window, ppm = data$ppm,
                      xmat = data$xmat,
                      bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
                      showPeaks = TRUE, show.mask.bounds = FALSE)
                        
            pexp <- expand_protofeature(p, xmat, ppm, half.window)
            driver <- pexp$driver
            
            wind <- pexp$specRegion.inds
            x <- wind
            specRegion = pexp$specRegion
            
            simplePlot(xmat[,wind], ppm[wind])
            simplePlot(specRegion, ppm[wind])
            
            mean.spec <- colMeans(specRegion)
            
            fits <- lapply(1:nrow(specRegion), function(m){
              # m <- m + 1
              fit <- fit_leastSquares(specRegion[m, ] %>% c, mean.spec, plots = FALSE)
              fit$fit
              # fit$plot
            })
 