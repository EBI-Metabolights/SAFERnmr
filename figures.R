devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/std/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'

feature <- readRDS(paste0(tmpdir, '/feature.final.RDS')) %>% expand_features

fse.result <- readRDS(paste0(tmpdir, '/fse.result.RDS'))
xmat <- fse.result$xmat
ppm <- fse.result$ppm

drivers <- lapply(1:length(feature$sfe), function(x){
  
  feature$position[x, 
                   feature$sfe[[x]]$feat$driver.relative]
}) %>% unlist %>% ppm[.]

ind <- vectInds(2.68, drivers)
drivers[ind]

a.drivers <- lapply(fse.result$storm_features, function(s){
  # s <- fse.result$storm_features[[1]]
  driver <- s$peak
  if(is.null(driver)){
    0
  } else {
    driver
  }
}) %>% unlist


vectInds(vectInds(2.68, ppm), a.drivers) %>% a.drivers[.] %>% ppm[.]
a.ind <- vectInds(vectInds(2.68, ppm), a.drivers)

# feature$stack[ind, ] %>% trim_sides %>% simplePlot

cols <- feature$position[ind,] %>% trim_sides(out = "inds")
feat.span <- span(cols)
reg <- feature$position[ind, cols] %>% expand_window(within = c(1,length(ppm)), by = feat.span)

xmat[, reg] %>% stackplot(xvect = ppm[reg], vshift = 3)
s=fse.result$storm_features[[a.ind]]
driver <- s$peak
which(reg %in% s$peak)

      specRegion <- xmat[, reg] %>% scale_between
      scale.to <- specRegion %>% Rfast::colMaxs(value = TRUE)
      
      
      sr <- list(
        cor = cor(specRegion, xmat[, driver]),
        cov = cor(specRegion, xmat[, driver]) %>% scale_between# %>% fit_batman(scale.to, exclude.lowest = 0.5) %>% .[["feat.fit"]]
      )

      g <- simplePlot(specRegion, ppm[reg], n_xticks = 4)
      # g <- stackplot(specRegion, ppm[reg], vshift = 1, hshift = 0)
      
      # # Make the plot
      #   g <- g + geom_path(data = data.frame(vals = sr$cov,
      #                                        ppms = ppm[reg]),
      #                      mapping = aes(x = ppms, y = vals),
      #                      colour = "black", linewidth = .5,
      #                      lineend = "round",
      #                      linejoin = "round",
      #                      linemitre = "5")

      # Plot the stocsy line over the current plot
        
        
        g <- g + new_scale_color() +
          geom_line(data = data.frame(refvals = sr$cov,
                                      refppms = ppm[reg],
                                      corr = sr$cor), 
                           mapping = aes(x = refppms, y = refvals, colour = corr),
                    linewidth = 1.25,
                    lineend = "round",
                    linejoin = "round",
                    linemitre = "5") +
          scale_colour_gradientn(colours = matlab.like2(10),
                                 limits = c(-1, 1))
          
      # Plot corr boundaries as vert lines  
        g <- g +
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")

        g

  plot_storm_refRegions(xmat = xmat, s = s, ppm = ppm, calcStocsy = F) +
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")
 