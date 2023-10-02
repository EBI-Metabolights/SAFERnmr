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
        cov = cor(specRegion, xmat[, driver]) %>% scale_between # %>% fit_batman(scale.to, exclude.lowest = 0.5) %>% .[["feat.fit"]]
      )

      specRegion[s$subset, ]  # pull subset
      specRegion[sortorder, ] # reorder rows low to high

      # Stackplot
      
        g <- specRegion %>%
          stackplot(ppm[reg], vshift = .3, hshift = 0)
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_stackplot.pdf', width = 12, height = 6)
            g
        dev.off()
      
            
      # Sorted stackplot with highlighted subset
      
        sortorder <- specRegion %>% rowSums(na.rm = TRUE) %>% order
        ss.sorted <- order(sortorder) %>% .[s$subset]
          
        g <- specRegion[sortorder,] %>%
          stackplot(ppm[reg], vshift = .3, hshift = 0, highlight.spec = ss.sorted)
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_stackplot_selected.pdf', width = 12, height = 6)
            g
        dev.off()
      
      
      # Plot the stocsy line over the simple overlay plot
        
        g <- simplePlot(specRegion, ppm[reg], n_xticks = 4)
        
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
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_stocsy.pdf', width = 12, height = 6)
          g
        dev.off()
        
        
        s$finalRegion <- reg
        sr <- list(
          cor = cor(specRegion[s$subset, ], xmat[s$subset, driver]),
          cov = cor(specRegion[s$subset, ], xmat[s$subset, driver]) %>% scale_between # %>% fit_batman(scale.to, exclude.lowest = 0.5) %>% .[["feat.fit"]]
        )
        sr$cor[sr$cor < 0.8] <- NA
        sr$cov[sr$cov < 0.8] <- NA
        s$corr <- sr$cor
        s$covar <- sr$cov
        devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
        g <- plot_storm_refRegions(xmat = xmat, s = s, ppm = ppm, calcStocsy = F, xlim = range(reg) %>% ppm[.]) +
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_storm.pdf', width = 12, height = 6)
          g
        dev.off()




 