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
          stackplot(ppm[reg], vshift = .3, hshift = 0)+
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_stackplot.pdf', width = 12, height = 6)
            g
        dev.off()
      
            
      # Sorted stackplot with highlighted subset
      
        sortorder.ss <- specRegion[s$subset, ] %>% rowSums(na.rm = TRUE) %>% order
        non.ss <- which(!(1:nrow(specRegion) %in% s$subset))
        reorder <- c(s$subset[sortorder.ss], non.ss) # size-sorted subset, then rest of matrix (unsorted)
        
        devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
        g <- specRegion[reorder,] %>%
          stackplot(ppm[reg], vshift = .3, hshift = 0, highlight.spec = 1:length(s$subset))+
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")
        
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
        
        s=fse.result$storm_features[[a.ind]]
        g <- plot_storm_refRegions(xmat = xmat, s = s, ppm = ppm, calcStocsy = F, xlim = range(reg) %>% ppm[.]) +
          geom_vline(xintercept = ppm[s$peak], linetype = 2, col = "black")
        
        pdf(file = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/citrate_peak_storm.pdf', width = 12, height = 6)
          g
        dev.off()

  # Citrate ref
  
    ldp <- readRDS(paste0(tmpdir,"/lib.data.processed.RDS"))
    ref.names <- c("Citrate", "AMP", "Norepinephrine", "Niacin", "L-Phenylalanine", "L-glutamine")
    
    lapply(ref.names, function(r){
      
      ref_pdf(ref.name = r, 
              ldp = ldp, 
              ppm = ppm, 
              where = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling')
      
    })

  ref_pdf <- function(ref.name, ldp, ppm, where = '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling'){
    message(ref.name)
    ref.num <- ldp %>% lapply(function(x) x$compound.name) %>% unlist %>% stringr::str_equal(ref.name) %>% which %>% .[1]
    ref <- ldp[[ref.num]] %>% expand_ref(ppm)
    ref$mapped$data[is.na(ref$mapped$data)] <- 0
    pdf(file = paste0(where,'/',ref.name,'_ref.pdf'), width = 12, height = 6)
     
      simplePlot(ref$mapped$data, 
                 xvect = ref$mapped$ppm, 
                 xdir = 'reverse', 
                 linecolor = 'blue', 
                 opacity = 0.5) %>% print
    dev.off()
    return(NULL)
  }
  