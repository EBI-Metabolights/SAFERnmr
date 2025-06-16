grid_plot_sats <- function(sats, xmat, ppm,
                               plotLoc = '.',
                               filename = 'grid_sats',
                               titles = 'number' # 'score', 'number_score'
                               ){

  message('Generating plots...')
  sats <- sat.list
  
  plots <- pblapply(sats, 
                    function(s){
        # s <- sats[[1]]
      plot_protofeature(p = data.frame(driver = s$peak),
                  half.window = half.window, ppm = data$ppm,
                  xmat = xmat[s$subset,],
                  bgplot = 'overlay', line.shape = 'covar', line.color = "corr",
                  showPeaks = FALSE, ref.mask = s$ref.idx, show.mask.bounds = TRUE)
    })
  
  # Add titles and scores to the individual plots, do formatting ####
  
      if (titles == 'number'){
        plots <- lapply(1:length(plots), function(x){
          
          plots[[x]] +
            ggtitle(x) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  plot.title = element_text(size=10))
          
        })
      }
  
      if (titles == 'score'){
        plots <- lapply(1:length(plots), function(x){
          
          plots[[x]] +
            ggtitle(paste0(", score = ",round(sfs$score[x],2)
                           )) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  plot.title = element_text(size=10))
          
        })
      }

      if (titles == 'number_score'){
        plots <- lapply(1:length(plots), function(x){
          
          plots[[x]] +
            ggtitle(paste0(x, ", score = ",round(sfs$score[x],2)
                           )) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  plot.title = element_text(size=10))
          
        })
      }

  # Print the plots to pdf
  
    message('Printing plots to file: ', paste0(plotLoc,'/',filename,"_",titles,".pdf"),'...')
    dim <- 3*round(sqrt(length(sats)))
    pdf(file = paste0(plotLoc,'/',filename,"_",titles,".pdf"),   # The directory you want to save the file in
        width = dim, # The width of the plot in inches
        height = dim)
    
      gridExtra::grid.arrange(grobs = plots)
    
    dev.off()

    
}
