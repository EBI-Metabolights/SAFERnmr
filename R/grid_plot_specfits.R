# grid_plot_specfits ############################################################################
#' 
#' Plot a grid of ref features fit to spectra. 
#'
#'
#' @param sfs long format spec feature information (includes backfit information)
#' @param feature features object for study
#' @param xmat spectral matrix for study
#' @param refmat reference spectral matrix for study (spectra on rows)
#' @param plotLoc plot location, ends with "/"
#' @param filename filename, excluding ".pdf"
#' @param titles options for data to display on title of each subplot: 
#' - 'number' (just the plot number, for tracking)
#' - 'score' (just the score)
#' - 'number_score' (both)
#' 
#' @return grid plot written to pdf specified
#' 
#' @importFrom magrittr %>%
#' @importFrom pbapply pblapply
#' @importFrom ggplot2 ggplot ggtitle theme element_blank element_text
#' @importFrom gridExtra grid.arrange
#' 
#' @export
grid_plot_specfits <- function(sfs, feature, xmat, refmat,
                               plotLoc = './',
                               filename = 'grid_specfits',
                               titles = 'number' # 'score', 'number_score'
                               ){

  message('Generating plots...')
  sfs <- arrange(sfs, score)
  
  plots <- pblapply(1:nrow(sfs), 
                    function(x){
        # x <- 27
        # res <- opt_specFit(sfs[x, ], feature, xmat, refmat)
        # plot_fit(list(feat.fit = res$feat,
        #               spec.fit = res$spec), type = 'auc')
        
        # fit <- res$sf %>% rebuild_specFit(feature, xmat, refmat)
        
        fit <- rebuild_specFit(sfs[x, ], feature, xmat, refmat)
          
        # # Zoom to feature
        #   
        #   bottom <- min((c(fit$spec, fit$feat)), na.rm = TRUE)
        #   ys <- range(fit$feat - bottom, na.rm = TRUE)
        #   adj <- c(-0.1, 0.1)
        #   new.lims <- ys + ys*adj
        
        plot_fit(list(feat.fit = fit$feat,
                      spec.fit = fit$spec), 
                 type = 'auc')
        
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
    dim <- 3*round(sqrt(nrow(sfs)))
    pdf(file = paste0(plotLoc,'/',filename,"_",titles,".pdf"),   # The directory you want to save the file in
        width = dim, # The width of the plot in inches
        height = dim)
    
      gridExtra::grid.arrange(grobs = plots)
    
    dev.off()

    
}
