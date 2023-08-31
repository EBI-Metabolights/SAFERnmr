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
  plots <- pblapply(1:nrow(sfs), function(x){
      # sf <- as.data.frame(sfs[1,])
      sf <- sfs[x, ]
      # print(x)
      # Get the feature model
        
        f.matched.pts <- sf$feat.start:sf$feat.end
        feat.gaps <- feature$position[sf$feat, f.matched.pts] %>% is.na
    
      # Get the ref segment (rf)
      
        rf <- sf$ref.start:sf$ref.end %>% refmat[sf$ref,.]
        
        # NA-fill the feature gaps
            
          rf[feat.gaps] <- NA
          
      # Get the spectrum data
          
          spec.reg <- c(sf$spec.start,sf$spec.end) %>% fillbetween()#%>% expand_window(within = c(1,ncol(xmat)))
            lag <- sf$spec.start - min(spec.reg) # store this for now
          
          spec.segment <- xmat[sf$ss.spec, spec.reg]

      # Put the rf into a vector the size of spec.segment    
      
          ref.segment <- rep(NA, length(spec.segment))
          ref.segment[1:length(rf) + lag] <- rf
                         
      # Fit to spec using recorded backfit params
        
        ref.segment <- as.numeric(sf$fit.intercept) + (ref.segment * as.numeric(sf$fit.scale))

      # Plot
    
        plot_fit(list(feat.fit = ref.segment, 
                      spec.fit = spec.segment))
    
    })
  
  # Add titles and bffs to the individual plots, do formatting ####
  
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
            ggtitle(paste0(", bff.tot = ",round(sfs$score[x],2)
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
