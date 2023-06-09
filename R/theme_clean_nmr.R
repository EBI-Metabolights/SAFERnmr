#' ggplot theme for viewing nmr spectra
#' 
#' all white background, few things 
#'
#' @return changes theme on current plot 
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom scales breaks_pretty
#' 
#' @import shiny
#' @import plotly
#' 
#' @export
theme_clean_nmr <- function()
{

        # Reverse axis ####
          # ggplot2::scale_x_reverse(breaks = scales::breaks_pretty())

      # Make simple ####
          ggplot2::theme_bw() %+replace%
          ggplot2::theme(axis.text = element_text(colour = "black",size = 12),
                        legend.position = "none",
                        axis.text.y = ggplot2::element_blank(),
                        # axis.title.x = element_text(size = 16,vjust = -0.5),
                        axis.title.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank(),
                        # axis.title = ggplot2::element_text(size = 12,vjust = 0.5),
                        panel.border =  ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        # panel.grid.major = ggplot2::element_line(color = "gray",
                        #                                         size = 0.1,
                        #                                         linetype = 1),
                        panel.grid.major = ggplot2::element_blank())
  

}

# theme_clean_nmr <- function(){ 
#     font <- "Georgia"   #assign font family up front
#     
#     theme_minimal() %+replace%    #replace elements we want to change
#     
#     theme(
#       
#       #grid elements
#       panel.grid.major = element_blank(),    #strip major gridlines
#       panel.grid.minor = element_blank(),    #strip minor gridlines
#       axis.ticks = element_blank(),          #strip axis ticks
#       
#       #since theme_minimal() already strips axis lines, 
#       #we don't need to do that again
#       
#       #text elements
#       plot.title = element_text(             #title
#                    family = font,            #set font family
#                    size = 20,                #set font size
#                    face = 'bold',            #bold typeface
#                    hjust = 0,                #left align
#                    vjust = 2),               #raise slightly
#       
#       plot.subtitle = element_text(          #subtitle
#                    family = font,            #font family
#                    size = 14),               #font size
#       
#       plot.caption = element_text(           #caption
#                    family = font,            #font family
#                    size = 9,                 #font size
#                    hjust = 1),               #right align
#       
#       axis.title = element_text(             #axis titles
#                    family = font,            #font family
#                    size = 10),               #font size
#       
#       axis.text = element_text(              #axis text
#                    family = font,            #axis famuly
#                    size = 9),                #font size
#       
#       axis.text.x = element_text(            #margin for axis text
#                     margin=margin(5, b = 10))
#       
#       #since the legend often requires manual tweaking 
#       #based on plot content, don't define it here
#     )
# }