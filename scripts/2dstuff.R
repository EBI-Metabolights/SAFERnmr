 devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
 # browse_evidence('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1700653867')
 # plot_spec(hs[500, ], 1:length(hs[50,]))
 
 hs <- readRDS('/Users/mjudge/Downloads/hsqc.RDS')
 data <- hs[-1, ]
 ppm <- hs[1, ] 
 
 # scattermore::scattermoreplot(x = 1:length(data), y = sort(data))

 heatmap_raster <- function(data, 
                            ppm1 = NULL, ppm2 = NULL, 
                            signal.cutoff = 1E9, contours = 10, 
                            transparency = 0.1, point.size = .5,
                            reverse = TRUE)
 {
   
   if (is.null(ppm1)){
     ppm1 <- 1:ncol(data)
   }
   if (is.null(ppm2)){
     ppm2 <- 1:nrow(data)
   }
   if (reverse){
     xrange <- range(ppm1) %>% rev
     # yrange <- range(ppm2) %>% rev
   }
   
   subs <- ind2subR(1:length(data), nrow(data))
   pos <- abs(data) >= signal.cutoff
   subs <- subs[pos, ]
   subs$val <- data[pos]
   
   intrange <- range(subs$val, na.rm = TRUE)
   col.levels <- seq(intrange[1], intrange[2], length.out = contours)
   col.levels[1] <- 0
   ntimes = 1:(length(col.levels)-1)
   # ntimes = round(ntimes ^ 1.5 )
   
   subs$level <- cut(subs$val, breaks = col.levels, labels = ntimes, include.lowest = TRUE)
   
   subs <- subs %>% 
     dplyr::group_by(level) %>% 
     mutate(group_number = group_indices()) %>% 
     dplyr::slice(rep(1:n(), each = level[1]))
   
   scattermore::scattermoreplot(x = ppm1[subs$cols],
                                y = ppm2[subs$rows],
                                xlab = 'ppm',
                                ylab = '',
                                size = dim(data),
                                cex = point.size,
                                # ylim = yrange,
                                xlim = xrange,
                                col = alpha('red', alpha = transparency),
                                yaxt="n",
                                xaxs = "i", 
                                yaxs = "i")
     
 }

 heatmap_raster(data, 
                ppm1 = ppm, 
                signal.cutoff = 1E9, contours = 10, 
                transparency = 0.1, point.size = .5)

