#' Plot features onto ref region using selected data from select_evidence_refmet()
#' Accessory plotting function for browse_evidence(), makes EST-style lines that 
#' allow user to browse the feature-level evidence for a given reference (instead
#' of sample-specific). Answers the question: what is the coverage for different
#' parts of the reference spectrum?
#'
#' @param reg currently plotted ppm region of reference
#' @param metab.evidence selected metabolite evidence for given samples
#' @param features.c compressed feature obj.
#' @param ppm ppm vector
#' @param plt.pars plot parameters (main one is plt.pars$pixels)
#' 
#' @return scattermore plot with two layers (lines and blocks)
#' @importFrom magrittr %>%
#' @importFrom scattermore scattermoreplot
#' 
#' @export 
feature_est_plot <- function(reg, metab.evidence, features.c, ppm, plt.pars){

 # Reset data ####

        reg <- sort(reg)
        rfs <- metab.evidence$match.info.ss
        rngs <- rbind(ppm[rfs$ref.start], ppm[rfs$ref.end])
        ft.inds <- rfs$feat %>% unique
        
    # Read in feature data and get feature positions for those matches ####
        
        feature <- features.c %>% expand_features(rfs$feat)
        f.seq <- as.numeric(factor(rfs$feat, levels = rfs$feat %>% unique))
        f.pos <- feature$position[f.seq, , drop = F] # duplicate the rows as necessary (each feature can have multiple fits)

      # Move the feature so the feat.start == ref.start
        f.starts <- sub2indR(rows = 1:nrow(f.pos), 
                             cols = rfs$feat.start, 
                             m = nrow(f.pos))
        f.pos <- f.pos - (f.pos[f.starts] - rfs$ref.start)
        
          pos <- f.pos
          pos.lines <- apply(pos, 1, function(r){
            # r <- pos[1,]
            inds <- trim_sides(r, out = "inds")
            r[inds] <- inds %>% range %>% r[.] %>% fillbetween
            return(r)
          }) %>% t
          
    # Convert to points for scattermore ####
        ft <- ind2subR(1:length(pos),nrow(pos))
        
        # boxes
        ft.pts <- data.frame(x = ppm[c(pos)],
                             y = ft$rows) %>% filter(!is.na(x)) #%>% filter(x>=reg[1] & x <= reg[2])
        # lines
        ft.lines <- data.frame(x = ppm[c(pos.lines)],
                               y = ft$rows) %>% filter(!is.na(x)) # %>% filter(x>=reg[1] & x <= reg[2])
        
        ################
    # Plot the points (raster) ####  
          
          xrange <- reg %>% rev
          # pdf(file = 'test.plot.pdf')
          scattermore::scattermoreplot(
                                        x = ft.lines$x,
                                        y = ft.lines$y,
                                        xlab = 'ppm',
                                        ylab = '',
                                        size = plt.pars$pixels,
                                        cex = .1,
                                        xlim = xrange,
                                        col = 'black',
                                        yaxt="n",
                                        xaxs = "i", 
                                        yaxs = "i"
                                      ) 
          par(new=TRUE)
          scattermore::scattermoreplot(
                                        x = ft.pts$x,
                                        y = ft.pts$y,
                                        xlab = 'ppm',
                                        ylab = '',
                                        size = plt.pars$pixels,
                                        cex = 2,
                                        xlim = xrange,
                                        col = 'blue',
                                        yaxt="n",
                                        xaxs = "i", 
                                        yaxs = "i"
                                      )           
        # dev.off()
        
}

