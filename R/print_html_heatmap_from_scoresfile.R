#' print_html_heatmap_from_scoresfile
#'
#' Accessory function to generate and print a heatmap in the browse_evidence style
#' (using drawHeatmap()) and print it to a self-contained html file.
#' Note: requires pandoc to be installed if running on Linux systems. If it is not
#' installed, the file will be written in a different format (YAML metadata e.g. 
#' for Jekyll??) that browsers won't interpret. If this is the case, you can run
#' check_heatmap_html_output() on the zip file for the results
#' 
#'
#' @param pars a list containing necessary parameters for the function
#' @return RDS file containing match scores between library reference spectra and the best subset spectrum score for each.
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom htmlwidgets saveWidget
#'
#'
print_html_heatmap_from_scoresfile <- function(tmpdir){
  
        
        scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
        ss.ref.mat <- scores$ss.ref.mat %>% t
        
      # Plot the matrix as an HCA'd heatmap
        ann.cmpds <- tryCatch({
          Rfast::colsums(ss.ref.mat) > 0
        }, error = function(cond){
          TRUE
        })
        
        
        tryCatch(
          {
              ss.ref.mat <- ss.ref.mat[, ann.cmpds, drop = FALSE] %>% t
              
            # Always cluster samples
              sample.order <- hclust(dist(t(ss.ref.mat)))$order
              ss.ref.mat <- ss.ref.mat[, sample.order,drop = FALSE]
              
            # Only cluster refs if there are > 1
              if (nrow(ss.ref.mat) == 1){
              
                ref.order <- 1
                
              } else {
              
                ref.order <- hclust(dist(ss.ref.mat))$order
                
              }
              
                heatmap.file <- paste0(tmpdir,'/match_scores_sample_x_compound.html') %>% stringr::str_replace_all('//', '/')
                drawHeatmap(ss.ref.mat[ref.order,,drop=FALSE], dropRowNames = FALSE) %>% 
                  htmlwidgets::saveWidget(file = heatmap.file, selfcontained = TRUE)
            

          },
          error = function(cond){
            message('\n\tScore matrix plotting error..')
          }
        )

        # Clean up any plotly js dependency files (not needed)
          
          js.dependency.dir <- paste0(tmpdir, '/match_scores_sample_x_compound_files')
          if (dir.exists(js.dependency.dir)){
            unlink(js.dependency.dir, recursive = TRUE)
          }
          
      return(sum(ann.cmpds))
}