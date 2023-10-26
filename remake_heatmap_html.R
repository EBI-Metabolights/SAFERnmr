devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698286479' # .6

browse_evidence(tmpdir)



  pars <- yaml::yaml.load_file(paste0(tmpdir,'/params.yaml'), eval.expr = TRUE)
  pars$dirs$temp <- tmpdir
  # Turn this into a nonduplicate ss-ref matrix  ####
  
      scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
      ss.ref.mat <- scores$ss.ref.mat %>% t

      # Plot the matrix as an HCA'd heatmap ####
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
              
                drawHeatmap(ss.ref.mat[ref.order,,drop=FALSE], dropRowNames = FALSE) %>% 
                  htmlwidgets::saveWidget(file = paste0(tmpdir,"/match_scores_sample_x_compound.html"), selfcontained = TRUE)
            
          },
          error = function(cond){
            message('\n\tScore matrix plotting error..')
          }
        )

      # Clean up any plotly js dependency files (not needed) ####
          
          js.dependency.dir <- paste0(pars$dirs$temp, '/match_scores_sample_x_compound_files')
          if (dir.exists(js.dependency.dir)){
            unlink(js.dependency.dir, recursive = TRUE)
          }
          
