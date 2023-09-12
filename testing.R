# This is my working script for testing of SAFER
# 
# 

devtools::document('/Users/mjudge/Documents/GitHub/SAFER')

# First, we would like to see a random subset of ref-features -> spectra, at each bff score level, for a given study ####

  # Choose study

    study <- 'MTBLS1'
    tmpdir <- '/Users/mjudge/Downloads/nfs 3/production/odonovan/nmr_staging/pipeline_tests/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'
    # study <- 'MTBLS424'
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS'
    # study <- 'MTBLS430'
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/MTBLS430_1r_noesygppr1d.comp_spectralMatrix.RDS'
    # study <- 'MTBLS395'
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/MTBLS395_1r_cpmgpr1d.comp_spectralMatrix.RDS'

  # Import the study results ####
  
    results.dir <- paste0(tmpdir,'/')
    run_params <- yaml::yaml.load_file(paste0(tmpdir, '/params.yaml'), eval.expr = TRUE)
    
    # Read in match data
    
        backfit.results <- readRDS(paste0(results.dir,"smrf.RDS"))
          backfits <- backfit.results$backfits
          match.info <- backfit.results$match.info

    # Read in library data
    
        # lib.data.processed <- readRDS(paste0(results.dir, "lib.data.processed.RDS"))
        refmat <- readRDS(paste0(tmpdir, "/temp_data_matching/ref.mat.RDS")) %>% cstack_expandRows
          # refmat <- refmat %>% apply(2, function(x) x/sum(x, na.rm = T))
          # refmat <- refmat %>% t
          # refmat.c <- refmat.c %>% compress_stack
        
    # Read in spectral matrix data
    
        fse.result <- readRDS(paste0(results.dir, "fse.result.RDS")) 
          xmat <- fse.result$xmat
          ppm <- fse.result$ppm
          rm(fse.result)

    # Read in the features 
        feature <- readRDS(paste0(results.dir, "feature.final.RDS")) %>% expand_features
          # features.c <- feature %>% compress_features

    # Read in scores matrix 
      scores <- readRDS(paste0(results.dir,"scores.RDS"))
        scores.mat <- scores$ss.ref.mat 

      rfs.used <- scores$rfs.used

  # Pull a random set of backfits for a given range of scores  ####
  
    # Which scores? For each match, pick the best rf-ss fit for each score. This 
    # will allow us to examine each match in its best context while vastly reducing
    # numbers.
      
      # List out the backfit index with the match index and the score
      
        specfits <-
          mclapply(1:length(backfit.results$backfits), function(x)
            {
            
              # Combine the match info and backfit info into one (merged) table
              
                mi <- match.info[x, ]
                scores <- backfit.results$backfits[[x]]
                # scores <- scores[scores$bffs.tot > 0.75,]
                # scores <- scores[scores$rmse < 0.1,]
                # scores <- scores[scores$rmse.biased > 0.1,]
                best.scores <- c(which.max(scores$fit.fsa),
                                 which.max(scores$fit.rval)) %>% unique
                scores <- scores[best.scores, ]
                
                # how is it possible to get a match whose best score is 0?
                
                if (nrow(scores) == 0){ return(NULL)}
              
                # Add any vals from match info we may want for plotting
                
                  scores$feat <- mi$feat
                  scores$ref <- mi$ref
                  scores$match.rmse <- mi$rmse
                  scores$match.rval <- mi$rval
                  scores$feat.start <- mi$feat.start
                  scores$feat.end <- mi$feat.end
                  scores$ref.start <- mi$ref.start
                  scores$ref.end <- mi$ref.end
                  scores$match.id <- mi$id

              scores
              
            }, mc.cores = 4
          ) %>% rbindlist
  
  # Take a random subset of those (to keep compute down; will be sampling anyways) ####
    
      sf.copy <- specfits
      # specfits <- sf.copy
      specfits <- specfits %>% slice_sample(n = 5000, replace = FALSE)

       
  # Plotting ####
  
      # Compute combined score (exploration): ####
      
        specfits <- as.data.frame(specfits)
        sfs <- specfits
        f.models <- feature$position

        # Recalculate the fits
          
          sfs <- mclapply(1:nrow(specfits), function(x){
          
              # x <- 10
              # print(x)
                tryCatch(
                  {
                    sf <- sfs[x, ]
                    
                    res <- opt_specFit(sf, feature, xmat, refmat)
                    
                    res$sf
                    # res$fit %>% plot_fit(type = 'auc')
                  }, 
                  error = function(cond)
                  {
                    sf <- sfs[x, ]
                    sf$fit.fsa <- Inf
                    sf$fit.rval <- Inf
                    sf$fit.intercept <- Inf
                    sf$fit.scale <- Inf
                    sf$spec.start <- Inf
                    sf$spec.end <- Inf
                  }
                )
                
              # plot_fit(list(feat.fit = res$feat,
              #               spec.fit = res$spec), type = 'auc')
          
          }, mc.cores = 4
          ) %>% do.call(rbind,.)

        # scoreType = 'frac.accounted.spec'
        scoreType = 'rval_x_frac.accounted.spec'
        # sfs$score <- sfs$fraction.spec.accounted
        sfs$score <- sfs$fit.rval * sfs$fraction.spec.accounted
        sfs.copy <- sfs
        sfs <- sfs[!is.infinite(sfs$score), ]
        
        scattermoreplot(x = 1:nrow(sfs), 
                        y = (sfs$score) %>% sort)
          abline(h = 0, col = 'red')
        
        sfs <- sfs[sfs$score > 0, ]
        
      # Bin the data by score ####
      
        bins <- c(0, .2, .4, .6, .8, 1)
        
        specfits.sliced <- sfs %>% 
          mutate(bin = cut(score, breaks= bins)) %>%
          group_by(bin) %>%
          slice_sample(n = 25, replace = FALSE) %>% 
          arrange(score) %>% 
          as.data.frame
          
                  
      # Plot ####
        message('--- Plotting numbers only for ', study,' ---')
        grid_plot_specfits(specfits.sliced, feature, xmat, refmat, plotLoc = '/Users/mjudge/Desktop', 
                           filename = paste0(study, '_grid_specfits_', scoreType), titles = 'number')
        
        message('--- Plotting number and scores for ', study,' ---')
        grid_plot_specfits(specfits.sliced, feature, xmat, refmat, plotLoc = '/Users/mjudge/Desktop', 
                           filename = paste0(study, '_grid_specfits_', scoreType), titles = 'number_score')
        
        rm(backfit.results)
        rm(backfits)
        rm(rfs.used)
        
        save.image(file = paste0('/Users/mjudge/Desktop/',study, '_specfits_', scoreType, '.RData'))
      
        
        