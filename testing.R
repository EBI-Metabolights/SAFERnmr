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
      
        
########### parameter sensitivity testing #########

# Write out a table of all pars tested ####
  # Locate files  ####
    parfile.loc <- '/Users/mjudge/Desktop/param_templates_sensitivity_testing/'
    
    # parsets <- parfile.loc %>% dir 
     parsets <- c("sloppy", "less.sloppy", "standard", "less.tight" , "tight")
  
  # Get the pars from each set-study, and format as df  ####
    pars <- lapply(parsets, function(parset){
      
      this.loc <- paste0(parfile.loc,parset)
      files <- this.loc %>% dir(pattern = '\\.yaml$', full.names = TRUE)
      
      lapply(files, function(f){
        
        study.name <- stringr::str_extract(f, pattern = "(?<=/)MTBLS\\d+")
        pars <- yaml::yaml.load_file(f, eval.expr = TRUE)
        df <- pars %>% rrapply::rrapply(how = 'melt')
        df$pars <- parset
        df$study.name <- study.name
        df
        
      }) %>% do.call(rbind,.)
    
    }) %>% do.call(rbind,.)
    
  # Clean up the df ####
    pars$run <- paste0(pars$study.name, '-', pars$pars)
      par.names <- pars[, 1:3]
      par.names[is.na(par.names)] <- ''
      par.names <- tidyr::unite(par.names, name, sep = '_')
      pars$name <- par.names$name
    
    pars <- subset(pars, 
                   select = -c(L1, L2, L3, run))
    
    pars <- tidyr::pivot_wider(pars, id_cols = c("study.name","pars"), names_from = name, values_from = value)

  # Write to csv ####
    pars <- apply(pars,2,as.character)
    write.csv(pars, file = paste0(parfile.loc,'par.table.csv'), na = '', row.names = FALSE)
    
# Extract results ####        

  # Loop through run results to get summary stats ####
    test.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/'
    parsets <- c("sloppy", "std", "tight")
    
    results <- 
      lapply(parsets, function(parset){
        
        files <- dir(paste0(test.dir, parset), full.names = TRUE)
        runs <- files %>% stringr::str_detect(".zip") %>% "!"(.) %>% files[.]
        
        lapply(runs, function(run){
          
          study.name <- stringr::str_extract(run, pattern = "(?<=/)MTBLS\\d+")
          
          message(parset, ' - ', study.name)
          
          feature <- readRDS(paste0(run, '/feature.final.RDS'))
            n.features <- feature$sfe %>% length
            rm(feature)
            
          match.info <- tryCatch({
              matches <- readRDS(paste0(run, '/matches.RDS'))
  
              # if any invalid matches for the feature
                matches <- matches[!is_nullish(matches)]
              # per feature, was NA returned, or was ('matches', 'peak.quality')?
                nomatch <- (lapply(matches, length) %>% unlist) == 1 
                matches <- matches[!nomatch]
              
              # For each feature: 
                # Check if there was an error message
                  errors <- matches[names(matches) %in% c('call', 'message')]
                # Check if both of the expected fields are present
                  matches <- matches[names(matches) %in% c('matches', 'peak.quality')]
        
                matches.split <- split(matches, names(matches))
                rm(matches)
                
              match.info <- rbindlist(matches.split$matches)
              
                list(
                  n.matches = nrow(match.info),
                  backfits = NA
                )
              
              
          }, warning = function(cond){
            
            match.info <- readRDS(paste0(run, '/smrf.RDS'))
            return(
                list(
                  n.matches = nrow(match.info$match.info),
                  backfits = lapply(match.info$backfits, nrow) %>% unlist %>% sum
                )
            )
              
          }, error = function(cond){
            
            match.info <- readRDS(paste0(run, '/smrf.RDS'))
            return(
                list(
                  n.matches = nrow(match.info$match.info),
                  backfits = lapply(match.info$backfits, nrow) %>% unlist %>% sum
                )
            )
            
          })
          
          scores <- tryCatch({
              scores <- readRDS(paste0(run, '/scores.RDS'))
              mat <- scores$ss.ref.mat
                n.compounds <- mat %>% Rfast::rowMaxs(., value = TRUE) %>% ">"(.,.5) %>% sum
                n.compound.spec <- mat %>% ">"(.,.5) %>% sum
                max.score <- max(mat)
                  
                list(
                  compounds = n.compounds,
                  comp.specs = n.compound.spec,
                  max.score = max.score
                )
              
                
            }, warning = function(cond){
              
                list(
                  compounds = NA,
                  comp.specs = NA,
                  max.score = NA
                )
  
            }, error = function(cond){
              
                list(
                  compounds = NA,
                  comp.specs = NA,
                  max.score = NA
                )
  
            })
          
            df <- data.frame(study = study.name,
                             pars = parset,
                             features = n.features,
                             matches = match.info$n.matches,
                             backfits = match.info$backfits,
                             compounds = scores$compounds,
                             max.score = scores$max.score,
                             ms.pairs = scores$comp.specs)
            print(df)
            return(df)
            
        }) %>% do.call(rbind,.)
        
      }) %>% do.call(rbind,.)
     
# Get the number of matches (~ rval)  ####

    # same ones, just filtered?
    
# smrfs (backfits)
       
    # same ones, just filtered? 
    
# compound-spectra by score
        
    # same ones, just filtered?
        
# For each std run, randomly select n scores  ####
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/std/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'

study.name <- stringr::str_extract(tmpdir, pattern = "(?<=/)MTBLS\\d+")
out <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/compound_level_test/'
dir.create(out, showWarnings = FALSE)

scores <- scores <- readRDS(paste0(tmpdir, '/scores.RDS'))
scores.mat <- scores$ss.ref.mat
rownums <- 1:nrow(scores.mat)

score.summary <- lapply(rownums, function(m){
  
  # Get the top 5 scores for this metabolite
  
    met.scores <- scores.mat[m,]
      score.order <- order(met.scores, decreasing = TRUE)
      
      # df <- data.frame(x = 1, score = met.scores)
      # ggplot(df, aes(x = x, y = score)) +
      #   geom_violin(alpha = 0.5) +
      #   geom_point(position = position_jitter(seed = 1, width = 0.2)) +
      #   theme(legend.position = "none")  
      
    top.5.samples <- score.order[1:5]
    top.5.scores <- met.scores[top.5.samples]
    
  # Report the mean score
    mean.score <- mean(top.5.scores)
    
    list(mean.score = mean.score,
         top.5.samples = top.5.samples,
         top.5.scores = top.5.scores)
}) %>% unlist(recursive = FALSE)

score.summary <- split(score.summary, names(score.summary))

mean.scores <- score.summary$mean.score %>% unlist

df <- data.frame(mean.score = mean.scores,
                 ref.num = 1:nrow(scores.mat),
                 compound.name = rownames(scores.mat))

# Bin the data by score ####

  bins <- c(0, .25, .5, .75)
  max.hits <- 5
  
  scores.sliced <- df %>% 
    arrange(desc(mean.score)) %>%
    mutate(bin = cut(mean.score, breaks= bins)) %>%
    group_by(bin) %>% 
    slice(                                       
            1:min( 
                    c(  length(mean.score),  
                        max.hits  
                      )
                  )                      
    ) %>% 
    as.data.frame %>% 
    apply(2, as.character)

    # scores.sliced <- df %>% 
    # slice_sample(n = 5, replace = FALSE) %>% 
    # arrange(mean.score) %>% 
    # as.data.frame %>% apply(2, as.character)

  # Write to csv ####
    write.csv(scores.sliced, file = paste0(out,'_',study.name, '_scores.sliced.csv'), na = '', row.names = FALSE)

# Just use top 10 from the study
  
  df <- df %>% arrange(desc(mean.score))
  head(df, 20)
  
# Viewing 
  devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
  browse_evidence(results.dir = tmpdir)

  
         