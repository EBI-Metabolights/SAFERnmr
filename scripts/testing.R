# This is my working script for testing of SAFER
# 
# 

devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655
# tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1697747670'
# tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698347163'
browse_evidence(tmpdir)

# backfit.results <- readRDS(paste0(tmpdir, "/smrf.RDS"))
# scores <- readRDS(paste0(tmpdir,"/scores.RDS"))

pars <- yaml::yaml.load_file(paste0(tmpdir,'/params.yaml'), eval.expr = TRUE)
pars$dirs$temp <- tmpdir
pars$debug$enabled <- TRUE
pars$par$ncores <- 4
pars$debug$all.outputs <- TRUE
pars$matching$filtering$max.backfits <- 1E5

# Accessory ####

data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'
index_studies <- function(data.dir, exclude = NULL){
  
  # Extract out run summary data ####
  
      data <- dir(data.dir)
        dontuse <- grepl(paste('.zip$','.csv$','.pdf$','.xlsx$',exclude %>% paste(collapse = "|"), sep = "|"), data)
        unzipped <- !dontuse
      data.unzipped <- paste0(data.dir,data[unzipped])


    df <- lapply(data.unzipped, function(x){
        dat <- read.csv(paste0(x, '/run.summary.csv'))
        names(dat) <- names(dat) %>% stringr::str_replace_all('\\.', '_')
        dat$local_id <- dat$run_id
        dat
    }) %>% bind_rows
    
    # df$total_time[df$total_time > 10]  <- df$total_time[df$total_time > 10] / 60
    # df$matching_elapsed_time[df$matching_elapsed_time > 10]  <- df$matching_elapsed_time[df$matching_elapsed_time > 10] / 60
    df$local_path <- data.unzipped
    df <- df[,!(names(df) == 'X')]
    df$start <- df$run_id %>% as.POSIXct(origin = "1970-01-01")
    df
}
unzip_studies <- function(data.dir, exclude = NULL){
  

    # Unzip any non-unzipped files ####
      # data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'
      data <- dir(data.dir)
        excluded <- which(grepl(exclude %>% paste(collapse = "|"), data))
        zipped <- grepl('.zip', data)
        not.zip <- which(!zipped)
        zipped <- which(zipped)
          zip.names <- data[zipped] %>% stringr::str_remove_all('.zip')
        already.unzipped <- (zip.names %in% data[not.zip]) %>% zipped[.]
        dont.unzip <- c(already.unzipped, excluded, not.zip)
        do.unzip <- data[-dont.unzip]
          
      if (length(do.unzip) > 0){
        fail <- lapply(do.unzip, function(d){
          # d <- data[-dont.unzip] %>% .[1]
          message('Unzipping ', d, ' ...')
          
          # New ####
            tryCatch({
              run.id <- d %>% stringr::str_remove('.zip$')
              out <- paste0(data.dir, '/',run.id) %>% stringr::str_replace_all('//','/')
              unzip(paste0(data.dir, '/',d) %>% stringr::str_replace_all('//','/'), junkpaths = FALSE, exdir = out)
              check_heatmap_html_output(out)
              return(0)
              
            }, error = function(cond){
              1
            })
        })
        
      } else {
        
        fail <- list('There are no zipped studies which need to be unzipped.')
      }

  if (is.numeric(fail)){
    return(fail %>% unlist)
  } else {
    return("nothing to unzip")
  }
  
}

    unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
    run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))

########### Just look at one ######
  
  browse_evidence(res.dir, select.sample = select.sample)
    
########### Random subset of smrf fits ####
# First, we would like to see a random subset of ref-features -> spectra, at each bff score level, for a given study ####

  # Choose study

    study <- 'MTBLS1'
    tmpdir <- '/Users/mjudge/Downloads/nfs 3/production/odonovan/nmr_staging/pipeline_tests/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS'ø
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
      
        
########### parameter sensitivity testing (pre-summ file) #########

tmpdir <- "/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/tight/MTBLS1_nmrML_pulProg_missing_spectralMatrix_tight"
# tmpdir <- "/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/sloppy/MTBLS1_nmrML_pulProg_missing_spectralMatrix_sloppy"
# tmpdir <- "/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/std/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS"
# tmpdir <- "/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/tight/MTBLS1_nmrML_pulProg_missing_spectralMatrix_tight"
        
# Get the number of features
  test.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/'
  parsets <- c("sloppy", "std", "tight")
  
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
   
# Get the number of matches (~ rval)

    # same ones, just filtered?
    
# smrfs (backfits)
       
    # same ones, just filtered? 
    
# compound-spectra by score
        
    # same ones, just filtered?
        
        
# For each std run, randomly select n scores 
        

    
    
devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/pars_sens_2/MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS_rescore.RDS'
pars <- yaml::yaml.load_file(paste0(tmpdir,'/params.yaml'), eval.expr = TRUE)
pars$dirs$temp <- tmpdir
pars$debug$enabled <- TRUE
pars$par$ncores <- 4
pars$debug$all.outputs <- TRUE


# Map the results from MTBLS1 maf file using chebis
cmpds <- c('CHEBI:1148','CHEBI:60647','CHEBI:24898','CHEBI:25017','CHEBI:27266','CHEBI:60645','CHEBI:24898','CHEBI:35932','CHEBI:16530','CHEBI:16236','CHEBI:20067','unknown','unknown','CHEBI:30860','unknown','CHEBI:78320','unknown','CHEBI:50129','CHEBI:16449','CHEBI:64390','CHEBI:74911','CHEBI:25017','CHEBI:25094','CHEBI:18257','CHEBI:18211','CHEBI:15366','CHEBI:15366','CHEBI:18257','CHEBI:24898','CHEBI:17533','unknown','CHEBI:21547','unknown','unknown','CHEBI:17533','CHEBI:28300','CHEBI:28300','CHEBI:28300','CHEBI:30772','CHEBI:15347','CHEBI:15344','CHEBI:74903','CHEBI:64390','CHEBI:30744','CHEBI:32816','CHEBI:18261','CHEBI:6650','CHEBI:15741','CHEBI:20067','CHEBI:30915','CHEBI:18237','CHEBI:28300','CHEBI:30769','CHEBI:30769','CHEBI:17170','CHEBI:15611','CHEBI:17724','CHEBI:16628','CHEBI:18139','CHEBI:17724','unknown','unknown','CHEBI:16810','CHEBI:16919','CHEBI:16737','CHEBI:27389','CHEBI:16919','CHEBI:17497','CHEBI:16704','CHEBI:16737','CHEBI:16704','CHEBI:78320','CHEBI:16704','CHEBI:16704','CHEBI:6650','CHEBI:64399','CHEBI:64399','unknown','unknown','CHEBI:15676','CHEBI:16335','CHEBI:18012','CHEBI:27838','CHEBI:27410','unknown','CHEBI:27838','unknown','unknown','unknown','CHEBI:32980','CHEBI:104011','CHEBI:27637','CHEBI:18186','CHEBI:27637','unknown','unknown','CHEBI:27570','CHEBI:32980','CHEBI:43355','CHEBI:64414','CHEBI:18186','CHEBI:43355','CHEBI:28044','CHEBI:18089','CHEBI:18089','CHEBI:27410','CHEBI:18123','CHEBI:18123','CHEBI:15940','unknown','CHEBI:30751','CHEBI:64399','CHEBI:18123','CHEBI:64399','CHEBI:18123')
cmpds <- tolower(cmpds)
cmpds <- unique(cmpds)
data.chebis <- cmpds[!grepl('nknown',cmpds)] 

gissmo.cmpds <- readxl::read_xlsx('/Users/mjudge/Documents/ftp_ebi/gissmo/gissmo_bmrb2chebi.xlsx')
gissmo.chebis.full <- gissmo.cmpds$database_identifier %>% tolower
gissmo.chebis <- gissmo.chebis.full %>% tolower %>% unique %>% na.omit


chebis.matched <- which(gissmo.chebis %in% data.chebis) %>% gissmo.chebis[.]

lapply(gissmo.cmpds$metabolite_identification, function(x){
  rownames(scores.mat)
})

matches <- 
  lapply(chebis.matched, function(x){
    
    match.rows <- gissmo.chebis.full %in% 
    names <- c(gissmo.cmpds$metabolite_identification[match.rows],
               gissmo.cmpds$`Compound Name`[match.rows]) %>% unique
    
    match.in.dataset <- 
    scores <- 
    data.frame(names = paste(names, collapse = ', '),
               chebi = x)
    
  }) %>% do.call(rbind,.)

scores <- readRDS(paste0(tmpdir, '/scores.RDS'))
scores.mat <- scores$ss.ref.mat
rownums <- 1:nrow(scores.mat)

score.summary <- lapply(rownums, function(m){
  
  # Get the top 5 scores for this metabolite
  
    met.scores <- scores.mat[m,]
      score.order <- order(met.scores, decreasing = TRUE)
      
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

# For each of the matched chebis, report the highest mean score



############ Get a list of metabolites for each sample type ###########

  # tmpdir
    tmpdir <- '/Users/mjudge/Documents/ftp_ebi/study_metabolites/'
  
  # Urine
  
  # Blood serum
  
  # Cells

# Try reading in HMDB spectrum ####
  # source('/Users/mjudge/Documents/GitHub/MARIANA_setup_chron/R/readJCAMPDX.R')
  # jd <- readJCAMPDX('/Users/mjudge/Downloads/HMDB0000159_142590_predicted_H_500.jdx', 'HMDB0000159_142590_predicted_H_500')
  # source('~/Documents/GitHub/MARIANA_setup_chron/R/readnmrML.R')
  # jd <- readnmrML('/Users/mjudge/Downloads/500_1H_for_JSviewer_HMDB0000159.nmrML', 'HMDB0000159_142590_predicted_H_500')
   
    
  # write a lib data function for HMDB files
  
# Try reading in processed NPC data ####
  
  # write a lib data function for 1r files

# Backfit cutoff sensitivity gradient ####

    # For a given study
      devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
      tmpdir<- '/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/1697627674'
      pars <- yaml::yaml.load_file(paste0(tmpdir,'/params.yaml'), eval.expr = TRUE)
      pars$dirs$temp <- tmpdir
      pars$debug$enabled <- TRUE
      pars$par$ncores <- 4
      pars$debug$all.outputs <- TRUE
      pars$matching$filtering$max.backfits <- 1E5
      bf.limit <- pars$matching$filtering$max.backfits
      
    # Progressively lower the backfit limit and rescore
    
      # Get the features obj
      
        feature.c <- readRDS(paste0(tmpdir,"/feature.final.RDS"))
      
      # Get the backfits
        
        backfit.results <- readRDS(paste0(tmpdir, "/smrf.RDS"))
        
      # Generate bf files at each bf.limit
      
      bf.limits <- c(1E3, 1E4, 1E5, 1E6, 1E7, 1E8)
        
      # report <- lapply(bf.limits, function(bf.limit){
      selections <- lapply(bf.limits, function(bf.limit){
        
          # Reset data ####
          
              match.info <- backfit.results$match.info
              backfits <- backfit.results$backfits
    
          # Take all the backfits for the maximum bf.limit and calculate their contributions ####
        
              ss.lengths <- feature.c$sfe %>% lapply(function(x){
                x$feat$ss %>% length
              }) %>% unlist
              
              
            
          # Sort matches by rval, then take the top n until # estimated backfits < limit ####
          
              sort.order <- order(match.info$rval, decreasing = TRUE)
              
                match.info <- match.info[sort.order,]
                backfits <- backfits[sort.order]
                
                contribution <- ss.lengths[match.info$feat]
                    est.bfs <- sum(contribution)
                  cum.bfs <- cumsum(contribution)
              
                  cutoff <- max(which(cum.bfs < bf.limit))
                    keep <- 1:cutoff
                      match.info <- match.info[keep, ]
                      backfits <- backfits[keep]
                      
                      lost <- length(sort.order)-length(keep)
              
              
              message('\n\t', lost, ' matches were jettisoned (', round(lost/length(sort.order)*100),'%)')
              message('\n\tThe new effective match rval cutoff is ', min(match.info$rval), '.')
              # message('Saving backfits...\n\n\n')
              
              # bfr <- list(match.info = match.info,
              #             backfits = backfits,
              #             all.succeeded = backfit.results$all.succeeded)
              
              # saveRDS(bfr, paste0(tmpdir,"/smrf_",
              #                     bf.limit %>% formatC(format = "e", digits = 0) %>% 
              #                       toupper %>% stringr::str_remove_all(pattern = "\\+"),
              #                                 ".RDS"))
              
          return(sort.order[keep])
          # return(
          #        data.frame(bf.limit = bf.limit,
          #                   rval.effective = min(match.info$rval),
          #                   n.matches = nrow(match.info),
          #                   n.matches.lost = lost,
          #                   pct.lost = lost/length(sort.order)
          #                   )
          #      )
            
      # }) %>% do.call(rbind,.)              
      })
      
      # plot(report$bf.limit, report$rval.effective, log = 'x', pch = 19, col = 'blue')
      
      lapply(1:length(bf.limits), function(x){
        message(' ---------------------------------------------------------------------------------------------------- ')
        message(bf.limits[x])
        score_matches(pars, 
                      selection = selections[[x]], 
                      alt.name = bf.limits[x] %>% 
                        formatC(format = "e", digits = 0) %>% toupper %>% stringr::str_remove_all(pattern = "\\+")
                      )
        message(' ---------------------------------------------------------------------------------------------------- ')
        message('  ')
      })
      
      
      # Compare to what it would be if selected based on top scoring backfits (if every fit was actually computed)
    
        # rescore()
    
# 424 Tests for parameter sensitivity ####
        
      data <- c('/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/1697627634',
            '/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/1697627674',
            '/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/1697627731',
            '/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/1697627789')
  
      df <- lapply(data, function(x){
        read.csv(paste0(x, '/run.summary.csv'))
      }) %>% do.call(rbind,.)
       
      write.csv(df, paste0('/Users/mjudge/Documents/ftp_ebi/pipeline_tests/MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS_std/run.summary.csv'))

  # Correlations between params and outcomes
      
      fse.pars <- c('protofeatures','did_not.converge_SATs','empty_subset_SATs',
                    'peak_contains.NULL_SATs','reference_degenerated_SATs',
                    'subset_degenerated_SATs','succeeded_SATs','n_feats','med_ss_length')
      match.pars <- c('matches_init','matching_elapsed_time','filtered_matches',
                      'used_matches','match_r_effective')
      bf.pars <- c('backfits','max_score','n_compounds','n_best_bfs','total_time')
      
      test <- fse.pars
      against <- 'storm_r'
      
      test <- c(match.pars, bf.pars)
      against <- 'match_r'
    
      
      layout(matrix(seq(length(test)),nrow=1))
      
      Map(function(x,y) 
        plot(df[c(x,y)], type = "b"), 
          x = names(df[against]), 
          y = names(df[test])
        )
      
      plot(x= df$match_r, y=df$filtered_matches/df$n_feats)
    
    
      
# MWB Tests ####

            
  # Convert downloads to SM
  
      source("./R/readnmrML.R")
      source("./R/readBrukerNMR.R")
      source("./R/readJCAMPDX.R")  
      
      downloads <- '/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies'
      
      fnames <- dir(downloads)
      fnames <- fnames[!grepl('*.zip', fnames)]
      '/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies/NMR raw data/1D'
      
files <- list(
      # Design of Experiments for Maximizing T cell endpoints
      ST001476 = c('/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies/ST001476/Tcell_data/data/batch1', '/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies/ST001476/Tcell_data/data/batch2'),
                  
      # NMR metabolomics on Salmonella enterica
      ST001308 = c('/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies/ST001308/For_upload/raw_processed_NMR_data/pellet_data/nmrpipe/GG_slm_4_jan_800_prj_v2-2/ft_1h1d','/Users/mjudge/Documents/metabolomics_workbench/downloaded_studies/ST001308/For_upload/raw_processed_NMR_data/media_data/nmrpipe/GG_slm_21_dec_media/ft_1h1d')
)
      
      
      filepaths.spectra <- 
      readBrukerNMR(filepaths, rnames)
      
      
      
# Simplified Gradient Tests for parameter sensitivity (19-20OCT) ####
    
    # functionalized indexing ####
    
      devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
      unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
      run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))
      
        # run.idx$write_time %>% as.POSIXct() %>% as.numeric %>% sort %>% .[16] # %>% plot
        # run.idx <- run.idx[(run.idx$write_time %>% as.POSIXct() %>% as.numeric) <= 1698006395, ]
        # run.idx$run_id %>% cat
        
      # run.idx <- run.idx[run.idx$run_id %in% c('1697740712','1697740820','1697740956','1697741076','1697747670','1697751222','1697751875','1697753401','1697793288','1697793355','1697793655','1697793709','1697823986','1697836051','1697842288','1697865804'),]
      run.idx$start <- run.idx$run_id %>% as.POSIXct(origin = "1970-01-01")
      run.idx <- run.idx[run.idx$start > "2023-10-25 22:14:39 EDT", ]
      # run.idx <- run.idx[run.idx$max_backfits == 5E7,]
      
    # Correlations between params and outcomes ####
      
      studies <- run.idx$study %>% unique
      colors <- RColorBrewer::brewer.pal(length(studies), 'Set2')
        
        fse.pars <- c('protofeatures','did_not_converge_SATs','empty_subset_SATs',
                      'peak_contains_NULL_SATs','reference_degenerated_SATs',
                      'subset_degenerated_SATs','succeeded_SATs','n_feats','med_ss_length')
        match.pars <- c('matches_init','matching_elapsed_time','filtered_matches',
                        'used_matches','match_r_effective')
        bf.pars <- c('backfits','max_score','n_compounds','n_best_bfs','total_time')
        
        
        against <- 'storm_r'
        
        test <- c(fse.pars, match.pars, bf.pars)
     
       plots <- lapply(1:length(test), function(x){
       
          # Basic plot
            this.par <- test[x]
            dft <- data.frame(study = run.idx$study,
                              x = run.idx[,against],
                              y = run.idx[,this.par])
          
            
            ggplot(data=dft,
                   aes(x=x, y=y, colour=study)) + 
              geom_line() +
              labs(y = this.par, x = against, title = "") +
              theme_minimal()
       })
       
        pdf(file = paste0(data.dir,against,'_vs_all.pdf'), 
            width = 15, height = 15)
       
          gridExtra::grid.arrange(grobs = plots)
          
        dev.off()
        
        
      
# Backfit cutoff sensitivity (independent runs) ####
  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
  unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
  run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))

  run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))
  run.idx <- run.idx %>% tail(6)
  data <- run.idx$local_path
  
  test <- c('filtered_matches', 'used_matches','match_r_effective','backfits','max_score','n_compounds','total_time','corr')
  
  # pars <- data %>% lapply(function(run){
  #   pars <- yaml::yaml.load_file(paste0(run,'/params.yaml'), eval.expr = TRUE)
  #   run.idx <- pars %>% unlist %>% as.list %>% as.data.frame
  #   run.idx$local.file <- run
  #   run.idx
  # }) %>% bind_rows
  
  mats <- data %>% pblapply(function(run){
    # run <- data[1]
    scores <- readRDS(paste0(run,'/scores.RDS'))
    scores$ss.ref.mat
  })

  key <- run.idx$max_backfits %>% which.max
  key.mat <- mats[[key]]
  # key.mat[key.mat == 0] <- NA
  
  run.idx$corr <- lapply(mats, function(mat){
    # mat[mat == 0] <- NA
    # cor( c(key.mat), c(mat), use = 'pairwise.complete.obs')
    cor( c(key.mat), c(mat), use = 'pairwise.complete.obs')
  }) %>% unlist

  run.idx$r2 <- df$corr ^ 2
  
  runs.info <- df
  layout(matrix(seq(length(test)),nrow=1))

  Map(function(x,y)
    {
      plot(df[c(x,y)], type = "b", log = 'x', ylim = c(0, max(df[y])))
    },
    x = names(df['max_backfits']),
    y = names(df[test])
    )
  par(mfrow=c(1,1))
  plot(x = df$max_backfits, y = df$r2, 
       type = "b", log = 'x', ylim = c(0, max(df$r2)),
       ylab = paste0('r^2 vs. ', max(df$max_backfits)),
       xlab = 'max_backfits')
  
  
  # N compounds is ~ coverage of metabolites ####
  
    bf.number <- 1
    # Loop through spectra and count % coverage ####
    # - allow repeats
     
    
    run <- data[bf.number]
      message('Reading ', run, '...')
      fse.result <- readRDS(paste0(run, '/fse.result.RDS'))
        xmat <- fse.result$xmat
        ppm <- fse.result$ppm
        
        xcov <- matrix(0, nrow(xmat), ncol(xmat))
        
      backfits <- readRDS(paste0(run, '/smrf.RDS'))
        bf.rows <- backfits$backfits %>% rbindlist
        
      # tally each bf
        message('Computing coverage matrix ...')
        for (i in 1:nrow(bf.rows)){
          cols <- bf.rows$spec.start[i]:bf.rows$spec.end[i]
          xcov[bf.rows$ss.spec[i], cols] <- xcov[bf.rows$ss.spec[i], cols] + 1
        }
        
        xmat.numel <- xcov %>% dim %>% prod
        pts.covered <- sum(xcov > 0)
        
        # What about coverage/signal?
        
        # Saturation
        # cor(c(xmat), c(xcov))
        
        coverage.x <- pts.covered/xmat.numel
        max.depth <- max(xcov)
        
        # scattermore::scattermoreplot(x = seq_along(xcov), xlab = 'pts',
        #                              y = sort(xcov), ylab = 'depth')
        
        coverage.mean <- Rfast::colmeans(xcov)
        spectral.mean <- Rfast::colmeans(xmat)
        
        # cor(c(coverage.mean), c(spectral.mean))
        
        message('Plotting coverage spectrum ...')
        library(ggplot2)

        df <- data.frame(y1 = spectral.mean,
                         y2 = coverage.mean,
                         x = ppm)
        
        max.ratio <- max(df$y2)/max(df$y1)

        p <- ggplot(df, aes(x = x)) +
          geom_line(aes(y = y1*max.ratio, color = "y1"), linewidth = 1) +
          geom_line(aes(y = y2, color = "y2"), linewidth = 1) +
          scale_color_manual(values = c("y1" = "red", "y2" = "blue"), labels = c('Spectral Intensity','Coverage Depth')) +
          labs(y = "xcov", x = "1H Chemical Shift ∂ (ppm)") +
          theme_minimal() +
          theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
                legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
          scale_y_continuous(
            name = "Mean Depth",
            sec.axis = sec_axis(~./max.ratio, name = "Spectral Intensity")
          ) + 
          scale_x_reverse() #+ 
          #theme(legend.position = "inside")
        
        # plotly::ggplotly(p)
        
        plot.loc <- run %>% strsplit('/') %>% .[[1]] %>% rev %>% .[-1] %>% rev %>% paste(collapse = '/')
        pdf(file = paste0(plot.loc, '/',
                          runs.info$study[1], 
                          '_bf.',runs.info$max_backfits[bf.number] %>% 
                            formatC(format = "e", digits = 0) %>% toupper %>% 
                            stringr::str_remove_all(pattern = "\\+"),
                          '_spectral.coverage.pdf'),
            width = 8, height = 4)
          p
        dev.off()

       
    # Loop through ref spectra and count % coverage ####
    
      run <- data[1]
      fse.result <- readRDS(paste0(run, '/fse.result.RDS'))
        xmat <- fse.result$xmat
        ppm <- fse.result$ppm
        
        xcov <- matrix(0, nrow(xmat), ncol(xmat))
        
      backfits <- readRDS(paste0(run, '/smrf.RDS'))
        bf.rows <- backfits$backfits %>% rbindlist


    mi <- backfits$match.info
    
    ldp <- readRDS(paste0(run, '/lib.data.processed.RDS'))
    ldp <- ldp %>% lapply(function(x) x %>% expand_ref(ppm))
    
    ref.counts <- matrix(0, nrow = length(ldp), ncol = length(ppm))
    
    mi$counts <- lapply(1:nrow(mi), function(x){
      return(backfits$backfits[[x]] %>% nrow)
    })
    
    for (i in 1:nrow(mi)){
      cols <- mi$ref.start[i]:mi$ref.end[i]
      ref.counts[mi$ref[i], cols] <- ref.counts[mi$ref[i], cols] + 1
    }
    
    # Raw coverage
      sum(ref.counts > 0)/(ref.counts%>%dim%>%prod)
    
    # Since ref spectra are mostly empty, however, it makes sense to calculate
    # % useful signal covered:
    
      refmat <- lapply(ldp, function(x) x$mapped$data) %>% do.call(rbind,.)
        refmat <- refmat / Rfast::rowsums(refmat, na.rm = T)
        refmat[is.na(refmat)] <- 0
        
      ref.counts.normed <- (ref.counts > 0) * refmat
      mean(ref.counts.normed)/mean(refmat)
      
      coverage.mean.ref <- Rfast::colmeans(ref.counts)
      
      plot_spec(coverage.mean.ref, ppm)
    
    
    
    
    
    
  
      
# Top n for each study ####

      devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
      unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
      run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))
      
      # run.idx <- run.idx[run.idx$start > "2023-10-25 22:14:39 EDT", ]
      
      # write.csv(run.idx, file = '/Users/mjudge/Documents/ftp_ebi/param_templates_sensitivity_testing/run.idx.csv')
      
      study <- which(run.idx$run_id %in% '1698347163') # MTBLS1
      # study <- which(run.idx$run_id %in% '1698353488') # MTBLS395
      # study <- which(run.idx$run_id %in% '1698599776') # MTBLS424
      # study <- which(run.idx$run_id %in% '1698408936') # MTBLS430
      results.dir <- run.idx$local_path[study] 

      # Pick the top 10 metabolites
        n <- 5
        scores <- readRDS(paste0(results.dir, '/scores.RDS'))
        
        mean.scores <- lapply(1:nrow(scores$ss.ref.mat), function(x){
          scores$ss.ref.mat[x, ] %>% sort(decreasing = TRUE) %>% .[1:n] %>% mean
        }) %>% unlist
      
        top.n <- order(mean.scores, decreasing = TRUE) %>% .[1:n]
        # rownames(scores$ss.ref.mat) %>% .[top10]
        browse_evidence(results.dir, select.rows = top.n)
        
        
      browse_evidence('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1701255992')
      
# ####

     ## Kick data out for chenomx
     
      
      
     # MTBLS1
     
      run <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1697747670'
      fse.result <- readRDS(paste0(run, '/fse.result.RDS'))
        xmat <- fse.result$xmat
        ppm <- fse.result$ppm

       
# Local test run 1MAR2024 ####
# did fits get fixed?
  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
        
  browse_evidence('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698599655') # MTBLS1 - good
  # browse_evidence('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1697747670') # MTBLS1 .8 - good
        
  # pipeline(params_loc = '/Users/mjudge/Documents/20240229_test_local.yaml') # CE data
  pipeline(params_loc = '/Users/mjudge/Documents/20240229_test_local.yaml')
      # -> browser call before sfe in tina
      #   -> select feature of interest:
      #     while in browser, do the following to see if we have the same number of features:
      #     tmpdir <- '/Users/mjudge/Documents/ftp_ebi/local_outputs/1709318713'
      #     results.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1701255992/'
      #     features.c <- readRDS(paste0(results.dir, "feature.final.RDS"))
      #     -> Loop through and compare the feature ranges
                # currentrngs <- lapply(1:nrow(feature$position), function(x){
                #   # get the range
                #    feature$position[x, ] %>% range(na.rm = TRUE)
                # }) %>% do.call(rbind,.)
                
                f.int <- features.c %>% expand_features(657)
                prevrng <- f.int$position %>% range(na.rm = TRUE)
                which.min(abs(prevrng[1] - currentrngs[, 1]))
                which.min(abs(prevrng[2] - currentrngs[, 2]))
      #     <- [not going to work; different data]
      #     -> 
            features <- readRDS(paste0(tmpdir, "/feature.final.RDS"))    
            length(features$sfe)
                
        
# CAF file ####

  scores <- readRDS('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1699016942/scores.RDS')
  ss.ref.mat <- scores$ss.ref.mat %>% t
    ss.ref.mat[ss.ref.mat >= 0.5] <- 1
    ss.ref.mat[ss.ref.mat < 0.5] <- 0
  heatmap(ss.ref.mat)
    
# MTBLS1 comparison to Chenomx for specific samples ####
  
  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

  res.dir <- run.idx$local_path[nrow(run.idx)]
  scores <- readRDS(paste0(res.dir, '/scores.RDS'))
  samples <- colnames(scores$ss.ref.mat)
    select.sample <- grepl(pattern = 'ADG10003u_007', x = samples) %>% which
  browse_evidence(res.dir, select.sample = select.sample)

  