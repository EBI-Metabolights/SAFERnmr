# Generate param files
  devtools::document('/Users/mjudge/Documents/GitHub/SAFER')
      
      template <- '/Users/mjudge/Downloads/backfit_limit_gradient_3/MTBLS1_bf_1E08_params.yaml'
      pars <- yaml::yaml.load_file(template, eval.expr = TRUE)
      
      mat.loc <- '/nfs/production/odonovan/nmr_staging/spectral_matrices/'
      end.loc <- '/nfs/production/odonovan/nmr_staging/safer_runs/'
      # out.dir <- '/Users/mjudge/Documents/ftp_ebi/param_templates_sensitivity_testing/try_4/'
      # out.dir <- '/Users/mjudge/Documents/ftp_ebi/param_templates_sensitivity_testing/backfit_limit_gradient_1/'
      out.dir <- '/Users/mjudge/Documents/ftp_ebi/param_templates_sensitivity_testing/storm_gradient_1/'
      
      studies <- list(
                      list(
                          name = 'MTBLS1',
                          spectral.matrix = 'MTBLS1_1r_noesypr1d_spectralMatrix.RDS',
                          library = '/nfs/production/odonovan/nmr_staging/gissmo_ref/data.list_700MHz.RDS',
                          spec.freq = 700
                          ),
                      list(
                          name = 'MTBLS395',
                          spectral.matrix = 'MTBLS395_1r_cpmgpr1d.comp_spectralMatrix.RDS',
                          library = '/nfs/production/odonovan/nmr_staging/gissmo_ref/data.list_600MHz.RDS',
                          spec.freq = 600
                          ),
                      list(
                          name = 'MTBLS424',
                          spectral.matrix = 'MTBLS424_1r_cpmgpr1d_spectralMatrix.RDS',
                          library = '/nfs/production/odonovan/nmr_staging/gissmo_ref/data.list_600MHz.RDS',
                          spec.freq = 600
                          ),
                      list(
                          name = 'MTBLS430',
                          spectral.matrix = 'MTBLS430_1r_noesygppr1d.comp_spectralMatrix.RDS',
                          library = '/nfs/production/odonovan/nmr_staging/gissmo_ref/data.list_600MHz.RDS',
                          spec.freq = 600
                          )
      )
      
      # Make parsets for each study ####
      
      
      edit_pars_studies <- function(pars, studies)
      {
        pars.studies <- lapply(studies, function(s){
          
          pars$dirs$temp <- paste0(end.loc,s$spectral.matrix)
          # pars$dirs$temp <- end.loc
          pars$study$id <- s$name
          pars$study$spectrometer.frequency <- s$spec.freq
          pars$files$spectral.matrix <- paste0(mat.loc,s$spectral.matrix)
          pars$files$lib.data <- s$library
          return(pars)
        })
        
        return(pars.studies)
        
      }
      
      
      parset <- pars
          parset$corrpockets$rcutoff <- 0.5
          parset$corrpockets$noise.percentile <- 0.99
          parset$corrpockets$half.window <- 0.06
          
          parset$storm$correlation.r.cutoff <- 0.6
          
          parset$tina$min.subset <- 5
          
          parset$matching$r.thresh <- .8
          parset$matching$filtering$ppm.tol <- .1
          parset$matching$filtering$max.backfits <- 5E7
          parset$par$ncores <- 48
          
      pars.6 <- parset
      pars.7 <- parset
        pars.7$storm$correlation.r.cutoff <- 0.7
      pars.8 <- parset
        pars.8$storm$correlation.r.cutoff <- 0.8
      pars.9 <- parset
        pars.9$storm$correlation.r.cutoff <- 0.9

      par.templates <- list(pars.6, pars.7, pars.8, pars.9)
      
      # # Backfit gradient
      # parset <- pars.7
      # pars.3 <- parset
      #   pars.3$matching$filtering$max.backfits <- 1E3
      # pars.4 <- parset
      #   pars.4$matching$filtering$max.backfits <- 1E4
      # pars.5 <- parset
      #   pars.5$matching$filtering$max.backfits <- 1E5
      # pars.6 <- parset
      #   pars.6$matching$filtering$max.backfits <- 1E6
      # pars.7 <- parset
      #   pars.7$matching$filtering$max.backfits <- 1E7
      # pars.8 <- parset
      #   pars.8$matching$filtering$max.backfits <- 1E8
      # 
      # par.templates <- list(pars.3,pars.4,pars.5,pars.6,pars.7,pars.8)
           
    # Apply par changes for different studies 
      par.list <- lapply(par.templates, function(pt){
         
        # edit_pars_studies(pt, studies[1]) # MTBLS1
        edit_pars_studies(pt, studies) # all studies
        
      })
        
      par.list <- par.list %>% unlist(recursive = F)
      
      dir.create(out.dir)
      par.df <- lapply(par.list, function(ps){
        # ps <- par.list[[1]]
        fname <- paste0(ps$study$id, '_r', ps$storm$correlation.r.cutoff, '_params.yaml')
        yaml::write_yaml(ps, paste0(out.dir, fname))
        
        # yaml::write_yaml(ps, paste0(out.dir, ps$study$id, '_bf_', 
        #                             ps$matching$filtering$max.backfits %>% 
        #                               formatC(format = "e", digits = 0) %>% 
        #                               toupper %>% stringr::str_remove_all(pattern = "\\+"), 
        #                             '_params.yaml'))
        message(fname,' - written')
        return(ps %>% unlist %>% as.list %>% as.data.frame)
      }) %>% bind_rows
      
      # lapply(par.list, function(x) x %>% unlist %>% as.list %>% as.data.frame) %>% bind_rows
      names(par.df) <- names(par.df) %>% stringr::str_replace_all('\\.', '_')
      
      dir_pop <- function(d, times=1){d %>% strsplit('/') %>% .[[1]] %>% rev %>% .[-(1:times)] %>% rev %>% paste(collapse = '/')}
      write.csv(par.df, file = paste0(out.dir %>% dir_pop, '/par.table.csv'))


  # Read base set
  
# Read.logs ####

  summ_table_fromLogs <- function(logfiles){
    
    log.list <- lapply(logfiles, function(logfile){
      print(logfile)
      # Read the text file into a character vector
      text <- readLines(logfile)
      start <- which(grepl('local_id', text))
      end <- which(grepl('run_id', text))
  
      # Remove excessive white spaces and split the lines into key-value pairs
      df <- lapply(text[start:end], function(x) {
        # x <- text[1]
        x <- x %>% stringr::str_split('"') %>% .[[1]] %>% lapply(trimws)
        key <- x[[1]]
        value <- x[[2]]
        data.frame(key = key,
                   value = value)
      }) %>% do.call(rbind,.)
      
      # Convert to data frame
      tdf <- as.data.frame(t(df$value))
      colnames(tdf) <- df$key
      
      tdf  
      
    })
    
    log.list %>% dplyr::bind_rows()
    
    
  }

  catsums <- summ_table_fromLogs(dir('/Users/mjudge/Downloads/october_20', full.names = TRUE)) %>% arrange(desc(study))

  write.csv(catsums, paste0('/Users/mjudge/Downloads/october_20', '/summary.csv'))
    

