# Function for annotation matching

# Give metabolite name
# Get the 10 top-scoring browser samples

devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1709742511' # fits mtbls1

study <- 'MTBLS1'
smrf <- readRDS(paste0(tmpdir, "/smrf.RDS"))
head(smrf$match.info, 100) %>% .[c(1,2,3,8),c(1,2,4,8,10)] %>% rbind(tail(smrf$match.info, 1) %>% .[,c(1,2,4,8,10)]) %>% write.csv('matches.1709742511.csv', col.names = TRUE)

cmpd = 'Citrate'

# browse_evidence(tmpdir, clusterSamples = F)
# 
# 
# # Get the sample names
# 
#   fse.result <- readRDS(paste0(tmpdir, "/fse.result.RDS"))
#     xmat <- fse.result$xmat
#     ppm <- fse.result$ppm
#   
#   samples <- rownames(xmat)
  
# Choose best author-scored samples for the metabolite
  maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS1/download/4ZWHUHHlKR?file=m_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
  tmpdir <- '/Users/mjudge/Edison_Lab@UGA Dropbox/Michael Judge/MJ_UGA_Root/Scheduling/safer_manuscript/data/study_metabolites/'
  study.dir <- paste0(tmpdir,study%>%tolower,'/')
  dir.create(study.dir)
  # download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
  maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                         fill = TRUE, header = TRUE, 
                         sep = '\t', quote = "")
  maf.data$metabolite_identification <- maf.data$metabolite_identification %>% tolower

  sample.dir <- paste0(study.dir,'sample_data')
  
  # Find the sample name which best represents the metabolite
  
# Cross-reference the author list with the Chenomx library ####
  # Read in the Chenomx library info ####
    
    # Acquired from T. Payne (EBI); mapped from chenomx_compounds.txt
    chenomx.cmpds <- read.table('/Users/mjudge/Edison_Lab@UGA Dropbox/Michael Judge/MJ_UGA_Root/Scheduling/safer_manuscript/data/chenomx/chenomx2chebi.csv', header = TRUE, sep = ',', fill = TRUE)
    
    chenomx <- data.frame(chebi = chenomx.cmpds[,'database_identifier'],
                          name = chenomx.cmpds[,'metabolite_identification'] %>% tolower) %>%
      filter(name != 'details.xml' & name != '' & chebi != '') %>%
      group_by(chebi) %>%
      summarise(names = paste(name, collapse = " | "))
    
  # Read in the author annotations ####
   
  annot <- read.table('/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/safer_manuscript/data/annotations_MTBLS1_auth_safer_chenomx.scored.csv', header = TRUE, sep = ',', fill = TRUE)
  author <- data.frame(chebi = maf.data$database_identifier,
                       name = maf.data$metabolite_identification %>% tolower) %>% 
    distinct %>% 
    filter(chebi != 'unknown' & chebi != '') %>% 
    group_by(chebi) %>%
    summarise(names = paste(name, collapse = " | ")) %>% as.data.frame
    
  # Intersect the author list with the chenomx library ####
    cna <- match(author$chebi, chenomx$chebi)
    message(sum(!is.na(cna)), '/', nrow(author), ' unique author annotations were present in the Chenomx library.')
    
    author$chenomx.name <- chenomx$names[cna]

  # ####
  best_samples <- function(cmpd, maf.data, key = NULL, top = 1, pause = 2){
  
      
    # Identify the sample cols (quant data cols)
      maf.samples <- maf.data %>% colnames %>% .[19:ncol(maf.data)]
      message('Pulling best ', top, ' sample(s) for ', cmpd, '...')
      
    # Get the data for the rows
      if (stringr::str_detect(cmpd, 'CHEBI')){
        # Use chebi
        cmpd.chebi <- cmpd
        
      } else {
        # Find metabolite in annotation table (names chebi-matched)
          key[is.na(key)] <- ''
          cmpd.chebi <- key$chebi[key$name == cmpd %>% tolower] %>% .[1]
      }
      
      
      quants <- maf.data[maf.data$database_identifier == cmpd.chebi, maf.samples] %>% colMeans(na.rm = TRUE)
        # plot(sort(quants))
        
    # Find the sample with the best average concentration
      q.order <- order(quants, decreasing = TRUE)
      sample.ids <- maf.samples[q.order[1:top]] # %>% paste(collapse = ' | ')
      
    # Pull sample from the web
    
      sample.dir <- sample.dir %>% paste0('/') %>% str_replace("//", "/") # handle trailing / or not
      dir.create(sample.dir,showWarnings = F)
      
      downloads <- lapply(sample.ids, function(sample.id){
        file.name <- paste0(sample.id,'.zip')
        sample.url <- paste0('https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/',study,'/FILES/',file.name)
        if (!(file.name %in% dir(sample.dir, '*.zip'))){
          Sys.sleep(pause)
          download.file(sample.url, destfile = paste0(sample.dir,'/',file.name))
        }
        unzip(paste0(sample.dir,file.name), exdir = sample.dir)
        return(sample.id)
      }) %>% unlist
      
      r.files <- lapply(sample.ids, function(sample.id){
        sample.data <- paste0(sample.dir,sample.id)
        dir(sample.data, pattern = '1r$', recursive = TRUE, full.names = TRUE) %>% paste(collapse = ' | ')
      }) %>% do.call(rbind,.)

      message('The following samples are in ', sample.dir, ':')
      message('\t',downloads %>% paste(collapse = '\n\t'))
      # message('If 1r files were found for the samples, their paths are listed in the output. To use in Chenomx:')
      # message('\t Copy the 1r path.')
      # message('\t Chenomx > Profiler > ["No" to sample spectrum] > \n
      #         Cmd-O > Cmd-Shift-G > [Right click, paste, hit Enter, hit Enter]')
      # system2("open", c("-a", "Chenomx\\ Profiler", shQuote(r.files[[1]])))
      # system2("open", c("-a", "Chenomx\\ NMR\\ Suite", shQuote(r.files[[1]])))
      df <- data.frame(samples = downloads %>% unlist)
      df <- cbind(df, r.files)
      return(df)
      
  }
  
  best <- best_samples(cmpd = 'Citrate', key = author, top = 1, pause = 2)
  
  # Loop through each compound
      
    sample.tab <- pblapply(1:nrow(author), function(x){
      cmpd.row <- author[x, ]
      bs <- best_samples(cmpd = cmpd.row$chebi, maf.data = maf.data, top = 1, pause = 2)
      cmpd.row$sample <- bs$samples[1]
      cmpd.row$spectrum <- bs$r.files[1]
      cmpd.row
    }) %>% do.call(rbind,.)

  
# Pick [random] samples

# do chenomx annotation

# do SAFER annotation

# 