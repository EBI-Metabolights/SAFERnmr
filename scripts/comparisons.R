# Function for annotation matching

# Give metabolite name
# Get the 10 top-scoring browser samples

devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655
tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1709742511' # fits mtbls1

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
  study.dir <- paste0(tmpdir,'mtbls1/')
  dir.create(study.dir)
  # download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
  maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                         fill = TRUE, header = TRUE, 
                         sep = '\t', quote = "")

  sample.dir <- paste0(study.dir,'sample_data')
  
  # Find the sample name which best represents the metabolite
  
  best_sample <- function(cmpd){
  
        # Find metabolite in annotation table (names chebi-matched)
      annot <- read.table('/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/safer_manuscript/data/annotations_MTBLS1_auth_safer_chenomx.scored.csv', header = TRUE, sep = ',', fill = TRUE)
        annot[is.na(annot)] <- ''
        
    # Identify the sample cols (quant data cols)
      maf.samples <- maf.data %>% colnames %>% .[19:ncol(maf.data)]
        
    # Get the data for the rows
      cmpd.chebi <- annot$chebi[annot$auth_name == cmpd %>% tolower] %>% .[1]
      quants <- maf.data[maf.data$database_identifier == cmpd.chebi, maf.samples] %>% colMeans(na.rm = TRUE)

    # Find the sample with the best average concentration
      sample.id <- which.max(quants) %>% maf.samples[.]
      
    # Record in index which sample was used
      
      
      
    # Pull sample from the web
      
      file.name <- paste0(sample.id,'.zip')
      dir.create(sample.dir,showWarnings = F)
      sample.url <- paste0('https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS1/FILES/',file.name)
      download.file(sample.url, destfile = paste0(sample.dir,file.name))
      unzip(paste0(sample.dir,file.name))
      # system2("open", c("-a", "Chenomx\\ Profiler", paste0(sample.dir,sample.id)))
      # system2("open", c("-a", "Chenomx\\ NMR\\ Suite", paste0(sample.dir,sample.id)))
      
      message(sample.id)
    
  }
      
      
      
      
# Pick [random] samples

# do chenomx annotation

# do SAFER annotation

# 