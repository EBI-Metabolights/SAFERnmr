
# Get the author-submitted annotations for mtbls1

    maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS1/download/4ZWHUHHlKR?file=m_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
    study <- 'MTBLS1'
    tmpdir <- '/Users/mjudge/Edison_Lab@UGA Dropbox/Michael Judge/MJ_UGA_Root/Scheduling/safer_manuscript/data/study_metabolites/'
    study.dir <- paste0(tmpdir,study%>%tolower,'/')
    dir.create(study.dir)
    # download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
    maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                           fill = TRUE, header = TRUE, 
                           sep = '\t', quote = "")
    maf.data$metabolite_identification <- maf.data$metabolite_identification %>% tolower
    
    maf.chebis <- maf.data$database_identifier %>% unique
      maf.chebis <- maf.chebis[-which(grepl('unknown', maf.chebis))] 

# What are the gissmo names that match to the author-annotated names (via common ChEBI) ####
  
      lib.data.700 <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_700MHz.RDS')
      gissmo.cmpds <- readxl::read_xlsx('/Users/mjudge/Documents/ftp_ebi/gissmo/gissmo2chebi_2024.xlsx')
      
      source('/Users/mjudge/Documents/GitHub/MARIANA_setup_chron/R/add_chebiIDs.R') # on "no-zip" branch
      lib.data.700 <- add_chebiIDs(lib.data = lib.data.700, key = gissmo.cmpds)
      
      in.author.list <- lapply(lib.data.700, function(x) x$chebi.id %in% maf.chebis) %>% unlist
      
      lib.data.700.slim <- lib.data.700[in.author.list]
      saveRDS(lib.data.700.slim, file = '/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_700MHz_mtbls1_annots.RDS')
      
      gissmo <- data.frame(name = lib.data.700.slim %>% lapply(function(x) x$compound.name) %>% unlist,
                           chebi = lib.data.700.slim %>% lapply(function(x) x$chebi) %>% unlist)

      
  