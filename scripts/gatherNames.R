devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')


# For each study, note the metabolites id'd ####
    
# tmpdir <- '/Users/mjudge/Documents/ftp_ebi/study_metabolites/'
tmpdir <- '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/safer_manuscript/data/study_metabolites/'
    
    # Load MTBLS1 ####
      
      maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS1/download/4ZWHUHHlKR?file=m_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
      study.dir <- paste0(tmpdir,'mtbls1/')
      dir.create(study.dir)
      download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
      maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                             fill = TRUE, header = TRUE, 
                             sep = '\t', quote = "")
      compounds <- data.frame(db.id = maf.data$database_identifier,
                              name = maf.data$metabolite_identification) %>% distinct
      compounds$study <- 'MTBLS1'
      
      cmpd.list <- list(compounds)
      
    # Load MTBLS395 ####
    
      maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS395/download/92353320-750c-40e8-a9af-031af8d7f8f4?file=m_MTBLS395_amiflorenceii_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
      study.dir <- paste0(tmpdir,'mtbls395/')
      dir.create(study.dir)
      download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
      maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                             fill = TRUE, header = TRUE,
                             sep = '\t')
      compounds <- data.frame(db.id = maf.data$database_identifier,
                              name = maf.data$metabolite_identification) %>% distinct
      compounds$study <- 'MTBLS395'
      
      # Get Pantelis' annotations 
        pt.ids <- readxl::read_excel(paste0(tmpdir,'pt/List_of_metabolites.xlsx'),col_names = FALSE)
        cmpd.names <- lapply(1:nrow(pt.ids), function(x) pt.ids[x,1] %>% as.character) %>% unlist
        pt.cmpds <- lapply(1:nrow(pt.ids), function(x) {
          chebis <- pt.ids[x,-1] %>% as.character()
          
          # If any of the chebis for the pt annotation match, just use that
          matched <- chebis %in% compounds$db.id
          if (any(matched)){
            # Already present, use that one
            chebis <- chebis[matched]
          } else {
            # If no match, use the first chebi (not L or R; NA if no chebi)
            chebis <- chebis[1]
          }
          
          data.frame(db.id = chebis,
                     name = cmpd.names[x],
                     study = 'MTBLS395')

        }) %>% do.call(rbind,.)
        
        
        compounds <- rbind(compounds, pt.cmpds)
        compounds <- compounds[!duplicated(compounds$db.id),]
        
      cmpd.list[[2]] <- compounds
      
    # Load MTBLS424 ####
    
      maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS424/download/3ca6dc8a-0f60-421a-b0af-0d49bd12420c?file=m_MTBLS424_breast_cancer_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
      study.dir <- paste0(tmpdir,'mtbls424/')
      dir.create(study.dir)
      download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
      maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                             fill = TRUE, header = TRUE,
                             sep = '\t')
      compounds <- data.frame(db.id = maf.data$database_identifier,
                              name = maf.data$metabolite_identification) %>% distinct
      compounds$study <- 'MTBLS424'
      cmpd.list[[3]] <- compounds
      
    # Load MTBLS430 ####
    
      maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS430/download/9d17bf2d-18b7-41ce-9446-ee600aef4f5a?file=m_MTBLS430_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
      study.dir <- paste0(tmpdir,'mtbls430/')
      dir.create(study.dir)
      download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
      maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                             fill = TRUE, header = TRUE,
                             sep = '\t')
      compounds <- data.frame(db.id = maf.data$database_identifier,
                              name = maf.data$metabolite_identification) %>% distinct
      compounds$study <- 'MTBLS430'
      cmpd.list[[4]] <- compounds
      
        
    
    
      
    # Combine into single df ####
      cmpd.list <- do.call(rbind,cmpd.list) %>% distinct
      cmpd.list <- cmpd.list %>% filter(db.id != "" & db.id != "unknown")
      # cmpd.list$name[!duplicated(cmpd.list$db.id)] 
      saveRDS(cmpd.list, file = paste0(tmpdir,'all.metabolites.RDS'))
      
    # Make venn diagram ####
    # Load library
    library(VennDiagram)
     
    # Prepare palette:
    library(RColorBrewer)
    myCol <- brewer.pal(length(unique(cmpd.list$study)), "Pastel2")
    
    venn.diagram(
            x = cmpd.list$db.id %>% split(cmpd.list$study),
            category.names = unique(cmpd.list$study),
            filename = paste0(tmpdir,'cmpds_venn_diagramm.png'),
            output=TRUE,
            
            # Output features
            imagetype="png" ,
            height = 480*2 , 
            width = 480*2 , 
            resolution = 300,
            compression = "lzw",
            
            # Circles
            lwd = 2,
            lty = 'blank',
            fill = myCol,
            
            # Numbers
            cex = .6,
            fontface = "bold",
            fontfamily = "sans",
            
            # Set names
            cat.cex = 0.6,
            cat.fontface = "bold",
            cat.default.pos = "outer",
            # cat.pos = c(-27, 27, 135),
            # cat.dist = c(0.055, 0.055, 0.085),
            cat.fontfamily = "sans",
            # rotation = 1
    )
        
          
          

  # Which compounds are in the 600 and 700 MHz Gissmo data? #####
  
    lib.data.600 <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_700MHz.RDS')
    lib.data.700 <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_600MHz.RDS')
    
    gissmo.cmpds <- readxl::read_xlsx('/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/safer_manuscript/data/gissmo/gissmo_bmrb2chebi.xlsx')
    
    source('/Users/mjudge/Documents/GitHub/MARIANA_setup_chron/R/add_chebiIDs.R')
    lib.data.600 <- add_chebiIDs(lib.data = lib.data.600, key = gissmo.cmpds)
      
    
    gissmo.chebis.full <- gissmo.cmpds$database_identifier %>% tolower
    gissmo.chebis <- gissmo.chebis.full %>% toupper %>% unique %>% na.omit
    
    # cmpds <- lapply(lib.data.600, function(x) x$compound.name)
    
    # data.chebis <- cmpds[!grepl('nknown',cmpds)] 
    data.chebis <- lib.data.600 %>% lapply(function(x) x$chebi.id) %>% unlist
    data.chebis <- cmpd.list$db.id %>% unique

    chebis.needed <- (!(data.chebis %in% gissmo.chebis)) %>% data.chebis[.]
    chebis.matched <- (gissmo.chebis %in% data.chebis) %>% gissmo.chebis[.]
    
    !(data.chebis %in% gissmo.chebis)
    
  # Match to author reported metabolites ####
  
    intersect(cmpd.list$db.id[cmpd.list$study=='MTBLS1'], chebis.matched) %>% length
    
    
  # Match to specific result ####
    # 
    # md <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/par_sens_oct_3/std/'
    # tmpdir <- paste0(md,'MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS_std/')
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/1697097481.99655/'
    # scores <- readRDS(paste0(tmpdir, 'scores.RDS'))
    # 
    # scores.mat <- scores$ss.ref.mat
    # data.cmpds <- scores.mat %>% rowSums %>% order(decreasing = TRUE) %>% scores.mat[., ] %>% rownames %>% unique
    # dma <- which(data.cmpds %in% cmpd.list$name)
    #   dma.score <- scores.mat %>% rowSums %>% .[dma]
    #   dma <- data.frame(name)
    #   plot(x = )
    # author.ann.found.in.top.10 <- which(data.cmpds %in% cmpd.list[cmpd.list$study == 'MTBLS424',]$db.id)
    # 
    # lapply(gissmo.cmpds$metabolite_identification, function(x){
    #   rownames(scores.mat)
    # })
    # 
    # matches <- 
    #   lapply(chebis.matched, function(x){
    #     
    #     match.rows <- gissmo.chebis.full %in% 
    #     names <- c(gissmo.cmpds$metabolite_identification[match.rows],
    #                gissmo.cmpds$`Compound Name`[match.rows]) %>% unique
    #     
    #     match.in.dataset <- 
    #     scores <- 
    #     data.frame(names = paste(names, collapse = ', '),
    #                chebi = x)
    #     
    #   }) %>% do.call(rbind,.)
    # 
    # scores <- readRDS(paste0(tmpdir, '/scores.RDS'))
    # scores.mat <- scores$ss.ref.mat
    # rownums <- 1:nrow(scores.mat)

