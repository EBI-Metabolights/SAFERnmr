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
    # lib.data.600 <- add_chebiIDs(lib.data = lib.data.600, key = gissmo.cmpds)
    lib.data.700 <- add_chebiIDs(lib.data = lib.data.700, key = gissmo.cmpds)
      
    
    gissmo.chebis.full <- gissmo.cmpds$database_identifier %>% tolower
    gissmo.chebis <- gissmo.chebis.full %>% toupper %>% unique %>% na.omit
    
    # cmpds <- lapply(lib.data.600, function(x) x$compound.name)
    
    # data.chebis <- cmpds[!grepl('nknown',cmpds)] 
    # data.chebis <- lib.data.600 %>% lapply(function(x) x$chebi.id) %>% unlist
    data.chebis <- lib.data.700 %>% lapply(function(x) x$chebi.id) %>% unlist
    # data.chebis <- cmpd.list$db.id %>% unique

    chebis.needed <- (!(data.chebis %in% gissmo.chebis)) %>% data.chebis[.]
    chebis.matched <- (gissmo.chebis %in% data.chebis) %>% gissmo.chebis[.]
    
    matched.names <- (data.chebis %in% gissmo.chebis) %>% lib.data.700[.] %>% lapply(function(x) x$compound.name) %>% unlist %>% unique
    
  # Match to author reported metabolites ####
  
    # cmpd.list$db.id[cmpd.list$study=='MTBLS1'] %>% intersect(chebis.matched)
    
    
  # Match to specific result ####
    
    md <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/par_sens_oct_3/std/'
    tmpdir <- paste0(md,'MTBLS1_nmrML_pulProg_missing_spectralMatrix.RDS_std/')
    # tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/1697097481.99655/'
    scores <- readRDS(paste0(tmpdir, 'scores.RDS'))

    # What are the gissmo names that match to the author-annotated names (via common ChEBI)
    
    scores.mat <- scores$ss.ref.mat
    # * the rownames have added whitespace
    data.cmpds <- scores.mat %>% rownames %>% tolower %>% stringr::str_remove_all(pattern = '  ') %>% unique
    cmpd.list$name <- cmpd.list$name %>% tolower
    
    dma <- which(data.cmpds %in% cmpd.list$name)
      dma.score <- scores.mat %>% rowSums %>% .[dma]
      dma <- data.frame(name)
      
    author.ann.found.in.top.10 <- which(data.cmpds %in% cmpd.list[cmpd.list$study == 'MTBLS424',]$db.id)

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

## Which compounds had evidence in MTBLS1? ####

    # What are the gissmo names that match to the author-annotated names (via common ChEBI) ####
        
        lib.data.700 <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/data.list_600MHz.RDS')
        gissmo.cmpds <- readxl::read_xlsx('/Users/mjudge/Documents/ftp_ebi/gissmo/gissmo2chebi_2024.xlsx')
        
        source('/Users/mjudge/Documents/GitHub/MARIANA_setup_chron/R/add_chebiIDs.R')
        lib.data.700 <- add_chebiIDs(lib.data = lib.data.700, key = gissmo.cmpds)
        gissmo <- data.frame(name = lib.data.700 %>% lapply(function(x) x$compound.name) %>% unlist,
                             chebi = lib.data.700 %>% lapply(function(x) x$chebi) %>% unlist)

    # Get ChEBIs for the run ####
    
      # md <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs/par_sens_oct_3/std/'
      tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1698599655/'
      scores <- readRDS(paste0(tmpdir, 'scores.RDS'))
    
      scores.mat <- scores$ss.ref.mat
      
      score.names <- scores.mat %>% rowSums %>% ">"(.,0) %>% 
        scores.mat[.,] %>% rownames %>% stringr::str_remove_all(pattern = '  ')
      # (gissmo.names %in% score.names) %>% sum
      
      score.chebis <- score.names %>% lapply(function(x) (gissmo$name %in% x) %>% which %>% gissmo$chebi[.] %>% unique) %>% unlist
      
      safer.annots <- data.frame(chebi = score.chebis,
                                 name = score.names)
      
    # Get ChEBIs for the author annotations ####
    
      tmpdir <- '/Users/mjudge/Dropbox (Edison_Lab@UGA)/MJ_UGA_Root/Scheduling/safer_manuscript/data/study_metabolites/'
    
      # maf.link <- 'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS1/download/4ZWHUHHlKR?file=m_MTBLS1_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
      study.dir <- paste0(tmpdir,'mtbls1/')
      # dir.create(study.dir, warn = FALSE)
      # download.file(maf.link, destfile = paste0(study.dir,'maf.tsv'))
      maf.data <- read.table(paste0(study.dir,'maf.tsv'), 
                             fill = TRUE, header = TRUE, 
                             sep = '\t', quote = "")
      
      compounds <- data.frame(chebi = maf.data$database_identifier,
                              name = maf.data$metabolite_identification)# %>% distinct
      compounds <- compounds %>% filter(chebi != '' & chebi != 'unknown')
      compounds$study <- 'MTBLS1'
      auth.annots <- compounds
    
    # Intersect them
      
      auth.n.safer <- auth.annots$chebi %in% safer.annots$chebi
      safer.n.auth <- safer.annots$chebi %in% auth.annots$chebi
        auth.n.safer %>% auth.annots$name[.] %>% unique
    
      combined <- data.frame(chebi = safer.annots$chebi %>% c(., auth.annots$chebi) %>% unique)
      
        combined$author <- lapply(combined$chebi, function(x) 
          {
            (auth.annots$chebi == x) %>% auth.annots$name[.] %>% unique %>% paste(collapse=' | ')
          }
        ) %>% unlist

        combined$safer <- lapply(combined$chebi, function(x) 
          {
            (safer.annots$chebi == x) %>% safer.annots$name[.] %>% unique %>% paste(collapse=' | ')
          }
        ) %>% unlist

              
    # Print the author annotation name and the name of the run compound name
      
      matched.names <- combined %>% filter(author != "" & safer != "")
      
    
    