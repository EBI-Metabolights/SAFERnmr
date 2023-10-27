#' score_matches function
#'
#' This function computes match scores between subset spectra and library reference spectra.
#' Uses pair_score_summation() to actually compute the scores
#' It builds a ss-ref matrix and looks for any compounds that match known annotations.
#'
#' @param pars a list containing necessary parameters for the function
#' @return RDS file containing match scores between library reference spectra and the best subset spectrum score for each.
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom htmlwidgets saveWidget
#'
#'
score_matches <- function(pars, selection=NULL, alt.name = ''){
  
  # Banner ####
  message('-------------------------------------------------------')
  message('-------------------  Match Scoring  -------------------')
  message('-------------------------------------------------------')
  
  
##################################################################################################################
 # Params and data: ####
    # pars <- yaml::yaml.load_file("../data/params.yaml", eval.expr=TRUE)
    
        tmpdir <- pars$dirs$temp

    message("Loading data from files...")

    fse.result <- readRDS(paste0(tmpdir,"/fse.result.RDS"))
          xmat <- fse.result$xmat
          ppm <- fse.result$ppm
    feature <- readRDS(paste0(tmpdir,"/feature.final.RDS")) %>% expand_features()

    # match.info <- readRDS(paste0(tmpdir,"/match.info.RDS"))

    lib.data.processed <- readRDS(paste0(tmpdir, "/lib.data.processed.RDS"))
    backfit.results <- readRDS(paste0(tmpdir, "/smrf.RDS"))
      match.info <- backfit.results$match.info
      backfits <- backfit.results$backfits
      if (!is.null(selection)){
        match.info <- match.info[selection]
        backfits <- backfits[selection]
      }
      
    ######################### Build match matrix  #############################    
    # Get the processed library data:
      
      if (!dir.exists(paste0(tmpdir, "/temp_data_matching"))){
         unzip(paste0(tmpdir, "/temp_data_matching.zip"), 
              exdir = paste0(tmpdir, "/temp_data_matching"),
              junkpaths = TRUE)
      }
      
      refmat <- readRDS(paste0(tmpdir, "/temp_data_matching/ref.mat.RDS")) %>% cstack_expandRows()


      cmpd.names <- lapply(lib.data.processed, function(x) x$compound.name) %>% do.call(rbind,.)
      # write(cmpd.names,"/Users/mjudge/Documents/ftp_ebi/gissmo/gissmo.cmpd.names.txt", sep = '\t')
      
    # For all subset spectrum - reference pairs, record the % of reference spectrum matched # ####
      message('Building match pair list (indexing matches)...')
      
      # Partition the ss.ref.pair operation by refs matched
      
        # match.info <- match.info[order(match.info$ref), ]
        # ncores <- pars$par$ncores
        # nrefs <- length(unique(match.info$ref))
        # chunk.size <- max(1, nrefs / pars$par$ncores)
        # f.grp <- ceiling((1:nrefs) / chunk.size)
        # 

      ss.ref.pairs <- pblapply(1:nrow(match.info), function(m) 
        {
          # m <- 1
          # Get data for this match
            # print(m)
            bf <- backfits[[m]]
            mi <- match.info[m,]
          
          # relevant backfit fields are ref.feature-specific - not spectrum specific
          # all fits in a backfit obj share the same ref.region, but differ in feature scores and ss
          # so just use the first one to get that info, and luckily we already extracted feature scores 
          # We do need to loop out the ss.specs.
          
            # Fast way (should be fine) ####
              
              # pct.ref <- sum(mi$ref.start:mi$ref.end %>% refmat[mi$ref, .], na.rm = T)
                # no need to sum the whole spectrum again; already normed to 1.
                # this is now done during backfitting.
              
            # return slimmed df (expanded this score to all ss.spec x rf combinations) 
            # - this is just for scoring - needs very little data
              data.frame(match = mi$id, # match # = backfit #
                         ref = mi$ref,
                         feat = mi$feat,
                         feat.start = mi$feat.start,
                         feat.end = mi$feat.end,
                         ref.start = mi$ref.start,
                         ref.end = mi$ref.end,
                         ss.spec = bf$ss.spec,
                         fit.fsa = bf$fit.fsa,
                         fit.rval = bf$fit.rval,
                         match.rval = mi$rval,
                         # rmse = bf$rmse,
                         # rmse.biased = bf$rmse.biased,
                         pct.ref = bf$pct.ref)
              
        }) %>% do.call(rbind,.)

      saveRDS(ss.ref.pairs, paste0(tmpdir, "/ss.ref.pairs.RDS"))
      # ss.ref.pairs <- readRDS(paste0(tmpdir, "/ss.ref.pairs.RDS"))
      
  # Turn this into a nonduplicate ss-ref matrix  ####
  
      # Compute scores in function
        
        scores <- pair_score_summation(pars, refmat) #  refs on rows
        
        ss.ref.pair.scores <- scores$ss.ref.pair.scores
        rfs.used <- scores$rfs.used
        
          # Put into matrix
            
              # scores <- ss.ref.pair.scores$score.rmseb
              vals <- ss.ref.pair.scores$score.fsaxrval
              # scattermore::scattermoreplot(seq_along(vals), sort(vals))
              # hist(scores, breaks = 1000)
              
              ss.ref.mat <- matrix(0, nrow = nrow(xmat), ncol = nrow(refmat))
              linds <- sub2indR(rows = ss.ref.pair.scores$ss.spec, 
                                cols = ss.ref.pair.scores$ref, 
                                m = nrow(xmat))
              
              ss.ref.mat[linds] <- vals
    
          # Add colnames (compounds) to scores matrix 
              
              colnames(ss.ref.mat) <- cmpd.names
              rownames(ss.ref.mat) <- rownames(xmat)
 
        scores$ss.ref.mat <- ss.ref.mat %>% t
        saveRDS(scores, paste0(tmpdir, "/scores",alt.name,".RDS"))
        # scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
        # ss.ref.mat <- scores$ss.ref.mat %>% t
        unlink(paste0(tmpdir, "/ss.ref.pairs.RDS"))
        
        
      # Plot the matrix as an HCA'd heatmap
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

        # Clean up any plotly js dependency files (not needed)
          
          js.dependency.dir <- paste0(pars$dirs$temp, '/match_scores_sample_x_compound_files')
          if (dir.exists(js.dependency.dir)){
            unlink(js.dependency.dir, recursive = TRUE)
          }
          # scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
          # rfs.used <- scores$rfs.used
          # ss.ref.mat <- scores$ss.ref.mat
        
        # Make a caf file ####
        
          # Construct values under sample names ####
          
            ss.ref.mat <- scores$ss.ref.mat %>% t
          
            caf.cmpds <- tryCatch({
              Rfast::colMaxs(ss.ref.mat,value = TRUE) > 0.5
            }, error = function(cond){
              FALSE
            })
            
              
            caf.mat <- ss.ref.mat[, caf.cmpds, drop = FALSE] %>% t

            new.cols <- tryCatch(
              {
                split.names <- colnames(caf.mat) %>% lapply(function(x) x %>% strsplit(split = '\\|') %>% lapply(trimws) %>% unlist)
                
                named.names <- lapply(split.names, function(x){
                  
                  msgs <- NULL
                  
                  if(length(x) != 5){
                    stop(c(msgs,'not enough elements extracted from rownames of xmat. Expected format: "file.number | study.id | sample.name | data.format | pulprog". See spectral matrix conversion.'))
                  }
                  
                  df <- data.frame(processed.file.name = x[1],
                             study.name = x[2],
                             sample.name = x[3],
                             data.format = x[4],
                             pulprog = x[5])
                  
                  if(df$study.name != pars$study$id){
                    msgs <- c(msgs,'study name extracted from xmat rownames does not match pars$study$id')
                  }
                  
                  if(length(msgs) > 0){
                    stop(msgs %>% paste(collapse = '. Also: \n'))
                  } else {
                    df
                  }
                  
                }) %>% do.call(rbind,.)
              },
              error = function(cond){
                message('\nSample names were not parsed successfully. Probable reason(s):')
                message(cond)
                
                df <- as.data.frame(colnames(caf.mat))
                names(df) <- 'sample.name'
                return(df)
              }
            )
            
            colnames(caf.mat) <- new.cols$sample.name
            caf.df <- caf.mat %>% as.data.frame(row.names = 1:nrow(caf.mat))

            caf.df <- cbind(data.frame(metabolite_identification = rownames(caf.mat)), caf.df)
            
          
          # Add chebi_ids and alt_ids (if available) ####
          
            
            chebi.ids <- tryCatch(
              {
                lapply(lib.data.processed[caf.cmpds], function(x) {
                  if (!is.null(x$chebi.id)){
                    x$chebi.id
                  } else {
                    NA
                  }
                  
                }) %>% unlist

              }, error = function(cond){
                return(NA)
              }
            )
            
            caf.df <- cbind(data.frame(database_identifier = chebi.ids), caf.df)
              
          
          # Write the file
            
            if (!any(caf.cmpds)){
                message('\nNo compounds were annotated with evidence score > 0.5 : .caf file will be empty!\n')
            }

            runid <- tmpdir %>% strsplit('/') %>% .[[1]] %>% rev %>% .[1]
            write.table(caf.df, 
                        file = paste0(tmpdir, '/',runid,'.caf.tsv'), 
                        row.names = FALSE, col.names = TRUE, sep = '\t')
            
          
######################################################################################################
  # ####       
  message('----------------------------------------------------------------')
  message('-------------------  Match Scoring Completed -------------------')
  message('----------------------------------------------------------------')
        
  return(
    data.frame(
      score.metric = 'fsaxrval',
      max.score = max(ss.ref.mat),
      n.compounds = sum(ann.cmpds),
      n.best.bfs = rfs.used$fsaxrval %>% unlist %>% length
    )
  )
}
         
         