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
#'
#'
score_matches <- function(pars){
  
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
              
              ss.ref.mat.nd <- matrix(0, nrow = nrow(xmat), ncol = nrow(refmat))
              linds <- sub2indR(rows = ss.ref.pair.scores$ss.spec, 
                                cols = ss.ref.pair.scores$ref, 
                                m = nrow(xmat))
              
              ss.ref.mat.nd[linds] <- vals
    
          # Add colnames (compounds) to scores matrix 
              
              colnames(ss.ref.mat.nd) <- cmpd.names
 
        scores$ss.ref.mat <- ss.ref.mat.nd %>% t
        saveRDS(scores, paste0(tmpdir, "/scores.RDS"))
        # scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
        # ss.ref.mat.nd <- scores$ss.ref.mat %>% t
        unlink(paste0(tmpdir, "/ss.ref.pairs.RDS"))
        
      # Plot the matrix as an HCA'd heatmap
        ann.cmpds <- tryCatch({
          Rfast::colsums(ss.ref.mat.nd) > 0
        }, error = function(cond){
          TRUE
        })
        
        tryCatch(
          {
            pdf(file = paste0(tmpdir,"/match_scores_sample_x_compound.pdf"),   # The directory you want to save the file in
                width = 8, # The width of the plot in inches
                height = 8) # The height of the plot in inches
              ss.ref.mat.nd <- ss.ref.mat.nd[, ann.cmpds]
              heatmap(ss.ref.mat.nd, scale = 'none') # Rowv = NA, Colv = NA,
              
            dev.off()
            
          },
          error = function(cond){
            message('\n\tScore matrix plotting error..')
          }
        )

          
          # scores <- readRDS(paste0(tmpdir, "/scores.RDS"))
          # rfs.used <- scores$rfs.used
          # ss.ref.mat.nd <- scores$ss.ref.mat
        
######################################################################################################
  # Given scores matrix, do any compounds match known annotations? ####
        # scores.mat <- ss.ref.mat.nd
        # scores.mat <- ss.ref.mat
        
          # lib.data.processed[[47]]$mapped$data %>% t %>% trim_sides %>% simplePlot(linecolor = 'black') + ggtitle(paste0(cmpd.names[47]))
          # scattermore::scattermoreplot(1:length(scores.mat), sort(scores.mat))
          # abline(h = mean(scores.mat))

        # rep.ref <- apply(scores.mat, 2, max)
        
        # See if annotations are in maf file list of annotations

                  # matches <- data.frame(gissmo.result = cmpd.names,
                  #                       best.refscore = rep.ref)
                  # matches <- matches[order(matches$best.refscore, decreasing = T),]
                  # 
            # saveRDS(matches, paste0(tmpdir, "/matches_scored_named.RDS"))
            # matches <- readRDS(paste0(tmpdir, "/matches_scored_named.RDS"))
  # ####       
  message('----------------------------------------------------------------')
  message('-------------------  Match Scoring Completed -------------------')
  message('----------------------------------------------------------------')
        
  return(
    data.frame(
      score.metric = 'fsaxrval',
      max.score = max(ss.ref.mat.nd),
      n.compounds = sum(ann.cmpds),
      n.best.bfs = rfs.used$fsaxrval %>% unlist %>% length
    )
  )
}
         
         