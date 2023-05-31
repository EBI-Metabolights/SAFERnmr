#' score.matches function
#'
#' This function computes match scores between subset spectra and library reference spectra.
#' Uses pair.score.summation() to actually compute the scores
#' It builds a ss-ref matrix and looks for any compounds that match known annotations.
#'
#' @param pars a list containing necessary parameters for the function
#' @return RDS file containing match scores between library reference spectra and the best subset spectrum score for each.
#' @export
#' 
#' @importFrom magrittr %>%
#'
#'
score.matches <- function(pars){
  
  message('-------------------------------------------------------')
  message('-------------------  Match Scoring  -------------------')
  message('-------------------------------------------------------')
  
  
##################################################################################################################
 # Params and data: ####
    # pars <- yaml::yaml.load_file("../data/params.yaml", eval.expr=TRUE)
    
        tmpdir <- pars$dirs$temp
        this.run <- paste0(tmpdir)

    message("Loading data from files...")

    fse.result <- readRDS(paste0(this.run,"/fse.result.RDS"))
          xmat <- fse.result$xmat
          ppm <- fse.result$ppm
    clusters <- readRDS(paste0(this.run,"/cluster.final.RDS"))
    feature <- readRDS(paste0(this.run,"/feature.final.RDS"))
    lib.info <- readRDS(paste0(this.run,"/lib.info.RDS"))
    
    match.info <- readRDS(paste0(this.run,"/match.info.RDS"))
    peak.qualities <- readRDS(paste0(this.run,"/peak.qualities.RDS"))
    fits.feature <- readRDS(paste0(this.run,"/fits.RDS"))
    
    lib.data.processed <- readRDS(paste0(this.run, "/lib.data.processed.RDS"))
    backfits <- readRDS(paste0(this.run, "/backfits.RDS"))
    
    
    ######################### Build match matrix  #############################    
    # Get the processed library data:
      
      refmat <- lapply(lib.data.processed, function(x) x$mapped$data) %>% do.call(rbind,.)
        
      # refmat <- readRDS(paste0(this.run,"/temp_data_matching/ref.mat.RDS")) %>% t
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
      
      ss.ref.pairs <- lapply(backfits, function(bf) 
        {
          # bf <- backfits[[1]]
          
          # relevant backfit fields are ref.feature-specific - not spectrum specific
          # all fits in a backfit obj share the same ref.region, but differ in bff scores and ss
          # so just use the first one to get that info, and luckily we already extracted bffs. 
          # We do need to loop out the ss.specs.
          
            # Fast way (should be fine) ####
              fit <- bf$fits[[1]]
              pct.ref <- sum(fit$ref.region %>% as.numeric %>% fillbetween %>% refmat[fit$ref, .], na.rm = T)
                # no need to sum the whole spectrum again; already normed to 1.
              
            # return df (expanded this score to all ss.spec x rf combinations) 
              data.frame(match = fit$match, # match # = backfit #
                         ref.start = fit$ref.region[1],
                         ref.end = fit$ref.region[2],
                         ref = fit$ref,
                         # ref.start = fit$ref.region$ref.start,
                         # ref.end = fit$ref.region$ref.end,
                         feat = fit$feat,
                         ss.spec = lapply(bf$fits, function(f) f$ss.spec) %>% unlist,
                         bff.res = bf$bffs.res,
                         bff.tot = bf$bffs.tot,
                         pct.ref = pct.ref)
              
        }) %>% do.call(rbind,.)

      message('Exporting match pair data for scoring ...')
      saveRDS(ss.ref.pairs, paste0(this.run, "/ss.ref.pairs.RDS"))
      # ss.ref.pairs <- readRDS(paste0(this.run, "/ss.ref.pairs.RDS"))
      
  # Turn this into a nonduplicate ss-ref matrix  ####
  
      # Compute scores in separate R instance (cleaner, for parallel): "./pair.score.summation.R"
        
        pair.score.summation(pars, refmat, ss.ref.pairs)
        
        ss.ref.pair.scores <- readRDS(paste0(this.run, "/ss.ref.pair.scores.RDS"))
        rfs.used <- readRDS(paste0(this.run, "/rfs.used.RDS"))
        
      
      # Put into matrix
        
          scores <- ss.ref.pair.scores$score.tot
          # scattermore::scattermoreplot(seq_along(scores), sort(scores))
          # hist(scores, breaks = 1000)
          
          ss.ref.mat.nd <- matrix(0, nrow = nrow(xmat), ncol = max(ss.ref.pairs$ref))
          linds <- sub2indR(rows = ss.ref.pair.scores$ss.spec, 
                            cols = ss.ref.pair.scores$ref, 
                            m = nrow(xmat))
          
          ss.ref.mat.nd[linds] <- scores
      
      # Update the rfs.used list to reflect these matrix positions (naming them)
          
          rfs.used$score.mat.coords$mat.pos <- linds
          
          colnames(ss.ref.mat.nd) <- cmpd.names[1:max(unique(ss.ref.pair.scores$ref))]
          # ss.ref.mat.nd <- ss.ref.mat.nd[, colSums(ss.ref.mat.nd) != 0]
          # ss.ref.mat.nd <- t(ss.ref.mat.nd)
        
          pdf(file = paste0(tmpdir,"/match_scores_sample_x_compound.pdf"),   # The directory you want to save the file in
              width = 8, # The width of the plot in inches
              height = 8) # The height of the plot in inches

            heatmap(ss.ref.mat.nd, scale = 'none') # Rowv = NA, Colv = NA,
            
          dev.off()

        
######################################################################################################
  # Given scores matrix, do any compounds match known annotations?
        # scores.mat <- ss.ref.mat.nd
        # scores.mat <- ss.ref.mat
        
          # lib.data.processed[[47]]$mapped$data %>% t %>% trim.sides %>% simplePlot(linecolor = 'black') + ggtitle(paste0(cmpd.names[47]))
          # scattermore::scattermoreplot(1:length(scores.mat), sort(scores.mat))
          # abline(h = mean(scores.mat))

        # rep.ref <- apply(scores.mat, 2, max)
        
        # See if annotations are in maf file list of annotations ####

                  # matches <- data.frame(gissmo.result = cmpd.names,
                  #                       best.refscore = rep.ref)
                  # matches <- matches[order(matches$best.refscore, decreasing = T),]
                  # 
            # saveRDS(matches, paste0(this.run, "/matches_scored_named.RDS"))
            saveRDS(ss.ref.mat.nd, paste0(this.run, "/ss.ref.sumScores.RDS"))
            saveRDS(rfs.used, paste0(this.run, "/rfs.used.RDS"))
            
            # matches <- readRDS(paste0(this.run, "/matches_scored_named.RDS"))
            
  message('----------------------------------------------------------------')
  message('-------------------  Match Scoring Completed -------------------')
  message('----------------------------------------------------------------')
        
}
         
         