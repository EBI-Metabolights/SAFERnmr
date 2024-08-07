# Alternative scoring 7AUG2024
# 
# The idea here is to incorporate the overall signature shape into scoring. 
# There are many features (spec-features) which degenerately map ref (r) data 
# to dataset spectrum (s) data. 
# Each of these has:
# - rval: feat vs ref
# - fsa: ref vs spec
# - fra: ref explained
# 
# Each backfit has
# - intensity in r
# - intensity in s
# 
# We want to plot these against each other
# can we evaluate lines of best fit to decide on scaling?
# 
# can we use backfit depth to affect the scaling and positions?
# 
# If we get the scaling down, do we also get position?
# 

# Read in some data pre-scoring:

    tmpdir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1722972385'
    pars<-list(par=list(ncores = 4))
    selection <- NULL; alt.name <- ''
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
            
##############                  
                  
        f.nums <- rf.specFits$feat %>% unique
        f.models <- features.c %>% expand_features(f.nums) %>% .[["position"]]
        
        
        # Unlist all the ref feats into spectrum-fit ref feats, and expand their spectrum positions: ####  
          fit.feats <- lapply(1:nrow(rf.specFits), function(x) {
            # x <- 1
            # Get the feature model
              sf <- rf.specFits[x, ]
              
              # Starting from feature and match info:
                feat.model <- f.models[which(sf$feat==f.nums),]
                f.matched.pts <- sf$feat.start:sf$feat.end
                feat.gaps <- feat.model[ f.matched.pts ] %>% is.na
          
            # Get the ref segment (rf)
            
              # Starting from refmat and match info and feature model:
              
                # rf <- ld$mapped$data %>% scale_between(0,1) %>% .[ sf$ref.start:sf$ref.end ]
                rf <- ld$mapped$data %>% .[ sf$ref.start:sf$ref.end ]
              
              # NA-fill the feature gaps
                  
                rf[feat.gaps] <- NA

            # Calculate the fit feature profile 
                
                # plot_fit(list(feat.fit = rf * sf$fit.scale + sf$fit.intercept,
                #               spec.fit = xmat[sf$ss.spec, sf$spec.start:sf$spec.end]),
                #               type = 'auc')
              
                feat.fit <- rf * sf$fit.scale + sf$fit.intercept
                
                return(feat.fit)
          })
            
          
          fit.positions <- lapply(1:nrow(rf.specFits), function(x) 
            {
                rff <- rf.specFits[x,]
              # Get the positions in the xrow
                pos <- rff$spec.start:rff$spec.end
                pos[is.na(fit.feats[[x]])] <- NA
                return(pos)
            }) 
##############            
            
        }) %>% do.call(rbind,.)

      scores <- pair_scores(pars, refmat)jfffjjj
      
      