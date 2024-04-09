devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

tmp <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1712450825'

# browse_evidence(tmp)

## Setup ####
  
  # for the selectizer:
  sort.choices <- c("Scores - cluster", "Search for compound...") # "Compound names - cluster", 'Compound names - alphabetical', 
  
  results.dir <- tmp
  select.rows = NULL
## Params ###############################################################################
  
  message('Evidence Viewer: reading data...')
  # Locate us ####
  
    # Handle whether or not user adds /
      if (stringr::str_sub(results.dir, start= -1) != "/"){
        results.dir <- paste0(results.dir, '/')
      }

  pars <- yaml::yaml.load_file(paste0(results.dir,'params.yaml'), eval.expr = TRUE)
  
  # Apply same validation function as in the pipeline (so pars are interpreted the same throughout)
    par.list <- pars %>% valid_pars(quiet = TRUE)
    pars <- par.list$pars
  
  study <- pars$study$id
  ppm.tolerance = pars$matching$filtering$ppm.tol
  hshift = 0
  

##########################     Setup/Read Data    ####################################

    # Read in spectral matrix data
        fse.result <- readRDS(paste0(results.dir, "fse.result.RDS"))
          xmat <- fse.result$xmat
          ppm <- fse.result$ppm
          rm(fse.result)

    # Read in the features 
    
        features.c <- readRDS(paste0(results.dir, "feature.final.RDS"))

    # Read in scores matrix 
      scores <- readRDS(paste0(results.dir,"scores.RDS"))
        scores.matrix <- scores$ss.ref.mat
        colnames(scores.matrix) <- 1:ncol(scores.matrix)
        
          if (is.null(select.rows)) {
            select.rows <- 1:nrow(scores.matrix)
          }
        
    # Read in match data
        backfit.results <- readRDS(paste0(results.dir,"smrf.RDS"))
          backfits <- backfit.results$backfits
          match.info <- backfit.results$match.info
          

    # Read in library data
        lib.data.processed <- readRDS(paste0(results.dir, "lib.data.processed.RDS"))
                
          # ppm isn't needed anymore; using spectral matrix ppm.
            
            lib.data.processed <- lib.data.processed %>% lapply(function(x) {x$data <- NULL; x$ppm <- NULL; return(x)})
            
          # add compound names as column in match.info
            compound.names.refmat <- lapply(lib.data.processed, function(x) x$compound.name) %>% unlist
            match.info$ref.name <- match.info$ref %>% compound.names.refmat[.]
          
      rfs.used <- scores$rfs.used
      
##########################     Filter    ####################################          
          
      # Remove all compounds without any scores > 0.5

        keeprefs <-  intersect(which(apply(scores.matrix, 1, max) > 0), select.rows)
        keepsamples <-  apply(scores.matrix, 2, max) > 0
        refs.used <- keeprefs
        compound.names.refmat <- compound.names.refmat[refs.used]
        samples.used <- which(keepsamples)
        
        scores.matrix <- scores.matrix[keeprefs,,drop=F]
        scores.matrix <- scores.matrix[,keepsamples,drop=F]
        lib.data.processed <- lib.data.processed[keeprefs]
        
        fits.keep <- match.info$ref %in% refs.used
        match.info <- match.info[fits.keep, ]
        backfits <- backfits[fits.keep]

##########################     Cluster    ####################################      

      
      if (nrow(scores.matrix) > 1){
        clust.refs <- T
        clust.samples <- T
      } else {
        clust.refs <- F
        clust.samples <- T
      }

      # hclust: column subsetting doesn't work when names are used. ### #
      ref.order <- 1:nrow(scores.matrix)
      sample.order <- 1:ncol(scores.matrix)
      
      if (clust.refs){
          ref.order <- hclust(dist(scores.matrix))$order
      }
      if (clust.samples){
          sample.order <- hclust(dist(t(scores.matrix)))$order
      }
      
      refs <- data.frame(number = seq_along(refs.used), # this is the initial row number (before sort). lib info matches this. 
                         id = refs.used %>% as.numeric, # this is the ref number in the full library (also match.info$ref)
                         name = compound.names.refmat) # name
        refs <- refs[ref.order, ]
        refs$row.mat <- 1:nrow(refs)         # this is the row number in mat

      samples <- data.frame(number = seq_along(samples.used),              # number is column number upon sort
                            id = samples.used %>% as.numeric,              # id is the sample number upon import
                            name = 1:ncol(scores.matrix) %>% as.character) # just use column number for now
        samples <- samples[sample.order, ]
        samples$col.mat <- 1:nrow(samples)                    

      scores.matrix<- scores.matrix[ref.order, ,drop =F]
      scores.matrix<- scores.matrix[, sample.order, drop =F]
      mat <- scores.matrix
      
      

# Plot the saved result using browse_evidence() code ####
# select evidence 
# devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
selectedRow <- which(refs$name %in% 'L-glutamine')
region <- c(2.0, 2.9)
selectedCols <- 35#:24
ss.spec <- selectedCols %>% samples$id[.]
metab.evidence <- select_evidence_refmet(ref = selectedRow %>% refs[.,],
                                         sample = selectedCols %>% samples[.,],
                                         # Big objects to subset using ref.ind:
                                           features.c = features.c,
                                           match.info = match.info,
                                           backfits = backfits,
                                           rfs.used = rfs.used, # new object with inds of rfs which contributed to scores
                                           lib.data.processed = lib.data.processed,
                                         # # Spectral data:
                                         #   xmat, 
                                           ppm = ppm)


# Is the correct spectrum being plotted?

rf.fits <- metab.evidence$rf.fits   
in.range <- TRUE

bfs <- list(fit.feats = rf.fits$fit.feats[in.range, , drop = F],
            fit.positions = rf.fits$fit.positions[in.range, , drop = F],
            fit.xrow = rf.fits$fit.xrow[in.range],
            pass.fit = T)

# look at one fit
  x <- 3
  mi <- metab.evidence$match.info.ss[x, ]
  bfs$fit.feats <- bfs$fit.feats[x, ,drop = FALSE]
  bfs$fit.positions <- bfs$fit.positions[x, ,drop = FALSE]
  bfs$fit.xrow <- bfs$fit.xrow[x]

  simplePlot(rbind(bfs$fit.feats, 
                   xmat[bfs$fit.xrow, bfs$fit.positions]))

plt.pars <- list(vshift = 1, 
                 pixels = c(512, 512), # inc.res
                 pointsize = 0, 
                 interpolate = T, 
                 # exp.by = 0.05,
                 xlim = c(2.0,2.6))

# feature_est_plot(reg = plt.pars$xlim,
#                  metab.evidence,
#                  features.c,
#                  ppm,
#                  plt.pars)

# devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
 fastStack.withFeatures(xmat, ppm, raster = T, bfs = bfs, plt.pars, res.ratio = .1)
 
 
# Look at feature match to ref ####
 
  # we need something that has the same feature (657), 
  # and the same reference (236 (when using unsorted refs)), 
  # and the same ss.spec (5)
 
  # f <- 591
  f <- 657
  feature <- features.c %>% expand_features(f)
  
  ref <- selectedRow %>% refs[.,]
  ld <- lib.data.processed[[ref$number]] %>% expand_ref(ppm)
  
  mi.rows <- which(match.info$feat == f)
  mi <- match.info[mi.rows[1], ]
  
  feat <- feature$stack[mi$feat.start:mi$feat.end]
    simplePlot(feat)
  ref <- ld$mapped$data[mi$ref.start:mi$ref.end]
    simplePlot(ref)
  
  # simplePlot(rbind(ref, feat))
  # simplePlot(rbind(ref, mi$fit.scale * (feat - mi$fit.intercept)))
  # simplePlot(rbind(ref, mi$fit.scale*feat + mi$fit.intercept))
  # plot_spec(ref, ppm = ppm[mi$ref.start:mi$ref.end])
  # plot_spec(feat, ppm = feature$position %>% ppm[.])
  
  plot_fit(list(feat.fit = ref,
                spec.fit = mi$fit.scale*feat[mi$feat.start:mi$feat.end] + mi$fit.intercept),
           type = 'auc')
  # fit_batman(ref, feat, exclude.lowest = 0.5) %>% plot_fit(type = 'auc')
  # fit_leastSquares(ref, feat, plots = F, scale.v2 = T) %>% plot_fit(type = 'auc')
  
# Troubleshoot the individual backfit ####
  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
  refmat.c <- readRDS('/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/1700653867/temp_data_matching/ref.mat.RDS')

  ss.spec <- 5
  mi.rows <- which(match.info$feat == f)
  mi <- match.info[mi.rows[1], ]
  
  backfit.results <- backfit_rfs3(match.info = mi,
                                  feature.c = features.c, # has sfe data
                                  xmat = xmat,
                                  refmat.c = refmat.c,
                                  ncores = 1)
  
  mi <- backfit.results$match.info
  feature <- features.c %>% expand_features()
  ref.reg <- refmat.c %>% cstack_selectRows(mi$ref) %>% cstack_expandRows %>% .[, mi$ref.start:mi$ref.end, drop = FALSE]
  
  bfres <- backfit.results$backfits[[1]] # correspond with row #s of match info
  
    bf <- bfres[which(bfres$ss.spec == ss.spec), ]
    
    bfs <- list(fit.feats = bf$fit.scale*ref.reg + bf$fit.intercept,
                fit.positions = matrix(bf$spec.start:bf$spec.end, nrow = 1),#feature$position[mi$feat, mi$feat.start:mi$feat.end,drop = FALSE],
                fit.xrow = ss.spec)

    fastStack.withFeatures(xmat = xmat, ppm = ppm, raster = TRUE, bfs = bfs, plt.pars = plt.pars)
    
# Check that the newly computed backfits are correctly reconstructed by select_evidence ####

    
selectedRow <- which(refs$name %in% 'Glutaric-acid')
  ss.spec <- 5
  mi.rows <- which(match.info$feat == f)
  mi <- match.info[mi.rows[1], ]
  
  backfit.results <- backfit_rfs3(match.info = mi,
                                  feature.c = features.c, # has sfe data
                                  xmat = xmat,
                                  refmat.c = refmat.c,
                                  ncores = 1)

region <- c(1.676, 2.273)
selectedCols <- 59#:24
metab.evidence <- select_evidence_refmet(ref = selectedRow %>% refs[.,],
                                         sample = selectedCols %>% samples[.,],
                                         # Big objects to subset using ref.ind:
                                           features.c = features.c,
                                           match.info = mi,
                                           backfits = backfit.results$backfits,
                                           rfs.used = rfs.used, # new object with inds of rfs which contributed to scores
                                           lib.data.processed = lib.data.processed,
                                         # # Spectral data:
                                         #   xmat, 
                                           ppm = ppm)


rf.fits <- metab.evidence$rf.fits   
in.range <- TRUE

bfs <- list(fit.feats = rf.fits$fit.feats[in.range, , drop = F],
            fit.positions = rf.fits$fit.positions[in.range, , drop = F],
            fit.xrow = rf.fits$fit.xrow[in.range],
            pass.fit = T)

# look at one fit
  x <- 1 # there is only 1 fit now
  mi <- metab.evidence$match.info.ss[x, ]
  bfs$fit.feats <- bfs$fit.feats[x, ,drop = FALSE]
  bfs$fit.positions <- bfs$fit.positions[x, ,drop = FALSE]
  bfs$fit.xrow <- bfs$fit.xrow[x]

  simplePlot(rbind(bfs$fit.feats, 
                   xmat[bfs$fit.xrow, bfs$fit.positions]))

plt.pars <- list(vshift = 1, 
                 pixels = c(512, 512), # inc.res
                 pointsize = 0, 
                 interpolate = T, 
                 # exp.by = 0.05,
                 xlim = c(2.1,2.2))

# feature_est_plot(reg = plt.pars$xlim,
#                  metab.evidence,
#                  features.c,
#                  ppm,
#                  plt.pars)

# devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
 fastStack.withFeatures(xmat, ppm, raster = T, bfs = bfs, plt.pars, res.ratio = .1)
 
