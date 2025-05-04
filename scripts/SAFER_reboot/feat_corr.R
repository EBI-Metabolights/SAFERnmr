
# Correlate features across samples # TRUE STOCSY
# - Relatively quantify the features using sfe fits (spec or fit AUC?)
# - use complete pairwise correlations to build network of features
# - there will be a lot of true or partial (adjacent) duplicate features
#   - these will cover the same spectral region of the spectra
#   - also the same regions in the ref
#   - there needs to be a good way of estimating the relative quantification for 
#     each combination of features
#   - imagine the reference spectrum getting annotated by features (or ref features)
#     - each spectrum-reference pair accumulates quantifications from 
#       matched features
#     - 
# - remember: you could be looking for two separate things:
#   1) which features should go together? cluster based on STOCSY
#     - then, is the signature for a feature cluster represented in the ref?
#     - in this case, feature fits are clustered by STOCSY then compiled onto ref spectrum
#       
#   2)  how well does the overall signature fit?
#     - i.e. when looking at (spec-)feature matches to ref, where each gets a quant
#     vote, how well do the overall signatures match?
#     - in this case, backfits are compiled onto ref spectrum, then assessed for
#       STOCSY 

# Simplified Gradient Tests for parameter sensitivity (19-20OCT) ####
    
    # functionalized indexing ####
    
      devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
      unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
      run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923'))
      
        # run.idx$write_time %>% as.POSIXct() %>% as.numeric %>% sort %>% .[16] # %>% plot
        # run.idx <- run.idx[(run.idx$write_time %>% as.POSIXct() %>% as.numeric) <= 1698006395, ]
        # run.idx$run_id %>% cat
        
      # run.idx <- run.idx[run.idx$run_id %in% c('1697740712','1697740820','1697740956','1697741076','1697747670','1697751222','1697751875','1697753401','1697793288','1697793355','1697793655','1697793709','1697823986','1697836051','1697842288','1697865804'),]
      run.idx$start <- run.idx$run_id %>% as.POSIXct(origin = "1970-01-01")
      run.idx <- run.idx[run.idx$start >= "2023-10-26 BST", ]
      run.idx <- run.idx[run.idx$max_backfits == 5E7,]

  # Load Study Data: ####
      run <- run.idx$local_path[1]
      message('Reading ', run, '...')
      fse.result <- readRDS(paste0(run, '/fse.result.RDS'))
        xmat <- fse.result$xmat
        ppm <- fse.result$ppm
        
      feature <- readRDS(paste0(run, '/feature.final.RDS')) %>% expand_features()

  # RQ Features: ####
  
    rq.big <- 
      mclapply(1:length(feature$sfe), function(f){
        
        # f <- 100

        # Make the information for reconstructing fits
          
          f.profile <- feature$stack[f,] %>% trim_sides
          f.pos <- feature$position[f,] %>% trim_sides
          fit.info <- feature$sfe[[f]]$fits %>% do.call(rbind,.) %>% cbind(ss = feature$sfe[[f]]$feat$ss)
          
          rqs <- lapply(1:nrow(fit.info), function(x){
            # x <- 1
            fi <- fit.info[x, ]
            lag <- feature$sfe[[f]]$lags[x]
            f.pos.ss <- (f.pos%>%complete_indsVect) + lag
            spec <- xmat[fi$ss, f.pos.ss]
              # simplePlot(spec, xvect = ppm[f.pos.ss])
              
            f.fit <- list(feat.fit = f.profile*fi$ratio + fi$intercept,
                          spec.fit = spec)
            
            # plot_fit(fit = f.fit, ppm = ppm[f.pos.ss], type = 'auc')
            # x <- x + 1
            return(data.frame(feat = f,
                              ss = fi$ss,
                              auc.feat = f.fit$feat.fit %>% sum(na.rm = TRUE),
                              auc.spec = f.fit$spec.fit %>% sum(na.rm = TRUE))
                   )
          }) %>% do.call(rbind,.)
          
      }, mc.cores = 5) %>% do.call(rbind,.)
      
      # scattermore::scattermoreplot(x = rq.big$auc.feat, xlab = 'auc.feat',
      #                              y = rq.big$auc.spec, ylab = 'auc.spec')
      
    # Build the quant matrix
      
      mat <- matrix(NA, nrow = nrow(xmat), ncol = length(feature$sfe))
      
      linds <- sub2indR(rows = rq.big$ss, cols = rq.big$feat, m = nrow(xmat))
      mat[linds] <- rq.big$auc.feat
      
      # mat[is.na(mat)] <- 0; heatmap(mat, scale = 'col')
      
    r <- cor(mat, use = "pairwise.complete.obs")
    
    #  r[is.na(r)] <- 0; heatmap(r, scale = 'none')

    r[is.na(r)] <- 0
    r.thresh <- 0.75
    scattermore::scattermoreplot(x = 1:length(r), y = sort(r))
    abline(h = r.thresh, col = 'red')
    r[r<r.thresh] <- 0
    
    # Try clustering: ####
    
          t1<-Sys.time()
          message('Doing OPTICS clustering on features using correlations...')
            max.eps = 50
            minPts = 2
            eps.stepsize = .01
            opres <- dbscan::optics(1-r, eps = 5, minPts = 2) # kdtree only available for euclidean dist
          print(Sys.time()-t1)
          
          message('Evaluating cluster number at a range of eps_cl values...')
            eps.vals <- seq(0, max.eps/minPts, eps.stepsize)
            num.clusters <- pblapply(eps.vals, function(x) {
              dbres <- dbscan::extractDBSCAN(opres, eps_cl = x)
              run.labels(dbres$cluster) %>% max
              # run.labels(zeroes) + 
            }) %>% unlist
            
            eps.val <- num.clusters %>% which.max %>% eps.vals[.]
          
          # Plot
          # pdf(file = paste0(plot.loc,'/optics.opt.epsvalue.pdf'),   # The directory you want to save the file in
          #   width = 3, # The width of the plot in inches
          #   height = 3)
          
            plot(eps.vals, num.clusters)
            abline(v = eps.val, col = 'red')
            abline(h = max(num.clusters), col = 'red')
          # dev.off()
          
            
          message('Using maximum cluster number: ', max(num.clusters), ' (eps_cl = ', eps.val, ')')
            dbres <- dbscan::extractDBSCAN(opres, eps_cl = eps.val)
            clusters <- run.labels(dbres$cluster)
            zeroes <- clusters == 0
            
          message('OPTICS Clustering Results: ')
          message(sum(zeroes), ' features (', round(sum(zeroes)/length(zeroes), 2)*100,'%) were excluded by clustering. The other ', sum(!zeroes), ' features (', round(sum(!zeroes)/length(zeroes), 2)*100,'%) fell into ', max(clusters),' shape clusters...')
          
          cluster.sizes <- lapply(unique(clusters[clusters>0]), function(x) which(clusters == x) %>% length) %>% unlist
          
          # pdf(file = paste0(plot.loc,'/optics.clusters.sizeDistribution.pdf'),   # The directory you want to save the file in
          #   width = 3, # The width of the plot in inches
          #   height = 3)

            # hist(cluster.sizes, breaks = 100, xlab = 'Cluster Size', ylab = 'Frequency')
            plot(sort(cluster.sizes), xlab = 'Cluster', ylab = 'Cluster Size', log = 'y')

          # dev.off()    
      
            
    #      
        