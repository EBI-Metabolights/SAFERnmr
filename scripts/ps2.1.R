## Functions for scoring2

  # Setup (from Line 230 in ps2.R):
  # refspec needs to be gotten. 
      # refspec <- lapply(lib.data.processed[ref], function(x)
      #   {
      #     x$mapped$data.compressed %>% expand_stacklist(which.stacks = 'data') %>% .[[1]]
      #   }
      # ) %>% do.call(rbind, .)

  # Select relevant backfits and matches ####

    res <- lapply(1:nrow(comb), function(combn){
      combn <- 1
      rp.rows <- which(ref.pairs$ss.spec == ss.spec)
      rp <- ref.pairs[rp.rows,]
      
      }
      res <- signature_similarity_multifeature(rp, cutoff.r2 = 0.9)
      
      model <- mblm(y ~ x, data = data)
                summary(model)
                plot(data$x, data$y)
                abline(model, col = "orange")
      
      list(stats = res$stats,
           exclude = res$excluded)
    })
    
      
      

signature_similarity_multifeature <- function(rp, sd_tol = 2, cutoff.r2 = 0.9){
  # The goal of this function is to gauge the signature-level similarity between
  # backfit features and reference spectra. 

  # Functions: 

  df <- data.frame(x = rp$pct.ref %>% scale_between(0,1),
                   y = (rp$pct.ref * rp$fit.scale) %>% scale_between(0,1))
  
  if (nrow(df %>% unique) > 2){
    # res <- itfit(df, cutoff.r2 = 0.9)#, show = TRUE)
    model <- MASS::lqs(y ~ x, data = df)
      model$residuals
    
    return(list(excluded = res$excluded,
            stats = res$stats))
  } else {
    
    return(list(excluded = res$excluded,
            stats = data.frame(r2 = 1,
                               corr = 1,
                               slope = NA,
                               intercept = NA)))
  }
  
}

plot_features <- function(rp, backfits, match.info, selectrows = NULL, plt.rng = NULL){
  
    get.intensities <- function(rf.specFits, feat.models, refspec){
    # Unlist all the ref feats into spectrum-fit ref feats, and expand their spectrum positions: ####  
    fit.feats <- lapply(1:nrow(rf.specFits), function(x) {
      # x <- 1
      # Get the feature model
        sf <- rf.specFits[x, ]
        
        # Starting from feature and match info:
          feat.model <- feat.models[which(sf$feat==f.nums),]
          f.matched.pts <- sf$feat.start:sf$feat.end
          feat.gaps <- feat.model[ f.matched.pts ] %>% is.na
    
      # Get the ref segment (rf)
      
        # Starting from refmat and match info and feature model:
        
          # rf <- ld$mapped$data %>% scale_between(0,1) %>% .[ sf$ref.start:sf$ref.end ]
          rf <- refspec %>% .[ sf$ref.start:sf$ref.end ]
        
        # NA-fill the feature gaps
            
          rf[feat.gaps] <- NA
    
      # Calculate the fit feature profile 
          
          feat.fit <- rf * sf$fit.scale + sf$fit.intercept
          
          # plot_fit(list(feat.fit = feat.fit,
          #               spec.fit = xmat[sf$ss.spec, sf$spec.start:sf$spec.end]),
          #               type = 'auc')
          
          return(feat.fit)
    })
    return(fit.feats)
    }
    get.positions <- function(rf.specFits, fit.feats){
    fit.positions <- lapply(1:nrow(rf.specFits), function(x) 
    {
        rff <- rf.specFits[x,]
      # Get the positions in the xrow
        pos <- rff$spec.start:rff$spec.end
        pos[is.na(fit.feats[[x]])] <- NA
        return(pos)
    }) 
    return(fit.positions)         
    }  
    feature.matrices <- function(rf.specFits, fit.feats, fit.positions){
    
    # Make matrices from those lists, add to bestfits list object ####
    
    maxlen.ff <- lapply(fit.feats, length) %>% unlist %>% max
    
    ffint <- lapply(1:length(fit.feats), function(i){
    
              c( 
                 fit.feats[[i]], 
                 rep(NA, maxlen.ff-length( fit.feats[[i]] ) )
                )
    
            }) %>% do.call(rbind,.)
    
    ffpos <- lapply(1:length(fit.positions), function(i){
    
              c( 
                 fit.positions[[i]], 
                 rep(NA, maxlen.ff-length( fit.positions[[i]] ) )
                )
    
            }) %>% do.call(rbind,.)
    
    fit.xrow <- rf.specFits$ss.spec
    
    return(list(fit.feats = ffint,
              fit.positions = ffpos,
              fit.xrow = fit.xrow))
    }

## Plot the features ####
  # Get all backfits relevant to this pair
    if(!is.null(selectrows)){
      rp <- rp[selectrows, ]
    }
  
    rfs.selection <- rp$match %>% unique
  
  # Put their fit data into a df
    rf.specFits <- lapply(rfs.selection, function(x) 
    {
      # x <- rfs.selection[1]
      rf <- backfits[[x]]
      rf <- rf[rf$ss.spec == ss.spec,]
      rf$feat <- match.info$feat[x]
      rf$feat.start <- match.info$feat.start[x]
      rf$feat.end <- match.info$feat.end[x]
      rf$ref.start <- match.info$ref.start[x]
      rf$ref.end <- match.info$ref.end[x]
      rf$id <- match.info$id[x]
      rf
    }) %>% do.call(rbind, .)
    
  # Plot them
    
    f.nums <- rf.specFits$feat %>% unique # note: rp$feat indexes for feat.models in chonk 
    # why do these features not line up with the matched regions
      
    fit.feats <- get.intensities(rf.specFits, feat.models, refspec)
    fit.positions <- get.positions(rf.specFits, fit.feats)
    f.mats <- feature.matrices(rf.specFits, fit.feats, fit.positions)
    
    if (is.null(plt.rng)){
      plt.rng <- apply(f.mats$fit.positions, 1, function(x) range(x, na.rm = TRUE)) %>% range %>% ppm[.]
    }
    
    plt.pars <- list(vshift = .1, 
                 pixels = c(1080, 1080), # inc.res
                 pointsize = 0, 
                 interpolate = T, 
                 # exp.by = 0.05,
                 xlim = plt.rng # this should never arrive here as NULL - screen out beforehand
                 )
    
    fastStack.withFeatures(xmat, ppm, raster = T, bfs = f.mats, plt.pars, 
                       res.ratio = .1)
    

}

lqs_features <- function(sd_tol = 1, show = FALSE)
    
    df <- data.frame(x = rp$pct.ref,# %>% scale_between(0,1),
                     y = (rp$pct.ref * rp$fit.scale))# %>% scale_between(0,1))
    
    model <- MASS::lqs(y ~ x, data = df)
      absr <- abs(model$residuals)
      rsd <- sd(absr)
      outliers <- absr > rsd * sd_tol
      r2 <- 1-(sum(model$residuals^2) / sum((df$y - mean(df$y))^2))
      
      if (show){
        plot(df$x, df$y)
        abline(model, col = "orange")
        # points(x = df$x, y = abs(model$residuals) %>% scale_between, col = 'blue')
        points(x = df$x[outliers], y = df$y[outliers], col = 'blue', pch = 19)
        abline(b = model$coefficients[2], a = model$coefficients[1] + rsd, col = 'blue', lty=2)
        abline(b = model$coefficients[2], a = model$coefficients[1] - rsd, col = 'blue', lty=2)
      }
      
      list(coefficients = model$coefficients,
           outliers = outliers,
           residuals = model$residuals,
           sd.residuals = rsd,
           r2 = r2)

    return(list(stats = res$stats,
                exclude = res$excluded))
}
      
    # model$coefficients
    # model$crit
    # model$fitted.values
    # model$residuals
    # model$model
    
    # model <- mblm(y ~ x, data = df)
    #           summary(model)
    #           plot(df$x, df$y)
    #           abline(model, col = "orange")
    #
    # model$coefficients
    # model$residuals
    # model$fitted.values
    # model$slopes
    # model$intercepts
    # model$model
    # model$effects          
    
    # res <- itfit(df, cutoff.r2 = 0.9)#, show = TRUE)
    
    # rp %>% plot_features(backfits, match.info, plt.rng = c(2.5,2.75))
    # rp[-res$excluded,] %>% plot_features(backfits, match.info, plt.rng = c(2.5,2.75))
      # stdv = sd(model$residuals)
      # cutoff <- stdv*sds
      # within_sds <- abs(model$residuals) <= cutoff
      # 
    # Plot feature quants
      
  

# Method: 
  fits <- lapply(1:nrow(comb), function(i){
    #i <- i+1
    #i <- 55
    ss.spec <- comb$ss.spec[i]
    rp <- ref.pairs[ref.pairs$ss.spec == ss.spec,]
    
    # Get the Pearson's r for the feature quants
      corr <- cor(rp$pct.ref, (rp$pct.ref * rp$fit.scale))
      
    # Fit the unit line to the feature quants
      fit = fit_unit_line(rp, sds = 2, show = TRUE)
    
    # fit = fit_unit_line(rp, sds = 2, show = TRUE) # debugging/examples
    # plot_features(rp, backfits, match.info),
    #               selectrows = res$within_sds) # debugging/examples
    
    return(list(pair = comb[i,],
                fit = res,
                corr = corr))
  })
  cors <- fits %>% lapply(function(x) x$corr) %>% unlist #%>% sort %>% plot
  which(cors<0)
  
  
  fits <- lapply(1:nrow(comb), function(i){
    #i <- i+1
    ss.spec <- comb$ss.spec[i]
    rp <- ref.pairs[ref.pairs$ss.spec == ss.spec,]
    
    # Want to know for each combination of features, how much of the ref is covered
    #   Assort them along the ref
    #   use linear assignment?
    #     issue: how do you assign multiple? no dist mat
    #     for nonoverlapping features
    #       find potential edges
    #     
    #     amount of ref accounted
    #     quality of fit
    #     relative intensities matched
    
    return(list(pair = comb[i,],
                fit = res,
                corr = corr))
  })
  
  
}



  # Calculate the signature-wide correlation for the pair:
  #   Note: calc these all in one, then loop through by pair to corr.
    # feature.in.ref <- rp$pct.ref
    # feature.in.spec <- (rp$pct.ref * rp$fit.scale) #/max(rp$fit.scale)
    
    # df <- data.frame(y = feature.in.spec, x = feature.in.ref)
    
    # plot(df$x, xlab = 'feature in ref', df$y, ylab = 'feature in spec')
    # 
    # model <- lm(I(y - x) ~ 1, data = df)
    #   intercept <- coef(model)
    #   abline(a = intercept, b = 1, col = "blue")
    #   ms <- summary(model)
    # 
    #   abline(model, col="blue")
      
    
    # Pick the points closest to the line
      # for each spectral point, pick the points closest to the line
      
      # dummy.ref <- rep(0, length(refspec))
    
    # Make a matrix housing the positions of all the features for the spec
      
      # mat <- matrix(data = 0, nrow = length(rp.rows), ncol = length(refspec))
      
      # ref.pairs$feat
    
}

  fit_unit_line <- function(rp, sds, show = FALSE){
    
    df <- data.frame(x = rp$pct.ref %>% scale_between(0,1),
                     y = (rp$pct.ref * rp$fit.scale) %>% scale_between(0,1))
    
    # model <- lm(df$y ~ df$x, data = df)
      # abline(model, col="blue")
      # ms <- summary(model)
      stdv <- sd(residuals)
      cutoff <- stdv*sds
      within_sds <- abs(residuals) <= cutoff
        # these are the rows of rp that would give a good whole-signature fit
      
    # Plot feature quants
    
      if(show){
      plot(df$x, xlab = 'feature in ref', df$y, ylab = 'feature in spec')
        abline(a = 0, b = 1, col = "blue")
        residuals <- df$x - df$y
        
        points(df$x[within_sds], df$y[within_sds], pch=19, col="blue")
        abline(a = -cutoff, b = 1, col = "blue", lty = 2)
        abline(a = cutoff, b = 1, col = "blue", lty = 2)
        
      }
  
      return(list(residuals = residuals,
                  sd = stdv,
                  within_sds = within_sds))
}
  
  itfit <- function(df, cutoff.r2 = 0.9, show = FALSE){
      df.copy <- df
      n <- nrow(df)
      df$id = 1:nrow(df)
      r2 <- 0
      
      record <- data.frame(high_leverage_points = NA,
                           r2 = NA,
                           m = NA,
                           b = NA,
                           iteration = NA,
                           corr = NA)
      j <- 1
      while(r2 < cutoff.r2 & 
            nrow(df %>% unique) > 1 & # don't want to end up with a single point (or duplicates thereof)
            j < n){
        
          model <- lm(y ~ x, data = df) # recompute model
          r2 <- 1 - sum((model$residuals)^2)/sum((df$y-mean(df$y))^2)
          leverage <- hatvalues(model)
          threshold <- 2 * (length(coef(model)) / nrow(df))
          hlp <- which(leverage > threshold)
  
          if (any(hlp)){
            
            # plot(df$x, xlab = 'feature in ref', df$y, ylab = 'feature in spec')
            #   abline(a = model$coefficients[1], b = model$coefficients[2], col = "blue")
            #   points(df$x[hlp], df$y[hlp], pch=19, col="blue")
            
            record <- rbind(record, 
                            data.frame(high_leverage_points = df$id[hlp],
                                       r2 = r2,
                                       m = model$coefficients[2] %>% as.numeric,
                                       b = model$coefficients[1] %>% as.numeric,
                                       iteration = j,
                                       corr = cor(df$x, df$y)
                                       )
                            )
            df <- df[-hlp,]
          }
          j <- j + 1
  
      }
      
      record <- record[-1,] # remove NA init row
      
      if (nrow(record)>0){
        iters <- max(record$iteration)-1
        hlps <- record$high_leverage_points[record$iteration <= iters]
      } else {
        iters <- 1
        hlps <- NULL
      }
      
      corr <- cor(df$x[-hlps], df$y[-hlps])
      
      if (show){
        
        plot(df.copy$x, xlab = 'feature in ref', df.copy$y, ylab = 'feature in spec')
          abline(a = model$coefficients[1], b = model$coefficients[2], col = "blue")
          if (nrow(record)>0){
            points(df.copy$x[hlps], df.copy$y[hlps], pch=19, col="blue")
          }
          
          text(x = .6, y=.7, paste0('Slope = ', model$coefficients[2] %>% round(4), 
                                    '\nIntercept = ', model$coefficients[1] %>% round(4),
                                    '\nR^2 = ', r2 %>% round(4),
                                    '\nCorr = ', corr %>% round(4)))
      }
      
      res <- data.frame(r2 = r2,
                 corr = corr,
                 slope = model$coefficients[2],
                 intercept = model$coefficients[1])
      
      return(
              list(model = model,
                    stats = res,
                    record = record,
                    excluded = hlps)
                   )
}
