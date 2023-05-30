#' Filter features based on specified criteria
#' - null features
#' - features outside of ppm range
#' - features with no runs >= min.runlength
#' - features derived from < min.subset spectra
#' - features with strong baseline effect 
#' 
#' @param feature A feature object to be filtered
#' @param ppm A vector of ppm values for each point in the spectra
#' @param ppm.range A vector specifying the range of ppm values to include
#' @param min.runlength The minimum length of a run of consecutive data points to include
#' @param min.subset The minimum number of spectra a feature must be present in to be included
#' @param prom.ratio The maximum ratio of peak prominence to intensity range for a feature to be considered monotonic
#' @param give A character string indicating whether to return the filtered feature object or the filter as a logical vector
#' 
#' @return If \code{give = "features"}, a feature object with the filtered features. If \code{give = "filter"}, a logical vector indicating which features passed the filter.
#' 
#' @importFrom magrittr %>%
#'
#' 
#' @export
filterFeatures <- function(feature, ppm, ppm.range, min.runlength = 3, min.subset = 5, prom.ratio = 0.3, give = "filter", max.features = 10000){
  
  
  # No empty features, please
    message("Filtering out null features")
    nullfeatures <- apply(feature$position, 1, function(x) all(is.na(x)))
    
  # Filter out features with NULL driver
    message("Filtering out features with NULL driver...")
    null.driver <- is.na(feature$driver.relative)
    
  # In the selected ppm range?
    message("Filtering out features outside of (", ppm.range[1], ") - (", ppm.range[2], ") ppm...")
    ppmrngs <- apply(feature$position, 1, function(x) range(x, na.rm = TRUE) %>% ppm[.] %>% rev) %>% t
    
    # plot(ppmrngs[TRUE])
    inbounds <- apply(ppmrngs, 1, function (x) (x[1] > ppm.range[1] & x[2] < ppm.range[2]) &
                                            (x[2] > ppm.range[1] & x[2] < ppm.range[2])  ) %>% 
                                                        unlist
      
  # Filter out features with runs too short
    message("Filtering out features with no runs >= ", min.runlength, " points...")
    rl.pass <- apply(feature$stack, 1,
                   function(x) {run.lens <- x %>% is.na %>% "!"(.) %>% runs.labelBy.lengths
                                return((run.lens >= min.runlength) %>% any)})
    
  # Filter out features without enough subset size 
    message("Filtering out features derived from < ", min.subset, " spectra...")
    ss.pass <- feature$subset$sizes >= min.subset 
  
  # Filter out monotonic features
    message("Filtering out features with strong baseline effect (max true peak prominence < ", prom.ratio, "*intensity range")
    message("or strong correlation (r >= .99) between valley and peak intensities). Progress...")
    
    # Are there any peaks > 0.3 * feature intensity range? (Monotonic?) ####
      bl.effect <- pblapply(1:nrow(feature$stack), function(f)
      {
         # print(f)
         # f <- 98
         feature$stack[f, , drop = FALSE] %>% 
           trim.sides %>%
           detect.baseline.effect 
      })
      
      bl.effect.ul <- bl.effect %>% unlist %>% as.logical
      pass.prom <- bl.effect.ul[c(T,F,F)]
      pass.fit <- bl.effect.ul[c(F,T,F)]
      not.singlet <- bl.effect.ul[c(F,F,T)]
      
      not.monotonic <- pass.prom & pass.fit & not.singlet
      
      # filt <- !nullfeatures & inbounds & rl.pass & ss.pass
      # feature$stack[which(!not.monotonic),] %>% apply(1,scale.between) %>% t %>% trim.sides %>% simplePlot
      # which(
      #   which(
      #      not.monotonic) %in% c(1068,1102,1116,1130,1136,1144,1151,1163,1166,1172,1173,1192,1225,1227)
      #   )
      # )
      # feature$stack[c(1068,1102,1116,1130,1136,1144,1151,1163,1166,1172,1173,1192,1225,1227),] %>% apply(1,scale.between) %>% t %>% trim.sides %>% simplePlot
      # sum(!not.monotonic)
      
  # Build the filter set
    
    all.filts = list(not.null = !nullfeatures,
                     not.null.driver = !null.driver,
                     inbounds.ppm = inbounds,
                     rl.pass = rl.pass,
                     ss.pass = ss.pass, 
                     pass.prom = pass.prom,
                     not.singlet = not.singlet,
                     not.monotonic = pass.fit)
    
    filt <- !nullfeatures & !null.driver & inbounds & rl.pass & ss.pass & not.monotonic 
    message('Filtering complete. ', sum(filt), '/', length(filt), ' features passed filters.')
    
    all.filts$rand.subset = rep(T, length(filt))
    
  # After filtering, if number of features exceeds max.feature, limit them randomly:
      
      passed <- which(filt)
      
      if (length(passed) > max.features){

        subset <- runif(n = max.features, 
                        min = 0, 
                        max = length(passed)) %>% ceiling %>% unique %>% passed[.]
        
        filt[-subset] <- F
        all.filts$rand.subset[-subset] <- F
        
        number.excluded <- length(passed) - sum(filt)
        message('Excluding ', number.excluded, ' / ', length(passed), ' passing features (', (number.excluded/length(passed)*100) %>% round,' %) to satisfy limit of ', pars$debug$throttle_features, ' features ...')
      }
  
      
  if (give == "features"){
    # Go through feature object and apply filter
      feature$stack <- feature$stack[filt, ]
      feature$position <- feature$position[filt, ]
      feature$subset$ss.all <- feature$subset$ss.all[filt, ]
    return(feature)
  }
  if (give == "filter"){
    return(
      list(combined = filt,
                all = all.filts)
    )
  }

}