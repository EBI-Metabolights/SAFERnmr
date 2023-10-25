#' Remove singlet matches
#'
#' Opting for a much lighter, much faster peak characterization approach. 
#' Characterize peak maxima locations first for features and refs. Use 
#' compressed matrix storage for peak models. 
#' 
#' For a spectral signature (i.e. feature, or ref signature):
#' Determine matched points: those which were both 
#' - not NA in the signature
#' - contained within the matched region
#' 
#' How many of the true peak maxima in the signature exist within the matched points?
#'
#' In other words, how many resonances do each of the fits contain?
#'
#' @param pars A list of input parameters.
#' @return A list of filtered and processed matched peak information, including back-fits to the original spectra.
#' @import yaml
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist
#' 
#' @import pbapply
#' 
#' @export
filter_matches_singlets2 <- function(match.info, feature.c, refmat.c, ncores = 1){
  
  # Pre-compute all feature peaks ####
  
    feats <- feature.c$stack %>% cstack_expandRows
    
    feats <- lapply(1:nrow(feats), function(x) feats[x, ])
    
    f.models <- lapply(feats, function(f) which(!is.na(f)))
    
    feat.pks <- mclapply(feats, function(f){
      pk_maxs(f, mask = !is.na(f))
    }, mc.cores = ncores)
    
    
  # Pre-compute all ref peaks ####
  
    refs <- refmat.c %>% cstack_expandRows
    
    refs <- lapply(1:nrow(refs), function(x) refs[x, ])

    r.models <- lapply(refs, function(r) which(!is.na(r)))

    ref.pks <- mclapply(refs, function(r){
      pk_maxs(r, mask = !is.na(r))
    }, mc.cores = ncores)
    

  # Split up match.info into a list ####
  
    mi.list <- lapply(1:nrow(match.info), function(x) match.info[x,])

  # Loop over each mi.list element (match.info row) and find how many peaks were matched
  
    npeaks <- mclapply(mi.list, function(mi){
    # npeaks <- pblapply(mi.list, function(mi){
      tryCatch({
      # How many f.pks (after filtering for matched region and na-gaps) were included?
  
        f.pks <- feat.pks[[mi$feat]]
        f.model <- f.models[[mi$feat]]
        
        f.pts.matched <- f.model[f.model >= mi$feat.start &
                                 f.model <= mi$feat.end]
        
        npks.feat <- sum(f.pks %in% f.pts.matched)
        if (npks.feat %>% is.null){npks.feat <- 0}
      
      # Same for ref peaks?
      
        r.pks <- ref.pks[[mi$ref]]
        r.model <- r.models[[mi$ref]]
        
        r.pts.matched <- r.model[r.model >= mi$ref.start &
                                 r.model <= mi$ref.end]
        
        npks.ref <- sum(r.pks %in% r.pts.matched)
        if (npks.ref %>% is.null){npks.ref <- 0}

      # Put in df
      
        data.frame(feat = npks.feat,
                   ref = npks.ref)
        
      }, error = function(cond){
        data.frame(feat = 0,
                   ref = 0)
      }
      )
        
    }, mc.cores = ncores
    ) %>% rbindlist
    
  match.info$numpeaks.feat <- npeaks$feat
  match.info$numpeaks.ref <- npeaks$ref
  
  filt <- match.info$numpeaks.feat > 1 & match.info$numpeaks.ref > 1
    
  match.info <- match.info[filt, ]
  
  return(match.info)
}
