
# Run the pipeline on a specific region

region <- c(2.5007,2.7558)
region <- c(3.8693,3.9265)
devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

parent <- "/Users/mjudge/Downloads"
dirs <- c("results", "params", "matrices", "libraries")
args <- paste("-p", file.path(parent, "safer", dirs), collapse = " ")

system2("mkdir", args)
tmpdir <- file.path(parent, "safer")
params.file <- file.path(tmpdir,"results",'params.yaml')

pars <- pars <- yaml::yaml.load_file(params.file, eval.expr = TRUE)
pars$dirs$temp <- file.path(tmpdir,"results")
pars$files$spectral.matrix <- file.path(file.path(tmpdir, "matrices"),dir(file.path(tmpdir, "matrices")))
pars$files$lib.data <- '/Users/mjudge/Downloads/safer/libraries/data.list_700MHz_mtbls1_annots.RDS'
# pars$files$lib.data <- '/Users/mjudge/Downloads/safer/libraries/data.list_700MHz.RDS'
pars$corrpockets$only.region.between <- c(1,10)
 
# pipeline(params_obj = pars)

pars$dirs$temp <- '/Users/mjudge/Downloads/safer/results/1739774749'
pars$matching$ref.sig.SD.cutoff <- 0.01
match_features2refs_par_setup(pars)
match_features2refs_par_explicit(pars)

    # Weed out empty matches or those which failed
    matches <- readRDS(paste0(tmpdir, "/matches.RDS"))
      
      # if any invalid matches for the feature
        matches <- matches[!is_nullish(matches)]
      # per feature, was NA returned, or was ('matches', 'peak.quality')?
        nomatch <- (lapply(matches, length) %>% unlist) == 1 
        matches <- matches[!nomatch]
      
      # For each feature: 
        # Check if there was an error message
          errors <- matches[names(matches) %in% c('call', 'message')]
        # Check if both of the expected fields are present
          matches <- matches[names(matches) %in% c('matches', 'peak.quality')]

    # Format matches ####
      # At this point, matches have been thoroughly validated and can be rbinded.
        matches.split <- split(matches, names(matches))
        rm(matches)
        
      match.info <- rbindlist(matches.split$matches)

ldp <- readRDS('/Users/mjudge/Downloads/safer/results/1739774749/lib.data.processed.RDS')
match.info <- match.info[order(match.info$rval)]
refs <- ldp[match.info$ref]
lapply(refs, function(x) x$compound.name) %>% unlist


