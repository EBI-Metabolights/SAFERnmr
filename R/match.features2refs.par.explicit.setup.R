#' match.features2refs.par.explicit.setup function
#'
#' Description:
#' Description
#'
#' @param pars a list of input parameters
#'
#' @return None
#'
#' @export
match.features2refs.par.setup <- function(pars) {
    message("-------------------------------------------------------")
    message("-------------------     Matching    -------------------")
    message("-------------------------------------------------------")
    message("\n\n\n")

    ################ Read parameters file ##################

    # library(yaml)
    # pars <- yaml.load_file("./data/params.yaml", eval.expr=TRUE)

    tmpdir <- pars$dirs$temp
    this.run <- paste0(tmpdir)

    ##################################################################################################################
    # Read data and set up ####
    message("### Matching STORM features to GISSMO database")
    message("### MTJ 2023")
    message("")
    message("Loading data from files...\n\n\n")

    fse.result <- readRDS(paste0(this.run, "/fse.result.RDS"))
    clusters <- readRDS(paste0(this.run, "/clusters.RDS"))
    feature <- readRDS(paste0(this.run, "/feature.RDS"))

    ##################################################################################################################

    ## Ref data import ####
    message("Loading and processing reference spectrum data...\n")
    lib.info <- readRDS(paste0(this.run, "/lib.info.RDS"))

    # Import and process the spectra for this dataset ####

    lib.data <- readRDS(paste0(this.run, "/lib.data.RDS"))
    message(" - interpolating ref data to study ppm axis...\n\n")
    lib.data.processed <- prepRefs.for.dataset(lib.data,
        ppm.dataset = fse.result$ppm,
        ref.sig.SD.cutoff = pars$matching$ref.sig.SD.cutoff
    )
    rm(lib.data)
    saveRDS(lib.data.processed, paste0(this.run, "/lib.data.processed.RDS"))

    ##################################################################################################################

    # Setup ####

    # Put features in a matrix
    f.subset <- lapply(unique(clusters$cluster.labs), function(x) {
        (clusters$cluster.labs == x) %>%
            which() %>%
            .[1]
    }) %>%
        unlist()
        .[1:pars$debug$throttle_matches]
    nfeats <- length(f.subset)
    f.stack <- feature$stack[f.subset, , drop = F]
    f.position <- feature$position[f.subset, , drop = F]

    # Put ref spectra in a matrix
    ref.mat <- lapply(
        lib.data.processed,
        function(x) x$mapped$data
    ) %>%
        do.call(rbind, .)

    # Clean up heavy objects
    rm(lib.data.processed)

    # Pre-compute fts for features and refs ####
    message("Pre-computing fts for features and refs (will take a minute)...\n\n\n")
    # Pad size for ref needs to be max.length(features)

    pad.size <- f.stack %>% ncol() - 1

    # Ref mat (padded)
    ref.mat.p <- ref.mat %>% padmat(use = 0, col.by = pad.size)
    ref.mat.p[is.na(ref.mat.p)] <- 0

    # Make padded feature matrix
    # Matrix of zeros of size (refs)

    f.stack <- f.stack - matrixStats::rowMins(f.stack, na.rm = T)
    # f.stack[is.na(f.stack)] <- 0
    f.mat <- matrix(0, nrow = nfeats, ncol = ncol(ref.mat.p))

    # Get the positions of each feature in feature matrix, with rows = feature number
    coords <- which(!is.na(f.stack), arr.ind = T)

    # Linear index the feature vals in
    li <- sub2indR(rows = coords[, "row"], cols = coords[, "col"], m = nfeats)
    f.mat[li] <- f.stack[li]

    # Loop through feature matrix, compute Conj(fftw::fft())
    # f.mat <- apply(f.mat, 1, function(f) Conj(fftw::FFT(f)))
    f.mat <- apply(f.mat, 1, function(f) Conj(fftw::FFT(f))) %>% t()

    # Loop through spec matrix, compute fftw::fft()
    r.mat <- apply(ref.mat.p, 1, function(r) fftw::FFT(r)) %>% t()

    # Transpose matrices (to work with do, dopar) ####
    message("Transposing feature and reference matrices (takes a few seconds)...\n\n")

    r.mat <- t(r.mat)
    ref.mat <- t(ref.mat)

    dir.create(paste0(this.run, "/temp_data_matching"))

    # Ref data applies to all nodes:
    # saveRDS(r.mat, paste0(this.run,'/temp_data_matching/rmat.RDS'))
    # saveRDS(ref.mat, paste0(this.run,'/temp_data_matching/ref.mat.RDS'))

    # Feat data is split between nodes:


    chunk.size <- max(1, nfeats / pars$par$ncores)
    f.grp <- ceiling((1:nfeats) / chunk.size)
    split.scheme <- lapply(unique(f.grp), function(g) {
        list(
            f.inds = which(f.grp == g),
            f.subset = which(f.grp == g) %>% f.subset[.]
        )
    })

    f.mat <- t(f.mat)

    f.mat.split <- lapply(unique(f.grp), function(x) f.mat[, f.grp == x, drop = F])
    # rm(f.mat)

    f.stack <- t(f.stack)
    f.stack.split <- lapply(unique(f.grp), function(x) f.stack[, f.grp == x, drop = F])
    rm(f.stack)

    saveRDS(f.mat.split, paste0(this.run, "/temp_data_matching/f.mat.split.RDS"))
    saveRDS(f.stack.split, paste0(this.run, "/temp_data_matching/f.stack.split.RDS"))

    rm(lib.data)
    rm(ref.mat.p)
    rm(lib.info)
    rm(ref.mat.p)
    rm(fse.result)
    rm(feature)
    rm(f.position)
    rm(clusters)
    rm(coords)

    saveRDS(pad.size, paste0(this.run, "/temp_data_matching/pad.size.RDS"))
    saveRDS(split.scheme, paste0(this.run, "/temp_data_matching/split.scheme.RDS"))

    # x <- lapply(1:length(f.mat.split), function(x){
    #   message('Saving features chunk ', x, '...')
    #   saveRDS(f.mat.split[[x]], paste0(this.run, "/temp_data_matching/f.mat.",x,".RDS"))
    #   saveRDS(f.stack.split[[x]], paste0(this.run, "/temp_data_matching/f.stack.",x,".RDS"))
    # })
}
