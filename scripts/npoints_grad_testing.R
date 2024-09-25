devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')#1697793655

data.dir <- '/Users/mjudge/Documents/ftp_ebi/pipeline_runs_new/'

    unzip_studies(data.dir, exclude = c('1697751863','1697759923'))
    run.idx <- index_studies(data.dir, exclude = c('1697751863','1697759923','1713439191','1726497509', '1726501646'))
    run.idx <- run.idx %>% filter(write_time > '2024-09-23 17:00') %>% arrange(n_points)
    
browse_evidence(results.dir = run.idx$local_path[4], 
                select.compounds = 'Citrate',
                clusterSamples = F,select.samples = 1:10)

run.idx$n_points[4]


browse_evidence(results.dir = run.idx$local_path[5], 
                select.compounds = 'Citrate',
                clusterSamples = F,select.samples = 1:10)

run.idx$n_points[5]


# Issue with spectra being reversed.

xmat <- readRDS('/Users/mjudge/Documents/ftp_ebi/spectral_matrices/MTBLS1_1r_noesypr1d_spectralMatrix.RDS')
ppm <- xmat[1,]
xmat <- xmat[-1,]
npoints <- 64000
cores <- 2

simplePlot(xmat[1,], xvect = ppm)
rs <- resample_spectra(xmat, ppm, npoints = 10000, cores = pars$par$ncores)
simplePlot(rs$spectra[1,], rs$ppm)
