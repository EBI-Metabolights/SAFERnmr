
# Function to read a matrix from a .mat file
readMatFile <- function(file_path, variable_name) {
  
  # Load the .mat file
  mat_data <- R.matlab::readMat(file_path)

  # Extract the matrix using the specified variable name
  matrix_data <- mat_data[[variable_name]]

  return(matrix_data)
}

# Read data
  data <- readMatFile('/Users/mjudge/Downloads/Data.mat', 'X.data')
  ppm <- readMatFile('/Users/mjudge/Downloads/ppm_data.mat', 'ppm.data')

  data <- rbind(ppm, data)
  rownames(data) <- c('ppm', 1:(nrow(data)-1))
  
  saveRDS(data, '/Users/mjudge/Downloads/spectral.matrix.RDS')
  
# Read ref
  ref <- readMatFile('/Users/mjudge/Downloads/Fraction.mat', 'X.fraction')
  refppm <- readMatFile('/Users/mjudge/Downloads/ppm_fraction.mat', 'ppm.fraction')

  ref <- rbind(refppm, ref)
  
  ref.lib <- lapply(1:nrow(ref), function(x)
    {
      list(
        tag = as.character(x),
        ref.name = as.character(x),
        compound.name = paste0('fraction - ', x),
        ppm = refppm,
        data = ref[x, ],
        info.row = data.frame(
          Compound.Name = NA,
          Entry.ID = NA,
          Simulation.Name = NA,
          Field.Strength = NA,
          Data = NA,
          selector = NA,
          url = NA,
          localpath = NA,
          file = NA,
          id = NA,
          tag = NA,
          local = NA,
          localFile = NA,
          spectra.json = NA)
        )
    }
  )
  
  # lib.data <- readRDS('/Users/mjudge/Documents/ftp_ebi/gissmo/700MHz_tiny.RDS')

# Save in SAFER format
  saveRDS(ref, '/Users/mjudge/Downloads/ce1/fraction.matrix.RDS')
  saveRDS(ref.lib, '/Users/mjudge/Downloads/fractions.lib.data.RDS')
  
## 

  devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')
  pipeline('/Users/mjudge/Downloads/ce1/CE1.local.usefracs.params.yaml')
  pipeline('/Users/mjudge/Downloads/ce1/CE1.local.fraction.opt.params.yaml')
  pipeline('/Users/mjudge/Downloads/ce1/CE1.local.dataset.opt.params.yaml')
  
  
  