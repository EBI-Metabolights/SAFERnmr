# Setup Dirs

# devtools::document('/Users/mjudge/Documents/GitHub/SAFERnmr')

setup_local <- function(parent = "/Users/mjudge/Downloads",
                        safer.dir = "safer", study='mtbls1'){
  
  dirs <- c("results", "params", "matrices", "libraries")
  args <- paste("-p", file.path(parent, safer.dir, dirs), collapse = " ")
  system2("mkdir", args)
  
  tmpdir <- file.path(parent, safer.dir)
  params.file <- file.path(tmpdir,"results",paste0('params_',study,'.yaml'))
  pars <- pars <- yaml::yaml.load_file(params.file, eval.expr = TRUE)
  return(pars)
}
