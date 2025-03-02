# Setup Dirs

parent <- "/Users/mjudge/Downloads"
dirs <- c("results", "params", "matrices", "libraries")
args <- paste("-p", file.path(parent, "safer", dirs), collapse = " ")
system2("mkdir", args)

tmpdir <- file.path(parent, "safer")
params.file <- file.path(tmpdir,"results",'params.yaml')
pars <- pars <- yaml::yaml.load_file(params.file, eval.expr = TRUE)


