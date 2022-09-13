
pkgs <- c("batch", "stringr", "speaq", "devtools", "pbapply", "yaml")
for (pkg in pkgs) {
    suppressPackageStartupMessages(stopifnot(library(pkg, quietly=TRUE,logical.return=TRUE, character.only=TRUE)))
}


list_args <- parseCommandArgs(evaluate = FALSE)

# reasoning that we don't need to perform a save operation if we are
# just going to
# load the files on the very next function call.
list[target, picked_peaks] <- stat_decomp()

# do the annotation matching call here