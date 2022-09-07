# R Package for ICL NMR tool

A package 'wrapper' for the tools and functions written by michael judge and goncalo graca, as part of efforts to make the overall NMR pipeline functional in galaxy.

## Installation
You can install the package using the R CLI straight from this repository.

1 - run `install.package("devtools")` to install the devtools utility library which includes the `install_github` function.

2 - run `library(devtools)` to load the library into the session.

3 - run `install_github("EBI-Metabolights/icl_nmr_R")` to install the package.

4 - run `library(ImperialNMRTool)` to load the package into the session.

To quickly verify if the package is installed you can run `?indsubr` . This will bring up the manual page for this function, and if loaded means that the package itself was installed and loaded successfully. The manual pages for each function contain the information you will need to run them. More detailed instruction may be added to either this repository or the manual pages at a later date.
