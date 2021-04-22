#!/usr/bin/env Rscript

# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

# install to local user library path
install.packages("methods",    lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("optparse",   lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("stringr",    lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("R.matlab",   lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("Matrix",     lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("matrixcalc", lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("expm",       lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("numDeriv",   lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("plyr",       lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
install.packages("pracma",     lib = Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")

library(stringr)


.libPaths(rev(dir('~/R/x86_64-pc-linux-gnu-library', full.names=TRUE)))
cat( 'R package path: ' )
cat( .libPaths(), '\n', sep=':')


thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}


install.packages(str_c( dirname(thisFile()), '/ReDCM-0.2.1.tar.gz' ), repos = NULL, type="source")
