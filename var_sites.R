#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
nloci <- args[1] + 1

var_sites <- rnorm(nloci,mean=15,sd=5)
cat(var_sites)
