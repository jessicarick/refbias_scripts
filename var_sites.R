#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
nloci <- args[1]

var_sites <- round(rnorm(nloci,mean=30,sd=5))

var_sites[var_sites < 1] <- 1

cat(var_sites)
