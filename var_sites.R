#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
nloci <- args[1]
SNPs <- args[2]

var_sites <- round(rnorm(nloci,mean=150,sd=50))

var_sites[var_sites < 1] <- 1

cat(var_sites)
