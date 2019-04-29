#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
nloci <- args[1]

var_sites <- round(rnorm(nloci,mean=15,sd=5))

if (var_sites < 1) {
	var_sites <- 1
}

cat(var_sites)
