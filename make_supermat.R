#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]
miss <- args[2]
maf <- args[3]

library(phylotools)

pattern <- paste0("^gene[0-9]+_miss",miss,"_maf",maf,".[A-Z]+.noInv.phy")
infiles <- list.files(path=".",pattern=pattern)
supermat(infiles, 
         outfile=paste(prefix,".all.noInv.phy",sep=""),
         partition.file=paste(prefix,".all.partitions.txt",sep=""))

