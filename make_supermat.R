#!/bin/env Rscript

options(echo=FALSE)
args = commandArgs(trailingOnly=TRUE)
prefix <- args[1]

library(phylotools)

infiles <- list.files(path=".",pattern="^gene[0-9]+.[A-Z]+.noInv.phy")
supermat(infiles, 
         outfile=paste(prefix,".all.noInv.phy",sep=""),
         partition.file=paste(prefix,".all.partitions.txt",sep=""))
