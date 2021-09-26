#!/bin/R

library(ape)
library(phangorn)

args <- commandArgs(trailingOnly=TRUE)

trees <- read.tree(args[1])
rf <- RF.dist(trees,normalize=TRUE,rooted=FALSE)
cat(paste0("sim",args[2]," ",round(mean(rf)[1],5),"\n"))
