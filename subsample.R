#!/bin/env Rscript

options(echo=FALSE)
library(phylotools)
library(stringr)

# read in arguments from command line
args = commandArgs(trailingOnly=TRUE)
phy <- paste0(args[1],".phy")
nsnp <- args[2]
max <- args[3]

# calculate number to remove
diff <- nsnp-max

# randomly choose which to remove
omit <- sample(1:nsnp, diff, replace=F)

# read in phylip and remove randomly chosen snps
alignment <- read.phylip(phy,clean_name=FALSE)
align.tab <- str_split(alignment$seq.text,"",simplify=TRUE)
align.subsamp <- align.tab[,-omit]
phy.subsamp <- apply(align.subsamp, 1, paste0, collapse="")

# write subsampled phylip
final.phy <- data.frame(seq.text=phy.subsamp)
row.names(final.phy) <- alignment$seq.name
final.snp <- ReadSNP(final.snp)
WriteSNP(final.snp, file=paste0(args[1],".subsamp.phy"),format='phylip')
