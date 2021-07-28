#!/bin/env Rscript

#options(echo=FALSE)
library(ape)
library(phylotools)
library(phrynomics)
library(stringr)

# read in arguments from command line
args = commandArgs(trailingOnly=TRUE)
phy <- paste0(args[1],".phy")
nsnp <- as.numeric(args[2])
#max <- args[3]
max <- 5000

print(paste("nsnp is ",nsnp," and max is ",max))

# calculate number to remove
diff <- nsnp - max

# randomly choose which to remove
if (diff > 0) {
	omit <- sample(1:nsnp, diff, replace=F)
	} else {
	omit <- sample(1:nsnp, 1, replace=F)
}

# read in phylip and remove randomly chosen snps
alignment <- read.phylip(phy,clean_name=FALSE)
align.tab <- str_split(alignment$seq.text,"",simplify=TRUE)
align.subsamp <- align.tab[,-omit]
phy.subsamp <- apply(align.subsamp, 1, paste0, collapse="")

# write subsampled phylip
final.phy <- data.frame(seq.text=phy.subsamp)
row.names(final.phy) <- alignment$seq.name
final.snp <- ReadSNP(final.phy)
WriteSNP(final.snp, file=paste0(args[1],".subsamp.phy"),format='phylip')

