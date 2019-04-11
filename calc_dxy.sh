#!/bin/sh

PATH=$PATH:/home/bin/jrick/genomics_general

## script for calculating pairwise Dxy
input=$1 # vcf file to analyze
output=`echo $input | sed 's/\.vcf//'`

sim=$2
tree_height=$3
int=$4
taxa_ref=$5

python parseVCF.py -i $input --skipIndels | gzip > ${output}.geno.gz

python distMat.py -g ${output}.geno.gz -f phased --windType cat -o ${output}.dist

R CMD BATCH dxy_rscript.R
