#!/bin/sh

PATH=$PATH:/home/bin/jrick/genomics_general

input=$1 # vcf file to analyze
output=`echo $input | sed 's/\.vcf//'`

python parseVCF.py -i $input --skipIndels | gzip > ${output}.geno.gz

python distMat.py -g ${output}.geno.gz -f phased --windType cat -o ${output}.dist

dxy_rscript 
