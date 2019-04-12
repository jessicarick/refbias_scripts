#!/bin/sh

source /project/phylogenref/scripts/refbias_config.txt
PATH=$PATH:/home/jrick/bin/genomics_general:/home/jrick/bin/genomics_general/VCF_processing/

## script for calculating pairwise Dxy
input=$1 # vcf file to analyze
output=`echo $input | sed 's/\.vcf//'`

sim=$2
tree_height=$3
int=$4
taxa_ref=$5

parseVCF.py -i ${input} --skipIndels | gzip > ${output}.geno.gz

distMat.py -g ${output}.geno.gz -f phased --windType cat -o ${output}.dist

avg_dist=`grep 'sim_${taxa_ref}' ${output}.dist | cut -f 2- -d' ' | awk '{s+=$1}END{print s/NR}' RS=" "`

echo "$sim $tree_height $int $taxa_ref $avg_dist" >> ${output_dir}/$date-refdist.txt