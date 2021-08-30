#!/bin/sh

source /project/phylogenref/scripts/refbias_config_emp.txt

#module load python/3.6.3
#module load py-numpy/1.14.3-py36

## script for calculating pairwise Dxy
input=$1 # vcf file to analyze
output=`echo $input | sed 's/\.vcf//'`

sim=$2
tree_height=$3
int=$4
taxa_ref=$5
day=$6

python3 /home/jrick/bin/genomics_general/VCF_processing/parseVCF.py -i ${input} --skipIndels | gzip > ${output}.geno.gz

python3 /home/jrick/bin/genomics_general/distMat.py -g ${output}.geno.gz -f phased --windType cat -o ${output}.dist

avg_dist=`grep $taxa_ref ${output}.dist | cut -f 2- -d' ' | sed 's/ /\n/g' | awk '{s+=$1} END {print s/NR}'`

echo "$sim $tree_height $int $taxa_ref $avg_dist" >> ${output_dir}/${day}-emp-refdist.txt

