#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --time=8:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=gtrf
module load gcc miniconda3 r parallel raxml

file=$1
height=$2
date=$3
echo "working with $file"

base=`echo $file | sed 's/\.noInv\.phy//g'`
sim_int=`echo $base | cut -f 1 -d'-'`
mkdir $base

nsites=`head -n 1 $file | cut -f 2 -d' '`
seq 1 200 $nsites > $base/start
seq 200 200 $nsites > $base/end
echo $nsites >> $base/end
paste -d'-' $base/start $base/end | awk '{$1="DNA, part" FSi++ FS " = " FS $1;}1' > ${base}/partitions

python /project/phylogenref/scripts/sim_scripts/unconcatenate_phylip.py -P ${base}/genetrees_ $file ${base}/partitions

cd $base
ls *.phylip | parallel -j 8 "raxmlHPC-PTHREADS-AVX -s {} -n {.} -m ASC_GTRCAT -f a -x 12345 -N 100 --asc-corr=lewis -p 12345 -T 2 -V"
cat RAxML_*bestTree* > RAxML_all_genetrees.tre

Rscript /project/phylogenref/scripts/sim_scripts/calc_genetree_rf.R RAxML_all_genetrees.tre $sim_int >> /project/phylogenref/scripts/output/new/${date}_${height}_gt_rf.txt

echo "all done"
date
