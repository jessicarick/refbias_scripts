#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --mem=124G
#SBATCH --array=1-10

module load gcc
module load vcftools
module load python/3.6.3
module load py-numpy/1.14.3-py36

sim=${SLURM_ARRAY_TASK_ID}
base=OUTFILE.q40_s${sim}.INT

echo "working with $base"

vcftools --vcf ${base}.vcf --max-missing 0.5 --recode --out ${base}

python3 /home/jrick/bin/genomics_general/VCF_processing/parseVCF.py -i ${base}.recode.vcf --skipIndels --gtf flag=DP min=10 | gzip > ${base}.recode.geno.gz && python3 /home/jrick/bin/genomics_general/distMat.py -g ${base}.recode.geno.gz -f phased --windType sites -w 1000000 -T 16 -o ${base}.recode.dist

echo "done calculating dxy for $base"
