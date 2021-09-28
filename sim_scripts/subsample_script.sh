#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --job-name=slurm.refbias.subsample
#SBATCH --time=3-00:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --no-requeue

module load gcc miniconda3 raxml parallel vcftools r

mkdir /lscratch/SLURM_$SLURM_JOB_ID
cd /lscratch/SLURM_$SLURM_JOB_ID
cp -r /project/phylogenref/scripts/sim_scripts/ .

source sim_scripts/refbias_config_subsamp.txt

sim=$1
tree_height=$2
int=$3
taxa_ref=$4
day=$5
vcf=${output_dir}/${day}-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz

source activate new_env

cp $vcf .

echo "beginning parallel jobs per mac"
export sim
export day
export tree_height
export int
export QUAL

parallel --delay 1 --jobs 4  --line-buffer --env QUAL --env sim --env day --env tree_height --env int "bash sim_scripts/each_mac_subsamp.sh {}" ::: $mac_list ::: $miss_list

mkdir s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.subsamp.phylip_tree_files
mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}*.phy s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.subsamp.phylip_tree_files
mv RAxML* s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.subsamp.phylip_tree_files

echo "done with all subsampling processes for $tree_height sim${sim} ${int}; transferring files"

## transfer files to results directory
cd ..  ## to lscratch
tar czf SLURM_${SLURM_JOBID}_results.tgz SLURM_${SLURM_JOBID}/
cp SLURM_${SLURM_JOBID}_results.tgz  $RESULTS_DIR
tar zxf ${RESULTDIR}/SLURM_${SLURM_JOBID}_results.tgz

## clean up
if [ -e "${RESULTDIR}/SLURM_${SLURM_JOBID}_results.tgz" ]
then
   rm -rf SLURM_${SLURM_JOBID}/
   rm  SLURM_${SLURM_JOBID}_results.tgz
fi

echo "done transferring files; exiting now"

date

