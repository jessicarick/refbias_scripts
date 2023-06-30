#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --time=2-00:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --array=1-10
# #SBATCH --array=4
#SBATCH --mem=124G
#SBATCH --no-requeue

module load gcc
module load miniconda3
module load perl
module load vcftools
module load bwa
module load samtools
module load htslib
module load bcftools
module load openmpi
module load raxml
module load r
module load parallel
module load jdk

source activate new_env
PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source

#####################################
####################################
#####################################

source /project/phylogenref/scripts/sim_scripts/refbias_config_emp.txt
tree_height=$1
int=EXT
day=050823
sim=${SLURM_ARRAY_TASK_ID}

cd /gscratch/jrick/phylogenref/emp_tmp/${tree_height}/
mkdir /lscratch/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}/


if [ "$tree_height" == "lates" ]; then
	if [ "$int" == "EXT" ]; then
		taxa_ref="Lcal"
		ref_genome="${REF_PATH}/references/lcalcarifer_genome_v3_chromosomal.fa"
	else
		taxa_ref="Lmar"
		ref_genome="${REF_PATH}/references/lmariae_genome_Feb2018.fa"
	fi
elif [ "$tree_height" == "cichlids" ]; then
	if [ "$int" == "EXT" ]; then
		taxa_ref="Onil"
		ref_genome="${REF_PATH}/references/oreochromis.fasta"
	else
		taxa_ref="Pnye"
		ref_genome="${REF_PATH}/references/puncross.gapsEstimated.fasta"
	fi
else
	echo "reference genome unknown; exiting now"
	exit 1
fi
    
echo "starting empirical $tree_height analysis for $int, simulation $sim"
#declare fastq_list=`ls /project/phylogenref/data/${tree_height}/${int}_ref/*/*.fastq.gz | xargs -n 1 basename | sed 's/\.fastq\.gz//'`

##############################################
### Align to reference  ######################
### Make sure everything in indexed! #########
##############################################
if false; then #skip vcf creation
#export ref_genome

#parallel --env ref_genome "echo 'Fastq: {} \n' && echo 'Mapping reads for {} \n' && bwa mem -t 16 ${ref_genome} /project/phylogenref/data/lates/{}.fastq.gz > aln_{}.sam && echo 'Converting sam to bam for {}' && samtools view -b -S -o aln_{}.bam aln_{}.sam && echo 'Sorting and indexing bam files for {} \n' && samtools sort aln_{}.bam -o aln_{}.sorted.bam && samtools index -c aln_{}.sorted.bam" ::: $fastq_list

#rm -f *.sam
#rm -f *[0-9].bam

#fi #DEBUGGING
#sims=1
#for sim in `seq $sims`; do
	
	if [[ ! -s rand_ind_sim${sim}.txt ]]; then
		if [ "$tree_height" == "lates" ]; then
			ls /project/phylogenref/data/${tree_height}/${int}_ref/*.sorted.bam | grep -v 'L' | shuf -n 49 > rand_ind_sim${sim}.txt
			ls /project/phylogenref/data/${tree_height}/${int}_ref/aln_*.L*.sorted.bam >> rand_ind_sim${sim}.txt
		elif [ "$tree_height" == "cichlids" ]; then
			cat /project/phylogenref/data/${tree_height}/tropheine_bamfiles_highpct_100kReads | grep $int | grep -v 'SRR' | shuf -n 49 > rand_ind_sim${sim}.txt
			ls /project/phylogenref/data/${tree_height}/${int}_ref/bwa_assem*/aln_SRR*.sorted.bam >> rand_ind_sim${sim}.txt
		fi
	elif [ "$tree_height" == "lates" ]; then
                sed -i "s/INT/${int}/g" rand_ind_sim${sim}.txt
                sed -i "s/EXT/${int}/g" rand_ind_sim${sim}.txt
        elif [ "$tree_height" == "cichlids" ] && [ "$int" == "EXT" ]; then
                sed -i 's#INT_ref/bwa_assem_pundamilia#EXT_ref/bwa_assem_oreochromis#g' rand_ind_sim${sim}.txt
        elif [ "$tree_height" == "cichlids" ] && [ "$int" == "INT" ]; then
		sed -i 's#EXT_ref/bwa_assem_oreochromis#INT_ref/bwa_assem_pundamilia#g' rand_ind_sim${sim}.txt
	fi

	QUAL=40	
	if [[ ! -s OUTFILE_q${QUAL}_s${sim}.${int}_RN.bcf ]] && [[ ! -s OUTFILE.q${QUAL}_s${sim}.${int}.vcf ]]; then
		samtools mpileup -g -t DP,AD \
			--skip-indels \
    	 	-P ILLUMINA \
        	-q $QUAL \
    	 	-f ${ref_genome} \
   	    	-b rand_ind_sim${sim}.txt \
     	  	-o OUTFILE_q${QUAL}_s${sim}.${int}.bcf 

      		if [[ -s OUTFILE_q${QUAL}_s${sim}.${int}.bcf ]]
        		then 
          			echo "OUTFILE_q${QUAL}_s${sim}.${int}.bcf not empty; moving on"
        		else
          			echo "OUTFILE_q${QUAL}_s${sim}.${int}.bcf empty; something went wrong"
          			exit 1
   			fi

		if [ "$tree_height" == "lates" ]; then
			cat rand_ind_sim${sim}.txt | sed 's#/project/phylogenref/data/${tree_height}/${int}_ref/##g' > names     
		elif [ "$tree_height" == "cichlids" ]; then
			if [ "$int" == "INT" ]; then
				cat rand_ind_sim${sim}.txt | sed 's#/project/phylogenref/data/${tree_height}/${int}_ref/bwa_assem_pundamilia/##g' > names
			elif [ "$int" == "EXT" ]; then
				cat rand_ind_sim${sim}.txt | sed 's#/project/phylogenref/data/${tree_height}/${int}_ref/bwa_assem_oreochromis/##g' > names
			fi
		fi

		bcftools reheader -s names OUTFILE_q${QUAL}_s${sim}.${int}.bcf > OUTFILE_q${QUAL}_s${sim}.${int}_RN.bcf
	fi


##############################################
#### CREATING RAW VARIANTS FILE FROM BAMs ####
#### DO ONCE FOR EACH MAPPING QUALITY ########
#### FOR EACH GENOME #########################
##############################################
#if false; then
	if [[ ! -s OUTFILE.q${QUAL}_s${sim}.${int}.vcf ]]; then
		bcftools call -m \
                --variants-only \
                --format-fields GQ \
                --skip-variants indels \
                OUTFILE_q${QUAL}_s${sim}.${int}_RN.bcf | bcftools filter \
                        --set-GTs . \
                        --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                                --min-alleles 2 \
                                --max-alleles 2 \
                                --types snps \
                                --apply-filter "PASS" \
                                --output-type v \
                                --output-file OUTFILE.q${QUAL}_s${sim}.${int}.vcf
	fi

#fi                
#                if [[ "$QUAL" -eq 40 ]]; then
                    #echo "Calculating Dxy from calc_dxy_emp script."
                    #${REF_PATH}/sim_scripts/calc_dxy_emp.sh OUTFILE.q${QUAL}_s${sim}.${int}.vcf $sim $tree_height $int $taxa_ref $day
                    #echo "done calculating Dxy"
#                else
#		    		echo "QUAL is $QUAL; not calculating Dxy this time"
#                fi

## delete bcf file if vcf file is successfully created

	if [[ -f "OUTFILE.q${QUAL}_s${sim}.${int}.vcf" ]]; then
		echo "VCF file exists and is not empty; deleting BCF file now"
		rm -f OUTFILE_q${QUAL}_s${sim}.${int}_RN.bcf
		rm -f OUTFILE_q${QUAL}_s${sim}.${int}.bcf
	else
		"VCF file is empty; exiting now"
		exit 1
	fi
          
fi # skip vcf creation
  
#################################
#### FILTERING WITH VCFTOOLS ####
### AND INFERRING PHYLOGENIES ###
#################################

	echo "beginning parallel jobs per maf"
	export QUAL
	export sim
	export REF_PATH
	export output_dir
	export day
	export tree_height
	export int
	export miss_list
	export genes

	parallel --delay 1 --jobs 2  --line-buffer --env genes --env sim --env QUAL --env miss_list --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/sim_scripts/each_maf_emp_jun2022.sh {}" ::: $maf_list ::: $miss_list

#    mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#   	mv *OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#   	mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files

				
####Create a batch file with all normal trees in order by file type, and create a name file
#cat s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-emp-batch.trees
#ls s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-emp-tree.names

#rm -f *sim${sim}.${int}*

cd /gscratch/jrick/phylogenref/slurm_results/
tar -zcvf ${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}-${tree_height}-s${sim}-${int}.tar.gz /lscratch/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}/
rm -rf /lscratch/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}/

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
