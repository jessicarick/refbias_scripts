#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --array=1-2
#SBATCH --mem=124G
#SBATCH --no-requeue

module load gcc
module load miniconda3
module load perl
module load vcftools
module load bwa
module load samtools/1.6
module load htslib
module load bcftools
module load raxml
module load r
module load parallel

source activate new_env
PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source

#####################################
####################################
#####################################

source /project/phylogenref/scripts/refbias_config_emp.txt
tree_height="lates"
int=EXT
day=080721
sim=${SLURM_ARRAY_TASK_ID}

cd /gscratch/jrick/phylogenref/emp_tmp

if [ "$int" == "EXT" ]; then
	taxa_ref="Lcal"
	ref_genome="${REF_PATH}/references/lcalcarifer_genome_v3_chromosomal.fa"
else
	taxa_ref="Lmar"
	ref_genome="${REF_PATH}/references/lmariae_genome_Feb2018.fa"
fi

    
echo "starting empirical $tree_height analysis for $int, simulation $sim"
declare fastq_list=`ls /project/phylogenref/data/lates/*.fastq.gz | xargs -n 1 basename | sed 's/\.fastq\.gz//'`

##############################################
### Align to reference  ######################
### Make sure everything in indexed! #########
##############################################

#export ref_genome

#parallel --env ref_genome "echo 'Fastq: {} \n' && echo 'Mapping reads for {} \n' && bwa mem -t 16 ${ref_genome} /project/phylogenref/data/lates/{}.fastq.gz > aln_{}.sam && echo 'Converting sam to bam for {}' && samtools view -b -S -o aln_{}.bam aln_{}.sam && echo 'Sorting and indexing bam files for {} \n' && samtools sort aln_{}.bam -o aln_{}.sorted.bam && samtools index -c aln_{}.sorted.bam" ::: $fastq_list

#rm -f *.sam
#rm -f *[0-9].bam

#fi #DEBUGGING
#sims=1
#for sim in `seq $sims`; do
	
	if [[ ! -s OUTFILE_q${QUAL}_s${sim}_RN.bcf ]]; then
	ls /project/phylogenref/data/lates/*.sorted.bam | grep -v 'SRR' | shuf -n 30 > rand_ind_sim${sim}.txt
	echo "/project/phylogenref/data/lates/aln_SRR3140997.sorted.bam" >> rand_ind_sim${sim}.txt
	fi
	
	QUAL=40	
	if [[ ! -s OUTFILE_q${QUAL}_s${sim}_RN.bcf ]]; then
		samtools mpileup -g -t DP,AD \
			--skip-indels \
    	 	-P ILLUMINA \
        	-q $QUAL \
    	 	-f ${ref_genome} \
   	    	-b rand_ind_sim${sim}.txt \
     	  	-o OUTFILE_q${QUAL}_s${sim}.bcf 

      		if [[ -s OUTFILE_q${QUAL}_s${sim}.bcf ]]
        		then 
          			echo "OUTFILE_q${QUAL}_s${sim}.bcf not empty; moving on"
        		else
          			echo "OUTFILE_q${QUAL}_s${sim}.bcf empty; something went wrong"
          			exit 1
   			fi

		cat rand_ind_sim${sim}.txt > names     
		bcftools reheader -s names OUTFILE_q${QUAL}_s${sim}.bcf > OUTFILE_q${QUAL}_s${sim}_RN.bcf
	fi


##############################################
#### CREATING RAW VARIANTS FILE FROM BAMs ####
#### DO ONCE FOR EACH MAPPING QUALITY ########
#### FOR EACH GENOME #########################
##############################################
#if false; then
		bcftools call -m \
                --variants-only \
                --format-fields GQ \
                --skip-variants indels \
                OUTFILE_q${QUAL}_s${sim}_RN.bcf | bcftools filter \
                        --set-GTs . \
                        --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                                --min-alleles 2 \
                                --max-alleles 2 \
                                --types snps \
                                --apply-filter "PASS" \
                                --output-type v \
                                --output-file OUTFILE.s${sim}_q${QUAL}.vcf
#fi                
#                if [[ "$QUAL" -eq 40 ]]; then
#                    echo "Calculating Dxy from calc_dxy script."
#                    ${REF_PATH}/calc_dxy.sh OUTFILE.q${QUAL}.vcf $sim $tree_height $int $taxa_ref
#                    echo "done calculating Dxy"
#                else
#		    		echo "QUAL is $QUAL; not calculating Dxy this time"
#                fi

## delete bcf file if vcf file is successfully created

	if [[ -f "OUTFILE.s${sim}_q${QUAL}.vcf" ]]; then
		echo "VCF file exists and is not empty; deleting BCF file now"
		rm -f OUTFILE_q${QUAL}_s${sim}_RN.bcf
		rm -f OUTFILE_q${QUAL}_s${sim}.bcf
	else
		"VCF file is empty; exiting now"
		exit 1
	fi
                
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

	parallel --delay 5 --jobs 2  --line-buffer --env genes --env sim --env QUAL --env miss_list --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/each_maf_emp.sh {}" ::: $maf_list ::: $miss_list

#    mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#   	mv *OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#   	mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files

				
####Create a batch file with all normal trees in order by file type, and create a name file
#cat s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-emp-batch.trees
#ls s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-emp-tree.names

date

rm -f *sim${sim}*

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
