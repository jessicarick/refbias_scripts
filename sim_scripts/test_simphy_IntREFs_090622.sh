#!/bin/sh

cp -r /project/phylogenref/scripts/sim_scripts .

source sim_scripts/refbias_config.txt
source activate $conda_env
module load jdk

#####################################
########Simulate (START)#############
####################################
#####################################

sim=$1
tree_height=$2
taxa_ref=$(grep $tree_height ${output_dir}/100922-output-int-taxa_ref.csv | grep $sim | cut -f 3 -d',')
#taxa_ref=$(($(($RANDOM%$num_sp))+1))_0_0
int=INT
subsample=true
    
##############################################
### Simulate reads for each gene w/ TTR ######
###############################################
echo "starting analysis for $tree_height, sim $sim, $int"

if [ ! -f ${output_dir}/${day}-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz ]; then

nloci=$(($genes + 1))
var_sites=(`Rscript sim_scripts/var_sites.R $nloci $varsites`)
echo ${var_sites[*]} > ${output_dir}/${day}-varSites-${tree_height}-sim${sim}-${int}

export sim
export taxa_ref
export int
export tree_height

seq -w $genes | parallel --delay 2 --jobs 16 --env sim --env int --env taxa_ref --env tree_height "sim_scripts/run_treetoreads.sh {}"

echo "done with treetoreads runs!"

for gene in `seq -w $genes`
	do if [ "$gene" -eq "0001" ]; then
		echo ">gene${gene}" > ${reference_prefix}_sim${sim}_${int}.fa
		grep '^[ACGTN]' gene${gene}_sim${sim}/fasta_files/sim_${taxa_ref}.fasta >> ${reference_prefix}_sim${sim}_${int}.fa
	else
		echo ">gene$gene" >> ${reference_prefix}_sim${sim}_${int}.fa
		grep '^[ATGCN]' gene${gene}_sim${sim}/fasta_files/sim_${taxa_ref}.fasta >> ${reference_prefix}_sim${sim}_${int}.fa
	fi
done

echo "indexing simulated genome"
bwa index ${reference_prefix}_sim${sim}_${int}.fa

echo "calculating RF index of gene tree discordance"
cat gene*_sim${sim}/simtree.tre >> genetrees_TTR_sim${sim}.tre
Rscript sim_scripts/calc_genetree_rf.R genetrees_TTR_sim${sim}.tre $sim >> ${output_dir}/${day}-${tree_height}-${int}.TTR.gt_rf

##############################################
### Concatenate reads into single file########
########### for each species #################
##############################################
echo "beginning alignment of fastq reads"
fastq_list=0

while [ $(echo "$(expr length "$fastq_list")") -lt 2 ]
do
        i=$(seq -w $genes | shuf -n 1)
        if [ ! -f gene${i}_sim${sim}/fastq/sim_0_0_0/sim_0_0_0_1.fq.gz ]
                then
                        echo "File empty. try again"
                else
                        echo "File present, well done"
                        declare fastq_list=`ls gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1\.fq\.gz//'`
        fi
done


for fastq in $fastq_list		
        do mkdir -p fastq_reads/sim${sim}/fastq/${fastq}
        zcat gene*_sim${sim}/fastq/${fastq}/${fastq}_1.fq.gz | gzip > fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz
	zcat gene*_sim${sim}/fastq/${fastq}/${fastq}_2.fq.gz | gzip > fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz
done

##############################################
### Align to reference  ######################
### Make sure everything in indexed! #########
##############################################
#fi #debugging
export reference_prefix
export sim
export int

parallel --delay 2 --env int --env reference_prefix --env sim -j 16 "bash sim_scripts/run_bwa.sh {}" ::: $fastq_list
	
rm -f *.sam
rm -f *[0-9].bam

#fi #DEBUGGING
samtools mpileup -g -t DP,AD \
		--skip-indels \
    	 	-P ILLUMINA \
        	-q $QUAL \
    	 	-f ${reference_prefix}_sim${sim}_${int}.fa \
   	    	*.sorted.bam \
     	  	-o OUTFILE_q${QUAL}_${int}.bcf 
    
if [[ -s OUTFILE_q${QUAL}_${int}.bcf ]]
        then 
          	echo "OUTFILE_q${QUAL}_${int}.bcf not empty; moving on and removing bam files"
                rm -f *.sorted.bam 
        else
          	echo "OUTFILE_q${QUAL}_${int}.bcf empty; something went wrong"
          	exit 1
   	fi

ls gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//' > names     
bcftools reheader -s names OUTFILE_q${QUAL}_${int}.bcf > OUTFILE_q${QUAL}_${int}_RN.bcf

rm -rf gene*_sim${sim}/

##############################################
#### CREATING RAW VARIANTS FILE FROM BAMs ####
##############################################

bcftools call -m \
        --variants-only \
        --format-fields GQ \
        --skip-variants indels \
        OUTFILE_q${QUAL}_${int}_RN.bcf | bcftools filter \
                --set-GTs . \
                --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                        --min-alleles 2 \
                        --max-alleles 2 \
                        --types snps \
                        --apply-filter "PASS" \
                        --output-type v \
                        --output-file OUTFILE_s${sim}_q${QUAL}_${int}.vcf
                
gzip -c OUTFILE_s${sim}_q${QUAL}_${int}.vcf > ${output_dir}/${day}-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz

echo "Calculating Dxy from calc_dxy script."
        ./sim_scripts/calc_dxy.sh OUTFILE_s${sim}_q${QUAL}_${int}.vcf $sim $tree_height $int $taxa_ref $day
echo "done calculating Dxy"

rm -f OUTFILE_q${QUAL}.bcf
rm -f OUTFILE_q${QUAL}_RN.bcf
 
fi

###############################################
## STARTING SUBSAMPLING SCRIPT WITH THIS VCF ##
###############################################
## note: not tested yet!!
if false; then #debugging
if [ "$subsample" = true ]; then
        echo "spawning subsampling script for $sim $tree_height $int for the $day data"
        sbatch sim_scripts/subsample_script.sh $sim $tree_height $int $taxa_ref $day
fi
fi # debugging

#cp ${output_dir}/${day}-${tree_height}-OUTFILE_s${sim}_q${QUAL}.${int}.vcf OUTFILE_s${sim}_q${QUAL}.vcf 

#################################
#### FILTERING WITH VCFTOOLS ####
#### AND INFERRING PHYLOGENY ####
#################################

echo "beginning parallel jobs per mac and miss"

export sim
export int
export taxa_ref
export tree_height

# the "each_mac.sh" script uses 4 threads for each job
parallel --delay 2 --jobs 4  --line-buffer --env tree_height --env taxa_ref --env sim --env int "bash sim_scripts/each_mac_jun2022.sh {}" ::: $mac_list ::: $miss_list

# compile phylogenies
#mkdir s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
#mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}*.phy s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
#mv RAxML* s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files

mkdir -p sim${sim}/config_files
mkdir -p sim${sim}/ref.fasta_files
		
mv *_config sim${sim}/config_files
mv *.fa sim${sim}/ref.fasta_files
				
####Create a batch file with all normal trees in order by file type, and create a name file
#cat s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-batch.trees
#ls s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-tree.names

rm -f gene*.phy
rm -f gene*/fasta_files/*.fasta

date

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
