#!/bin/sh

source activate new_env
PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source


#####################################
######## Simulate (START)############
#####################################
#####################################

source /project/phylogenref/scripts/refbias_config.txt
sim=$1
tree_height=$2
taxa_ref=$(($(($RANDOM%$num_sp))+1))_0_0
int=INT

  
##############################################
### Simulate reads for each gene w/ TTR ######
###############################################
echo "starting analysis for $tree_height, sim $sim, $int"

nloci=$(($genes + 1))
var_sites=(`Rscript ${REF_PATH}/var_sites.R $nloci $varsites`)
echo ${var_sites[*]} >> ${output_dir}/${day}-varSites-${tree_height}-sim${sim}-${int}

#if false; then # DEBUGGING

export REF_PATH
export tree_height
export sim
export var_sites
export taxa_ref
export reference_prefix
export error
export int
export day

seq -w $genes | parallel --delay 2 --jobs 16  'snps=`head -n 1 ${REF_PATH}/output/new/${day}-varSites-${tree_height}-sim${sim}-${int} | tr " " "\n" | head -n {} | tail -n 1` && python2 ${REF_PATH}/write_config.py -treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees{}.trees -v `echo "$snps"` -ref $taxa_ref -path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa -o gene{}_sim${sim} -rate rat.matrix -g 5 -r 150 -f 500 -s 50 -c 20 -pre sim_ -errorfile $error > gene{}_sim${sim}_config && python2 /project/phylogenref/programs/TreeToReads/treetoreads.py gene{}_sim${sim}_config && grep -v "^>" ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa > ${reference_prefix}_gene{}.fa'

cat ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_sim${sim}.fa > ${reference_prefix}_sim${sim}.fa 


##############################################
### Concatenate reads into single file########
########### for each species #################
##############################################

echo "beginning processing fastq data"
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
                        #export "$fastq_list"
        fi
        ##export $fastq_list
        ##echo $fastq_list
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

bwa index ${reference_prefix}_sim${sim}.fa

export reference_prefix
export sim

parallel --env reference_prefix --env sim -j 16 --delay 2 "echo 'Mapping reads for {}' && bwa mem ${reference_prefix}_sim${sim}.fa fastq_reads/sim${sim}/fastq/{}/{}_${sim}_read1.fq.gz fastq_reads/sim${sim}/fastq/{}/{}_${sim}_read2.fq.gz > aln_{}.sam && echo 'Converting sam to bam for {}' && samtools view -b -S -o aln_{}.bam aln_{}.sam && echo 'Sorting and indexing bam files for {}' && samtools sort aln_{}.bam -o aln_{}.sorted.bam && samtools index -c aln_{}.sorted.bam" ::: $fastq_list
        
rm -f *.sam
rm -f *[0-9].bam

#fi #DEBUGGING
QUAL=40
samtools mpileup -g -t DP,AD \
                --skip-indels \
                -P ILLUMINA \
                -q $QUAL \
                -f ${reference_prefix}_sim${sim}.fa \
                *.sorted.bam \
                -o OUTFILE_q${QUAL}.bcf 
    
if [[ -s OUTFILE_q${QUAL}.bcf ]]
        then 
                echo "OUTFILE_q${QUAL}.bcf not empty; moving on"
        else
                echo "OUTFILE_q${QUAL}.bcf empty; something went wrong"
                exit 1
        fi

ls gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//' > names     
bcftools reheader -s names OUTFILE_q${QUAL}.bcf > OUTFILE_q${QUAL}_RN.bcf

rm -f gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz

##############################################
#### CREATING RAW VARIANTS FILE FROM BAMs ####
##############################################

bcftools call -m \
        --variants-only \
        --format-fields GQ \
        --skip-variants indels \
        OUTFILE_q${QUAL}_RN.bcf | bcftools filter \
                --set-GTs . \
                --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                        --min-alleles 2 \
                        --max-alleles 2 \
                        --types snps \
                        --apply-filter "PASS" \
                        --output-type v \
                        --output-file OUTFILE.s${sim}_q${QUAL}.vcf
                
        
echo "Calculating Dxy from calc_dxy script."
        ${REF_PATH}/calc_dxy.sh OUTFILE.s${sim}_q${QUAL}.vcf $sim $tree_height $int $taxa_ref $day
echo "done calculating Dxy"

rm -f OUTFILE_q${QUAL}.bcf
rm -f OUTFILE_q${QUAL}_RN.bcf
                
#################################
#### FILTERING WITH VCFTOOLS ####
#### AND INFERRING PHYLOGENY ####
#################################

echo "beginning parallel jobs per maf and miss to filter and infer phylogeny"

export QUAL
export sim
export REF_PATH
export output_dir
export day
export tree_height
export int
export miss_list
export genes

# the "each_maf.sh" script uses 8 threads for each job
parallel --delay 2 --jobs 4  --line-buffer --env genes --env sim --env QUAL --env miss_list --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/each_maf.sh {}" ::: $maf_list ::: $miss_list

# compile phylogenies
#mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files

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
