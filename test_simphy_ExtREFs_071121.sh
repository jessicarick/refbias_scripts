#!/bin/sh

source activate new_env
PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source


#####################################
########Simulate (START)#############
####################################
#####################################

source /project/phylogenref/scripts/refbias_config.txt
sim=$1
tree_height=$2
taxa_ref=0_0_0
int=EXT

#reference_prefix=lmariae_${tree_height}_${sim}_
    
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

seq -w $genes | parallel --delay 5 --jobs 32  'snps=`head -n 1 ${REF_PATH}/output/new/${day}-varSites-${tree_height}-sim${sim}-${int} | tr " " "\n" | head -n {} | tail -n 1` && python ${REF_PATH}/write_config.py -treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees{}.trees -v `echo "$snps"` -ref $taxa_ref -path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa -o gene{}_sim${sim} -rate rat.matrix -g 5 -r 150 -f 500 -s 50 -c 20 -pre sim_ -errorfile $error > gene{}_sim${sim}_config && python /project/phylogenref/programs/TreeToReads/treetoreads.py gene{}_sim${sim}_config && grep -v "^>" ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa > ${reference_prefix}_gene{}.fa'

#exit 0 # debugging

#seq -w $genes | parallel --delay 5 --jobs 32 --env REF_PATH --env tree_height --env sim --env var_sites --env taxa_ref --env reference_prefix --env error "index=$(echo {} | sed 's/^0*//') &&  echo 'index is $index' && echo 'varsites is `echo ${var_sites[$index]}`' && exit 0 && python ${REF_PATH}/write_config.py -treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees{}.trees -v `echo ${var_sites[$index]}` -ref $taxa_ref -path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa -o gene{}_sim${sim} -rate rat.matrix -g 5 -r 150 -f 500 -s 50 -c 20 -pre sim_ -errorfile $error > gene{}_sim${sim}_config && python /project/phylogenref/programs/TreeToReads/treetoreads.py gene{}_sim${sim}_config && grep -v '^>' ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa > ${reference_prefix}_gene{}.fa"

#echo ">lmariae_genome_Feb2018_${ref_length}" > ${reference_prefix}_sim${sim}.fa
cat ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_sim${sim}.fa > ${reference_prefix}_sim${sim}.fa 

#fi # DEBUGGING



##############################################
### Concatenate reads into single file########
########### for each species #################
##############################################
echo "beginning RAxML analyses"
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

parallel --env reference_prefix --env sim -j 32 --delay 5 "echo 'Mapping reads for {}' && bwa mem ${reference_prefix}_sim${sim}.fa fastq_reads/sim${sim}/fastq/{}/{}_${sim}_read1.fq.gz fastq_reads/sim${sim}/fastq/{}/{}_${sim}_read2.fq.gz > aln_{}.sam && echo 'Converting sam to bam for {}' && samtools view -b -S -o aln_{}.bam aln_{}.sam && echo 'Sorting and indexing bam files for {}' && samtools sort aln_{}.bam -o aln_{}.sorted.bam && samtools index -c aln_{}.sorted.bam" ::: $fastq_list
	
rm -f *.sam
rm -f *[0-9].bam

#fi #DEBUGGING

for QUAL in $qual_list
	do samtools mpileup -g -t DP,AD \
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
#### DO ONCE FOR EACH MAPPING QUALITY ########
#### FOR EACH GENOME #########################
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
                
                if [[ "$QUAL" -eq 40 ]]; then
                    echo "Calculating Dxy from calc_dxy script."
                    ${REF_PATH}/calc_dxy.sh OUTFILE.s${sim}_q${QUAL}.vcf $sim $tree_height $int $taxa_ref $day
                    echo "done calculating Dxy"
                else
		    echo "QUAL is $QUAL; not calculating Dxy this time"
                fi
                
#################################
#### FILTERING WITH VCFTOOLS ####
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

	parallel --delay 5 --jobs 4  --line-buffer --env genes --env sim --env QUAL --env miss_list --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/each_maf.sh {}" ::: $maf_list ::: $miss_list


#        for maf in $maf_list
#        	do for miss in $miss_list
#                	do vcftools --vcf OUTFILE.s${sim}_q${QUAL}.vcf \
#                        	--out OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf} \
#                        	--remove-filtered-all \
#                        	--maf $(printf "$maf") \
#                        	--max-missing $(printf "$miss") \
#                        	--recode \
#                        	--recode-INFO-all \
#                        	--minDP 1 # keep the same?
#

#######################################
#### CONVERTING VCF TO PHYLIP FILE ####
#### and removing invariant sites #####
#### making one phylip per gene #######
#######################################
#export miss
#export maf
#export QUAL
#export sim
#export REF_PATH
#export output_dir
#export day
#export tree_height
#export int
#
#seq -w ${genes} | parallel --delay 5 --jobs 8 --env sim --env QUAL --env miss --env maf --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/filter_gene.sh {}"

#         for gene in `seq -w ${genes}`
#                 do vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf \
#                       --chr gene${gene} \
#                       --recode \
#                       --out gene${gene}.REF
#                       
#                 python ${REF_PATH}/vcf2phylip.py -i gene${gene}.REF.recode.vcf \
#                     	-o gene${gene}.REF.phy
#                 
#		sites=`head -n 1 gene${gene}.REF.phy | awk '{print $2}'`
#		if [ "$sites" -eq "0" ]; then
#			echo "no sites left"
#			cat gene${gene}.REF.phy > gene${gene}.REF.noInv.phy
#		else
#                	printf "library(ape)\nlibrary(phrynomics)\nReadSNP('gene${gene}.REF.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='gene${gene}.REF.noInv.phy',format='phylip')" > Rscript.R
#                 	R --vanilla --no-save < Rscript.R    	
#                fi
#
#                 nsnps=`cat gene${gene}.REF.noInv.phy | head -n 1 | awk '{print $2}'`
#                 echo "gene${gene},s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-${tree_height}-${int}
#         done
         
         ## combine into one supermatrix
         
#         Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF
#         rm -f gene*.REF.noInv.phy
         
#######################################
#### Running RaxML w/ Ref #############
#######################################
              
#         sites_ref=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.phy | head -n 1 | awk '{print $2}' )
 			        
#         if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out ]]
# 		        then
#                 rm -f RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out
#           	else
#            			echo "everything is good; running raxml"
#         fi
         	     	  
#         echo "running raxml on concatenated SNPs"
#        raxmlHPC-PTHREADS-AVX -T 1 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466 

#######################################
#### Running RaxML w/o Ref ############
#######################################
#### CONVERTING VCF TO PHYLIP FILE ####
#### and removing invariant sites #####
#### making one phylip per gene #######
#######################################

#        for gene in `seq -w ${genes}`
#                do vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf \
#                      --chr gene${gene} \
#                      --recode \
#                      --out gene${gene}.NOREF
#                      
#                python ${REF_PATH}/vcf2phylip.py -i gene${gene}.NOREF.recode.vcf \
#                    	-o gene${gene}.NOREF.phy
#                
#                sites=`cat gene${gene}.NOREF.phy | head -n 1 | awk '{print $2}'`
#		if [ "$sites" -eq "0" ]; then
#			cat gene${gene}.NOREF.phy > gene${gene}.NOREF.noInv.phy
#		else
#			printf "library(ape)\nlibrary(phrynomics)\nReadSNP('gene${gene}.NOREF.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='gene${gene}.NOREF.noInv.phy',format='phylip')" > Rscript.R
#                	R --vanilla --no-save < Rscript.R    	
#                fi
#                nsnps=`cat gene${gene}.NOREF.noInv.phy | head -n 1 | awk '{print $2}'`
#                
#		echo "gene${gene},s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-${tree_height}-${int}
#        done
#        
#        ## combine into one supermatrix
#
#        Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF
#        rm -f gene*.NOREF.noInv.phy   
#        
#        
#      	sites_noref=$(head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.all.noInv.phy | awk '{print $2}' )
#
#        if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out ]]
# 			        then
#  		            rm -f RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out
#          		else
#           		    echo "everything is good"
#        fi
#        		
#      	echo "running raxml on concatenated SNPs"
#      	raxmlHPC-PTHREADS-AVX -T 16 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.all.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 323 -N 100 -p 476
#                  
#         
#        echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf},${sites_ref},${sites_noref}" >> ${output_dir}/${day}-filteredSites-${tree_height}-${int}

#######################################
#### Subsample & Run RAxML ############
#######################################

#           Rscript ${REF_PATH}/subsample.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv $sites_ref $maxSNP
# 
# 	sites_samp=$(head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy | awk '{print $2}')
#           
#           echo "running raxml on subsampled phylip"
#               	raxmlHPC-PTHREADS-AVX -T 16 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_samp}.REF.${int}.subsamp.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466 

        mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
       	mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
       	mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
#
#                done
#         done

    done

mkdir -p sim${sim}/config_files
mkdir -p sim${sim}/ref.fasta_files
		
mv *_config sim${sim}/config_files
mv *.fa sim${sim}/ref.fasta_files
				
####Create a batch file with all normal trees in order by file type, and create a name file
cat s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-tree.names

## ###Create a batch file with all subsampled trees in order by file type, and create a name file
## cat s${sim}*q*miss*/*bipartitions*subsamp* >> ${output_dir}/${day}-${tree_height}-subsamp-batch.trees
## ls s${sim}*q*miss*/*bipartitions*subsamp* >> ${output_dir}/${day}-${tree_height}-subsamp-tree.names

rm -f gene*.phy
rm -f gene*/fasta_files/*.fasta

date

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
