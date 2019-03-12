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

taxa_ref=$(($(($RANDOM%$num_sp))+1))_0_$(($RANDOM%$num_ind))
    
##############################################
### Simulate reads for each gene w/ TTR ######
##############################################

	for i in `seq $genes`;
	 do python ${REF_PATH}/write_config.py -treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees${i}.trees -v $var_sites -ref $taxa_ref -path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_${i}.fa -o gene${i}_sim${sim} -rate rat.matrix -g 5 -r 150 -f 400 -s 20 -c 13 -pre sim_ -errorfile $error > gene${i}_sim${sim}_config
		python /project/phylogenref/programs/TreeToReads/treetoreads.py gene${i}_sim${sim}_config
		
		cat ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_${i}.fa | sed 's/^\(>lmariae_genome_Feb2018_10000\)*//' > ${reference_prefix}_gene${i}.fa
		
	done

		cat ${reference_prefix}_gene*.fa | sed '1 i\>lmariae_genome_Feb2018' > ${reference_prefix}_sim${sim}.fa
		

##############################################
### Concatenate reads into single file########
########### for each species #################
##############################################
fastq_list=0

while [ $(echo "$(expr length "$fastq_list")") -lt 2 ]
do
        i=$(shuf -i 1-$gene_length -n 1)
        if [ ! -f gene${i}_sim${sim}/fastq/sim*/sim_0_0_0*_1.fq.gz ]
                then
                        echo "File empty. try again"
                else
                        echo "File present, well done"
                        declare fastq_list=`ls gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//'`
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

	for fastq in $fastq_list
      	do echo "Fastq: $fastq \n"
        
       	sam=aln_${fastq}.sam
       	bam=aln_${fastq}.bam
       	sorted=aln_${fastq}.sorted.bam

      	 echo "Mapping reads for $fastq \n"
      	 bwa index ${reference_prefix}_sim${sim}.fa
       	 #bwa index fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz
       	 #bwa index fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz
         bwa mem -t 16 ${reference_prefix}_sim${sim}.fa fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz > $sam

       	 echo "Converting sam to bam for $fastq \n"
       	 samtools view -b -S -o $bam $sam
        
       	 echo "Sorting and indexing bam files for $fastq \n"
       	 samtools sort $bam -o $sorted
         samtools index -c $sorted
	done	

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

		rm *sam

		ls gene${i}_sim${sim}/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//' > names     
		bcftools reheader -s names OUTFILE_q${QUAL}.bcf > OUTFILE_q${QUAL}_RN.bcf

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
                        --include $(printf "QUAL>40")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                                --min-alleles 2 \
                                --max-alleles 2 \
                                --types snps \
                                --apply-filter "PASS" \
                                --output-type v \
                                --output-file OUTFILE.s${sim}_q${QUAL}.vcf

#################################
#### FILTERING WITH VCFTOOLS ####
#################################

        for maf in $maf_list
        do for miss in $miss_list
                do vcftools --vcf OUTFILE.s${sim}_q${QUAL}.vcf \
                        --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf} \
                        --remove-filtered-all \
                        --maf $(printf "$maf") \
                        --max-missing $(printf "$miss") \
                        --recode \
                        --recode-INFO-all \
                        --minDP 1 # keep the same?

#######################################
#### CONVERTING VCF TO PHYLIP FILE ####
#######################################

                python ${REF_PATH}/vcf2phylip.py -i OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf \
                        -o OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy \
                        
                rm OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf
     
#######################################
#### Removing invariant sites  ########
#######################################           

		            printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy', fileFormat='phy', extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy',format='phylip')" > Rscript.R
                
                R --vanilla --no-save < Rscript.R
                
#######################################
#### Running RaxML w/ Ref #############
#######################################
              
              let sites_ref=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy | head -n 1 | cut -f 2 -d' ')
			        
               if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.INT.filtered.out ]]
 		  			      then
  						      rm RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.INT.filtered.out
          			  else
           				  echo "everything is good"
        	     fi
        		
              raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.INT.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466
        
#######################################
#### Running RaxML w/o Ref ############
##### First delete ref from .phy ######
#######################################
###Here is an issue...awk ###

			numtaxa=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy | head -n1 | awk '{print $1;}')
			newtaxa=$(($numtaxa - 1))
			
			cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy | sed '/^'sim_${taxa_ref}'/ d' | sed "1s/$numtaxa/$newtaxa/" > OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy
			
			 printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy', fileFormat='phy', extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy',format='phylip')" > Rscript.R
                
                R --vanilla --no-save < Rscript.R
                
      sites_noref=$(head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy | cut -f 2 -d' ')

        	if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.INT.filtered.out ]]
 		  			then
  						rm RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.INT.filtered.out
          			else
           				echo "everything is good"
        	fi
        		
        		raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.INT.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 323 -N 100 -p 476
        
        
				mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.INT.ref-${taxa_ref}.phylip_tree_files
        		mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.INT.ref-${taxa_ref}.phylip_tree_files
        		mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.INT.ref-${taxa_ref}.phylip_tree_files

   		
                done
         done

    done

rm *.bam
   
mkdir -p sim${sim}/config_files
mkdir -p sim${sim}/ref.fasta_files
		
mv *_config sim${sim}/config_files
mv *.fa sim${sim}/ref.fasta_files
				
###Create a batch file with all trees in order by file type, and create a name file
cat s${sim}*q*miss*/*bestTree* >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}**q*miss*/*bestTree* >> ${output_dir}/${day}-${tree_height}-tree.names

/project/phylogenref/scripts/astral_script_012819.sh $sim $tree_height $taxa_ref INT

