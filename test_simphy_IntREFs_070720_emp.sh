#!/bin/sh

source activate new_env
PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source

#####################################
########Simulate (START)#############
####################################
#####################################

source /project/phylogenref/scripts/refbias_config.txt
tree_height="lates"
int=INT
if [ "$int" == "EXT" ]; then
	taxa_ref="Lcal"
	ref_genome="${REF_PATH}/references/lcalcarifer_genome_v3_chromosomal.fa"
else
	taxa_ref="Lmar"
	ref_genome="${REF_PATH}/references/lmariae_genome_Feb2018.fa"
fi

#reference_prefix=lmariae_${tree_height}_${sim}_
    
##############################################
### Simulate reads for each gene w/ TTR ######
###############################################
echo "starting analysis for $int"

declare fastq_list=`ls /project/phylogenref/data/fastqs/lates/*.fastq.gz | xargs -n 1 basename | sed 's/\.fastq\.gz//'`

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
      bwa mem -t 16 ${ref_genome} \
	      	/project/phylogenref/data/fastqs/${fastq}/${fastq}.fastq.gz > $sam	
	  
      echo "Converting sam to bam for $fastq \n"
      samtools view -b -S -o $bam $sam
        
      echo "Sorting and indexing bam files for $fastq \n"
      samtools sort $bam -o $sorted
      samtools index -c $sorted
done

rm -f *.sam
rm -f *[0-9].bam

#fi #DEBUGGING

for sim in `seq $sims`; do
	ls *.sorted.bam | grep -v 'SRR' | shuf -n 100 > rand_ind
	echo "SRR3140997.sorted.bam" >> rand_ind

	for QUAL in $qual_list
		do samtools mpileup -g -t DP,AD \
			--skip-indels \
    	 	-P ILLUMINA \
        	-q $QUAL \
    	 	-f ${ref_genome} \
   	    	-l rand_ind \
     	  	-o OUTFILE_q${QUAL}.bcf 
    
      		if [[ -s OUTFILE_q${QUAL}.bcf ]]
        		then 
          			echo "OUTFILE_q${QUAL}.bcf not empty; moving on"
        		else
          			echo "OUTFILE_q${QUAL}.bcf empty; something went wrong"
          			exit 1
   			fi

		echo "$fastq_list" > names     
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
                        --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                                --min-alleles 2 \
                                --max-alleles 2 \
                                --types snps \
                                --apply-filter "PASS" \
                                --output-type v \
                                --output-file OUTFILE.s${sim}_q${QUAL}.vcf
                
                if [[ "$QUAL" -eq 40 ]]; then
                    echo "Calculating Dxy from calc_dxy script."
                    ${REF_PATH}/calc_dxy.sh OUTFILE.q${QUAL}.vcf $sim $tree_height $int $taxa_ref
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

	parallel --delay 5 --jobs 4  --line-buffer --env genes --env sim --env QUAL --env miss_list --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/each_maf_emp.sh {}" ::: $maf_list ::: $miss_list

    mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
   	mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files
   	mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}-${taxa_ref}.phylip_tree_files

    done # closing qual loop
done # closing sims loop
				
####Create a batch file with all normal trees in order by file type, and create a name file
cat s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}*q*miss*/*bipartitions.*filtered* >> ${output_dir}/${day}-${tree_height}-tree.names

date

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
