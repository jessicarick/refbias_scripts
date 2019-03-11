#! /bin/sh

#############################################################################
## This script is intended to be added into the refbias_Fullsim.sh script ##
#############################################################################

PATH=$PATH:/project/phylogenref/programs/ASTRAL/Astral
source /project/phylogenref/scripts/refbias_config.txt
sim=$1
tree_height=$2
taxa_ref=$3

#######################################
##### Setting parameter values ########
#######################################

# FOR DEBUGGING ################ -- for actual run: source /project/phylogenref/scripts/refbias_config.txt
#sim=1
#i=1
#REF_PATH=/project/phylogenref/scripts/slurm_results/SLURM_2079583
#reference=/project/phylogenref/scripts/references/lmariae_genome_Feb2018_concat.fa
#taxa_ref=0_0_0
#tree_height=500000
#reference_prefix=lmariae_500000
#ref_length=500000
#gene_length=1000
#genes=$(expr $ref_length / $gene_length)
#var_sites=5
#num_sp=100
#num_ind=1
#tree_height_list='500000 2000000 10000000'
#qual_list='0 20 40 60'
#maf_list='0 0.01 0.05'
#miss_list='0 0.50 0.75 1.0'
###############################

module load gcc
module load samtools
module load perl
module load vcftools
module load bwa
module load bcftools
module load raxml
module load r
module load jdk

# FOR ACTUAL RUN
#sim_fastq=`ls $REF_PATH/gene${i}_sim${sim}/fastq/sim*/*1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//'`

#######################################
######## Taking fastq files from ######
#### simulation, aligning to genomes ##
#######################################

mkdir astral/
touch astral/mutations.txt
echo "gene;num_SNPs;num_nonInv" >> mutations.txt

for i in `seq 100 110`;

	do sim_fastq=`ls gene${i}_sim${sim}/fastq/sim*/*1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//'`

mkdir astral/gene${i}_sim${sim}
cd astral/gene${i}_sim${sim}/

for fastq in $sim_fastq
	do echo "Fastq: $fastq \n"
        
    	sam=aln_${fastq}.sam
    	bam=aln_${fastq}.bam
    	sorted=aln_${fastq}.sorted.bam

    	echo "Mapping reads for $fastq \n"
#    	bwa index $REF_PATH/fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz
#    	bwa index $REF_PATH/fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz
    	bwa index ../../sim${sim}/ref.fasta_files/${reference_prefix}_sim${sim}.fa
	bwa mem -t 16 ../../sim${sim}/ref.fasta_files/${reference_prefix}_sim${sim}.fa \
    		../../gene${i}_sim${sim}/fastq/${fastq}/${fastq}_1.fq.gz \
    		../../gene${i}_sim${sim}/fastq/${fastq}/${fastq}_2.fq.gz > $sam

    echo "Converting sam to bam for $fastq \n"
    samtools view -b -S -o $bam $sam
       
    echo "Sorting and indexing bam files for $fastq \n"
    samtools sort $bam -o $sorted
    samtools index -c $sorted
done	


##############################################
#### CREATING RAW VARIANTS FILE FROM BAMs ####
#### DO ONCE FOR EACH MAPPING QUALITY ########
#### FOR EACH GENOME #########################
##############################################

echo "Calling variants on gene${i}_sim${sim} \n"

samtools mpileup -g -t DP,AD \
    	--skip-indels \
    	-P ILLUMINA \
    	-f ../../sim${sim}/ref.fasta_files/${reference_prefix}_sim${sim}.fa \
   	aln_*.sorted.bam \
     	-o OUTFILE.bcf 
     
	rm -f *.sam
#	rm *.bam
	rm -f *.bam.csi

	ls aln_*.sorted.bam | xargs -n 1 basename | sed 's/\.sorted\.bam//' | sed 's/aln_//' > names     

#	bcftools reheader -s OUTFILE.bcf > OUTFILE_RN.bcf

for QUAL in $qual_list
	do echo "quality $QUAL"
        bcftools call -m \
                --variants-only \
                --format-fields GQ \
                --skip-variants indels \
                OUTFILE.bcf | bcftools filter \
                        --set-GTs . \
                        --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
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
        		do echo "maf $maf; miss $miss"
	                vcftools --vcf OUTFILE.s${sim}_q${QUAL}.vcf \
        	                --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf} \
                	        --remove-filtered-all \
                        	--maf $(printf "$maf") \
                        	--max-missing $(printf "$miss") \
                        	--recode \
                        	--minDP 1 # keep the same?

#######################################
#### CONVERTING VCF TO PHYLIP FILE ####
#######################################
			echo "converting vcf to phy"
                
			python /project/phylogenref/scripts/vcf2phylip.py -i OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf \
                        	-o OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy \
                        
               		## rm OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf
     
#######################################
#### Removing invariant sites  ########
#######################################           
			echo "removing invariant sites"
			printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy', fileFormat='phy', extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy',format='phylip')" > Rscript_astral.R
                
                	R --vanilla --no-save < Rscript_astral.R
                
#######################################
#### Running RaxML w/ Ref #############
#######################################

			echo "running RAxML to create gene trees"
		
#			sites_ref=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy | head -n 1 | cut -f 2 -d' ')
			if [[ -s RAxML_info.OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.filtered.out ]]
 		  		then
  					rm -f RAxML_info.OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.filtered.out
          			else
           				echo "everything is good"
        		fi
        		
        		raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy -n OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.filtered.out -j --no-bfgs --silent -m ASC_GTRCAT -V --asc-corr=lewis -f a -x 223 -N 100 -p 466

#######################################
## move gene trees to common folder ###
#######################################

			mkdir ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/
        		cp OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/
        		cp RAxML* ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/

#######################################
#### Running RaxML w/o Ref ############
##### First delete ref from .phy ######
#######################################

			numtaxa=$(cat ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy | head -n1 | awk '{print $1;}')
			newtaxa=$(($numtaxa - 1))
			
			cat ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy | sed '/^'sim_$taxa_ref'/ d' | sed "1s/$numtaxa/$newtaxa/" > OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy
        		
			echo "removing invariant sites"
                	printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy', fileFormat='phy', extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy',format='phylip')" > Rscript_astral.R

                	R --vanilla --no-save < Rscript_astral.R

#			sites_noref=$(head -n 1 OUTFILE_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy | cut -f 2 -d' ')

        		if [[ -s RAxML_info.OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.filtered.out ]]
 		  		then
  					rm -f RAxML_info.OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.filtered.out
          			else
           				echo "everything is good"
        		fi
        		
        		raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy -n OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.filtered.out -j --no-bfgs --silent -m ASC_GTRCAT -V --asc-corr=lewis -f a -x 323 -N 100 -p 476 		

#######################################
## move gene trees to common folder ###
#######################################

			mkdir ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.gene_tree_files
        		cp OUTFILE_gene${i}_s${sim}_q${QUAL}_miss${miss}_maf${maf}*noInv*.phy ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.gene_tree_files/
        		cp RAxML* ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.gene_tree_files/

#######################################
#### run ASTRAL on the gene trees #####
#######################################
## need to now close all of the for loops!

			done # closes miss
		done # closes maf
	done # closes QUAL
cd ../../

mutations=$(grep -v '^#' astral/gene${i}_sim${sim}/OUTFILE_s${sim}_q0_miss0_maf0.recode.vcf | wc -l)
nonInv=$(head -n 1 astral/gene${i}_sim${sim}/OUTFILE_gene${i}_s${sim}_q0_miss0_maf0.noInv.REF.phy | cut -f 2 -d' ')
echo "height${tree_height}_sim${sim}_gene${i};$mutations;$nonInv" >> astral/mutations.txt

done # closes genes

# now, have to re-open miss, maf, qual
echo "Beginning ASTRAL analyses... wish us luck!"
cd astral/
mkdir species_trees/
mkdir astral_logs/

for QUAL in $qual_list
	do for maf in $maf_list
		do for miss in $miss_list
			do echo "running ASTRAL on gene trees for sim${sim} q${QUAL} maf${maf} miss${miss}"

			cd sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.gene_tree_files/
			echo "running ASTRAL on gene trees with reference"
		
			for tree in `ls RAxML*bestTree*.out`
				do cat $tree >> s${sim}_q${QUAL}_miss${miss}_maf${maf}.gene_tree_in.tree
			done
		
			java -jar /project/phylogenref/programs/ASTRAL/Astral/astral.5.6.1.jar -i s${sim}_q${QUAL}_miss${miss}_maf${maf}.gene_tree_in.tree -o s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.species_tree_out.tre 2>s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.EXT.species_tree_out.log 

## if we wanted to use bootstrap trees!
			#java -jar astral.5.6.1.jar -i in.tree.best -b bs_paths -r 100 -o out.tre 2>out.log ## if you have bootstrap gene trees

			mv *species_tree_out.tre ../species_trees/
			mv *species_tree_out.log ../astral_logs/
		
			cd ../sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.gene_tree_files
	
			echo "running ASTRAL on gene trees without reference"
		
			for tree in `ls RAxML*bestTree*.out`
				do cat $tree >> s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.gene_tree_in.tree
			done
		
			java -jar /project/phylogenref/programs/ASTRAL/Astral/astral.5.6.1.jar -i s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.gene_tree_in.tree -o s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.species_tree_out.tre 2>s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.EXT.species_tree_out.log

			#java -jar astral.5.6.1.jar -i in.tree.best -b bs_paths -r 100 -o out.tre 2>out.log ## if you have bootstrap gene trees

			mv *species_tree_out.tre ../species_trees/
			mv *species_tree_out.log ../astral_logs/
		cd ../

#######################################
#### close all of the for loops #######
#######################################

       		done # closes $miss
         done # closes $maf
    done # closes $QUAL

day=`date +%m%d%y`
cat species_trees/*.tre > /project/phylogenref/scripts/sim_results/${day}-${tree_height}-astral.trees
ls species_trees/*.tre > /project/phylogenref/scripts/sim_results/${day}-${tree_height}-astral.names
