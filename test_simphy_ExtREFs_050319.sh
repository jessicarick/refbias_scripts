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
##############################################
nloci=$(($genes + 1))
var_sites=(`Rscript ${REF_PATH}/var_sites.R $nloci`)
echo ${var_sites[*]} >> /project/phylogenref/scripts/output/${day}-varSites-${tree_height}-sim${sim}-${int}

export REF_PATH
export tree_height
export sim
export var_sites
export taxa_ref
export reference_prefix
export error
  
seq -w $genes | parallel --delay 5 --env REF_PATH --env tree_height --env sim --env var_sites --env taxa_ref --env reference_prefix --env error "index=$(echo {} | sed 's/^0*//') &&  python ${REF_PATH}/write_config.py -treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees{}.trees -v `echo ${var_sites[$index]}` -ref $taxa_ref -path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa -o gene{}_sim${sim} -rate rat.matrix -g 5 -r 150 -f 400 -s 20 -c 13 -pre sim_ -errorfile $error > gene{}_sim${sim}_config && python /project/phylogenref/programs/TreeToReads/treetoreads.py gene{}_sim${sim}_config && grep -v '^>' ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_{}.fa > ${reference_prefix}_gene{}.fa"

echo ">lmariae_genome_Feb2018_${ref_length}" > ${reference_prefix}_sim${sim}.fa
cat ${reference_prefix}_gene*.fa >> ${reference_prefix}_sim${sim}.fa 
		
if [ ! -d astral/ ]; then
    mkdir astral/
fi

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

for fastq in $fastq_list
      do echo "Fastq: $fastq \n"
        
      sam=aln_${fastq}.sam
      bam=aln_${fastq}.bam
      sorted=aln_${fastq}.sorted.bam

      echo "Mapping reads for $fastq \n"
      bwa index ${reference_prefix}_sim${sim}.fa
       #bwa index fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz
       #bwa index fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz
      bwa mem -t 32 ${reference_prefix}_sim${sim}.fa fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read1.fq.gz fastq_reads/sim${sim}/fastq/${fastq}/${fastq}_${sim}_read2.fq.gz > $sam

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

	rm -f *sam

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
                        --include $(printf "QUAL>$QUAL")$(printf "&&")$(printf "FMT/GQ>10") | bcftools view \
                                --min-alleles 2 \
                                --max-alleles 2 \
                                --types snps \
                                --apply-filter "PASS" \
                                --output-type v \
                                --output-file OUTFILE.s${sim}_q${QUAL}.vcf
                
                if [[ "$QUAL" -eq 0 ]]; then
                    echo "QUAL=0; Calculating Dxy from calc_dxy script."
                    ${REF_PATH}/calc_dxy.sh OUTFILE.s${sim}_q${QUAL}.vcf $sim $tree_height $int $taxa_ref
                    echo "done calculating Dxy"
                else
		                echo "QUAL is $QUAL; not calculating Dxy this time"
		            fi
                
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
                    	-o OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy 

##################################    
#### CREATING PARTITIONS FILE ####
##################################
                  
		phy=OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy

                echo '' > sites.txt
                for locus in `seq -w $genes`
                      	do name=`echo gene${locus}`
                	sites=`vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf --chr lmariae_genome_Feb2018_${ref_length} --from-bp $((1000*$((10#${locus}-1))+1)) --to-bp $((1000*10#${locus})) --recode --stdout | grep -v  '^#' | wc -l`
                	echo "${name},${sites}" >> sites.txt
                done
                  
                echo '' > phy.tmp
                echo '' > partitions.txt
                echo '' > phynames
		for line in `cat sites.txt`;
                  	do gene=`echo "$line" | cut -f 1 -d','`
                      	nsnp=`echo $line | cut -f 2 -d','`
                      	if (( $nsnp > 0 )); then
                        	prev=`tail -n 1 partitions.txt | cut  -f 2 -d'-'`
                        	begin=$((prev+1))
                        	end=$((begin + nsnp))
                        	echo "DNA,  ${gene}=${begin}-${end}" >> partitions.txt
                        	echo $gene >> phynames
        
                      		cat $phy | tail -n +2 | awk '{print $2}' | cut -c${start}-${end} > seq.tmp
				if [ -f phy.tmp ]; then
						paste phy.tmp seq.tmp > phy2.tmp
						mv phy2.tmp phy.tmp
					else
						cp seq.tmp phy.tmp
				fi
				rm -f seq.tmp
                	fi
                done
		cat $phy | tail -n +2 | cut -f 1 -d' ' > indnames
		paste indnames phy.tmp > phy2.tmp
		mv phy2.tmp phy.tmp
		#locnames=`cat phynames | tr "\n" "\t"`
		#sed -i "1s/^/${locnames}\n/" phy.tmp
                  
                #rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf
     
#######################################
#### Removing invariant sites  ########
#######################################           

		printf "library(ape)\nlibrary(phrynomics)\nphynames<-scan('phynames',character())\nphy<-read.table('phy.tmp',row.names=1)\ncolnames(phy) <- phynames\nReadSNP(phy)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy',format='phylip')\nnsnps <- GetNumberOfSitesForLocus(snps)\nwrite(nsnps, file='nsnps_per_loc_ref', sep=' ', ncolumns=length(nsnps))\nwrite(names(nsnps), file='nsnps_locus_names_ref', sep=' ', ncolumns=length(nsnps))" > Rscript.R
                
                R --vanilla --no-save < Rscript.R
               
		echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy" >> /project/phylogenref/scripts/output/${day}-SNPs-${tree_height}-sim${sim}-${int}

                cat nsnps_locus_names_ref | tr ' ' '\n' > locnames.tmp
                cat nsnps_per_loc_ref | tr ' ' '\n' | awk 'p{print $0-p}{p=$0}' > locnsnps.tmp	

                paste locnames.tmp locnsnps.tmp >> /project/phylogenref/scripts/output/${day}-SNPs-${tree_height}-sim${sim}-${int}
                  
 
#######################################
#### Running RaxML w/ Ref #############
#######################################
              
              	sites_ref=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy | head -n 1 | cut -f 2 -d' ')
			        
               	if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out ]]
		        then
                        	rm -f RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out
          		else
           			echo "everything is good; running raxml"
        	fi
        	     	  
        	echo "creating partition for gene trees"
        	python ${REF_PATH}/snp_part_raxml.py -snpfile nsnps_per_loc_ref -namesfile nsnps_locus_names_ref > raxml_partitions.txt
        		      
        	echo "running raxml to create gene tree alignments"
        	raxmlHPC-PTHREADS-AVX -T 2 -f s -q raxml_partitions.txt -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy -m ASC_GTRGAMMA --asc-corr=lewis -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.partitions.out
        		      
        	echo "running raxml on concatenated SNPs"
              	raxmlHPC-PTHREADS-AVX -T 16 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466 

#               echo "running raxml to create gene trees"
		echo "not creating gene trees right now"
#		export int
#               ls OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy.locus*.phy | parallel --delay 10 --ungroup --env int 'echo "starting raxml on {}" && raxmlHPC-PTHREADS-AVX -T 2 -s {} -n {.}.REF.${int}.filtered.out -j --no-bfgs --silent -m ASC_GTRCAT -V --asc-corr=lewis -f d -F -p 466'

#                if [ ! -d  astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.gene_tree_files/ ]; then
#			mkdir astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.gene_tree_files
#		fi
        		
#        	mv RAxML*locus*.out astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.gene_tree_files/
                  
                echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy" >> /project/phylogenref/scripts/output/${day}-SNPs-${tree_height}-sim${sim}-${int}
                  
		rm -f locnames.tmp
		rm -f locsnps.tmp
                   
#######################################
#### Running RaxML w/o Ref ############
##### First delete ref from .phy ######
#######################################

		numtaxa=$(cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy | head -n1 | awk '{print $1;}')
		newtaxa=$(($numtaxa - 1))
			
		cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy | sed '/^'sim_${taxa_ref}'/ d' | sed "1s/$numtaxa/$newtaxa/" > OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy
			

##################################
#### CREATING PARTITIONS FILE ####
##################################

                phy=OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy
                
		#echo '' > sites.txt
                #for locus in `seq -w $genes`
                #    do name=`echo gene${locus}`
                #    sites=`vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf --chr lmariae_genome_Feb2018_${ref_length} --from-bp $((1000*$((10#${locus}-1))+1)) --to-bp $((1000*10#${locus})) --recode --stdout | grep -v  '^#' | wc -l`
                #    echo "${name},${sites}" >> sites.txt
                #done

                echo '' > phy.tmp
                echo '' > partitions.txt
                echo '' > phynames
                for line in `cat sites.txt`;
                      do gene=`echo "$line" | cut -f 1 -d','`
                      nsnp=`echo $line | cut -f 2 -d','`
                      if (( $nsnp > 0 )); then
                        prev=`tail -n 1 partitions.txt | cut  -f 2 -d'-'`
                        begin=$((prev+1))
                        end=$((begin + nsnp))
                        echo "DNA,  ${gene}=${begin}-${end}" >> partitions.txt
                        echo $gene >> phynames

                        cat $phy | tail -n +2 | awk '{print $2}' | cut -c${start}-${end} > seq.tmp
                        if [ -f phy.tmp ]; then
                                paste phy.tmp seq.tmp > phy2.tmp
                                mv phy2.tmp phy.tmp
                        else
                                cp seq.tmp phy.tmp
                        fi
                        rm -f seq.tmp
                      fi
                done
                cat $phy | tail -n +2 | cut -f 1 -d' ' > indnames
                paste indnames phy.tmp > phy2.tmp
                mv phy2.tmp phy.tmp

		printf "library(ape)\nlibrary(phrynomics)\nphynames<-scan('phynames',character())\nphy<-read.table('phy.tmp',row.names=1)\ncolnames(phy) <- phynames\nReadSNP(phy)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.noInv.phy',format='phylip')\nnsnps <- GetNumberOfSitesForLocus(snps)\nwrite(nsnps, file='nsnps_per_loc_noref', sep=' ', ncolumns=length(nsnps))\nwrite(names(nsnps), file='nsnps_locus_names_noref', sep=' ', ncolumns=length(nsnps))" > Rscript.R
              
                R --vanilla --no-save < Rscript.R
                
      		sites_noref=$(head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.noInv.phy | cut -f 2 -d' ')

        	if [[ -s RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out ]]
 			then
  		            rm -f RAxML_info.OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out
          		else
           		    echo "everything is good"
        	fi
        		
        	echo "creating partition for gene trees"
        	python ${REF_PATH}/snp_part_raxml.py -snpfile nsnps_per_loc_noref -namesfile nsnps_locus_names_noref > raxml_partitions.txt
       		      
        	echo "running raxml to create gene tree alignments"
        	raxmlHPC-PTHREADS-AVX -T 2 -f s -q raxml_partitions.txt -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy -m ASC_GTRGAMMA --asc-corr=lewis -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.NOREF.${int}.partitions.out

        	echo "running raxml on concatenated SNPs"
        	raxmlHPC-PTHREADS-AVX -T 16 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_noref}.NOREF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 323 -N 100 -p 476
                  
#                echo "running raxml to create gene trees"
		echo "not creating gene trees right now"  
#		export int
#		ls OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.NOREF.phy.locus*.phy | parallel --jobs 16 --delay 10 --ungroup --env int 'echo "starting raxml on {}" && raxmlHPC-PTHREADS-AVX -T 2 -s {} -n {.}.NOREF.${int}.filtered.out -j --no-bfgs --silent -m ASC_GTRCAT -V --asc-corr=lewis -f d -p 466'

#                if [ ! -d  astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.${int}.gene_tree_files/ ]; then
#			mkdir astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.${int}.gene_tree_files
#		fi
        		
#        	mv RAxML*locus*NOREF*.out astral/${tree_height}_sim${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.${int}.gene_tree_files/
                  
                echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.NOREF.phy" >> /project/phylogenref/scripts/output/${day}-SNPs-${tree_height}-sim${sim}-${int}

		wait
                  
		cat nsnps_locus_names_noref | tr ' ' '\n' > locnames.tmp
		cat nsnps_per_loc_noref | tr ' ' '\n' > locnsnps.tmp
		paste locnames.tmp locnsnps.tmp >> /project/phylogenref/scripts/output/${day}-SNPs-${tree_height}-sim${sim}-${int}
		rm -f locnames.tmp
		rm -f locnsnps.tmp
                  
                mkdir s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}.noref-${taxa_ref}.phylip_tree_files
        	  mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}.noref-${taxa_ref}.phylip_tree_files
        	  mv RAxML* s${sim}_q${QUAL}_miss${miss}_maf${maf}.${int}.noref-${taxa_ref}.phylip_tree_files

                done
         done

    done

rm -f *.bam

mkdir -p sim${sim}/config_files
mkdir -p sim${sim}/ref.fasta_files
		
mv *_config sim${sim}/config_files
mv *.fa sim${sim}/ref.fasta_files
				
###Create a batch file with all trees in order by file type, and create a name file
cat s${sim}*q*miss*/*bipartitions* >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}*q*miss*/*bipartitions* >> ${output_dir}/${day}-${tree_height}-tree.names

### Create a loop to run ASTRAL

date

echo "done with all processes for $tree_height sim${sim} ${int}; exiting now"
date
