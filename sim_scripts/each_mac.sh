#!/bin/bash

source sim_scripts/refbias_config.txt

mac=$1
miss=$2
sim=$3
int=$4

echo "working with mac $mac, miss $miss"
vcftools --vcf OUTFILE.s${sim}_q${QUAL}.vcf \
         --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac} \
         --remove-filtered-all \
         --mac $(printf "$mac") \
         --max-missing-count $(printf "$miss") \
	 	 --min-alleles 2 \
		 --max-alleles 2 \
         --recode \
         --recode-INFO-all \
         --minDP 5 # keep the same?

#######################################
#### CONVERTING VCF TO PHYLIP FILE ####
#### and removing invariant sites #####
#### making one phylip per gene #######
#######################################
export miss
export mac
export sim
export int

seq -w ${genes} | parallel --jobs 4 --env sim --env miss --env mac --env int "bash sim_scripts/filter_gene.sh {}"

## combine into one supermatrix
mkdir miss${miss}_mac${mac}_gene_phylips/
mv gene*miss${miss}_mac${mac}.REF.noInv.phy miss${miss}_mac${mac}_gene_phylips/
cd miss${miss}_mac${mac}_gene_phylips/
find . -type f -size -617w -name "gene*.noInv.phy" -delete #delete files with no sites
python ../sim_scripts/make_supermat.py ../OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy
cd ../

rm -rf gene*/

#######################################
#### Running RaxML w/ Ref #############
#######################################

sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy | head -n 1 | awk '{print $2}'`

echo "running raxml on concatenated SNPs"
raxmlHPC-PTHREADS-AVX -T 4 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy \
	-n OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out \
	-m ASC_GTRCAT -V \
	--asc-corr=lewis \
	-f a \
	-x 223 \
	-N 100 \
	-p 466

echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac},${sites_ref}" >> ${output_dir}/${day}-filteredSites-${tree_height}-${int}

# move phylogenies
if [ ! -d  s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files ]; then
	mkdir s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
fi

mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
mv RAxML*OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files

####Add trees to batch file with all normal trees in order by file type, and add name to file
cat s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/RAxML*bipartitions.*OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-tree.names

