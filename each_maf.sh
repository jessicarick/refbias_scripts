#!/bin/bash

source refbias_config.txt

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

seq -w ${genes} | parallel --delay 1 --jobs 4 --env sim --env QUAL --env miss --env mac --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/filter_gene.sh {}"

## combine into one supermatrix
Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF $miss $mac
#rm -f gene*.REF.noInv.phy

#######################################
#### Running RaxML w/ Ref #############
#######################################

    sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy | head -n 1 | awk '{print $2}'`

    echo "running raxml on concatenated SNPs"
    raxmlHPC-PTHREADS-AVX -T 4 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466

echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac},${sites_ref}" >> ${output_dir}/${day}-filteredSites-${tree_height}-${int}
	fi

# move phylogenies
if [ ! -d  s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files ]; then
	mkdir s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
fi

mv OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.phy s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
mv RAxML*OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files

####Add trees to batch file with all normal trees in order by file type, and add name to file
cat s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/RAxML*bipartitions.*OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-tree.names

