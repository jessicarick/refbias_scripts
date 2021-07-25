#!/bin/bash
maf=$1
miss=$2

echo "working with maf $maf, miss $miss"
vcftools --vcf OUTFILE.s${sim}_q${QUAL}.vcf \
         --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf} \
         --remove-filtered-all \
         --maf $(printf "$maf") \
         --max-missing $(printf "$miss") \
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
	export maf
	export QUAL
	export sim
	export REF_PATH
	export output_dir
	export day
	export tree_height
	export int

	seq -w ${genes} | parallel --delay 5 --jobs 8 --env sim --env QUAL --env miss --env maf --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/filter_gene.sh {}"

	## combine into one supermatrix

	Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF $miss $maf
	#rm -f gene*.REF.noInv.phy

#######################################
#### Running RaxML w/ Ref #############
#######################################

    sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.phy | head -n 1 | awk '{print $2}'`

    echo "running raxml on concatenated SNPs"
    raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_ref}.REF.${int}.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466

 echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf},${sites_ref}" >> ${output_dir}/${day}-filteredSites-${tree_height}-${int}
