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

	seq -w ${genes} | parallel --delay 5 --jobs 8 --env sim --env QUAL --env miss --env maf --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/filter_gene_subsamp.sh {}"

	## combine into one supermatrix

	Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF $miss $maf
	#rm -f gene*.REF.noInv.phy

#######################################
#### Running RaxML w/ Ref #############
#######################################

    sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.phy | head -n 1 | awk '{print $2}'`
    if [ "$sites_ref" -eq "$maxSNP" ]; then
    for rep in `seq 10`
	do Rscript ${REF_PATH}/subsample.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv $sites_ref $maxSNP
   	sites_samp=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy | awk '{print $2}'`
	echo "OUTFILE_${tree_height}_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}_rep${rep},${sites_noref},${sites_samp}" >> /project/phylogenref/scripts/output/new/${day}-subsampSNPs-all

	echo "running raxml on subsampled phylip for rep $rep"
    	raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466
     done
     else
	rep=1
	Rscript ${REF_PATH}/subsample.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv $sites_ref $maxSNP
        sites_samp=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy | awk '{print $2}'`
        echo "OUTFILE_${tree_height}_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}_rep${rep},${sites_noref},${sites_samp}" >> /project/phylogenref/scripts/output/new/${day}-subsampSNPs-all

        echo "running raxml on subsampled phylip for rep $rep"
        raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.all.noInv.subsamp.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466
     fi
