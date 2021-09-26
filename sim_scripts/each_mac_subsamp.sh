#!/bin/bash
mac=$1
miss=$2

source sim_scripts/refbias_config_subsamp.txt

echo "working with mac $mac, miss $miss"
vcftools --gzvcf ${day}-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz \
         --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac} \
         --remove-filtered-all \
         --mac $(printf "$mac") \
         --max-missing $(printf "$miss") \
	 --min-alleles 2 \
	 --max-alleles 2 \
         --recode \
         --recode-INFO-all \
         --minDP 5 # keep the same?

if false; then
	#######################################
	#### CONVERTING VCF TO PHYLIP FILE ####
	#### and removing invariant sites #####
	#### making one phylip per gene #######
	#######################################
	export miss
	export mac
	export QUAL
	export sim
	export REF_PATH
	export output_dir
	export day
	export tree_height
	export int

	seq -w ${genes} | parallel --delay 5 --jobs 4 --env sim --env QUAL --env miss --env mac --env REF_PATH --env output_dir --env day --env tree_height --env int "bash ${REF_PATH}/filter_gene_subsamp.sh {}"

	## combine into one supermatrix

	Rscript ${REF_PATH}/make_supermat.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF $miss $mac
	#rm -f gene*.REF.noInv.phy
fi

module unload py-numpy
module unload py-scipy
python sim_scripts/vcf2phylip.py3 -i OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.recode.vcf -o OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.phy

python sim_scripts/ascbias.py -p OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.phy -o OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.phy
rm -f *.felsenstein
rm -f *.stamatakis

# change tab in phylip to spaces -- required for r script to work!
sed -i 's/\t/ /g' OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.phy

#######################################
#### Running RaxML w/ Ref #############
#######################################

    sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.phy | head -n 1 | awk '{print $2}'`
    if [ "$sites_ref" -gt "$maxSNP" ]; then
    for rep in `seq 10`
	do Rscript sim_scripts/subsample.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv
   	sites_samp=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.subsamp.phy | awk '{print $2}'`
	echo "OUTFILE_${tree_height}_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}_rep${rep},${sites_noref},${sites_samp}" >> /project/phylogenref/scripts/output/new/${day}-subsampSNPs-all

	echo "running raxml on subsampled phylip for rep $rep"
    	raxmlHPC-PTHREADS-AVX -T 4 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.subsamp.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466

	ls RAxML*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out >> ${REF_PATH}/output/new/${day}-${tree_height}-subsamp-tree.names
        cat RAxML*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out >> ${REF_PATH}/output/new/${day}-${tree_height}-subsamp-tree.trees
     done
     else
	rep=1
	Rscript sim_scripts/subsample.R OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv $sites_ref $maxSNP
        sites_samp=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.all.noInv.subsamp.phy | awk '{print $2}'`
        echo "OUTFILE_${tree_height}_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}_rep${rep},${sites_noref},${sites_samp}" >> /project/phylogenref/scripts/output/new/${day}-subsampSNPs-all

        echo "running raxml on subsampled phylip for rep $rep"
        raxmlHPC-PTHREADS-AVX -T 4 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.subsamp.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466

	ls RAxML*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out >> ${REF_PATH}/output/new/${day}-subsamp-tree.names
	cat RAxML*bipartitions.OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_sites${sites_samp}.REF.${int}.rep${rep}.subsamp.out >> ${REF_PATH}/output/new/${day}-subsamp-tree.trees
     fi
