#!/bin/bash

source refbias_config.txt

gene=$1

echo "filtering for gene $gene, miss $miss, mac $mac"
vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.recode.vcf \
	--chr gene${gene} \
	--recode \
	--out gene${gene}_miss${miss}_mac${mac}.REF

python2 ${REF_PATH}/vcf2phylip.py \
	-i gene${gene}_miss${miss}_mac${mac}.REF.recode.vcf \
	-o gene${gene}_miss${miss}_mac${mac}.REF.phy  

#sites=`head -n 1 gene${gene}_miss${miss}_mac${mac}.REF.phy | awk '{print $2}'` # doesn't work for some reason

if [ "$sites" -eq "0" ]; then 
	echo "no sites left"
	cat gene${gene}_miss${miss}_mac${mac}.REF.phy > gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy
else 
	#printf "library(ape)\nlibrary(phrynomics)\nReadSNP('gene${gene}_miss${miss}_mac${mac}.REF.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy',format='phylip')" > Rscript_${gene}_miss${miss}_mac${mac}.R
	
	#R --vanilla --no-save < Rscript_${gene}_miss${miss}_mac${mac}.R
	#rm -f Rscript_${gene}_miss${miss}_mac${mac}.R

	python ${REF_PATH}/ascbias.py -p gene${gene}_miss${miss}_mac${mac}.REF.phy -o gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy
	rm -f gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy.felsenstein
	rm -f gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy.stamatakis
fi 

nsnps=`cat gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy | head -n 1 | awk '{print $2}'`

echo "gene${gene},s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-${tree_height}-${int}
