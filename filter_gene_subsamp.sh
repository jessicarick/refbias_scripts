#!/bin/bash
gene=$1

echo "filtering for gene $gene, miss $miss, maf $maf"
vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.recode.vcf \
	--chr gene${gene} \
	--recode \
	--out gene${gene}_miss${miss}_maf${maf}.REF

python ${REF_PATH}/vcf2phylip.py \
	-i gene${gene}_miss${miss}_maf${maf}.REF.recode.vcf \
	-o gene${gene}_miss${miss}_maf${maf}.REF.phy  

sites=`head -n 1 gene${gene}_miss${miss}_maf${maf}.REF.phy | awk '{print $2}'`

if [ "$sites" -eq "0" ]; then 
	echo "no sites left"
	cat gene${gene}_miss${miss}_maf${maf}.REF.phy > gene${gene}_miss${miss}_maf${maf}.REF.noInv.phy
else 
	printf "library(ape)\nlibrary(phrynomics)\nReadSNP('gene${gene}_miss${miss}_maf${maf}.REF.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='gene${gene}_miss${miss}_maf${maf}.REF.noInv.phy',format='phylip')" > Rscript_${gene}_miss${miss}_maf${maf}.R
	
	R --vanilla --no-save < Rscript_${gene}_miss${miss}_maf${maf}.R
	rm -f Rscript_${gene}_miss${miss}_maf${maf}.R
fi 

#nsnps=`cat gene${gene}_miss${miss}_maf${maf}.REF.noInv.phy | head -n 1 | awk '{print $2}'`
#echo "gene${gene},s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-${tree_height}-${int}
