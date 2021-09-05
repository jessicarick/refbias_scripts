#!/bin/bash

source sim_scripts/refbias_config.txt
source activate $conda_env

gene=$1

echo "filtering for gene $gene, miss $miss, mac $mac"
vcftools --vcf OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.recode.vcf \
	--chr gene${gene} \
	--recode \
	--out gene${gene}_miss${miss}_mac${mac}.REF

python2 sim_scripts/vcf2phylip.py \
	-i gene${gene}_miss${miss}_mac${mac}.REF.recode.vcf \
	-o gene${gene}_miss${miss}_mac${mac}.REF.phy  

sites=`head -n 1 gene${gene}_miss${miss}_mac${mac}.REF.phy | awk '{print $2}'` 

if [ "$sites" -eq "0" ]; then 
	echo "no sites left"
	#cat gene${gene}_miss${miss}_mac${mac}.REF.phy > gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy
else 
	python sim_scripts/ascbias.py -p gene${gene}_miss${miss}_mac${mac}.REF.phy -o gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy
	rm -f gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy.felsenstein
	rm -f gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy.stamatakis
fi 

nsnps=`cat gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy | head -n 1 | awk '{print $2}'`

echo "gene${gene},s${sim}_q${QUAL}_miss${miss}_mac${mac}.${tree_height}.REF.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-${int}
