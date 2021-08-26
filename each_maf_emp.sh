#!/bin/bash
maf=$1
miss=$2

echo "working with maf $maf, miss $miss"

#if false; then 
vcftools --vcf OUTFILE.q${QUAL}_s${sim}.${int}.vcf \
         --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int} \
         --remove-filtered-all \
         --maf $(printf "$maf") \
         --max-missing $(printf "$miss") \
         --recode \
         --recode-INFO-all \
         --minDP 5 # keep the same?


#######################################
#### Running RaxML w/ Ref #############
#######################################
python2 ${REF_PATH}/vcf2phylip.py \
	-i OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.recode.vcf \
	-o OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy  
#fi

if [[ -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy ]]; then
	rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.recode.vcf
else
	echo "phylip file doesn't exist or is empty; exiting now"
	exit 2
fi

sites=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy | awk '{print $2}'`

#printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy',format='phylip')" > Rscript_miss${miss}_maf${maf}.R
	
#R --vanilla --no-save < Rscript_miss${miss}_maf${maf}.R
#rm -f Rscript_miss${miss}_maf${maf}.R 

python /home/jrick/bin/raxml_ascbias/ascbias.py -p OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy -o OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.noInv.phy

rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy.felsenstein
rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.phy.stamatakis

nsnps=$(head -n 1  OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.noInv.phy | awk '{print $2}')
echo "s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${tree_height}.${int}.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-emp

#sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.noInv.phy | head -n 1 | awk '{print $2}'`

echo "running raxml on concatenated SNPs"
raxmlHPC-PTHREADS-AVX -T 8 -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}_${int}.noInv.phy -n OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.emp.filtered.out -j -m ASC_GTRGAMMA --asc-corr=lewis -f a -x 223 -N 100 -p 466

####Write number of sites to output file
#echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf},${sites_ref}" >> ${output_dir}/${day}-filteredSites-emp-${tree_height}-${int}

####Add tree and name to output files
cat RAxML*bipartitions.*s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.emp.filtered.out >> ${output_dir}/${day}-${tree_height}-emp-batch.trees
ls RAxML*bipartitions.*s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.emp.filtered.out >> ${output_dir}/${day}-${tree_height}-emp-tree.names

####Move RAxML results to a single directory
if [ ! -d s${sim}_q${QUAL}.${int}-lates-emp.phylip_tree_files ]; then
	mkdir s${sim}_q${QUAL}.${int}-lates-emp.phylip_tree_files
fi
mv *OUTFILE_s${sim}_q${QUAL}_miss${miss}_maf${maf}*.phy s${sim}_q${QUAL}.${int}-lates-emp.phylip_tree_files
mv RAxML*s${sim}_q${QUAL}_miss${miss}_maf${maf}.REF.${int}.emp.filtered.out s${sim}_q${QUAL}.${int}-lates-emp.phylip_tree_files

