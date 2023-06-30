#!/bin/bash

source sim_scripts/refbias_config.txt
module unload py-numpy
module unload py-scipy
source activate $conda_env

mac=$1
miss=$2

vcf=${output_dir}/100922-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz

#sim=`echo $vcf | cut -f 2 -d'_' | sed 's/s//g'`
#int=`echo $vcf | cut -f 2 -d'.'`
#height=`echo $vcf | cut -f 2 -d'-'`

taxa_ref=`cat ${output_dir}/100922-all-raxml.names | grep $tree_height | grep $int | grep s$sim_ | grep miss0_mac0 | cut -f 3 -d'-' | cut -f 1 -d'.'`


#if false; then #debugging

echo "working with mac $mac, miss $miss"
vcftools --gzvcf ${output_dir}/100922-${tree_height}-OUTFILE_s${sim}_q${QUAL}_rawvar.${int}.vcf.gz \
         --out OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac} \
         --remove-filtered-all \
         --mac $mac \
         --max-missing $miss \
 	 --min-alleles 2 \
	 --max-alleles 2 \
         --recode \
         --recode-INFO-all \
         --minDP 5 # keep the same?

if [ ! -f OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}.recode.vcf ]; then
	echo "filtered vcf not created; exiting now"
	exit 1
fi

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
mkdir miss${miss}_mac${mac}_${int}_gene_phylips/
mv gene*miss${miss}_mac${mac}.REF.noInv.phy miss${miss}_mac${mac}_${int}_gene_phylips/
cd miss${miss}_mac${mac}_${int}_gene_phylips/

for file in gene*.noInv.phy
	do sites=`head -n 1 $file | cut -f 2 -d' '`
	if [ "$sites" -eq "0" ]; then
		rm -f $file
		echo "removing $file with no sites"
	fi
done

python ../sim_scripts/make_supermat.py ../OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}.REF.all.noInv.phy

###########################################################
## RUN ASTRAL ON GENE TREES AND CALCUALTE GT DISCORDANCE ##
###########################################################

#if false; then # skip astral

seq -w $genes | parallel --jobs 2 --env miss --env mac --delay 1 "echo {} && raxmlHPC-PTHREADS -T 2 -s gene{}_miss${miss}_mac${mac}.REF.noInv.phy -n gene{}_miss${miss}_mac${mac}.out -m ASC_GTRCAT -V --asc-corr=lewis -p 123"

cat RAxML_bestTree*gene*miss${miss}_mac${mac}.out >> genetrees_miss${miss}_mac${mac}.tre
rm -f RAxML_*gene*miss${miss}_mac${mac}*
rm -f gene*miss${miss}_mac${mac}.*.phy

java -jar ${PROGRAM_DIR}/ASTRAL/astral.5.6.1.jar -i genetrees_miss${miss}_mac${mac}.tre -o astral_sim${sim}_miss${miss}_mac${mac}.tre

cat astral_sim${sim}_miss${miss}_mac${mac}.tre >> ${output_dir}/${day}-${tree_height}-ASTRAL-batch.trees
echo "sim${sim}_miss${miss}_mac${mac}_${int}" >> ${output_dir}/${day}-${tree_height}-ASTRAL-tree.names

Rscript ../sim_scripts/calc_genetree_rf.R genetrees_miss${miss}_mac${mac}.tre $sim >> ${output_dir}/${day}-${tree_height}-${int}-post.gt_rf

#fi # skip astral
###########

cd ../

#rm -rf gene*/
#fi

#######################################
#### Running RaxML w/ Ref #############
#######################################
if false; then # skip RAxML
sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}.REF.all.noInv.phy | head -n 1 | awk '{print $2}'`

echo "running raxml on concatenated SNPs"
raxmlHPC-PTHREADS -T 4 -s OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}.REF.all.noInv.phy \
	-n OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out \
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

mv OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}.REF.all.noInv.phy s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files
mv RAxML*OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files

####Add trees to batch file with all normal trees in order by file type, and add name to file
cat s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/RAxML*bipartitions.*OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-batch.trees
ls s${sim}_q${QUAL}_miss${miss}_mac${mac}.${int}-${taxa_ref}.phylip_tree_files/*bipartitions.OUTFILE_s${sim}_q${QUAL}_${int}_miss${miss}_mac${mac}_sites${sites_ref}.REF.${int}.filtered.out >> ${output_dir}/${day}-${tree_height}-tree.names

fi # skip RAxML
