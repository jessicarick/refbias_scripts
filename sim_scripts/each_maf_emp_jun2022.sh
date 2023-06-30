#!/bin/bash
source /project/phylogenref/scripts/sim_scripts/refbias_config_emp.txt

mac=$1
miss=$2

echo "working with mac $mac, miss $miss"

mkdir /lscratch/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}/sim${sim}_miss${miss}_mac${mac}_${int}
cd /lscratch/${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}/sim${sim}_miss${miss}_mac${mac}_${int}

#if false; then 
vcftools --gzvcf /gscratch/jrick/phylogenref/emp_tmp/${tree_height}/reheader.OUTFILE.q${QUAL}_s${sim}.${int}.vcf.gz \
         --out OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int} \
	 --remove-filtered-all \
         --mac $(printf "$mac") \
         --max-missing $(printf "$miss") \
         --recode \
	 --recode-INFO-all \
         --minDP 5 # keep the same?

#######################################
#### Running ASTRAL ###################
#######################################
echo "creating gene trees for ASTRAL analysis"
mkdir partition_files_s${sim}_miss${miss}_mac${mac}/

vcf=OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.recode.vcf
nsnp=100000

echo "snp interval is ${nsnp}bp"

cat /gscratch/jrick/phylogenref/emp_tmp/${tree_height}/chrom_${int} > chr

for chrom in `cat chr`; do
        echo "working on chrom $chrom"
        #if [ ${vcf: -7} == ".vcf.gz" ]; then
        #        vcftools --gvcf $vcf --chr $chrom --recode --out tmp
        #else
        #        vcftools --vcf $vcf --chr $chrom --recode --out tmp
        #fi
        #grep '^#' $vcf > header
        #grep -v '^#' tmp.${chrom}.recode.vcf | split -l $nsnp --additional-suffix=.vcf -d - ${chrom}_ 

        #first=$(grep -v '^#' $vcf | grep $chrom | head -n 1 | cut -f 2)
        last=$(grep -P "^${chrom}\t" $vcf | tail -n 1 | cut -f 2)
	echo "max bp for $chrom is $last"

        i=1
        j=1
        while [[ $i -le $last ]]; do
        	start=$i
                end=$((i+99999))

                vcftools --vcf $vcf --chr `echo $chrom` --from-bp $start --to-bp $end --recode --out partition_files_s${sim}_miss${miss}_mac${mac}/${chrom}_${j}

		numvar=`grep -v -c '^#' partition_files_s${sim}_miss${miss}_mac${mac}/${chrom}_${j}.recode.vcf`
		if [[ $numvar -eq 0 ]]; then
			echo "no variants in window"
			rm -f partition_files_s${sim}_miss${miss}_mac${mac}/${chrom}_${j}.recode.vcf
		fi
		
		j=$((j+1))
		i=$((i+100000))
        done
done

cd partition_files_s${sim}_miss${miss}_mac${mac}/

ls *.recode.vcf | parallel -j 8 "echo 'working with {}' && python /project/phylogenref/scripts/sim_scripts/vcf2phylip.py3 -i {} -o {.}.phy && python /project/phylogenref/scripts/sim_scripts/ascbias.py -p {.}.phy -o {.}.noInv.phy"

for file in *.noInv.phy
        do sites=`head -n 1 $file | cut -f 2 -d' '`
        if [[ "$sites" -le "2" ]]; then
                rm -f $file
                echo "removing $file with fewer than 2 sites"
        fi
done

ls *.noInv.phy | parallel -j 8 "bash /project/phylogenref/scripts/sim_scripts/remove_missing_phy.sh {} && raxmlHPC-PTHREADS -s {.}.reduced.phy -m ASC_GTRCAT -n {.} -T 1 -V --asc-corr=lewis -p 12345"

cat RAxML_bestTree* >> ../genetrees_s${sim}_miss${miss}_mac${mac}.tre

cd ../
rm -rf partition_files_s${sim}_miss${miss}_mac${mac}/

echo "done creating gene trees, now running ASTRAL"
java -jar ${PROGRAM_DIR}/ASTRAL/astral.5.6.1.jar -i genetrees_s${sim}_miss${miss}_mac${mac}.tre -o astral_sim${sim}_miss${miss}_mac${mac} && cat astral_sim${sim}_miss${miss}_mac${mac} >> ${output_dir}/${day}-${tree_height}-emp-ASTRAL-batch.trees && echo "sim${sim}_miss${miss}_mac${mac}_${int}" >> ${output_dir}/${day}-${tree_height}-emp-ASTRAL-tree.names

Rscript ${REF_PATH}/sim_scripts/calc_genetree_rf.R genetrees_s${sim}_miss${miss}_mac${mac}.tre $sim >> ${output_dir}/${day}-${tree_height}-${int}-post.gt_rf
echo "done with ASTRAL run"

#######################################
#### Running RaxML w/ Ref #############
#######################################
#echo "not running genome-wide raxml inference"
#if false; then  # don't run raxml
echo "beginning concatenated RAxML run"
python2 ${REF_PATH}/sim_scripts/vcf2phylip.py \
	-i OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.recode.vcf \
	-o OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy  
#fi

if [[ -s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy ]]; then
	rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.recode.vcf
else
	echo "phylip file doesn't exist or is empty; exiting now"
	exit 2
fi

sites=`head -n 1 OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy | awk '{print $2}'`

#printf "library(ape)\nlibrary(phrynomics)\nReadSNP('OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.phy',fileFormat='phy',extralinestoskip=1)->fullSNPs\nRemoveInvariantSites(fullSNPs, chatty=TRUE)->fullSNPs_only\nsnps <- RemoveNonBinary(fullSNPs_only, chatty=TRUE)\nWriteSNP(snps, file='OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.phy',format='phylip')" > Rscript_miss${miss}_mac${mac}.R
	
#R --vanilla --no-save < Rscript_miss${miss}_mac${mac}.R
#rm -f Rscript_miss${miss}_mac${mac}.R 

python ${REF_PATH}/sim_scripts/ascbias.py \
	-p OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy \
	-o OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.noInv.phy

rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy.felsenstein
rm -f OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.phy.stamatakis

nsnps=$(head -n 1  OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.noInv.phy | awk '{print $2}')
echo "s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.${tree_height}.${int}.noInv,${nsnps}" >> ${output_dir}/${day}-SNPs-emp

#sites_ref=`cat OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.noInv.phy | head -n 1 | awk '{print $2}'`

echo "running raxml on concatenated SNPs"
raxmlHPC-PTHREADS -T 8 \
	-s OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}_${int}.noInv.phy \
	-n OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.${int}.emp.filtered.out \
	-m ASC_GTRCAT -V \
	--asc-corr=lewis \
	-f a \
	-x 223 \
	-N 100 \
	-p 466

####Write number of sites to output file
#echo "OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac},${sites_ref}" >> ${output_dir}/${day}-filteredSites-emp-${tree_height}-${int}

####Add tree and name to output files
cat RAxML*bipartitions.*s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.${int}.emp.filtered.out >> ${output_dir}/${day}-${tree_height}-emp-batch.trees
ls RAxML*bipartitions.*s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.${int}.emp.filtered.out >> ${output_dir}/${day}-${tree_height}-emp-tree.names

####Move RAxML results to a single directory
if [ ! -d s${sim}_q${QUAL}.${int}-lates-emp.phylip_tree_files ]; then
	mkdir s${sim}_q${QUAL}.${int}-${tree_height}-emp.phylip_tree_files
fi

mv *OUTFILE_s${sim}_q${QUAL}_miss${miss}_mac${mac}*.phy s${sim}_q${QUAL}.${int}-${tree_height}-emp.phylip_tree_files
mv RAxML*s${sim}_q${QUAL}_miss${miss}_mac${mac}.REF.${int}.emp.filtered.out s${sim}_q${QUAL}.${int}-${tree_height}-emp.phylip_tree_files/

#fi # don't run raxml
