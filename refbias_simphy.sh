#!/bin/sh

module load gcc
module load miniconda2/4.3.30
module load samtools
module load perl

source activate new_env

PATH=$PATH:/project/phylogenref/programs/art_bin_GreatSmokyMountains:/project/phylogenref/programs/TreeToReads:/project/phylogenref/programs/ASTRAL:/project/phylogenref/programs/SimPhy_1.0.2/bin:/project/phylogenref/programs/Seq-Gen-1.3.4/source

source refbias_config.txt
max=100000000

echo "reference: $reference"
echo "ref_length: $ref_length"
echo "gene_length: $gene_length"
echo "number of genes: $genes"
echo "max length: $max"

##############################################
#### CREATING FOLDERS FOR SIMULATIONS ########
##############################################

for ils_level in $ils_level_list
	do mkdir sims_${ils_level}
	cd sims_${ils_level}
	for sim in $num_sims
		do mkdir sim${sim}
		cd sim${sim}
		
##############################################
### Random subsample of full ref, no overlap##
##############################################

#		reference_prefix=lmariae_${tree_height}_${sim}_

		python ${REF_PATH}/rand_subsample.py -s 1 -e $max -n $genes -l $gene_length > ${reference_prefix}.random.txt;

		header=`grep '^>' $reference | sed 's/>//g'`
		echo "header is $header"
		#echo ">${header}_${ref_length}" > ${reference_prefix}.random_sim${sim}.fa

		for i in `seq -w $genes`;
			do range=`sed -n "${i}p" ${reference_prefix}.random.txt`;
			samtools faidx $reference $header:$range >> ${reference_prefix}.random_${i}.tmp;
			echo ">${header}_${ref_length}" > ${reference_prefix}.random_${i}.fa
			grep -v '^>' ${reference_prefix}.random_${i}.tmp | tr -d '\n' >> ${reference_prefix}.random_${i}.fa
			rm -f ${reference_prefix}.random_${i}.tmp
			echo ">gene${i}" >> ${reference_prefix}.random_sim${sim}.fa
			grep -v '^>' ${reference_prefix}.random_${i}.fa | tr -d '\n' >> ${reference_prefix}.random_sim${sim}.fa
			echo >> ${reference_prefix}.random_sim${sim}.fa
			
			###Counts Ns and stores as variable counts if needed

			counts=`awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/N/,"") "\t" NR}' ${reference_prefix}.random_${i}.fa | awk 'FNR == 3 {print $1}'`
			echo "replacing $counts N's with random AGCTs"

			#####Draws random values from ACGT of length counts
			res=`head -n $counts /dev/urandom | tr -dc ACGT | head -n $counts ; echo ''`
			
			k=0
			while [ $k -le $counts ]
				do sed -i "s/N/${res:$k:1}/" ${reference_prefix}.random_${i}.fa
				k=$(($k + 1))
			done

#			bases="ACGT"
#			randbase=${bases:$(( RANDOM % ${#bases} )):1}
#			sed -i "s/N/$randbase/g" ${reference_prefix}.random_${i}.fa
#			sed -i "s/N/$randbase/g" ${reference_prefix}.random_sim${sim}.fa
		done

##############################################
### For each gene, simulate MSC genealogy#####
##############################################
	rand_num=`echo $RANDOM`
	echo "random number seed: ${rand_num}"
	
	if [ "${tree_height}" == "NH" ]; then
	        echo "simulation without specifying tree height"
		simphy_lnx64 -rl f:$genes \
        	        -rg 1 \
                	-rs 1 \
	                -sl f:$num_sp \
        	        -sb f:$sp_rate \
                	-si f:$num_ind \
	                -sp f:50000 \
        	        -so f:10 \
                	-cs $rand_num \
        	        -o species_tree${sim}
		tree_height=${ils_level}
	else
		simphy_lnx64 -rl f:$genes \
			-rg 1 \
			-rs 1 \
			-sl f:$num_sp \
			-sb f:$sp_rate \
			-si f:$num_ind \
			-sp f:50000 \
			-so f:10 \
			-st f:$tree_height \
			-cs $rand_num \
			-o species_tree${sim}
	fi
		
#	for file in species_tree${sim}/1/g_trees0*; do cp $file `echo $file | sed 's/0//'`; done
#	rename g_trees0 g_trees species_tree${sim}/1/g_trees0*
	cat species_tree${sim}/1/s_tree.trees >> ${output_dir}/${day}-s_tree.tree
	echo "height${ils_level}_sim${sim}_s_tree.trees" >> ${output_dir}/${day}-s_tree.names

	conda deactivate

	perl /project/phylogenref/scripts/wrap_slurm_teton_refbias.pl $sim $ils_level
	
	cd ../
	done

cd ../
done

echo "sim ils_level int taxa_ref avg_dist" >> ${output_dir}/${day}-refdist.txt
#echo "gene,num_SNPs,num_noRef,num_nonInv" >> ${output_dir}/${day}-mutations.txt

