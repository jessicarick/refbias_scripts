i#!/bin/sh

module load gcc
module load miniconda2
module load samtools
module load perl

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

for tree_height in $tree_height_list
	do mkdir sims_$tree_height
	cd sims_$tree_height
	for sim in `seq $num_sims`
		do mkdir sim${sim}
		cd sim${sim}

##############################################
### Random subsample of full ref, no overlap##
##############################################

		python ${REF_PATH}/rand_subsample.py -s 1 -e $max -n $genes -l $gene_length > ${reference_prefix}.random.txt;

		header=`grep '^>' $reference | sed 's/>//g'`
		echo "header is $header"

		for i in `seq $genes`;
			do range=`sed -n "${i}p" ${reference_prefix}.random.txt`;
			samtools faidx $reference $header:$range >> ${reference_prefix}.random_${i}.tmp;
			echo ">${header}_${length}" > ${reference_prefix}.random_${i}.fa
			grep -v '^>' ${reference_prefix}.random_${i}.tmp | tr -d '\n' >> ${reference_prefix}.random_${i}.fa
			rm -f ${reference_prefix}.random_${i}.tmp
			grep -v '^>' ${reference_prefix}.random_${i}.fa | tr -d '\n' >> ${reference_prefix}.random_sim${sim}.fa
		done

##############################################
### For each gene, simulate MSC genealogy#####
##############################################
	rand_num=`echo $RANDOM`
	echo "random number seed: ${rand_num}"
	source activate new_env
	simphy_lnx64 -rl f:$genes -rg 1 -rs 1 -sl f:$num_sp -sb f:0.0000001 -si f:$num_ind -sp f:50000 -st f:$tree_height -so f:2 -cs $rand_num -o species_tree${sim}
		
	rename g_trees0 g_trees species_tree${sim}/1/g_trees0*
	cat species_tree${sim}/1/s_tree.trees >> ${output_dir}/${day}-${tree_height}-s_tree.tree
	echo "height${tree_height}_sim${sim}_s_tree.trees" >> ${output_dir}/${day}-${tree_height}-s_tree.names
	source deactivate
	perl /project/phylogenref/scripts/wrap_slurm_teton_refbias.pl $sim $tree_height
	
	cd ../
	done

cd ../
done


