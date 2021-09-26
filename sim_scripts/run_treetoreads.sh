#!/bin/bash

source sim_scripts/refbias_config.txt

gene=$1

snps=`head -n 1 ${REF_PATH}/output/new/${day}-varSites-${tree_height}-sim${sim}-${int} | tr " " "\n" | head -n ${gene} | tail -n 1` 

python2 sim_scripts/write_config.py \
	-treefile ${REF_PATH}/sims_${tree_height}/sim${sim}/species_tree${sim}/1/g_trees${gene}.trees \
	-v `echo "$snps"` \
	-ref 0_0_0 \
	-path ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_${gene}.fa \
	-o gene${gene}_sim${sim} \
	-rate rat.matrix \
	-g 5 \
	-r 150 \
	-f 400 \
	-s 50 \
	-cl 5 \
	-ch 20 \
	-cs 1 \
	-p HS25 \
	-pre sim_ \
	-errorfile $error > gene${gene}_sim${sim}_config 

python2 ${PROGRAM_DIR}/TTR2/TreeToReads/treetoreads.py gene${gene}_sim${sim}_config 

#echo ">gene${gene}" >> ${reference_prefix}_sim${sim}.fa
#grep '^[ACGTN]' gene${gene}_sim${sim}/fasta_files/sim_${taxa_ref}.fasta >> ${reference_prefix}_sim${sim}.fa

#grep -v "^>" ${REF_PATH}/sims_${tree_height}/sim${sim}/${reference_prefix}.random_${gene}.fa > ${reference_prefix}_gene${gene}.fa
