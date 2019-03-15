#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --job-name=slurm.astral

file=$1
base=`basename $file`
dir=/gscratch/jrick/astral/$base
cd $dir

#for file in `ls /project/phylogenref/scripts/slurm_results/SLURM_2468*.tgz`;
#	base=`basename --suffix='_results.tgz' $file`;
#	tar -xzf $file;
#	mv $file phylogenref_tarred/;
#	cd $base;
	sim=`ls -d sim* | sed 's/sim//g'`;
	tree_height=`cat gene1_sim*/analysis_configuration.cfg | grep 'treefile' | grep -o -E '[0-9][0-9]+'`;
	taxa_ref=`cat gene1_sim*/analysis_configuration.cfg | grep 'base_genome_name' | grep -o -E '[0-9]+_0_0'`
	int=`ls -d s*phylip_tree_files | head -n 1 | grep -o -E '[A-Z][A-Z][A-Z]'`
	echo "starting astral concatenation for sim: $sim height: $tree_height taxa_ref: $taxa_ref int=$int";
	
#	/project/phylogenref/scripts/astral_script_031119.sh $sim $tree_height $taxa_ref $int
	ls astral/*.gene_tree_files/*.tre | sed 's!.*/!!' >> /project/phylogenref/scripts/output/031419-${tree_height}-astral.names
	cat astral/*.gene_tree_files/*.tre >> /project/phylogenref/scripts/output/031419-${tree_height}-astral.trees
#done

