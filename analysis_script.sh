#!/bin/sh

#SBATCH --account=phylogenref
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --mail-type=all

module load gcc
module load r

day=$1

# concatenating raxml trees
if [ ! -f output/${day}-all-subsamp-raxml.trees ]; then
	echo "concatenating raxml trees..."
	for height in 500000 2000000 10000000;
#	for height in `cat refbias_config.txt | grep 'tree_height_list' | sed 's/tree_height_list=\'//' | sed 's/\'//' `;
		do trees=output/${day}-${height}-subsamp-batch.trees;
		names=output/${day}-${height}-subsamp-tree.names;
		cat $trees >> output/${day}-all-subsamp-raxml.trees;
		for name in `cat $names`;
			do echo "${height}-${name}" >> output/${day}-all-subsamp-raxml.names;
		done;
	done
else
	echo "raxml trees already concatenated; moving on"
fi

# concatenating astral trees
#if [ ! -f output/${day}-all-astral.trees ]; then
#	echo "concatenating astral trees..."
#	for height in `cat refbias_config.txt | grep 'tree_height_list' | sed 's/tree_height_list=\'//' | sed 's/\'//'`;
#		do trees=output/${day}-${height}-astral.trees;
#		names=output/${day}-${height}-astral.names;
#		cat $trees >> output/${day}-all-astral.trees;
#		cat $names >> output/${day}-all-astral.names;
#	done 
#else
#	echo "astral trees already concatenated; moving on"
#fi

echo "beginning R analysis of astral and raxml trees!"

Rscript analysis/analyze_refbias_trees_raxml.R \
	--ml.trees ${day}-s_tree.tree \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-subsamp-raxml.trees \
	--raxml.tree.names ${day}-all-subsamp-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--refdist ${day}-refdist.txt \
	--output ${day}-output 

#Rscript analysis/analyze_refbias_trees_astral.R \
#	--ml.trees ${day}-s_tree.tree \
#	--ml.tree.names ${day}-s_tree.names \
#	--raxml.trees ${day}-all-raxml.trees \
#	--raxml.tree.names ${day}-all-raxml.names \
#	--astral.trees ${day}-all-astral.trees \
#	--astral.tree.names ${day}-all-astral.names \
#	--output ${day}-output 
