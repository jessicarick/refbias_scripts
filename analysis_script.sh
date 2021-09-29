#!/bin/sh

##SBATCH --account=phylogenref
##SBATCH --nodes=1
##SBATCH --time=0-01:00:00
##SBATCH --mem=124G
##SBATCH --mail-type=all

module load gcc
module load r

day=$1

# concatenating raxml trees
if [ ! -f output/new/${day}-all-raxml.trees ]; then
	echo "concatenating raxml trees..."
	for height in SHORT MED LONG;
#	for height in `cat refbias_config.txt | grep 'tree_height_list' | sed 's/tree_height_list=\'//' | sed 's/\'//' `;
		do trees=output/new/${day}-${height}-batch.trees;
		names=output/new/${day}-${height}-tree.names;
		cat $trees >> output/new/${day}-all-raxml.trees;
		for name in `cat $names`;
			do echo "${height}-${name}" >> output/new/${day}-all-raxml.names;
		done;
		if [ -f output/new/${day}-${height}-subsamp-tree.trees ]; then
			echo "concatenating subsampled trees..."
			sub_trees=output/new/${day}-${height}-subsamp-tree.trees
			sub_names=output/new/${day}-${height}-subsamp-tree.names
			cat $sub_trees >> output/new/${day}-all-subsamp-raxml.trees
			for name in `cat $sub_names`
				do echo "${height}-${name}" >> output/new/${day}-all-subsamp-raxml.names
			done
		fi

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

#if false

exit
echo "beginning R analysis of raxml trees!"

#for i in `seq -w 00 09`
#do echo "working with batch ${i} of 09"
Rscript analysis/analyze_refbias_trees_raxml.R \
	--ml.trees ${day}-s_tree.tree \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-raxml.trees \
	--raxml.tree.names ${day}-all-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--refdist ${day}-refdist.txt \
	--output ${day}-output
#done

#Rscript analysis/analyze_refbias_trees_astral.R \
#	--ml.trees ${day}-s_tree.tree \
#	--ml.tree.names ${day}-s_tree.names \
#	--raxml.trees ${day}-all-raxml.trees \
#	--raxml.tree.names ${day}-all-raxml.names \
#	--astral.trees ${day}-all-astral.trees \
#	--astral.tree.names ${day}-all-astral.names \
#	--output ${day}-output 

#fi
