#!/bin/sh

day=$1

# concatenating raxml trees
if [ ! -f output/${day}-all-raxml.trees ]; then
	echo "concatenating raxml trees..."
	for height in `cat refbias_config.txt | grep 'tree_height_list' | sed 's/tree_height_list=\'//' | sed 's/\'//'`;
		do trees=output/${day}-${height}-batch.trees;
		names=output/${day}-${height}-tree.names;
		cat $trees >> output/${day}-all-raxml.trees;
		for name in `cat $names`;
			do echo "${height}-${name}" >> output/${day}-all-raxml.names;
		done;
	done
else
	echo "raxml trees already concatenated; moving on"
fi

# concatenating astral trees
if [ ! -f output/${day}-all-astral.trees ]; then
	echo "concatenating astral trees..."
	for height in `cat refbias_config.txt | grep 'tree_height_list' | sed 's/tree_height_list=\'//' | sed 's/\'//'`;
		do trees=output/${day}-${height}-astral.trees;
		names=output/${day}-${height}-astral.names;
		cat $trees >> output/${day}-all-astral.trees;
		cat $names >> output/${day}-all-astral.names;
	done 
else
	echo "astral trees already concatenated; moving on"
fi

echo "beginning R analysis of astral and raxml trees!"

Rscript analysis/analyze_refbias_trees_raxml.R \
	--ml.trees ${day}-s_tree.trees \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-raxml.trees \
	--raxml.tree.names ${day}-all-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--output ${day}-output 

Rscript analysis/analyze_refbias_trees_astral.R \
	--ml.trees ${day}-s_tree.trees \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-raxml.trees \
	--raxml.tree.names ${day}-all-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--output ${day}-output 