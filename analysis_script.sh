#!/bin/sh

day=$1

# concatenating raxml trees
echo "concatenating raxml trees..."
for height in 500000 2000000 10000000;
	do trees=output/${day}-${height}-batch.trees;
	names=output/${day}-${height}-tree.names;
	cat $trees >> output/${day}-all-raxml.trees;
	for name in `cat $names`;
		do echo "${height}-${name}" >> output/${day}-all-raxml.names;
	done;
done

# concatenating astral trees
echo "concatenating astral trees..."
for height in 500000 2000000 10000000;
	do trees=output/${day}-${height}-astral.trees;
	names=output/${day}-${height}-astral.names;
	cat $trees >> output/${day}-all-astral.trees;
	cat $names >> output/${day}-all-astral.names;
done 

echo "beginning R analysis of astral and raxml trees!"

Rscript analysis/analyze_refbias_trees.R \
	--ml.trees ${day}-s_tree.trees \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-raxml.trees \
	--raxml.tree.names ${day}-all-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--output ${day}-output \

