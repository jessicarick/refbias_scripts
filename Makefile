day=$1

manuscript.pdf: 

analysis: ./analysis/univariate_plots_refbias.R ./analysis/model_refbias_trees.R ./analysis/analyze_refbias_trees.R output/${day}-*

univariate_plots_refbias: analysis/model_refbias_trees.R analysis/analyze_refbias_trees.R output/${day}-output.csv
	Rscript analysis/model_refbias_trees.R \
	--results ${day}-output.csv

model_refbias_trees: analysis/analyze_refbias_trees.R output/${day}-output.csv
	Rscript analysis/analyze_refbias_trees.R \
	--results ${day}-output.csv

analyze_refbias_trees: ./output/${day}-*
	Rscript analysis/analyze_refbias_trees.R \
	--ml.trees ${day}-s_tree.trees \
	--ml.tree.names ${day}-s_tree.names \
	--raxml.trees ${day}-all-raxml.trees \
	--raxml.tree.names ${day}-all-raxml.names \
	--astral.trees ${day}-all-astral.trees \
	--astral.tree.names ${day}-all-astral.names \
	--output ${day}-output

simulations: refbias_config.txt refbias_simphy.sh wrap_slurm_teton_refbias.pl test_simphy_ExtREFs_031119.sh test_simphy_IntREFs_031119.sh vcf2phylip.py astral_script_031119.sh latesGBS.txt 
	source refbias_simphy.sh