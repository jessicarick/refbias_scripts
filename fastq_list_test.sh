#!/bin/sh

#SBATCH --job-name Ref_Bias_Sims

#SBATCH --account=phylogenref

#SBATCH -o stdoutRef_file
#SBATCH -e stderrRef_file

#SBATCH --mail-type=END
#SBATCH --mail-user=cbrock2@uwyo.edu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --time=0-00:00:10

length=500
fastq_list=0

while [ $(echo "$(expr length "$fastq_list")") -lt 2 ]
do
	i=$(shuf -i 1-$length -n 1)
	if [ ! -f gene1_sim1/fastq/sim*/sim_0_0_0*_1.fq.gz ]
		then
 			echo "File empty. try again"
		else	
			echo "File present, well done"
			declare fastq_list=`ls gene1_sim1/fastq/sim*/*_1.fq.gz | xargs -n 1 basename | sed 's/_1.fq.gz//'`
			##export "$fastq_list"
	fi
	##export $fastq_list
	##echo $fastq_list
done 

for fastq in $fastq_list
	do echo $fastq
done
