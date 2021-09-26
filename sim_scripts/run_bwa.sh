#!/bin/bash

source sim_scripts/refbias_config.txt

ind=$1

echo "Mapping reads for ${ind}"

bwa mem ${reference_prefix}_sim${sim}_${int}.fa fastq_reads/sim${sim}/fastq/${ind}/${ind}_${sim}_read1.fq.gz fastq_reads/sim${sim}/fastq/${ind}/${ind}_${sim}_read2.fq.gz > aln_${ind}.sam

echo "Converting sam to bam for ${ind}"
samtools view -b -S -o aln_${ind}.bam aln_${ind}.sam

echo "Sorting and indexing bam files for ${ind}"
samtools sort aln_${ind}.bam -o aln_${ind}.sorted.bam
samtools index -c aln_${ind}.sorted.bam
