#!/bin/bash


data=results_curated_genome
results=get_bedfiles

mkdir -p  $results

module load bedtools

bam1_curated=$data/nanopore_wt1_curated.bam
bam2_curated=$data/nanopore_wt2_curated.bam


bedtools bamtobed -i  $bam1_curated > $results/nanopore_wt1_curated_stranded.bed &
bedtools bamtobed -i  $bam2_curated > $results/nanopore_wt2_curated_stranded.bed
wait

