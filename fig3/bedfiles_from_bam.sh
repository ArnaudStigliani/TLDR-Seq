#!/bin/bash
module load bedtools


sequencing_run=$1 # lib1 or lib2 or R10
data=../mapping_preprocessing/results/pipeline_map/$sequencing_run/mapping

results=./results/bedfiles_from_bam/

echo $results
mkdir -p $results

bam=$data/trimmed_primary.bam


bedtools bamtobed -i  $bam | awk -v FS="_" -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}'  |  \
    awk -v OFS="\t" '$4=="s1" {print $1,$2,$3,$4,$5,$6,$7,$8} $4=="s2" {print $1,$2,$3,$4,$5,$6,$7,($8=="+")?"-":"+"}' |\
    awk -v OFS="\t" '{print $1,$2,$3,$4"_"$5"_"$6,$7,$8}' \
	> $results/$sequencing_run.bed


