#!/bin/bash

module load samtools

bam_in=$1

bam_in="path to basecalled bam file"
results=./results/get_pA_reads_R10

mkdir -p  $results


polyA_list=$results/reads_pA.txt
no_polyA_list=$results/no_polyA.list

mkdir -p  $results


samtools view $bam_in | awk -v OFS="\t"  '$1 !~ /^@/ {for (i=12; i<=NF; i++) if ($i ~ /^pt:i:/) print $1, $i}'  > $results/reads_pA_length.tsv
sed 's/\tpt:i:/\t/'   $results/reads_pA_length.tsv |  awk '$2 > 15 {print $1}' > $polyA_list
cat $polyA_list <(samtools view $bam_in | awk '{print $1}') | awk '{count[$0]++} END {for (line in count) if (count[line] == 1) print line}'  > $no_polyA_list



wait 
