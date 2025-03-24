#!/bin/bash


sequencing_run=$1
data=../../results/TLDRSeq/pipeline_map/$sequencing_run

results=../../results/TLDRSeq/demultiplex_read_names/$sequencing_run
echo $results
mkdir -p $results 
fastq=$data/trimmed.fastq.gz

zcat $fastq |  sed -n 1~4p | awk '{print $1}' | awk -v FS="_" -v OFS="\t" '{print $1, $2, $3}' | sed 's/^@//'  > $results/reads_specifications.tsv

a=$(awk '{print $2}' $results/reads_specifications.tsv | sort | uniq)

echo -e "$a" | xargs -I {} --max-procs=13  bash -c  "awk -v bc={} '\$2==bc{print \$3}' $results/reads_specifications.tsv > $results/{}_reads.txt"



    


