#!/bin/bash


arg=$1
echo $geno


data=../results/$geno/pipeline_map/
results=../results/$geno/demultiplex_read_names/
fastq=$data/trimming/trimmed.fastq.gz

mkdir -p $results


zcat $fastq |  sed -n 1~4p | awk '{print $1}' | awk -v FS="_" -v OFS="\t" '{print $1, $2, $3}' | sed 's/^@//'  > $results/reads_specifications.tsv

a=$(awk '{print $2}' $results/reads_specifications.tsv | sort | uniq)

echo -e "$a" | xargs -I {} --max-procs=13  bash -c  "awk -v bc={} '\$2==bc{print \$3}' $results/reads_specifications.tsv > $results/{}_reads.txt"



    


