#!/bin/bash
export TPMDIR=/maps/projects/sandelin_main/scratch

module load samtools



bam_genome=../mapping_preprocessing/results/pipeline_map/mapping/trimmed_primary.bam
results=./results/get_processed_bam_files_R10/
mkdir -p  $results


echo $results


bam_header=$results/bam.header


samtools view -H    $bam_genome > $bam_header
samtools view   $bam_genome |   awk -v FS="\t" -v OFS="\t" '{gsub("_","\t",$1)} ($2==0 || $2==16) {print $0}'  |\
    awk  -v OFS="\t"  -v FS="\t"  '{a=($1=="s2"?-1:1); b=($4==16?-1:1); $4=(a*b==1?0:16)}  {print $0}'   |\
    awk -v dir=$results  -v OFS="\t" '(!($3 in max) || ($7 > max[$3])) {max[$3] = $7; fn[$3]=$2 ;$1=""; $2=""; line[$3] = substr($0,3)} END{ for (i in line) {file=fn[i]".sam"; print line[i] > dir"/"file}}'


find $results -name "*.sam" | xargs -I {} --max-procs=21 bash -c 'a=$(basename {} .sam) ; b=$(dirname {})/$a".bam"; cat "'$bam_header'" {} | samtools view -Sb | samtools sort | tee $b | samtools index - $b".bai"'







