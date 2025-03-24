#!/bin/bash

module load samtools
module load minimap2

results_mapped=result_minimap2
results=results_curated
mkdir -p $results_mapped $results

genome_index="path of GRCh38.mmi" # path of human genome index
genome_bedgff="path of anno.bed" # path of splice junctions generated with the paftools

minimap2 -ax splice -un --junc-bed  $genome_bedgff -t 50  $genome_index ERR5880581.fastq.gz   | samtools view  -S -b  > $results_mapped/nanopore_wt1.bam  2> $results_mapped/nanopore_wt1.log 
minimap2 -ax splice -un --junc-bed  $genome_bedgff -t 50  $genome_index ERR5880582.fastq.gz   | samtools view  -S -b  > $results_mapped/nanopore_wt2.bam  2> $results_mapped/nanopore_wt2.log 


bam1=$results_mapped/nanopore_wt1.bam
bam1_H=$results_mapped/nanopore_wt1.headers
bam1_primmary=$results_mapped/nanopore_wt1_primary.bam

bam2=$results_mapped/nanopore_wt2.bam
bam2_H=$results_mapped/nanopore_wt2.headers
bam2_primmary=$results_mapped/nanopore_wt2_primary.bam


samtools view -H $bam1 > $bam1_H
samtools view $bam1  -F 256 -F 4 | awk -v OFS="\t" '{gsub("AS:i:","",$14)}{print $0}' | sort -k1,1 -k14,14nr | sort -k1,1 -u | awk -v OFS="\t" '{gsub("^","AS:i:",$14)} {print $0}' | cat $bam1_H  - | samtools view -S -b   | samtools sort -@ 45 | tee $bam1_primmary  | samtools index  - $bam1_primmary.bai  

samtools view -H $bam2 > $bam2_H
samtools view $bam2  -F 256 -F 4 | awk -v OFS="\t" '{gsub("AS:i:","",$14)}{print $0}' | sort -k1,1 -k14,14nr | sort -k1,1 -u | awk -v OFS="\t" '{gsub("^","AS:i:",$14)} {print $0}' | cat $bam2_H  - | samtools view -S -b   | samtools sort -@ 45 | tee $bam2_primmary  | samtools index  - $bam2_primmary.bai  


bam1_curated=$results/nanopore_wt1_curated.bam
bam2_curated=$results/nanopore_wt2_curated.bam

samtools view  $bam1 |  awk -v FS="\t" -v OFS="\t" '($2==0 || $2==16 || $1~"^@") {print $0}'  |\
    sort -k1,1 -k5,5nr |\
    sort -k1,1 -u |\
    cat $bam1_H - |\
    samtools view -S -b   - |\
    samtools sort -@100  >    $bam1_curated
samtools index $bam1_curated $bam1_curated.bai

samtools view  $bam2 |  awk -v FS="\t" -v OFS="\t" '($2==0 || $2==16 || $1~"^@") {print $0}'  |\
    sort -k1,1 -k5,5nr |\
    sort -k1,1 -u |\
    cat $bam2_H - |\
    samtools view -S -b   - |\
    samtools sort -@100  >    $bam2_curated
samtools index $bam2_curated $bam2_curated.bai
