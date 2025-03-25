#!/bin/bash


module load samtools

geno=$1

data=./results/$geno/pipeline_map/mapping
results=./results/$geno/get_primary_processed_bam
echo $results
mkdir -p $results 

bam_genome=$data/trimmed_primary.bam

bam_polyA_processed=$results/trimmed_primary_processed.bam
bam_header=$results/bam.header



polyA_list=./results/$geno/characterize_polyA/reads_with_polyA.list


bam_polyA_processed=$results/trimmed_primary_processed.bam


bam_header=$results/bam.header

samtools view -H    $bam_genome > $bam_header

samtools view  $bam_genome |  awk -v FS="\t" -v OFS="\t" '{gsub("_","\t",$1)} ($2==0 || $2==16) {print $0}'  |\
    awk  -v OFS="\t"  -v FS="\t"  '{a=($1=="s2"?-1:1); b=($4==16?-1:1); $4=(a*b==1?0:16)}  {print $0}'   |\
    sed   's/^s\S\+\t//1' |\
    sed   's/^b\S\+\t//1' |\
    sort -k1,1 -k5,5nr |\
    sort -k1,1 -u |\
    cat $bam_header - |\
    samtools view -S -b   - |\
    samtools sort -@100  >    $bam_polyA_processed
samtools index $bam_polyA_processed $bam_polyA_processed.bai


echo "subsampling"

bam_polyA_plus=$results/trimmed_primary_polyA_plus.bam
bam_polyA_minus=$results/trimmed_primary_polyA_minus.bam

samtools view -h -b -N   $polyA_list $bam_polyA_processed  >  $bam_polyA_plus


samtools index  $bam_polyA_plus  $bam_polyA_plus.bai

samtools view $bam_polyA_plus | awk '{print $1}'  > $results/mapped_with_polyA.list
samtools view $bam_polyA_processed | awk '{print $1}'  > $results/mapped_all.list
cat $results/mapped_with_polyA.list $results/mapped_all.list | sort | uniq -c | awk '$1==1{print $2}' > $results/mapped_wo_polyA.list


samtools view -h -b -N   $results/mapped_wo_polyA.list $bam_polyA_processed  | samtools sort -@40 >  $bam_polyA_minus
samtools index  $bam_polyA_minus  $bam_polyA_minus.bai
