#!/bin/bash


module load samtools

data=$1 ## ../mapping_preprocessing/pipeline_map/{lib1 or lib2}/mappping/trimmed_primary.bam, don't try R10, the right files don't exist.
lib=$(sed 's/.*lib\(.\).*/lib\1/' <<< $data)
results=./results/get_polyA_bam_files/$lib
mkdir -p $results 


polyA_list=./results/characterize_polyA/$lib/reads_with_polyA.list
polyA_list_all=./results/characterize_polyA/$lib/reads_with_polyA_all.list
polyA_list_0_15=./results/characterize_polyA/$lib/reads_with_polyA_0-15.list

## input
bam_genome=$data
## awk output
bam_polyA_renamed=$results/trimmed_primary_renamed.bam

## normal threshold
bam_polyA_plus=$results/trimmed_primary_polyA_plus.bam
bam_polyA_minus=$results/trimmed_primary_polyA_minus.bam
## 0-15 threshold
bam_polyA_plus_0_15=$results/trimmed_primary_polyA_plus_0_15.bam
## >0 threshold
bam_polyA_plus_all=$results/trimmed_primary_polyA_plus_all.bam
bam_polyA_minus_all=$results/trimmed_primary_polyA_minus_all.bam

bam_header=$results/bam.header


samtools view -H    $bam_genome > $bam_header
samtools view   $bam_genome |  awk -v FS="\t" -v OFS="\t" '{gsub("_","\t",$1)} ($2==0 || $2==16 || $1~"^@") {print $0}'  |\
    awk  -v OFS="\t"  -v FS="\t"  '$1!~"^@"{a=($1=="s2"?-1:1); b=($4==16?-1:1); $4=(a*b==1?0:16)}  {print $0}'   |\
    sed   's/^s\S\+\t//1' |\
    sed   's/^b\S\+\t//1' |\
    awk  -v OFS="\t" '(!($1 in max) || ($5 > max[$1])) {max[$1] = $5 ; line[$1] = $0} END{ for (i in line) {print line[i]}}' |\
    cat $bam_header - |\
    samtools view -S -b   - |\
    samtools sort -@100  >    $bam_polyA_renamed
samtools index $bam_polyA_renamed $bam_polyA_renamed.bai

echo "subsampling"

samtools view -h -b -N   $polyA_list $bam_polyA_renamed  >  $bam_polyA_plus
samtools index  $bam_polyA_plus  $bam_polyA_plus.bai

samtools view -h -b -N   $polyA_list_all $bam_polyA_renamed  >  $bam_polyA_plus_all
samtools index  $bam_polyA_plus_all  $bam_polyA_plus_all.bai

samtools view -h -b -N   $polyA_list_0_15 $bam_polyA_renamed  >  $bam_polyA_plus_0_15
samtools index  $bam_polyA_plus_0_15  $bam_polyA_plus_0_15.bai

## normal threshold
cat <(samtools view $bam_polyA_plus | awk '{print $1}') <(samtools view $bam_polyA_renamed | awk '{print $1}') | sort | uniq -c | awk '$1==1{print $2}' > $results/mapped_wo_polyA.list
## all threshold (no need to generate a 0-15 complement)
cat <(samtools view $bam_polyA_plus_all | awk '{print $1}') <(samtools view $bam_polyA_renamed | awk '{print $1}') | sort | uniq -c | awk '$1==1{print $2}' > $results/mapped_wo_polyA_all.list



samtools view -h -b -N   $results/mapped_wo_polyA.list $bam_polyA_renamed  | samtools sort -@100 >  $bam_polyA_minus
samtools index  $bam_polyA_minus  $bam_polyA_minus.bai
samtools view -h -b -N   $results/mapped_wo_polyA_all.list $bam_polyA_renamed  | samtools sort -@100 >  $bam_polyA_minus_all
samtools index  $bam_polyA_minus_all  $bam_polyA_minus_all.bai

