#!/bin/bash
module load samtools
module load cutadapt
module load minimap2

geno=$1
results=./results/$geno/pipeline_map/

echo $results
mkdir -p  $results


fastq_in=../data/all_pass_${geno}.fastq.gz

results_trimmed=$results/trimming
results_mapped=$results/mapping

mkdir -p $results $results_trimmed $results_mapped 


## adapters/barcodes

A1rev=AAGAAAGTTGTCGGTGTCTTTGTGCGGATGAACATAGGATAGCGATTCACGGTATGTC
A2_bc9=TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC

A2_bc9_rev=$(echo -n $A2_bc9 | tr '[ACGT]' '[TGCA]' | rev)
A1=$(echo -n $A1rev | tr '[ACGT]' '[TGCA]' | rev)

#################  trimming


if [ ! -s  $results_trimmed/trimmed.fastq.gz ]
then
    
    cutadapt  -g $A2_bc9...$A1rev  -o $results_trimmed/bc9_s1.fastq.gz  -j 10 -O 8 -m 1 --discard-untrimmed -x s1_bc9_  $fastq_in > $results_trimmed/bc9_s1.log &
    cutadapt  -g $A1...$A2_bc9_rev  -o $results_trimmed/bc9_s2.fastq.gz  -j 10 -O 8 -m 1 --discard-untrimmed -x s2_bc9_  $fastq_in > $results_trimmed/bc9_s2.log &
fi
wait

cat $results_trimmed/bc9*.gz > $results_trimmed/trimmed.fastq.gz

############ mapping
genome_index=../data/Mus_musculus.GRCm39.dna.primary_assembly.fa
genome_bedgff=../data/anno_junc_mm.bed

echo "minimap2 -ax splice --junc-bed -un $genome_bedgff -t 10  $genome_index $results_trimmed/trimmed.fastq.gz"

minimap2 -ax splice -un --junc-bed  $genome_bedgff -t 10  $genome_index $results_trimmed/trimmed.fastq.gz   | samtools view  -S -b  > $results_mapped/trimmed.bam  2> $results_mapped/trimmed.log 


############ get primary alignment genome

bam=$results_mapped/trimmed.bam
bam_H=$results_mapped/trimmed.headers
bam_primmary=$results_mapped/trimmed_primary.bam

samtools view -H $bam > $bam_H
samtools view $bam  -F 256 -F 4 | awk -v OFS="\t" '{gsub("AS:i:","",$14)}{print $0}' | sort -k1,1 -k14,14nr | sort -k1,1 -u | awk -v OFS="\t" '{gsub("^","AS:i:",$14)} {print $0}' | cat $bam_H  - | samtools view -S -b   | samtools sort -@ 10 | tee $bam_primmary  | samtools index  - $bam_primmary.bai  




exit 0
