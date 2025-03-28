#!/bin/bash


module load samtools

data=get_nanocounts
data_genome=results_curated

results_curated_transcriptome=results_curated_transcriptom
results_curated_genome=results_curated_genome

mkdir -p $results_curated_transcriptome $results_curated_genome

nanopore1_bam_transcriptome=$data/nanopore1_N10.bam
nanopore2_bam_transcriptome=$data/nanopore2_N10.bam

bam1_curated_transcriptome=$results_curated_transcriptome/nanopore_wt1_curated.bam
bam2_curated_transcriptome=$results_curated_transcriptome/nanopore_wt2_curated.bam


bam_H_transcriptome=$results_curated_transcriptome/nanopore_transcriptome.headers
samtools view -H $nanopore1_bam_transcriptome > $bam_H_transcriptome

samtools view $nanopore1_bam_transcriptome  -F 256 -F 4 |\
    sort -k1,1 -k5,5nr |\
    sort -k1,1 -u |\
    cat $bam_H_transcriptome - |\
    samtools view -S -b   - > $bam1_curated_transcriptome

samtools view $nanopore2_bam_transcriptome  -F 256 -F 4 |\
    sort -k1,1 -k5,5nr |\
    sort -k1,1 -u |\
    cat $bam_H_transcriptome - |\
    samtools view -S -b   - > $bam2_curated_transcriptome



samtools view $bam1_curated_transcriptome | awk -v OFS="\t" -v FS="[|\t]" '{print $1,$3}' > $results_curated_transcriptome/reads_transcript1.tsv
samtools view $bam2_curated_transcriptome | awk -v OFS="\t" -v FS="[|\t]" '{print $1,$3}' > $results_curated_transcriptome/reads_transcript2.tsv


#######
gencode_transcript_annotation=../../../shared_data/human_genome/gencode_v39_all_transcripts_for_TLDR.bed  

join <(sort -k2,2 $results_curated_transcriptome/reads_transcript1.tsv )  <(sort -k4,4 $gencode_transcript_annotation) -1 2 -2 4  | awk -v OFS="\t" '{print $2,$7}' > $results_curated_transcriptome/reads_strand1.tsv
join <(sort -k2,2 $results_curated_transcriptome/reads_transcript2.tsv )  <(sort -k4,4 $gencode_transcript_annotation) -1 2 -2 4  | awk -v OFS="\t" '{print $2,$7}' > $results_curated_transcriptome/reads_strand2.tsv


nanopore1_bam_genome=$data_genome/nanopore_wt1_curated.bam
nanopore2_bam_genome=$data_genome/nanopore_wt2_curated.bam

bam_H_genome=$results_curated_genome/nanopore_genome.headers
samtools view -H $nanopore1_bam_genome > $bam_H_genome


bam1_curated_genome=$results_curated_genome/nanopore_wt1_curated.bam
bam2_curated_genome=$results_curated_genome/nanopore_wt2_curated.bam


join <(sort -k1,1 $results_curated_transcriptome/reads_strand1.tsv )  <(samtools view $nanopore1_bam_genome |  sort -k1,1 ) -1 1 -2 1 |\
    sed 's/ /\t/g' |\
    sed 's/\t+\t[0-9][0-9]*\t/\t0\t/' |\
    sed 's/\t-\t[0-9][0-9]*\t/\t16\t/' |\
    cat $bam_H_genome - |\
    samtools view -S -b > $bam1_curated_genome &


join <(sort -k1,1 $results_curated_transcriptome/reads_strand2.tsv )  <(samtools view $nanopore2_bam_genome |  sort -k1,1 ) -1 1 -2 1 |\
    sed 's/ /\t/g' |\
    sed 's/\t+\t[0-9][0-9]*\t/\t0\t/' |\
    sed 's/\t-\t[0-9][0-9]*\t/\t16\t/' |\
    cat $bam_H_genome - |\
    samtools view -S -b > $bam2_curated_genome

wait


samtools sort $bam1_curated_genome | tee  $results_curated_genome/nanopore_wt1_curated_sorted.bam | samtools index -
samtools sort $bam2_curated_genome | tee  $results_curated_genome/nanopore_wt2_curated_sorted.bam | samtools index -

