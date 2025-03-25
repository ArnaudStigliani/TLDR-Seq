#!/bin/bash

results=./results/get_motifs_from_read_starts
mkdir -p $results

genome_fas=../data/GRCh38.primary_assembly.genome.fa
module load bedtools 


sliccage_out_dir=../data/sliccage_bedfiles 
slic_cage_1_bed=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_1.bed
slic_cage_2_bed=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_2.bed
slic_cage_3_bed=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_3.bed

nanopore1_name=../map_ONTnanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed 
nanopore2_name=../map_ONTnanopore/get_bedfiles/nanopore_wt2_curated_stranded.bed 


nanopore_all_bed=$results/nanopore_all_curated_stranded.bed
cat $nanopore1_bed $nanopore2_bed > $nanopore_all_bed

slic_cage_all_bed=$results/slic_cage_all_curated_stranded.bed
cat $slic_cage_3_bed $slic_cage_2_bed $slic_cage_1_bed > $slic_cage_all_bed


#### output tss files
nanopore_tss_bed=$results/nanopore_tss.bed
nanopore_tss_fas=$results/nanopore_tss.fas

slic_cage_tss_bed=$results/slic_cage_tss.bed
slic_cage_tss_fas=$results/slic_cage_tss.fas

tldr_tss_bed=$results/tldr_tss.bed
tldr_tss_fas=$results/tldr_tss.fas

tldrlib2_tss_bed=$results/tldrlib2_tss.bed
tldrlib2_tss_fas=$results/tldrlib2_tss.fas


cat  $tldr_lib1_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="+" && $2-49 >0 {print $1,$2-49, $2+51, ".", ".",$6} $1~/chr/ &&  $6=="-" && $3-49 > 0 {print $1,$3-51, $3+49, ".", ".", $6}' >  $tldr_tss_bed ## trim 1bp 
cat  $tldr_lib2_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="+" && $2-49 >0 {print $1,$2-49, $2+51, ".", ".",$6} $1~/chr/ &&  $6=="-" && $3-49 > 0 {print $1,$3-51, $3+49, ".", ".", $6}' >  $tldrlib2_tss_bed ## trim 1bp 
cat  $nanopore_all_bed | awk   -v OFS="\t" '$1~/chr/ && $6=="+" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="-" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $nanopore_tss_bed
cat  $slic_cage_all_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="+" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="-" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $slic_cage_tss_bed


wait

bedtools getfasta -fi $genome_fas -fo $slic_cage_tss_fas -bed $slic_cage_tss_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $nanopore_tss_fas -bed $nanopore_tss_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $tldr_tss_fas -bed $tldr_tss_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $tldrlib2_tss_fas -bed $tldrlib2_tss_bed -s 2> /dev/null &
wait
exit 0

 
