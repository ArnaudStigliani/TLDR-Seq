#!/bin/bash

output=./results/write_distance_tables/
mkdir -p  $output


sliccage_out_dir=../data/sliccage_bedfiles 
slic_cage_1_name=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_1.bed
slic_cage_2_name=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_2.bed
slic_cage_3_name=$sliccage_out_dir/hg38.SLICCAGE_CPH4_SiGFB_3.bed

gencodeall_tss_name="../data/gencode_v39_all_transcripts_for_TLDR.bed"
gencodesubset_tss_name="../data/gencode_v39_subset_selected_transcripts_for_TLDR.bed"


lib1_tss_name=./results/bedfiles_from_bam/lib1.bed
lib2_tss_name=./results/bedfiles_from_bam/lib2.bed
R10_lib3_tss_name=./results/bedfiles_from_bam/R10.bed




lib1_tss_trimmed_name=$output/lib1_trimmed_1bp.bed
lib2_tss_trimmed_name=$output/lib2_trimmed_1bp.bed
R10_lib3_tss_trimmed_name=$output/R10_lib3_trimmed_1bp.bed


nanopore1_name=../map_ONTnanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed 
nanopore2_name=../map_ONTnanopore/get_bedfiles/nanopore_wt2_curated_stranded.bed 


awk -v OFS="\t" '$6=="+" {print $1,$2+1, $3,$4,$5,$6} $6=="-" {print $1,$2, $3-1, $4, $5, $6}' $lib1_tss_name >  $lib1_tss_trimmed_name &
awk -v OFS="\t" '$6=="+" {print $1,$2+1, $3,$4,$5,$6} $6=="-" {print $1,$2, $3-1, $4, $5, $6}' $lib2_tss_name >  $lib2_tss_trimmed_name &
awk -v OFS="\t" '$6=="+" {print $1,$2+1, $3,$4,$5,$6} $6=="-" {print $1,$2, $3-1, $4, $5, $6}' $R10_lib3_tss_name >  $R10_lib3_tss_trimmed_name 

wait


# R10 
Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $R10_lib3_tss_trimmed_name -b "${slic_cage_1_name},${slic_cage_2_name},${slic_cage_3_name}" -a -o $output/R10_lib3_in_slic_a.tsv -u $output/slic_in_R10_lib3_a.tsv -r "50" &
Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $R10_lib3_tss_trimmed_name -b $gencodesubset_tss_name -a -t 1 -o $output/R10_lib3_in_gencode_a.tsv -u $output/gencode_in_R10_lib3_a.tsv -r "50" &

wait


Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib1_tss_trimmed_name -b $gencodesubset_tss_name -a -t 1 -o $output/lib1_in_gencode_a.tsv -u $output/gencode_in_lib1_a.tsv -r "50" &
Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib2_tss_trimmed_name -b $gencodesubset_tss_name -a -t 1 -o $output/lib2_in_gencode_a.tsv -u $output/gencode_in_lib2_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $gencodesubset_tss_name -b "${slic_cage_1_name},${slic_cage_2_name},${slic_cage_3_name}" -s 1 -a -o $output/gencode_in_slic_a.tsv -u $output/slic_in_gencode_a.tsv -r "50" &
Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib1_tss_trimmed_name -b "${slic_cage_1_name},${slic_cage_2_name},${slic_cage_3_name}" -a -o $output/lib1_in_slic_a.tsv -u $output/slic_in_lib1_a.tsv -r "50" &
Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $"${nanopore1_name},${nanopore2_name}" -b "${slic_cage_1_name},${slic_cage_2_name},${slic_cage_3_name}" -a -o $output/nanopore_in_slic_a.tsv -u $output/slic_in_nanopore_a.tsv -r "50" &

sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib2_tss_trimmed_name -b "${slic_cage_1_name},${slic_cage_2_name},${slic_cage_3_name}" -a -o $output/lib2_in_slic_a.tsv -u $output/slic_in_lib2_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $gencodesubset_tss_name -b "${nanopore1_name},${nanopore2_name}" -s 1 -a -o $output/gencode_in_nanopore_a.tsv -u $output/nanopore_in_gencode_a.tsv -r "50" &




wait
