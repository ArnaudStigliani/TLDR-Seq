#!/bin/bash

output=./results/write_distance_tables_TTS
mkdir -p  $output



### prepare TLDRSeq polyA+ bed files
lib1_polyA_reads=../fig1/results/characterize_polyA/lib1/reads_with_polyA.list
lib1_tss_name=./results/bedfiles_from_bam/lib1.bed
lib1_tss_name_polyA=$output/lib1_polyA_reads_stranded.bed

lib2_polyA_reads=../fig1/results/characterize_polyA/lib2/reads_with_polyA.list
lib2_tss_name=./results/bedfiles_from_bam/lib2.bed
lib2_tss_name_polyA=$output/lib2_polyA_reads_stranded.bed

R10_lib1_polyA_reads=../mapping_preprocessing/results/get_pA_reads/reads_pA.txt
R10_lib1_tss_name=./results/reviews/bedfile_from_bam/R10.bed
R10_lib1_tss_name_polyA=$output/R10_lib1_polyA_reads_stranded.bed

R10_lib3_tss_name=./results/bedfiles_from_bam/R10.bed

grep -f  $lib1_polyA_reads  $lib1_tss_name >  $lib1_tss_name_polyA
grep -f  $lib2_polyA_reads  $lib2_tss_name >  $lib2_tss_name_polyA
grep -f  $R10_lib1_polyA_reads  $R10_lib1_tss_name >  $R10_lib1_tss_name_polyA



gencodeall_tss_name="../data/gencode_v39_all_transcripts_for_TLDR.bed"
gencodesubset_tss_name="../data/gencode_v39_subset_selected_transcripts_for_TLDR.bed"


QS_dir=../data/QS_dir/
####  mapped reads for quantSeq in bed format (noPAP = polyA selection without artificial poly adenylation)
quantSeq_nopap1=$QS_dir/siGFP_noPAP_in_1.bed 
quantSeq_nopap2=$QS_dir/siGFP_noPAP_in_2.bed 
quantSeq_nopap3=$QS_dir/siGFP_noPAP_in_3.bed 
quantSeq_nopap4=$QS_dir/siGFP_noPAP_in_4.bed 
quantSeq_nopap5=$QS_dir/siGFP_noPAP_in_5.bed

###  (xPAP = polyA selection with artificial  poly  adenylationof all the transcripts before pA seslection)
quantSeq_xpap1=$QS_dir/siGFP_xPAP_in_1.bed 
quantSeq_xpap2=$QS_dir/siGFP_xPAP_in_2.bed 
quantSeq_xpap3=$QS_dir/siGFP_xPAP_in_3.bed 
quantSeq_xpap4=$QS_dir/siGFP_xPAP_in_4.bed 
quantSeq_xpap5=$QS_dir/siGFP_xPAP_in_5.bed 

nanopore1_name=../map_ONTnanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed 
nanopore2_name=../map_ONTnanopore/get_bedfiles/nanopore_wt2_curated_stranded.bed 
allnanopore="${nanopore1_name},${nanopore2_name}"

allQuantSeq="${quantSeq_nopap1},${quantSeq_nopap2},${quantSeq_nopap3},${quantSeq_nopap4},${quantSeq_nopap5}"
allQuantSeq_xpap="${quantSeq_xpap1},${quantSeq_xpap2},${quantSeq_xpap3},${quantSeq_xpap4},${quantSeq_xpap5}"



### R10 for supplementary ####
Rscript ../sequencing_july_2023/compare_TTSs.r -f $R10_lib1_tss_name_polyA -b $gencode_tss_name -a -t 1 -o $output/R10_lib1pA_in_gencode_a.tsv -u $output/gencode_in_R10_lib1pA_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TTSs.r -f $R10_lib1_tss_name -b $gencodesubset_tss_name -a -t 1 -o $output/R10_lib1TotRna_in_gencode_a.tsv -u $output/gencode_in_R10_lib1TotRna_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib1_tss_name_polyA -b $allQuantSeq -a -o $output/lib1pA_in_quantSeq_a.tsv -u $output/quantSeq_in_lib1pA_a.tsv -r "50" &
sleep 10




### the rest
## lib1 and lib2 vs gencode



Rscript compare_TTSs.r -f $lib1_tss_name_polyA -b $gencode_tss_name -a -t 1 -o $output/lib1_in_gencode_a.tsv -u $output/gencode_in_lib1_a.tsv -r "50" &
sleep 10

Rscript compare_TTSs.r -f $lib2_tss_name_polyA -b $gencode_tss_name -a -t 1 -o $output/lib2_in_gencode_a.tsv -u $output/gencode_in_lib2_a.tsv -r "50" &
sleep 10

Rscript compare_TTSs.r -f $lib1_tss_name -b $gencodesubset_tss_name -a -t 1 -o $output/lib1xpap_in_gencode_a.tsv -u $output/gencode_in_lib1xpap_a.tsv -r "50" &
sleep 10

Rscript compare_TTSs.r -f $lib2_tss_name -b $gencodesubset_tss_name -a -t 1 -o $output/lib2xpap_in_gencode_a.tsv -u $output/gencode_in_lib2xpap_a.tsv -r "50" &
sleep 10



## lib1 and lib2  vs quantSeq
wait 
Rscript compare_TTSs.r -f $gencode_tss_name -b $allQuantSeq -s 1 -a  -o $output/gencode_in_quantSeq_a.tsv -u $output/quantSeq_in_gencode_a.tsv -r "50" &
sleep 10

wait 
Rscript compare_TTSs.r -f $lib1_tss_name_polyA -b $allQuantSeq -a -o $output/lib1_in_quantSeq_a.tsv -u $output/quantSeq_in_lib1_a.tsv -r "50" &
sleep 10

Rscript compare_TTSs.r -f $lib2_tss_name_polyA -b $allQuantSeq -a -o $output/lib2_in_quantSeq_a.tsv -u $output/quantSeq_in_lib2_a.tsv -r "50" &
sleep 10

Rscript compare_TTSs.r -f $lib1_tss_name -b $allQuantSeq_xpap -a -o $output/lib1_in_quantSeqxpap_a.tsv -u $output/quantSeqxpap_in_lib1_a.tsv -r "50" &
sleep 10
wait
Rscript compare_TTSs.r -f $lib2_tss_name -b $allQuantSeq_xpap -a -o $output/lib2_in_quantSeqxpap_a.tsv -u $output/quantSeqxpap_in_lib2_a.tsv -r "50" &
sleep 10

## lib1 and  lib2 vs ONT nanopore
Rscript compare_TTSs.r -f $lib1_tss_name_polyA -b $allnanopore -a -o $output/lib1_in_nanopore_a.tsv -u $output/nanopore_in_lib1_a.tsv -r "50" &
sleep 10
Rscript compare_TTSs.r -f $lib2_tss_name_polyA -b $allnanopore -a -o $output/lib2_in_nanopore_a.tsv -u $output/nanopore_in_lib2_a.tsv -r "50" &
sleep 10


Rscript compare_TTSs.r -f $lib1_tss_name -b $allnanopore -a -o $output/lib1xpap_in_nanopore_a.tsv -u $output/nanopore_in_lib1xpap_a.tsv -r "50" &
sleep 10
Rscript compare_TTSs.r -f $lib2_tss_name -b $allnanopore -a -o $output/lib2xpap_in_nanopore_a.tsv -u $output/nanopore_in_lib2xpap_a.tsv -r "50" &
sleep 10

wait 

## ONT nanopore vs gencode
Rscript compare_TTSs.r -f $allnanopore -b $gencode_tss_name -a -t 1 -o $output/nanopore_in_gencode_a.tsv -u $output/gencode_in_nanopore_a.tsv -r "50" &
sleep 10
Rscript compare_TTSs.r -f $allnanopore -b $allQuantSeq -a -o $output/nanopore_in_quantSeq_a.tsv -u $output/quantSeq_in_nanopore_a.tsv -r "50" &
wait
sleep 10

## lib1 vs lib2
Rscript compare_TTSs.r -f $lib2_tss_name -b $lib1_tss_name -a -o $output/lib2_in_lib1_a.tsv -u $output/lib1_in_lib2_a.tsv -r "50" 
Rscript compare_TTSs.r -f $lib2_tss_name_polyA -b $lib1_tss_name_polyA -a -o $output/lib2_in_lib1_polyA_a.tsv -u $output/lib1_in_lib2_polyA_a.tsv -r "50" &


## quantSeq vs ONT nanopore
Rscript compare_TTSs.r -f $allnanopore -b $allQuantSeq_xpap -a -o $output/nanopore_in_quantSeqxpap_a.tsv -u $output/quantSeqxpap_in_nanopore_a.tsv -r "50" &
sleep 10
wait

## quantSeq vs gencode
Rscript compare_TTSs.r -f $gencodesubset_tss_name -b $allQuantSeq_xpap -s 1 -a  -o $output/gencode_in_quantSeqxpap_a.tsv -u $output/quantSeqxpap_in_gencode_a.tsv -r "50" &



exit 0

