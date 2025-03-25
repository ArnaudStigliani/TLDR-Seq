#!/bin/bash

results=./results/get_motifs_from_read_end
mkdir -p $results

genome_fas=../data/GRCh38.primary_assembly.genome.fa

module load bedtools 


### prepare TLDRSeq polyA+ bed files
lib1_polyA_reads=../fig1/results/characterize_polyA/lib1/reads_with_polyA.list
lib1_tss_name=./results/bedfiles_from_bam/lib1.bed
lib1_tss_name_polyA=$results/lib1_polyA_reads_stranded.bed

lib2_polyA_reads=../fig1/results/characterize_polyA/lib2/reads_with_polyA.list
lib2_tss_name=./results/bedfiles_from_bam/lib2.bed
lib2_tss_name_polyA=$results/lib2_polyA_reads_stranded.bed


grep -f  $lib1_polyA_reads  $tldr_lib1_bed >  $tldrpA_lib1_bed &
grep -f  $lib2_polyA_reads  $tldr_lib2_bed >  $tldrpA_lib2_bed &

wait 


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

nanopore_all_bed=$results/nanopore_all_curated_stranded.bed
cat $nanopore1_bed $nanopore2_bed > $nanopore_all_bed

quantSeq_nopap_all_bed=$results/quantSeq_nopap_all_curated_stranded.bed
cat $quantSeq_nopap1 $quantSeq_nopap2 $quantSeq_nopap3 $quantSeq_nopap4 $quantSeq_nopap5 > $quantSeq_nopap_all_bed &


quantSeq_xpap_all_bed=$results/quantSeq_xpap_all_curated_stranded.bed 
cat $quantSeq_xpap1 $quantSeq_xpap2 $quantSeq_xpap3 $quantSeq_xpap4 $quantSeq_xpap5 > $quantSeq_xpap_all_bed &


wait 

#### output tss files
nanopore_tts_bed=$results/nanopore_tts.bed
nanopore_tts_fas=$results/nanopore_tts.fas

quantSeq_nopap_tts_bed=$results/quantSeq_nopap_tts.bed
quantSeq_nopap_tts_fas=$results/quantSeq_nopap_tts.fas

quantSeq_xpap_tts_bed=$results/quantSeq_xpap_tts.bed
quantSeq_xpap_tts_fas=$results/quantSeq_xpap_tts.fas

tldr_lib1_tts_bed=$results/tldr_lib1_tts.bed
tldr_lib1_tts_fas=$results/tldr_lib1_tts.fas

tldr_lib2_tts_bed=$results/tldr_lib2_tts.bed
tldr_lib2_tts_fas=$results/tldr_lib2_tts.fas

tldrpA_lib1_tts_bed=$results/tldrpA_lib1_tts.bed
tldrpA_lib1_tts_fas=$results/tldrpA_lib1_tts.fas

tldrpA_lib2_tts_bed=$results/tldrpA_lib2_tts.bed
tldrpA_lib2_tts_fas=$results/tldrpA_lib2_tts.fas




cat  $tldr_lib1_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $tldr_lib1_tts_bed &
cat  $tldr_lib2_bed | awk   -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $tldr_lib2_tts_bed &

cat  $tldrpA_lib1_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $tldrpA_lib1_tts_bed &
cat  $tldrpA_lib2_bed | awk   -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $tldrpA_lib2_tts_bed &

cat  $nanopore_all_bed | awk   -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $nanopore_tts_bed &
cat  $quantSeq_nopap_all_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $quantSeq_nopap_tts_bed &
cat  $quantSeq_xpap_all_bed | awk  -v OFS="\t" '$1~/chr/ && $6=="-" && $2-50 >0 {print $1,$2-50, $2+50, ".", ".",$6} $1~/chr/ &&  $6=="+" && $3-50 > 0 {print $1,$3-50, $3+50, ".", ".", $6}' >  $quantSeq_xpap_tts_bed &

wait

bedtools getfasta -fi $genome_fas -fo $tldr_lib1_tts_fas -bed $tldr_lib1_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $tldr_lib2_tts_fas -bed $tldr_lib2_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $tldrpA_lib1_tts_fas -bed $tldrpA_lib1_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $tldrpA_lib2_tts_fas -bed $tldrpA_lib2_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $nanopore_tts_fas -bed $nanopore_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $quantSeq_nopap_tts_fas -bed $quantSeq_nopap_tts_bed -s 2> /dev/null &
bedtools getfasta -fi $genome_fas -fo $quantSeq_xpap_tts_fas -bed $quantSeq_xpap_tts_bed  -s 2> /dev/null &
wait 
