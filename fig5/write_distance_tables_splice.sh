#!/bin/bash

output=./results/write_distance_tables_splice

## acceptor and  donor splice sites
dir_intron_start=$output/intron_start 
dir_intron_end=$output/intron_end

mkdir -p  $output $dir_intron_start $dir_intron_end

in_dir=./results/get_splice_sites

lib1_tss_name=$in_dir/tldr_intron.bed
lib2_tss_name=$in_dir/tldr_lib2_intron.bed
lib3_tss_name=$in_dir/tldr_lib3_intron.bed

gencode_tss_name=$in_dir/gencode_intron.bed 
nanopore_name=$in_dir/nanopore_intron.bed

lib1_polyA_reads=../fig1/results/characterize_polyA/lib1/reads_with_polyA.list
lib2_polyA_reads=../fig1/results/characterize_polyA/lib2/reads_with_polyA.list
lib3_polyA_reads=../mapping_preprocessing/results/get_pA_reads/reads_pA.txt

### prepare TLDRSeq polyA+ bed files
lib1_tss_name_polyA=$output/tldr_intron_polyA.bed
grep -f  $lib1_polyA_reads  $lib1_tss_name >  $lib1_tss_name_polyA

lib2_tss_name_polyA=$output/tldr_lib2_intron_polyA.bed
grep -f  $lib2_polyA_reads  $lib2_tss_name >  $lib2_tss_name_polyA

lib3_tss_name_polyA=$output/tldr_lib3_intron_polyA.bed
grep -f  $lib3_polyA_reads  $lib3_tss_name >  $lib3_tss_name_polyA 


output=../../results/sequencing_july_2023/shared/reviews/write_distance_tables_splice
dir_intron_start=$output/intron_start
dir_intron_end=$output/intron_end




## ## intron end
### nanopore vs gencode

Rscript ../sequencing_july_2023/compare_TTSs.r -f $nanopore_name -b $gencode_tss_name -a -t 1 -o $dir_intron_end/nanopore_in_gencode_a.tsv -u $dir_intron_end/gencode_in_nanopore_a.tsv -r "50" &
sleep 10

### lib1 vs stuff


Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib1_tss_name_polyA -b $nanopore_name -a -o $dir_intron_end/lib1polyA_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib1polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib1_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_end/lib1_in_gencode_a.tsv -u $dir_intron_end/gencode_in_lib1_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib1_tss_name -b $nanopore_name -a -o $dir_intron_end/lib1_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib1_a.tsv -r "50" &
sleep 10 

## lib2

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib2_tss_name_polyA -b $nanopore_name -a -o $dir_intron_end/lib2polyA_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib2polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib2_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_end/lib2_in_gencode_a.tsv -u $dir_intron_end/gencode_in_lib2_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib2_tss_name -b $nanopore_name -a -o $dir_intron_end/lib2_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib2_a.tsv -r "50" &
sleep 10 


## lib3

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib3_tss_name_polyA -b $nanopore_name -a -o $dir_intron_end/lib3polyA_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib3polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib3_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_end/lib3_in_gencode_a.tsv -u $dir_intron_end/gencode_in_lib3_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TTSs.r -f $lib3_tss_name -b $nanopore_name -a -o $dir_intron_end/lib3_in_nanopore_a.tsv -u $dir_intron_end/nanopore_in_lib3_a.tsv -r "50" &
sleep 10 


wait



#### ###### intron start

## nanopore vs gencode

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $nanopore_name -b $gencode_tss_name -a -t 1 -o $dir_intron_start/nanopore_in_gencode_a.tsv -u $dir_intron_start/gencode_in_nanopore_a.tsv -r "50" &
sleep 10

### lib1


Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib1_tss_name_polyA -b $nanopore_name -a -o $dir_intron_start/lib1polyA_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib1polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib1_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_start/lib1_in_gencode_a.tsv -u $dir_intron_start/gencode_in_lib1_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib1_tss_name -b $nanopore_name -a -o $dir_intron_start/lib1_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib1_a.tsv -r "50" &
sleep 10 

## lib2


Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib2_tss_name_polyA -b $nanopore_name -a -o $dir_intron_start/lib2polyA_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib2polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib2_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_start/lib2_in_gencode_a.tsv -u $dir_intron_start/gencode_in_lib2_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib2_tss_name -b $nanopore_name -a -o $dir_intron_start/lib2_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib2_a.tsv -r "50" &
sleep 10 

## lib3


Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib3_tss_name_polyA -b $nanopore_name -a -o $dir_intron_start/lib3polyA_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib3polyA_a.tsv -r "50" &
sleep 10 

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib3_tss_name -b $gencode_tss_name -a -t 1 -o $dir_intron_start/lib3_in_gencode_a.tsv -u $dir_intron_start/gencode_in_lib3_a.tsv -r "50" &
sleep 10

Rscript ../sequencing_july_2023/compare_TSSs_v2.r -f $lib3_tss_name -b $nanopore_name -a -o $dir_intron_start/lib3_in_nanopore_a.tsv -u $dir_intron_start/nanopore_in_lib3_a.tsv -r "50" &
sleep 10 
