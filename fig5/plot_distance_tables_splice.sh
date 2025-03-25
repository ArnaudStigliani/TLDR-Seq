#!/bin/bash


## acceptor and  donor splice sites
data_start=./results/write_distance_tables_splice/intron_start 
data_end=./results/write_distance_tables_splice/intron_end


results=./results/plot_distance_splice


mkdir -p  $results

### pannels E AND F

##lib1
Rscript ../fig3/plot_distance_v2.r -p $data_start/lib1_in_gencode_a.tsv    -o $results/R9_Amp_in_gencode_a_donor/ -y 30000
Rscript ../fig3/plot_distance_v2.r -p $data_end/lib1_in_gencode_a.tsv    -o $results/R9_Amp_in_gencode_a_acceptor/ -y 30000


## lib2
Rscript ../fig3/plot_distance_v2.r -p $data_start/lib2_in_gencode_a.tsv   -o $results/R9_NoAmp_in_gencode_a_donor/ -y 30000
Rscript ../fig3/plot_distance_v2.r -p $data_end/lib2_in_gencode_a.tsv    -o $results/R9_NoAmp_in_gencode_a_acceptor/ -y 30000


## lib3 for sup fig
Rscript ../fig3/plot_distance_v2.r -p $data_start/lib3_in_gencode_a.tsv   -o $results/R10_Amp_in_gencode_a_donor/ -y 30000
Rscript ../fig3/plot_distance_v2.r -p $data_end/lib3_in_gencode_a.tsv    -o $results/R10_Amp_in_gencode_a_acceptor/ -y 30000


# genocode vs nanopore

Rscript ../fig3/plot_distance_v2.r -p $data_start/nanopore_in_gencode_a.tsv    -o $results/ONT_in_gencode_a_donor/ -y 30000
Rscript ../fig3/plot_distance_v2.r -p $data_end/nanopore_in_gencode_a.tsv    -o $results/ONT_in_gencode_a_acceptor/ -y 30000


## pannel A and C

Rscript ../fig3/plot_distance_v2.r -p $data_start/lib3_in_nanopore_a.tsv   -o $results/R10_Amp_in_ONT_a_donor/ -y 60000
Rscript ../fig3/plot_distance_v2.r -p $data_end/lib3_in_nanopore_a.tsv    -o $results/R10_Amp_in_ONT_a_acceptor/ -y 60000

## pannel B and D

## nanopore  in lib3 for Sup fig
Rscript ../fig3/plot_distance_v2.r -p $data_start/nanopore_in_lib3_a.tsv   -o $results/ONT_in_R10_Amp_a_donor/ -y 60000
Rscript ../fig3/plot_distance_v2.r -p $data_end/nanopore_in_lib3_a.tsv    -o $results/ONT_in_R10_Amp_a_acceptor/ -y 60000

# Sup fig S2
# same with pA+ only


Rscript ../fig3/plot_distance_v2.r -p $data_start/nanopore_in_lib3polyA_a.tsv    -o $results/ONT_in_R10_Amp_a_pA_donor/ 
Rscript ../fig3/plot_distance_v2.r -p $data_end/nanopore_in_lib3polyA_a.tsv    -o $results/ONT_in_R10_Amp_a_pA_acceptor/ 




exit 0
