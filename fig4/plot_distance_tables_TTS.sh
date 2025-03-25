#!/bin/bash

data=./results/write_distance_tables_TTS
results=./results/plot_distance_TTS


mkdir -p  $results



Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_gencode_a.tsv  -o $results/lib1_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_lib1_a.tsv  -o $results/gencode_in_lib1_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib1xpap_in_gencode_a.tsv  -o $results/lib1xpap_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_lib1xpap_a.tsv  -o $results/gencode_in_lib1xpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_gencode_a.tsv  -o $results/lib2_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_lib2_a.tsv  -o $results/gencode_in_lib2_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib2xpap_in_gencode_a.tsv  -o $results/lib2xpap_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_lib2xpap_a.tsv  -o $results/gencode_in_lib2xpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeq_in_gencode_a.tsv  -o $results/quantSeq_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_quantSeq_a.tsv  -o $results/gencode_in_quantSeq_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeq_in_lib1_a.tsv  -o $results/quantSeq_in_lib1_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_quantSeq_a.tsv  -o $results/lib1_in_quantSeq_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeq_in_lib2_a.tsv  -o $results/quantSeq_in_lib2_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_quantSeq_a.tsv  -o $results/lib2_in_quantSeq_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeqxpap_in_lib1_a.tsv  -o $results/quantSeqxpap_in_lib1_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_quantSeqxpap_a.tsv  -o $results/lib1_in_quantSeqxpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeqxpap_in_lib2_a.tsv  -o $results/quantSeqxpap_in_lib2_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_quantSeqxpap_a.tsv  -o $results/lib2_in_quantSeqxpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_nanopore_a.tsv  -o $results/lib1_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_lib1_a.tsv  -o $results/nanopore_in_lib1_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib1xpap_in_nanopore_a.tsv  -o $results/lib1xpap_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_lib1xpap_a.tsv  -o $results/nanopore_in_lib1xpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_nanopore_a.tsv  -o $results/lib2_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_lib2_a.tsv  -o $results/nanopore_in_lib2_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib2xpap_in_nanopore_a.tsv  -o $results/lib2xpap_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_lib2xpap_a.tsv  -o $results/nanopore_in_lib2xpap_a/

Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_nanopore_a.tsv  -o $results/gencode_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_gencode_a.tsv  -o $results/nanopore_in_gencode_a/

Rscript ./fig3/plot_distance_v2.r -p $data/gencode_protCoding_in_nanopore.tsv  -o $results/gencode_protCoding_in_nanopore/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_gencode_protCoding.tsv  -o $results/nanopore_in_gencode_protCoding/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeq_in_nanopore_a.tsv  -o $results/quantSeq_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_quantSeq_a.tsv  -o $results/nanopore_in_quantSeq_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_lib2_a.tsv  -o $results/lib1_in_lib2_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_lib1_a.tsv  -o $results/lib2_in_lib1_a/

Rscript ./fig3/plot_distance_v2.r -p $data/lib1_in_lib2_polyA_a.tsv  -o $results/lib1_in_lib2_polyA_a/
Rscript ./fig3/plot_distance_v2.r -p $data/lib2_in_lib1_polyA_a.tsv  -o $results/lib2_in_lib1_polyA_a/

Rscript ./fig3/plot_distance_v2.r -p $data/quantSeqxpap_in_nanopore_a.tsv  -o $results/quantSeqxpap_in_nanopore_a/
Rscript ./fig3/plot_distance_v2.r -p $data/nanopore_in_quantSeqxpap_a.tsv  -o $results/nanopore_in_quantSeqxpap_a/


Rscript ./fig3/plot_distance_v2.r -p $data/quantSeqxpap_in_gencode_a.tsv  -o $results/quantSeqxpap_in_gencode_a/
Rscript ./fig3/plot_distance_v2.r -p $data/gencode_in_quantSeqxpap_a.tsv  -o $results/gencode_in_quantSeqxpap_a/



#### R10 


Rscript ./fig3/plot_distance_v2.r -p $data/R10_lib1pA_in_gencode_a.tsv   -o $results/R10_lib1pA_in_gencode_a/

Rscript ./fig3/plot_distance_v2.r -p $data/R10_lib1TotRna_in_gencode_a.tsv   -o $results/R10_lib1TotRna_in_gencode_a/


Rscript ./fig3/plot_distance_v2.r -p $data/R10_lib3pA_in_quantSeq_noPAP_a.tsv   -o $results/R10_lib3pA_in_quantSeq_noPAP_a

Rscript ./fig3/plot_distance_v2.r -p $data/R10_lib3TotRna_in_quantSeq_PAP_a.tsv   -o $results/R10_lib3TotRna_in_quantSeq_PAP_a
