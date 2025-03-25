#!/bin/bash

Rscript charaterize_polyA.r ## to generate > 15nt pA list of all pA list (no threshold) for lib1 and lib2.

./get_polyA_bam_files.sh lib1
./get_polyA_bam_files.sh lib2
./get_processed_bam_files_R10.sh

Rscript annotate_reads.r # make fig 1D and other similar supplementary figures
Rscript make_fig_1.r ## make the est of fig1

exit 0
