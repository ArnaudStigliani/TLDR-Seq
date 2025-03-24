#!/bin/bash

./get_polyA_bam_files_R10.sh
Rscript reproducibility.r
Rscript reproducibility_R10.r

Rscript compare_lib1_lib3.r

