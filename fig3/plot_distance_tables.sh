#!/bin/bash

data=./results/write_distance_tables
results=./results/plot_distance


mkdir -p  $results

##### R10

Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/R10_lib3_in_gencode_a.tsv    -o $results/R10_lib3_in_gencode_a/
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/R10_lib3_in_slic_a.tsv    -o $results/R10_lib3_in_slic_a/


### remake fig 3 with ymax

Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/lib1_in_slic_a.tsv    -o $results/libAmp_in_slic_a -y 32000 &
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/lib2_in_slic_a.tsv    -o $results/lib2NoAmp_in_slic_a -y 32000 &
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/nanopore_in_slic_a.tsv    -o $results/nanopore_in_slic_a -y 32000  &

Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/lib1_in_gencode_a.tsv    -o $results/lib1Amp_in_gencode_a -y 12000 &
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/lib2_in_gencode_a.tsv    -o $results/lib2NoAmp_in_gencode_a -y 12000 &
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/nanopore_in_gencode_a.tsv    -o $results/nanopore_in_gencode_a -y 12000 &
Rscript ../sequencing_july_2023/plot_distance_v2.r -p $data/slic_in_gencode_a.tsv    -o $results/slic_in_gencode_a -y 12000 &
wait


exit 0

