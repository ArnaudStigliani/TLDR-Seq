#!/bin/bash


./bedfiles_from_bam.sh lib1
./bedfiles_from_bam.sh lib2
./bedfiles_from_bam.sh R10

./write_distance_tables.sh

Rscript ./plot_distance_v2.r

get_motif_from_read_starts.sh

Rscript ./plot_motif.r

exit 0
