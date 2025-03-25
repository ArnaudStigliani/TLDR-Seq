#!/bin/bash


./write_distance_tables_TTS.sh

plot_distance_tables_TTS.sh

get_motif_from_read_ens.sh

Rscript ./plot_motifs_TTSs.r
Rscript ./plot_motifs_TTSs_gencode.r



exit 0
