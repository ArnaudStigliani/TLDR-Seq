#!/bin/bash


### sequencing run is either lib1 (amplified) or lib1 (non amplified) or R10 (amplified and sequenced on R10 flowcell)

sequencing_run=$1
fast5_dir=$2


./pipeline_map.sh  $sequencing_run
./format_adapter_informations.sh $sequencing_run
./analyze_adapter_informations.r $sequencing_run
./demultiplex_read_names.sh $sequencing_run
./demultiplex_reads_from_fast5.sh $sequencing_run $fast5_dir
./get_polyA_tailfindR.r $sequencing_run ## Don't run if R10

./get_pA_reads_R10.sh ## run if R10


exit 0
