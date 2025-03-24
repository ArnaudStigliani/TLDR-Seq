#!/bin/bash

sequencing_run=$1
fast5_dir=$2


pipeline_map.sh  $sequencing_run
format_adapter_informations.sh $sequencing_run
analyze_adapter_informations.r $sequencing_run
demultiplex_read_names.sh $sequencing_run
demultiplex_reads_from_fast5.sh $sequencing_run $fast5_dir
get_polyA_tailfindR.r $sequencing_run
characterize_polyA.r $sequencing_run
get_polyA_bam_files.sh $sequencing_run


exit 0
