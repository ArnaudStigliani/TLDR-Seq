#!/bin/bash

sequencing_run=$1
fast5_dir=$2 ## add the path of the fast5 files



data=./results/demultiplex_read_names/$sequencing_run

results=./results/demultiplex_reads_from_fast5/$sequencing_run
echo $results
mkdir -p $results 
read_names_dir=$data

find $read_names_dir -name "bc*" | xargs --max-procs=13 -I {} bash -c "a=\$(basename {} _reads.txt) ; fast5_subset -i $fast5_dir -s $results/\$a -l {}"

