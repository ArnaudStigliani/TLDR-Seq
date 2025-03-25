#!/bin/bash

./pipeline_map.sh ../data/all_pass_mES-ZCCHC8.fastq.gz & 
./pipeline_map.sh ../data/all_pass_mES-WT1.fastq.gz &
./pipeline_map.sh ../data/all_pass_mES-WT2.fastq.gz   &
./pipeline_map.sh ../data/all_pass_mES-RBM7.fastq.gz   &

wait

./demultiplex_read_names.sh ../data/all_pass_mES-ZCCHC8.fastq.gz & 
./demultiplex_read_names.sh ../data/all_pass_mES-WT1.fastq.gz &
./demultiplex_read_names.sh ../data/all_pass_mES-WT2.fastq.gz   &
./demultiplex_read_names.sh ../data/all_pass_mES-RBM7.fastq.gz   &
