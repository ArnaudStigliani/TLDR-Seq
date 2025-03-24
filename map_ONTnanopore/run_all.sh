#!/bin/bash

./pipeline_map.sh
./restrand_nanopore_reads.sh
./get_bedfiles.sh

exit 0
