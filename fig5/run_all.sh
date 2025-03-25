#!/bin/bash

Rscript get_splice_sites.r
Rscript get_splice_sites_R10.r
./write_distance_tables_splice.sh
./plot_distance_tables_splice.sh

exit 0
