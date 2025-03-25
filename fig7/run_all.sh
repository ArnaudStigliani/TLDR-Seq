#!/bin/bash

./pipeline_map.sh mES-ZCCHC8 &
./pipeline_map.sh mES-WT     &
./pipeline_map.sh mES-RBM7   &

Rscript get_polyA_tailfindr.r


