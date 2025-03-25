#!/bin/bash

./pipeline_map.sh mES-ZCCHC8 
./pipeline_map.sh mES-WT     
./pipeline_map.sh mES-RBM7   

Rscript get_polyA_tailfindr.r

./get_processed_bam_files.sh mES-ZCCHC8
./get_processed_bam_files.sh mES-WT
./get_processed_bam_files.sh mES-RBM7

Rscript characterize_polyA.r mES-ZCCHC8
Rscript characterize_polyA.r mES-WT
Rscript characterize_polyA.r mES-RBM7

### make eiffel tower plots

Rscript eiffel_tower_plot_step1.r

Rscript eiffel_tower_plot_step2.r mES-WT1  
Rscript eiffel_tower_plot_step2.r mES-ZCCHC8  
Rscript eiffel_tower_plot_step2.r  mES-RBM7 dnase1 

Rscript eiffel_tower_plot_step3.r mES-WT1 
Rscript eiffel_tower_plot_step3.r mES-ZCCHC8 
Rscript eiffel_tower_plot_step3.r mES-RBM7 

## make the enhancer one

Rscript eiffel_tower_enhancer_plot_step1.r


Rscript eiffel_tower_enhancer_plot_step2.r mES-WT1 
Rscript eiffel_tower_enhancer_plot_step2.r mES-ZCCHC8  
Rscript eiffel_tower_enhancer_plot_step2.r mES-RBM7 dnase1 


Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-WT1 
Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-ZCCHC8 
Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-RBM7 

Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-WT1 
Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-ZCCHC8 
Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-RBM7 

Rscript eiffel_tower_enhancer_plot_step4.r 

### read statistics

get_summary_stats_reads.r
get_summary_stats_reads_both.r
get_summary_stats_reads_enhancers.r

exit 0
