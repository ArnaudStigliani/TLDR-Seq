#!/bin/bash

# Rscript eiffel_tower_enhancer_plot_step1.r


Rscript eiffel_tower_enhancer_plot_step2.r mES-WT1 &
Rscript eiffel_tower_enhancer_plot_step2.r mES-ZCCHC8 & 
Rscript eiffel_tower_enhancer_plot_step2.r mES-RBM7 dnase1 &
wait 


Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-WT1 &
Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-ZCCHC8 &
Rscript eiffel_tower_enhancer_plot_step3_pAp.r mES-RBM7 &

Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-WT1 &
Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-ZCCHC8 &
Rscript eiffel_tower_enhancer_plot_step3_pAm.r mES-RBM7 &
wait

# Rscript eiffel_tower_enhancer_plot_step4.r 
