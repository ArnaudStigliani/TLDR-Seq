#!/bin/bash

# Rscript eiffel_tower_plot_step1.r


# Rscript eiffel_tower_plot_step2.r mES-WT1  &
# Rscript eiffel_tower_plot_step2.r mES-ZCCHC8 & 
# Rscript eiffel_tower_plot_step2.r  mES-RBM7 dnase1 &
# wait 


Rscript eiffel_tower_plot_step3.r mES-WT1 &
Rscript eiffel_tower_plot_step3.r mES-ZCCHC8 &
Rscript eiffel_tower_plot_step3.r mES-RBM7 &
wait
