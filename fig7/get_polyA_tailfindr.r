rm(list=ls())
library(tidyverse)
library(parallel)
library(tailfindr)
library(rbokeh)


## WT

in_dir <- "../data/WT_fast5out"
out_dir <-  file.path("../../results/sequencing_run_before_promethion2", "mES-WT", "get_polyA_tailfindR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


A2_bc9 <- "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"

A2.c <- c(bc9 = A2_bc9)
A1 <- "GACATACCGTGAATCGCTATCCTATGTTCATCCGCACAAAGACACCGACAACTTTCTT"


df <- find_tails(fast5_dir = in_dir,
                 save_dir = out_dir,
                 csv_filename = "cDNA_tails.csv",
                 num_cores = 40,
                 basecall_group="Basecall_1D_000", 
                 save_plots=FALSE,
                 plot_debug_traces=FALSE,
                 dna_datatype = 'custom-cdna',
                 front_primer=A2.c,
                 end_primer=A1,
                 )

## RBM7


in_dir <- "../data/RBM7_fast5out"
out_dir <-  file.path("../../results/sequencing_run_before_promethion2", "mES-RBM7", "get_polyA_tailfindR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


A2_bc9 <- "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"

A2.c <- c(bc9 = A2_bc9)
A1 <- "GACATACCGTGAATCGCTATCCTATGTTCATCCGCACAAAGACACCGACAACTTTCTT"


df <- find_tails(fast5_dir = in_dir,
                 save_dir = out_dir,
                 csv_filename = "cDNA_tails.csv",
                 num_cores = 40,
                 basecall_group="Basecall_1D_000", 
                 save_plots=FALSE,
                 plot_debug_traces=FALSE,
                 dna_datatype = 'custom-cdna',
                 front_primer=A2.c,
                 end_primer=A1,
                 )

## ZCCHC8


in_dir <- "../data/ZCCHC8_fast5out"
out_dir <-  file.path("../../results/sequencing_run_before_promethion2", "mES-ZCCHC8", "get_polyA_tailfindR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


A2_bc9 <- "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"

A2.c <- c(bc9 = A2_bc9)
A1 <- "GACATACCGTGAATCGCTATCCTATGTTCATCCGCACAAAGACACCGACAACTTTCTT"


df <- find_tails(fast5_dir = in_dir,
                 save_dir = out_dir,
                 csv_filename = "cDNA_tails.csv",
                 num_cores = 40,
                 basecall_group="Basecall_1D_000", 
                 save_plots=FALSE,
                 plot_debug_traces=FALSE,
                 dna_datatype = 'custom-cdna',
                 front_primer=A2.c,
                 end_primer=A1,
                 )
