rm(list=ls())
library(tidyverse)
library(parallel)
library(tailfindr)
library(rbokeh)



arg = commandArgs(trailingOnly=TRUE)
### arg[1] <- "Name of the  sequencing run"

in_dir <- file.path("../../results/TLDRSeq/demultiplex_reads_from_fast5/", arg[1])

out_dir <- file.path("../../results/TLDRSeq/get_polyA_tailfindr",  arg[1])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)




A2_bc1 <- "TGGCTTGATCCCTCATCTTGTGAAGTTGTTTCGGTCGATTCCGTTTGTAGTCGTCTGT"
A2_bc2 <- "ACTGGTGCAGCTTTGAACATCTAGAGAGGGTACTATGTGCCTCAGCACCGTACAGCAA"
A2_bc3 <- "ATGGACTTTGGTAACTTCCTGCGTCACCCACACTTACTTCAGGACGTACCAGTAGAAG"
A2_bc4 <- "GTTGAATGAGCCTACTGGGTCCTCTTCTGAAGTTCCTGGGTCTTGAACCAGACTTGGT"
A2_bc5 <- "TGAGAGACAAGATTGTTCGTGGACGACAGACACCGTTCATCGACTTTCGGACGAAGAA"
A2_bc6 <- "AGATTCAGACCGTCTCATGCAAAGTTCTCAGTCTTCCTCCAGACAAGGCTTACGAAGC"
A2_bc7 <- "CAAGAGCTTTGACTAAGGAGCATGCCGATCCTTGTGGCTTCTAACTTCATGTCCCAGT"
A2_bc8 <- "TGGAAGATGAGACCCTGATCTACGGTTTGTCATACTCGTGTGCTCACCCTAGACACCT"
A2_bc9 <- "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"
A2_bc10 <- "GCTAGGTCAATCTCCTTCGGAAGTTACAGTCCGAGCCTCATGTGATCTACAACGTCAT"
A2_bc11 <- "CAGGTTACTCCTCCGTGAGTCTGAACCGAGATCCTACGAATGGAGTGTGTTGGGTAAC"
A2_bc12 <- "TCAATCAAGAAGGGAAAGCAAGGTCCTGGGAGCATCAGGTAGTAACAGGGAGGAAACA"
A2_bc13 <- "CATGTTCAACCAAGGCTTCTATGGTAGCTGACTGTCTTCCATACCGACAAGAAGGCAC"

A2.c <- c(bc1 = A2_bc1,
          bc2 = A2_bc2,
          bc3 = A2_bc3,
          bc4 = A2_bc4,
          bc5 = A2_bc5,
          bc6 = A2_bc6,
          bc7 = A2_bc7,
          bc8 = A2_bc8,
          bc9 = A2_bc9,
          bc10 = A2_bc10,
          bc11 = A2_bc11,
          bc12 = A2_bc12, 
          bc13 = A2_bc13)

A1 <- "GACATACCGTGAATCGCTATCCTATGTTCATCCGCACAAAGACACCGACAACTTTCTT"



in_dir.ls <- list.files(in_dir)

for (elt in in_dir.ls)
{
    full.fast5.path <- file.path(in_dir, elt)
    full.out.path <- file.path(out_dir, elt)
    dir.create(full.out.path, showWarnings=FALSE, recursive=TRUE)
    A2.sel <- A2.c[names(A2.c)==elt]
    print(A2.sel)
    df <- find_tails(fast5_dir = full.fast5.path,
                     save_dir = full.out.path,
                     csv_filename = "cDNA_tails.csv",
                     num_cores = 40,
                     basecall_group="Basecall_1D_000", 
                     save_plots=FALSE,
                     plot_debug_traces=FALSE,
                     dna_datatype = 'custom-cdna',
                     front_primer=A2.sel,
                     end_primer=A1, 
                     )
    break
}
