rm(list=ls())
library(tidyverse)



arg <- commandArgs(trailingOnly = TRUE)
geno <- arg[1]


in_dir <- file.path("./results", geno, "get_polyA_tailfindR")
out_dir <- file.path("./results", geno, "characterize_polyA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
bc_polyA.name <- file.path(in_dir, "cDNA_tails.csv")


bc_polyA.df <- read.table(bc_polyA.name, sep=",", header=TRUE)  %>%
    dplyr::select(-file_path)



#### to generate bam files (with get_bam_polyA.sh)
reads_with_polyA <- bc_polyA.df  %>%
    dplyr::filter(!(tail_length == 0 |  is.na(tail_length))) %>%
    dplyr::filter(tail_length > 15)

write.table(reads_with_polyA %>% dplyr::select(read_id), file.path(out_dir, "reads_with_polyA.list"),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(reads_with_polyA  %>% dplyr::select(read_id, read_type, tail_length), file.path(out_dir, "reads_with_polyA.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

