library(tidyverse)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(stringr)
library(data.table)
library(parallel)
library(plyranges)

out_dir <- "./results/get_splice_sites/"
dir.create(out_dir, showWarnings = FALSE, recursive=TRUE)


bam_tldr_lib3.name <- "./results/fig2/results/compare_lib1_lib3/R10.bam"
bam_tldr_lib3.bed12.name <- file.path(out_dir, "tldr_lib3.bed12")


command.line <- paste("module load bedtools;  bedtools bamtobed  -i", bam_tldr_lib3.name ," -bed12 > ", bam_tldr_lib3.bed12.name)
system(command.line) 
tldr_lib3.bed12.df <- read.table(bam_tldr_lib3.bed12.name, sep="\t") %>%
    dplyr::select(chr = V1, start.region = V2, strand=V6,  block.start = V12, block.size = V11, read.name =V4 )

   

#### TLDR lib3

tldr_lib3.splice.left <-  tldr_lib3.bed12.df$block.start %>%
    str_split(",")  %>%
    mclapply(as.numeric, mc.cores = 50) %>%
    mcmapply(function(x, y) x + y, . , tldr_lib3.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr_lib3.bed12.df$chr , .,
             tldr_lib3.bed12.df$strand,
             tldr_lib3.bed12.df$read.name,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-1], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") %>% 
    mutate(start = as.numeric(start), end = as.numeric(end) -1 ) 
tldr_lib3.splice.right <-  mcmapply(function(x, y) x + y,
                                  tldr_lib3.bed12.df$block.start %>% str_split(",") %>% mclapply(as.numeric,mc.cores = 50),
                                  tldr_lib3.bed12.df$block.size %>% str_split(",") %>% mclapply(as.numeric, mc.cores = 50),
                                  SIMPLIFY=FALSE, mc.cores = 50)  %>%
    mcmapply(function(x, y) x + y, . , tldr_lib3.bed12.df$start.region, SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr_lib3.bed12.df$chr , .,
             tldr_lib3.bed12.df$strand,
             tldr_lib3.bed12.df$read.name,             
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-length(x)], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") 

tldr_lib3.intron.bed <- cbind.data.frame(tldr_lib3.splice.left, tldr_lib3.splice.right) %>%
    setNames(paste0(names(.), c(".l", ".l", ".l", ".l", ".l", ".l", ".r", ".r", ".r", ".r", ".r", ".r"))) %>%
    dplyr::select(chr = chr.l, start = start.r, end = end.l, name = name.l, score = score.l, strand = strand.l ) %>%
    mutate(score=1) %>%
    mutate(start = start %>% as.numeric(), end  = end %>% as.numeric)

options(scipen=15)
write.table(tldr_lib3.intron.bed, file.path(out_dir, "tldr_lib3_intron.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


