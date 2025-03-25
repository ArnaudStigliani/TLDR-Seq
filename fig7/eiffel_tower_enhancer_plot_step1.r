rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(parallel)
library(ggpubr)
library(ggcorrplot)
library(parallel)
library(stringr)
library(genomation)


out_dir <- "./results/shared/eiffel_tower_plot_enhancer/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )


#### intersect beds with annotation to remove reads that overlap genode transcripts
reads.curated_RBM7.name <- "./results/mES-RBM7/get_primary_processed_bam/trimmed_primary_processed.bam"
reads.curated_ZCCHC8.name <- "./results/mES-ZCCHC8/get_primary_processed_bam/trimmed_primary_processed.bam"
DHS.name <- "../data/ENCFF048DWN_mm39.bed"


gencode_annot.name <- "../data/Mus_musculus.GRCm39.103.bed6"
reads.curated_noAnnot.name <- file.path(out_dir, "all_reads_noAnnot.bam")

command1 <- paste("module load bedtools; module load samtools; samtools merge",
                  reads.curated_ZCCHC8.name, reads.curated_RBM7.name,
                  " -o - | bedtools intersect -v -abam - ",
                  "-b",gencode_annot.name,
                  "| bedtools intersect -wa -abam - -b", DHS.name ,">",reads.curated_noAnnot.name )
system(command1)

reads.no.annot.list <- file.path(out_dir, "reads_no_annot.list")
command3 <- paste("module load samtools; samtools view", reads.curated_noAnnot.name,
                  "| awk '{print $1}' >", reads.no.annot.list)
system(command3)

read_no_annot.df  <- read.table(reads.no.annot.list)


####### concatenate reads
read_curated_bed.name <- file.path(out_dir, "reads_curated.bed")
      
command1 <- paste("module load bedtools; module load samtools; samtools merge",
                  reads.curated_ZCCHC8.name, reads.curated_RBM7.name,
                  " -o - | bedtools bamtobed -i - >", read_curated_bed.name)
system(command1)

    
####

format_bed <- function(bed, cov) 
{
    bed.formatted <- bed %>%
        mutate(start.tss = ifelse(strand=="+", start, end)) %>%
        mutate(seq=paste0(seqnames,":",  start.tss, ":", strand)) %>%
        dplyr::rename(chr=seqnames) %>% 
        dplyr::select( -end, -width, -score, -start) %>% 
        group_by(seq) %>%
        add_count() %>%
        ungroup %>%
        dplyr::filter(n >= cov )  %>%
        distinct(seq, .keep_all=TRUE) %>%
        dplyr::select(-seq) %>% 
        arrange(chr, start.tss)
    return(bed.formatted)
}


determine.windows <- function(df)
{
    ord <- order(df$n, decreasing=TRUE)
    ws <- 1000
    i <- 1
    all.inx <- rep(TRUE, length(ord))
    for (i in ord)
    {
        if(!all.inx[i])
        {
            next
        }
        pos1 <- df[i,]$start.tss - ws
        pos2 <- df[i,]$start.tss + ws
        a <- which(df$start.tss >= pos1 &  df$start.tss < pos2)
        all.inx[a] <- FALSE
        all.inx[i] <- TRUE
    }
    a <- df[all.inx, ]    
}


read_curated_bed.df <- readBed(read_curated_bed.name)  %>%
    as.data.frame %>%
    dplyr::filter(name %in% unlist(read_no_annot.df)) %>% 
    arrange(seqnames, start)


reads_formatted.df <- format_bed(read_curated_bed.df, 1)

    
reads_formatted_split.df <-  reads_formatted.df %>% split(f=.$chr) %>% .[grepl("chr[0-9XY]", names(.) )]
reads_formatted_filtered.df <- mclapply(reads_formatted_split.df, determine.windows, mc.cores=30) 
reads_formatted_filtered_rbind.df <- do.call(rbind.data.frame, reads_formatted_filtered.df)


reads_formatted_filtered_rbind_plus <- reads_formatted_filtered_rbind.df %>%
    dplyr::filter(strand=="+")
reads_formatted_filtered_rbind_minus <- reads_formatted_filtered_rbind.df %>%
    dplyr::filter(strand=="-")

###### write tss -2500:2500 bed files 
tss_window_plus.bed.name <- file.path(out_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(out_dir, "tss_minus_window.bed")

tss_window_plus.bed <- reads_formatted_filtered_rbind_plus %>%
    mutate(start = start.tss - 2500, end = start.tss + 2500 ) %>%
    dplyr::select(chr, start, end)
tss_window_minus.bed <- reads_formatted_filtered_rbind_minus %>%
    mutate(start = start.tss - 2500, end = start.tss + 2500 ) %>%
    dplyr::select(chr, start, end)
write.table(tss_window_plus.bed, tss_window_plus.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(tss_window_minus.bed, tss_window_minus.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
