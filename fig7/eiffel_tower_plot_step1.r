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


out_dir <- "../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )




### generate a bed file from bam

reads.curated_polyA_RBM7.name <- "../../results/sequencing_run_before_promethion2/mES-RBM7/get_primary_processed_bam/trimmed_primary_polyA_plus.bam"
reads.curated_polyA_ZCCHC8.name <- "../../results/sequencing_run_before_promethion2/mES-ZCCHC8/get_primary_processed_bam/trimmed_primary_polyA_plus.bam"
read_curated_bed.name <- file.path(out_dir, "trimmed_primary_polyA_plus.bed")

read_curated_polyA_RBM7_bed.name <- file.path(out_dir, "trimmed_primary_polyA_plus_RBM7.bed")
read_curated_polyA_ZCCHC8_bed.name <- file.path(out_dir, "trimmed_primary_polyA_plus_ZCCHC8.bed")

gencode_prot_coding_bed.name <- "../../../shared_data/Mus_musculus.GRCm39.103.protein_coding.bed6"
gencode_annot.name <- "../../../shared_data/Mus_musculus.GRCm39.103.bed6" ### choose this one for the end

command1 <- paste("module load bedtools;  bedtools bamtobed -i", reads.curated_polyA_RBM7.name,
                  "| bedtools intersect -a -  -b", gencode_prot_coding_bed.name, " -wa | uniq >",
                  read_curated_polyA_RBM7_bed.name)

command2 <- paste("module load bedtools;  bedtools bamtobed -i", reads.curated_polyA_ZCCHC8.name,
                  "| bedtools intersect -a -  -b", gencode_prot_coding_bed.name, " -wa | uniq >",
                  read_curated_polyA_ZCCHC8_bed.name)

system(paste(command1, " & \n", command2))


################################################################

   

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

## old code
## determine.windows <- function(df)
## {
##     w.size <- 3000
##     a <- lapply(0:(w.size-1), function(x, y) (y[["start.tss"]] + x) %/% w.size , df ) %>%
##         setNames(paste0("window.", 0:(w.size-1))) %>%
##         do.call(cbind.data.frame, .) %>%
##         cbind.data.frame(df, .)
##     for(i in 0:(w.size-1))
##     {
##         oldname <- paste0("window.",i)
##         a <- a %>%
##             mutate(w := get(!!oldname)) %>%
##             arrange(strand, chr, w, ifelse(strand=="+", start.tss,  desc(start.tss) )) %>% 
##             group_by(chr, w, strand) %>% 
##             dplyr::slice(1)
##     }
##     a  <- a %>% ungroup
##     return(a)
## }


### has to be run on each strand
determine.windows <- function(df)
{
    ord <- order(df$start.tss)
    ws <- 3000
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


read_curated_polyA_bed.df <- c(readBed(read_curated_polyA_RBM7_bed.name),
                               readBed(read_curated_polyA_ZCCHC8_bed.name)) %>%
    as.data.frame %>%
    arrange(seqnames, start)


reads_formatted.df <- format_bed(read_curated_polyA_bed.df, 10)



reads_formatted_split.df <-  reads_formatted.df %>% split(f=paste0(.$chr,";",.$strand)) %>% .[grepl("chr[0-9XY]", names(.) )]
reads_formatted_filtered.df <- mclapply(reads_formatted_split.df, determine.windows, mc.cores=60) 
reads_formatted_filtered_rbind.df <- do.call(rbind.data.frame, reads_formatted_filtered.df) 


reads_formatted_filtered_rbind_plus <- reads_formatted_filtered_rbind.df %>%
    dplyr::filter(strand=="+")
reads_formatted_filtered_rbind_minus <- reads_formatted_filtered_rbind.df %>%
    dplyr::filter(strand=="-")

###### write tss -1500:1500 bed files 
tss_window_plus_temp.bed.name <- file.path(out_dir, "tss_plus_window_temp.bed")
tss_window_minus_temp.bed.name  <- file.path(out_dir, "tss_minus_window_temp.bed")

tss_window_plus_temp.bed <- reads_formatted_filtered_rbind_plus %>%
    mutate(start = start.tss - 1500, end = start.tss  ) %>%
    dplyr::select(chr, start, end) %>%
    mutate(name = paste0(chr, ":",start, "-",end), score = ".", strand ="+")
tss_window_minus_temp.bed <- reads_formatted_filtered_rbind_minus %>%
    mutate(start = start.tss, end = start.tss + 1500 ) %>%
    dplyr::select(chr, start, end) %>% 
    mutate(name = paste0(chr, ":",start, "-",end), score = ".", strand ="-") 

write.table(tss_window_plus_temp.bed, tss_window_plus_temp.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(tss_window_minus_temp.bed, tss_window_minus_temp.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


### filter out protein codin TSS  on the antisense strand
### (actually, filter out stuff that overlaps with something in the same direction)

#### interlude
reads.curated_RBM7.name <- "../../results/sequencing_run_before_promethion2/mES-RBM7/get_primary_processed_bam/trimmed_primary_processed.bam"
reads.curated_ZCCHC8.name <- "../../results/sequencing_run_before_promethion2/mES-ZCCHC8/get_primary_processed_bam/trimmed_primary_processed.bam"
read_curated_bed.name <- file.path(out_dir, "trimmed_primary.bed")

command1 <- paste("module load bedtools; module load samtools; samtools merge",
                  reads.curated_ZCCHC8.name, reads.curated_RBM7.name,
                  " -o - | bedtools bamtobed -i - >", read_curated_bed.name)
system(command1)

#######

tss_window_plus_temp_filtered.bed.name <- file.path(out_dir, "tss_plus_window_temp_filtered.bed")
tss_window_minus_temp_filtered.bed.name  <- file.path(out_dir, "tss_minus_window_temp_filtered.bed")

command1 <- paste("module load bedtools; bedtools intersect  -S -a", tss_window_plus_temp.bed.name, "-b ",
                  gencode_annot.name, "| bedtools intersect -S -a - -b", read_curated_bed.name,
                  " |awk '{print $4}' | sort -u >", tss_window_plus_temp_filtered.bed.name )

command2 <- paste("module load bedtools; bedtools intersect  -S -a", tss_window_minus_temp.bed.name, "-b ",
                  gencode_annot.name, "| bedtools intersect -S -a - -b", read_curated_bed.name,
                  " |awk '{print $4}' | sort -u >", tss_window_minus_temp_filtered.bed.name )

system(paste(command1, " & \n", command2, " & wait"))


tss_window_plus.bed <-  read.table(tss_window_plus_temp.bed.name, sep="\t") %>%
    setNames(c("chr", "start", "end", "name", "score", "strand")) %>%
    mutate( end = end + 1500) %>%
    dplyr::filter(! name %in% read.table(tss_window_plus_temp_filtered.bed.name)$V1 ) %>% 
    dplyr::select(chr, start, end)


tss_window_minus.bed <-  read.table(tss_window_minus_temp.bed.name, sep="\t") %>%
    setNames(c("chr", "start", "end", "name", "score", "strand")) %>%
    mutate( start = start - 1500 ) %>%
    dplyr::filter(! name %in% read.table(tss_window_minus_temp_filtered.bed.name)$V1 ) %>% 
    dplyr::select(chr, start, end)

tss_window_plus.bed.name <- file.path(out_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(out_dir, "tss_minus_window.bed")

write.table(tss_window_plus.bed, tss_window_plus.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(tss_window_minus.bed, tss_window_minus.bed.name, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
