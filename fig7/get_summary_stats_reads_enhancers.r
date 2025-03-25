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
library(ggh4x)

in_dir <- "./results/shared/eiffel_tower_plot_enhancer/"

genos <- c("mES-WT1", "mES-RBM7", "mES-ZCCHC8")

for (geno in genos)
{
    out_dir <- file.path("./results/shared/summary_stats_enhancers", geno)
    eiffel_dir  <-  file.path(in_dir, geno)
    ##
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )
    ##
   
### bam files
    reads.curated_noAnnot.name <- file.path(in_dir, geno, "all_reads_noAnnot.bam")

    ## TSS
    tss_window_all.bed.name <- file.path(in_dir, "peaks_for_stats.bed")


###### get reads
    tss_window.bed12.name <-   file.path(out_dir, "tss_window.bed12")


    command.line.1 <- paste("module load samtools; module load bedtools;",
                            " samtools view -Sb ", reads.curated_noAnnot.name,
                            "| bedtools intersect -abam - -b", tss_window_all.bed.name, 
                            "| bedtools bamtobed -bed12 -i - >", 
                            tss_window.bed12.name)


    system(command.line.1)

    tss_window_wawb.bed12.name <-   file.path(out_dir, "tss_window_wawb.bed12")

    command.line.1 <- paste(" module load bedtools;",
                            "bedtools intersect -wa -wb -a",
                            tss_window.bed12.name," -b", tss_window_all.bed.name, ">", 
                            tss_window_wawb.bed12.name)


    system(command.line.1)

##### analyze data
}

for (geno in genos)
{
    ## geno <- "mES-RBM7"
    out_dir <- file.path("./results/shared/summary_stats_enhancers", geno)
    tss_window_wawb.bed12.name <-   file.path(out_dir, "tss_window_wawb.bed12")

    enhancers.df <- read.table(tss_window_wawb.bed12.name, sep="\t")

    enhancers.summarized.df  <-  enhancers.df %>%
        dplyr::select( nblock = V10) %>%
        mutate(read_type =  ifelse(nblock > 1, "spliced", "unspliced")) %>%
        group_by( read_type) %>%
        dplyr::summarize(n = n()) %>%
        mutate(nreads = sum(n)) %>%
        mutate(fraction = n/nreads, geno = geno) 
    write.table(enhancers.summarized.df, file.path(out_dir, "summary_table.tsv"), sep="\t", quote =FALSE, row.names=FALSE )

    enhancers.summarized_detail.df  <-  enhancers.df %>%
        dplyr::select(nintrons = V10) %>%
        mutate(nintrons = nintrons - 1 ) %>% 
        group_by(nintrons) %>%
        dplyr::summarize(n = n()) %>%
        mutate(nreads = sum(n)) %>%
        mutate(fraction = n/nreads, geno = geno) 
    write.table(enhancers.summarized_detail.df, file.path(out_dir, "summary_table_detail.tsv"), sep="\t", quote =FALSE, row.names=FALSE )
    
    enhancers.summarized_length.df  <-  enhancers.df %>%
        mutate(TSS=paste0(V13,":",V14,"-",V15), genomic.length = V3-V2) %>% 
        dplyr::select( blocks = V11, TSS, genomic.length) %>%
        mutate(length = blocks %>% str_split(",") %>% lapply(as.numeric) %>% sapply(sum, na.rm =TRUE)) %>% 
        group_by(TSS) %>%
        dplyr::summarize(average_length = mean(length), average_genomic_length = mean(genomic.length)) %>%
        ungroup %>%
        mutate(geno = geno)
        write.table(enhancers.summarized_length.df, file.path(out_dir, "summary_table_length.tsv"), sep="\t", quote =FALSE, row.names=FALSE )
    
}

#### plot tables

all.summary.df <- file.path(dirname(out_dir), genos, "summary_table.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames()

g <- ggplot(all.summary.df, aes(x=read_type, y = fraction, fill = geno)) +
    geom_bar(stat = "identity", position = "dodge") +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic() 
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types.pdf")
ggsave(out_eiffel, g, width=5, height=4)


all.summary_detail.df <- file.path(dirname(out_dir), genos, "summary_table_detail.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    dplyr::filter(nintrons < 11) %>%
    mutate(nexons = nintrons + 1) %>% 
    mutate(nexons = factor(nexons, c(1:11) %>% as.character )) %>%
    drop_na 

g <- ggplot(all.summary_detail.df , aes(x=nexons, y = fraction, fill = geno)) +
    geom_bar(stat = "identity", position = "dodge") +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic()
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types_detail.pdf")
ggsave(out_eiffel, g, width=5, height=4)



all.summary_length.df <- file.path(dirname(out_dir), genos, "summary_table_length.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    drop_na %>%
    dplyr::rename(average_read_length = average_length) %>%
    reshape2::melt() %>%
    dplyr::filter(!(variable =="average_genomic_length" & value > 10000))

g <- ggplot(all.summary_length.df , aes(colour = geno, x = value)) +
    geom_density() +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2")[1:3])+
    theme_classic() +
    facet_grid(vars(variable), scales = "free_y") 
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types_length.pdf")
ggsave(out_eiffel, g, width=6, height=4)

g <- ggplot(all.summary_length.df , aes(colour = geno, x = value, y=..ndensity..)) +
    geom_density() +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2")[1:3])+
    theme_classic() +
    facet_grid( vars(variable), scales = "free_y") 
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types_length_ndensity.pdf")
ggsave(out_eiffel, g, width=6, height=4)
