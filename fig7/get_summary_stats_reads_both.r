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

in_dir <- "./results/shared/"
enhancer.dir <- file.path(in_dir, "summary_stats_enhancers")
prompt.dir <- file.path(in_dir, "summary_stats")
genos <- c("mES-WT1", "mES-RBM7", "mES-ZCCHC8")
out_dir <- file.path(in_dir, "get_summary_stats_reads_both")
dir.create(out_dir, showWarnings=FALSE)

#### enhancers #####
all.summary.df <- file.path(enhancer.dir, genos, "summary_table.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    mutate(transcript_type = "enhancer") %>%
    rbind.data.frame(
        file.path(prompt.dir, genos, "summary_table.tsv") %>%
        lapply(read.table, sep="\t", header=TRUE) %>%
        do.call(rbind.data.frame, .))%>%
    remove_rownames


all.summary_detail.df <- file.path(enhancer.dir, genos, "summary_table_detail.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    dplyr::filter(nintrons < 11) %>%
    mutate(nexons = nintrons + 1) %>% 
    mutate(nexons = factor(nexons, c(1:11) %>% as.character )) %>%
    drop_na  %>%
    mutate(transcript_type = "enhancer") %>%
    rbind.data.frame(
        file.path(prompt.dir, genos, "summary_table_detail.tsv") %>%
        lapply(read.table, sep="\t", header=TRUE) %>%
        do.call(rbind.data.frame, .) %>%
        remove_rownames() %>%
        dplyr::filter(nintrons < 11) %>%
        mutate(nexons = nintrons + 1) %>% 
        mutate(nexons = factor(nexons, c(1:11) %>% as.character )) %>%
        drop_na ) 


all.summary_length.df <- file.path(enhancer.dir, genos, "summary_table_length.tsv") %>%
    lapply(read.table, sep="\t", header=TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    drop_na %>%
    dplyr::rename(average_read_length = average_length) %>%
    reshape2::melt() %>%
    dplyr::filter(!(variable =="average_genomic_length" & value > 10000)) %>%
    mutate(transcript_type = "enhancer") %>%    
    rbind.data.frame(file.path(prompt.dir, genos, "summary_table_length.tsv") %>%
                     lapply(read.table, sep="\t", header=TRUE) %>%
                     do.call(rbind.data.frame, .) %>%
                     remove_rownames() %>%
                     drop_na %>%
                     dplyr::rename(average_read_length = average_length) %>%
                     reshape2::melt() %>%
                     dplyr::filter(!(variable =="average_genomic_length" & value > 10000))) %>%
    remove_rownames()


#### plot stuff

g <- ggplot(all.summary.df, aes(x=read_type, y = fraction, fill = geno)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic() +
    facet_grid(vars(transcript_type))
out_eiffel <- file.path(out_dir,  "summary_read_types.pdf")
ggsave(out_eiffel, g, width=5, height=4)


g <- ggplot(all.summary_detail.df , aes(x=nexons, y = fraction, fill = geno)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic() +
    facet_grid(vars(transcript_type), scales = "free_x") ## otherwise, use faacet_grid2 with independand ="x"
out_eiffel <- file.path(out_dir,  "summary_read_types_detail.pdf")
ggsave(out_eiffel, g, width=5, height=4)

g <- ggplot(all.summary_length.df , aes(colour = transcript_type, x = value)) +
    geom_density() +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic() +
    facet_grid(vars(geno), vars(variable), scales = "free_x") 
out_eiffel <- file.path(out_dir,  "summary_read_types_length.pdf")
ggsave(out_eiffel, g, width=6, height=4)

g <- ggplot(all.summary_length.df , aes(colour = transcript_type, x = value, y=..ndensity..)) +
    geom_density() +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2"))+
    theme_classic() +
    facet_grid(vars(geno), vars(variable), scales = "free_x") 
out_eiffel <- file.path(out_dir,  "summary_read_types_length_ndensity.pdf")
ggsave(out_eiffel, g, width=6, height=4)




