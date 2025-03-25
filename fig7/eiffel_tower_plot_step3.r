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

arg = commandArgs(trailingOnly=TRUE)

geno <- arg[1]
geno <- "mES-WT1"

in_dir <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot/", geno)
out_dir <- in_dir

reads_bed.name <- file.path(in_dir, "trimmed_primary.bed")
reads_bed.range <- readBed(reads_bed.name)

## reads_bed2.name <- file.path(in_dir, "trimmed_primary.bed")
## reads_bed2.range <- readBed(reads_bed2.name)


tss_window_plus_forward.bdg.name <-   file.path(in_dir, "tss_plus_window_forward.bdg")
tss_window_plus_reverse.bdg.name <-   file.path(in_dir, "tss_plus_window_reverse.bdg")
tss_window_minus_forward.bdg.name <-   file.path(in_dir, "tss_minus_window_forward.bdg")
tss_window_minus_reverse.bdg.name <-   file.path(in_dir, "tss_minus_window_reverse.bdg")


tss_window_plus_forward_polyA.bdg.name <-   file.path(in_dir, "tss_plus_window_forward_polyA.bdg")
tss_window_plus_reverse_polyA.bdg.name <-   file.path(in_dir, "tss_plus_window_reverse_polyA.bdg")
tss_window_minus_forward_polyA.bdg.name <-   file.path(in_dir, "tss_minus_window_forward_polyA.bdg")
tss_window_minus_reverse_polyA.bdg.name <-   file.path(in_dir, "tss_minus_window_reverse_polyA.bdg")


#### all reads

norm.factor <- reads_bed.range %>%
    as.data.frame %>%
    summarize(n = sum(end - start)) %>%
    unlist %>%
    {1/.} %>%
    { . * 10^9}



eiffel_forward.df <- rbind.data.frame(
    read.table(tss_window_minus_forward.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_forward.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(cov = cov * norm.factor) %>% 
    mutate(log.cov=log10(1+cov) %>% ifelse( .  > 2, 2, .)) %>%
    mutate(direction="forward") 

eiffel_reverse.df <- rbind.data.frame(
    read.table(tss_window_minus_reverse.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_reverse.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(cov = cov * norm.factor) %>% 
    mutate(log.cov=log10(1+cov) %>% ifelse( .  > 1, 1, .)) %>%
    mutate(direction="reverse") 

levels.y <- eiffel_reverse.df  %>% dplyr::filter(position < 0) %>%
    dplyr::filter(log.cov > 0) %>% 
    arrange( desc(position)) %>%
    group_by(region) %>%
    slice(1)  %>%
    arrange(desc(position)) %>%
    .$region %>%
    c(eiffel_reverse.df  %>% dplyr::filter(position < 0) %>%
      group_by(region) %>% dplyr::summarize(m = max(log.cov)) %>% dplyr::filter(m ==0) %>% .$region) %>%
    rev %>%
    {cbind( paste0(., ";+"),  paste0(., ";-") )} %>%
    t %>%
    c
    


eiffel_all.df <- rbind.data.frame(eiffel_forward.df, eiffel_reverse.df) %>%
    mutate(region = region %>% paste0(ifelse(direction =="forward", ";+", ";-"))) %>%
    dplyr::filter( -position > -500)  %>% 
    mutate(region = region %>% factor(., levels = levels.y )) %>%
    mutate(log.cov = ifelse(direction =="reverse",-log.cov, log.cov))


##### pA + 


eiffel_forward_polyA.df <- rbind.data.frame(
    read.table(tss_window_minus_forward_polyA.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_forward_polyA.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(cov = cov * norm.factor) %>% 
    mutate(direction="forward") 

eiffel_reverse_polyA.df <- rbind.data.frame(
    read.table(tss_window_minus_reverse_polyA.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_reverse_polyA.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(cov = cov * norm.factor) %>% 
    mutate(direction="reverse") 



eiffel_all_polyA.df <- rbind.data.frame(eiffel_forward_polyA.df, eiffel_reverse_polyA.df) %>%
    dplyr::filter(-position > -500) %>% 
    mutate(log.cov = log10(1+cov) %>% ifelse( .  > 2, 2, .)) %>%
    mutate(log.cov = ifelse(log.cov >= 1 & direction =="reverse", 1, log.cov)) %>% 
    mutate(region = region %>% paste0(ifelse(direction =="forward", ";+", ";-"))) %>%
    mutate(region = region %>% factor(., levels = levels.y )) %>%
    mutate(log.cov = ifelse(direction =="reverse",-log.cov, log.cov))




write.table(data.frame(levels.y), file.path(out_dir, "list_tss_sorted.tsv"), quote =FALSE, row.names=FALSE, col.names=FALSE )

all.df <- rbind.data.frame(eiffel_all_polyA.df %>% mutate(polyA = "pA+", geno = geno),
                           eiffel_all.df %>% mutate(polyA = "pA-", geno = geno))
write.table(all.df, file.path(out_dir, "cov_to_plot.tsv"), sep = "\t", quote = FALSE, row.names=FALSE )

