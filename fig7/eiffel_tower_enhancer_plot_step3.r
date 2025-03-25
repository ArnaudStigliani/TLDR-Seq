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
## geno <- "mES-RBM7"

in_dir <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot_enhancer/", geno)
out_dir <- in_dir




tss_window_plus_forward.bdg.name <-   file.path(in_dir, "tss_plus_window_forward.bdg")
tss_window_plus_reverse.bdg.name <-   file.path(in_dir, "tss_plus_window_reverse.bdg")
tss_window_minus_forward.bdg.name <-   file.path(in_dir, "tss_minus_window_forward.bdg")
tss_window_minus_reverse.bdg.name <-   file.path(in_dir, "tss_minus_window_reverse.bdg")


eiffel_forward.df <- rbind.data.frame(
    read.table(tss_window_minus_forward.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_forward.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(log.cov=log10(1+cov) %>% ifelse( .  > 2, 2, .)) %>%
    mutate(direction="forward") 


filter_out.forward <- eiffel_forward.df %>%
    dplyr::filter(position > 0) %>%
    dplyr::filter(position < 100) %>%    
    dplyr::filter(log.cov > 0.3) %>% 
    group_by(region) %>%
    slice(1)  %>%
    .$region

eiffel_reverse.df <- rbind.data.frame(
    read.table(tss_window_minus_reverse.bdg.name, sep="\t") %>%  mutate(V2 = V2 * -1),
    read.table(tss_window_plus_reverse.bdg.name, sep="\t")) %>%     
    setNames(c("region", "position", "cov")) %>%
    mutate(log.cov=log10(1+cov) %>% ifelse( .  > 2, 2, .)) %>%
    mutate(direction="reverse") 

levels.y <- eiffel_reverse.df  %>% dplyr::filter(region %in% filter_out.forward, position < 0) %>%
    dplyr::filter(log.cov > 0.3) %>% 
    arrange( desc(position)) %>%
    group_by(region) %>%
    slice(1)  %>%
    arrange(desc(position)) %>%
    .$region %>%
    rev %>%
    {cbind( paste0(., ";+"),  paste0(., ";-") )} %>%
    t %>%
    c
    

offset.df <- eiffel_reverse.df  %>% dplyr::filter(region %in% filter_out.forward, position < 0) %>%
    dplyr::filter(log.cov > 0.3) %>% 
    arrange( desc(position)) %>%
    group_by(region) %>%
    slice(1) %>%
    mutate(offset = round(position/2) ) %>% 
    dplyr::select(group.var = region, offset)



eiffel_all.df <- rbind.data.frame(eiffel_forward.df, eiffel_reverse.df) %>%
    mutate(group.var = region) %>%
    mutate(region = region %>% paste0(ifelse(direction =="forward", ";+", ";-"))) %>%
    mutate(region = region %>% factor(., levels = levels.y )) %>%
    mutate(log.cov = ifelse(direction =="reverse",-log.cov, log.cov)) %>%
    dplyr::filter(group.var %in% offset.df$group.var)


write.table(data.frame(levels.y), file.path(out_dir, "list_tss_sorted.tsv"), quote =FALSE, row.names=FALSE, col.names=FALSE )

all.df <- eiffel_all.df
write.table(all.df, file.path(out_dir, "cov_to_plot.tsv"), sep = "\t", quote = FALSE, row.names=FALSE )


