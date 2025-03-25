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


in_dir <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot/")
out_dir <- in_dir

genos <- c("mES-WT1", "mES-RBM7", "mES-ZCCHC8")

cov.names <- file.path(in_dir, genos, "cov_to_plot.tsv")


cov1.df <- cov.names %>%
    mclapply(read.table, header=TRUE, sep="\t", mc.cores=3) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames %>% 
    mutate(polyA = polyA %>% factor(levels = c("pA-", "pA+"))) %>% 
    mutate( geno = factor(geno, levels = genos))



offset.df <- cov1.df %>%
    dplyr::filter(position < 0, direction =="reverse") %>%
    dplyr::filter(log.cov < 0) %>% 
    arrange( desc(position)) %>%
    group_by(region) %>%
    slice(1)  %>%
    arrange( desc(position)) %>%
    mutate(offset = round(position/2) )  %>%
    dplyr::filter(offset > -400) %>%
    mutate(group.var = region %>% str_replace(";.*","")) %>%
    .$group.var %>% 
    {cbind( paste0(., ";+"),  paste0(., ";-") )} %>%
    t %>%
    c



levels.y <- offset.df %>% rev

## levels.y <- read.table(file.path(in_dir, genos[2], "list_tss_sorted.tsv")) %>%
##     .$V1 %>%
##     .[. %in% offset.df]

cov.df <- cov1.df %>%
    dplyr::filter(region %in% levels.y) %>%
    mutate(region = factor(region, levels = levels.y))


my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="RdBu") %>% str_replace("F7F7F7", "FFFFFF"))(100) %>% rev
col_breaks <- c(seq(-1,-0.55, length=5),
                seq(-0.54,-0.35,length=20),
                seq(-0.34, 0,length=25),
                seq(0.3, 1, length=25),
                seq(1.1, 2,length=25))



g <- ggplot(cov.df , aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn(colours=my_palette, values=scales::rescale(col_breaks), limits = c(-1,2)) +
    theme_classic() +
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())  +
    facet_grid( vars(polyA), vars(geno))


out_eiffel <- file.path(out_dir, "eiffel.png")
ggsave(out_eiffel, g, width=12, height=8)
## out_eiffel <- file.path(out_dir, "eiffel.pdf")
## ggsave(out_eiffel, g, width=12, height=12)


### for reviewers
my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="RdBu") %>% str_replace("F7F7F7", "FFFFFF"))(100) %>% rev
col_breaks <- c(seq(-1,-0.10,length=5),
                seq(-0.09, 0,length=45),
                seq(0.3, 1, length=25),
                seq(1.1, 2,length=25))
### for reviewers
g <- ggplot(cov.df %>% dplyr::filter(geno =="mES-WT1") , aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn(colours=my_palette, values=scales::rescale(col_breaks), limits = c(-1,2)) +
    theme_classic() +
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())  +
    facet_grid( vars(polyA), vars(geno))
out_eiffel <- file.path(out_dir, "eiffel_referee_WT.png")
ggsave(out_eiffel, g, width=4, height=8)
out_eiffel <- file.path(out_dir, "eiffel_referee_WT.pdf")
ggsave(out_eiffel, g, width=4, height=8)




#### dnase1

dnase.tss_plus.bdg.name <- file.path(in_dir, "dnase1_tss_plus.bdg")
dnase.tss_minus.bdg.name <- file.path(in_dir, "dnase1_tss_minus.bdg")


eiffel_dnase1.df <- rbind.data.frame(
    read.table(dnase.tss_minus.bdg.name, sep="\t") %>%  mutate(strand= "minus"),
    read.table(dnase.tss_plus.bdg.name, sep="\t") %>%  mutate(strand= "plus")) %>% 
    setNames(c("chr", "a.start", "a.end", "chr.b","b.start","b.end", "cov", "strand")) %>%
    mutate(b.start = ifelse(a.start > b.start, a.start, b.start),
           b.end = ifelse(b.end > a.end, a.end,  b.end)) %>%
    mutate(region = paste0(chr, ":",a.start,"-", a.end, ";", strand)) %>%
    dplyr::select(region, b.start, b.end, cov) 
    
index <- mapply(rep, 1:nrow(eiffel_dnase1.df),
                eiffel_dnase1.df$b.end - eiffel_dnase1.df$b.start,
                SIMPLIFY=FALSE ) %>% unlist


eiffel_dnase1_complete.df  <- eiffel_dnase1.df[index,]  %>%
    group_by(region) %>%
    mutate(position = 1:n() - 1500) %>%
    separate(region, into = c("region.2","strand"), sep=";" ) %>%
    dplyr::rename(region=region.2) %>%
    mutate(position = ifelse(strand =="minus", -position, position)) %>%
    dplyr::filter(-position > -500) %>%
    mutate(log.cov = log10(1+cov) ) %>%
    mutate(log.cov = ifelse(log.cov > 1, 1, log.cov ))

my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="Greens"))(100)  
col_breaks <- c(seq(0, 0.4,length=50),
                seq(0.41, 1,length=50))
                

levels.y.2 = levels.y[grepl("\\+", levels.y)] %>% str_replace(";.*", "")


out_eiffel <- file.path(out_dir, "eiffel_dnase1.png")
g <- ggplot(eiffel_dnase1_complete.df %>%
            dplyr::filter(region %in% levels.y.2) %>% 
            mutate(region = factor(region, levels = levels.y.2)),
            aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn( colours= my_palette, values=scales::rescale(col_breaks), limits= c(0,1)) +
    theme_classic() +
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank()) 
ggsave(out_eiffel, g, width=3, height=4)


out_eiffel <- file.path(out_dir, "eiffel_dnase1.pdf")
ggsave(out_eiffel, g, width=4, height=6)

################# write TSS locations for 


tss_window_plus.bed.name <- file.path(in_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(in_dir, "tss_minus_window.bed")

tss_window_plus.bed <- read.table(tss_window_plus.bed.name, sep="\t") %>%
    mutate(group.var = paste0(V1,":", V2,"-", V3), strand = "+")
tss_window_minus.bed <- read.table(tss_window_minus.bed.name, sep="\t") %>%
        mutate(group.var = paste0(V1,":", V2,"-", V3), strand = "-")
tss_window.bed <- rbind.data.frame(tss_window_plus.bed, tss_window_minus.bed) %>%
    dplyr::filter(group.var %in% levels.y.2)


regions_plus.bed <- tss_window.bed %>%
    dplyr::filter(strand == "+") %>% 
    dplyr::select(-strand, -group.var)

regions_minus.bed <- tss_window.bed %>%
    dplyr::filter(strand == "-") %>%
    dplyr::select(-strand, -group.var)

peaks_for_stats_plus <- file.path(out_dir, "peaks_for_stats_plus.bed")
write.table(regions_plus.bed, peaks_for_stats_plus, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
peaks_for_stats_minus <- file.path(out_dir, "peaks_for_stats_minus.bed")
write.table(regions_minus.bed, peaks_for_stats_minus, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


#### IGV browser

RoI.df <- cov1.df %>%
    dplyr::filter(position < 0, direction =="reverse") %>%
    dplyr::filter(log.cov <= -1) %>% 
    arrange( desc(position)) %>%
    group_by(region) %>%
    slice(1)  %>%
    arrange( desc(position)) %>%
    mutate(offset = round(position/2) )  %>%
    dplyr::filter(offset > -400) %>%
    mutate(group.var = region %>% str_replace(";.*","")) %>%
    .$group.var %>% 
    {cbind( paste0(., ";+"),  paste0(., ";-") )} %>%
    t %>%
    c
