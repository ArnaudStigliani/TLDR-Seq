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

in_dir <- file.path("./results/shared/eiffel_tower_plot_enhancer/")
out_dir <- in_dir

genos <- c("mES-WT1", "mES-RBM7", "mES-ZCCHC8")

cov.names <- c(file.path(in_dir, genos, "cov_to_plot_pAp.tsv"),
               file.path(in_dir, genos, "cov_to_plot_pAm.tsv"))

#### norm factor
reads_bed.name <- file.path(in_dir,genos, "trimmed_primary.bed")


norm.factor <- mclapply( reads_bed.name,readBed, mc.cores = 3) %>%
    lapply(as.data.frame, mc.cores=3) %>% 
    mapply(mutate,., geno=genos,  SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame,.) %>%
    group_by(geno) %>% 
    summarize(n = sum(end - start)) %>%
    mutate(f = 1/n * 10^10 )


cov1.df <- cov.names %>%
    mclapply(read.table, header=TRUE, sep="\t", mc.cores=6) %>%
    mapply(mutate,., geno=genos, polyA = c("plus", "plus", "plus", "minus", "minus", "minus"), SIMPLIFY=FALSE) %>%  
    do.call(rbind.data.frame, .) %>%
    remove_rownames %>% 
    mutate( geno = factor(geno, levels = genos[c(2,1,3)]))

offset.df <- cov1.df  %>%
    dplyr::filter(direction=="reverse") %>%
    dplyr::filter( position <= 0, position > - 600) %>%
    group_by(region, position, direction, group.var, geno) %>%
    dplyr::summarize(cov = sum(cov)) %>%
    mutate(log.cov = ifelse(direction =="reverse", -log(1+cov, 10), log(1+cov, 10))) %>%
    arrange( group.var, geno,  desc(position), desc(abs(log.cov))) %>%
    mutate(index = ifelse(position == 0, "a", NA)) %>%
    group_by(group.var, geno) %>%
    mutate(cov0 = first(log.cov)) %>%
    ## dplyr::filter(group.var =="chr17:10649281-10653281") %>%
    ## dplyr::filter(group.var =="chr10:110688799-110692799") %>%
    mutate(index = ifelse(position != 0 & log.cov == (log.cov %>% min) & log.cov < cov0, "c",
                   ifelse(position != 0 & log.cov < cov0, "b",
                          index))) %>% 
    ungroup %>% 
    drop_na %>%
    distinct(group.var, geno, index, .keep_all=TRUE) %>% 
    {left_join(reshape2::dcast(., group.var + geno ~ index, value.var ="log.cov"),
               reshape2::dcast(., group.var + geno ~ index, value.var ="position"), by=c("group.var","geno"))} %>%
    drop_na() %>%
    dplyr::filter(2* abs(a.x) < abs(c.x)) %>% 
    arrange(geno) %>%
    group_by(group.var) %>%
    slice(1)  %>%
    arrange( desc(b.y)) %>%
    mutate(offset = round(b.y/2, 0)) %>%
    dplyr::select(group.var, offset) %>%
    dplyr::filter(abs(offset) > 50, abs(offset) < 300)


levels.y <- offset.df$group.var %>%
    rev %>%
    {cbind( paste0(., ";+"),  paste0(., ";-") )} %>%
    t %>%
    c

    
cov.df <- cov1.df  %>%
    ## dplyr::filter(group.var =="chr10:110688299-110693299") %>%
    left_join(norm.factor %>% dplyr::select(geno, f)) %>%
    mutate(log.cov =  ifelse(direction =="reverse", -log(1+(cov *f), 10), log(1+(cov * f) , 10))) %>% 
    mutate(region = factor(region, levels = levels.y)) %>%
    mutate( geno = factor(geno, levels = genos), polyA=factor(polyA, levels = c( "minus","plus"))) %>%
    drop_na %>%
    mutate(log.cov = ifelse(abs(log.cov) >= 1, sign(log.cov) * 1, log.cov)) %>%
    left_join(offset.df) %>%
    mutate(position = position - offset)  %>%
    dplyr::filter(position >-2000, position < 2000)


##### for browser
RoI <- cov.df %>% dplyr::filter(geno =="mES-RBM7") %>%
    group_by(offset, group.var) %>%
    dplyr::summarize(min.cov = min(log.cov), max.cov = max(log.cov)) %>%
    dplyr::filter(offset < -150) %>%
    arrange(min.cov) 
    


    

########################### dnase1 ##################
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
    ## mutate(test = paste0(chr, ":",a.start,"-", a.end)) %>%
    dplyr::select(region, b.start, b.end, cov) 
    
index <- mapply(rep, 1:nrow(eiffel_dnase1.df),
                eiffel_dnase1.df$b.end - eiffel_dnase1.df$b.start,
                SIMPLIFY=FALSE ) %>% unlist




eiffel_dnase1_complete.df  <- eiffel_dnase1.df[index,]  %>%
    mutate(group.var = region %>% str_replace(";.*","")) %>%
    group_by(region) %>%
    mutate(position = 1:n() - 2500) %>%
    ungroup %>% 
    dplyr::filter(group.var %in% unique(cov.df$group.var) ) %>%
    separate(region, into = c("region.2","strand"), sep=";" ) %>%
    dplyr::rename(region=region.2) %>%
    mutate(position = ifelse(strand =="minus", -position, position)) %>%
    mutate(log.cov = log10(1+cov) ) %>%
    mutate(log.cov = ifelse(log.cov > 0.5, 0.5, log.cov )) %>%
    left_join(offset.df) %>%
    dplyr::filter(abs(offset) < 300) %>% 
    mutate(position = position - offset) %>%
    dplyr::filter(position > - 2000, position < 2000) 


my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="Greens"))(100)  
col_breaks <- c(seq(0, 0.2,length=50),
                seq(0.21, 0.5,length=50))
                
levels.y.2 = levels.y[grepl("\\+", levels.y)] %>% str_replace(";.*", "")


out_eiffel <- file.path(out_dir, "eiffel_dnase1_enhancers.png")
g <- ggplot(eiffel_dnase1_complete.df %>%
            mutate(region = factor(region, levels = levels.y.2)),
            aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn( colours= my_palette, values=scales::rescale(col_breaks), limits= c(0,0.5)) +
    theme_classic() +
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank()) 
ggsave(out_eiffel, g, width=5, height=7)




#################################### enhancers ################################

my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="RdBu") %>% str_replace("F7F7F7", "FFFFFF"))(40)  %>% rev
col_breaks <- c(seq(-1,-0.1,length=20),
                ## seq(-0.69, 0.69 ,length=20),
                seq(0.1, 1,length=20))

out_eiffel <- file.path(out_dir, "eiffel.png")
g <- ggplot(cov.df %>%
          dplyr::filter(group.var %in% eiffel_dnase1_complete.df$group.var), aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn(colours=my_palette, values=scales::rescale(col_breaks), limits = c(-1,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())  +
    facet_grid( vars(polyA), vars(geno))
ggsave(out_eiffel, g, width=6, height=4)



### for reviewers

my_palette <- colorRampPalette(RColorBrewer::brewer.pal(n=5,name="RdBu") %>% str_replace("F7F7F7", "FFFFFF"))(40)  %>% rev
col_breaks <- c(seq(-1,-0.1,length=20),
                ## seq(-0.69, 0.69 ,length=20),
                seq(0.1, 1,length=20))

out_eiffel <- file.path(out_dir, "eiffel_enhancer_referee.png")
g <- ggplot(cov.df             %>%
            dplyr::filter(geno =="mES-WT1") %>% 
            dplyr::filter(group.var %in% eiffel_dnase1_complete.df$group.var), aes(x=position, y=region, fill=log.cov)) +
    geom_tile() +
    scale_fill_gradientn(colours=my_palette, values=scales::rescale(col_breaks), limits = c(-1,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())  +
    facet_grid( vars(polyA), vars(geno))
ggsave(out_eiffel, g, width=4, height=8)
out_eiffel <- file.path(out_dir, "eiffel_enhancer_referee.pdf")
ggsave(out_eiffel, g, width=4, height=8)


## count the number of enhancers in WT
cov.df %>%
    mutate(region == region %>% as.character) %>% 
    dplyr::filter(geno =="mES-WT1") %>%
    dplyr::filter(cov !=0) %>%
    distinct(group.var, region) %>%
    group_by(group.var) %>% 
    add_count %>%
    dplyr::filter(n == 1)  %>%
    ungroup %>%
    dplyr::filter(n == 2) %>%
    .$group.var %>% unique %>% length
    

## for IGV
list.enhancer <- cov.df             %>%
    dplyr::filter(geno =="mES-WT1") %>% 
    dplyr::filter(group.var %in% eiffel_dnase1_complete.df$group.var) %>%
    arrange( log.cov, region) %>%
    distinct(region, .keep_all=TRUE) %>%
    arrange(desc(cov))





################# write TSS locations for 

tss_window_plus.bed.name <- file.path(in_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(in_dir, "tss_minus_window.bed")

tss_window_plus.bed <- read.table(tss_window_plus.bed.name, sep="\t") %>%
    mutate(group.var = paste0(V1,":", V2,"-", V3), strand = "+")
tss_window_minus.bed <- read.table(tss_window_minus.bed.name, sep="\t") %>%
        mutate(group.var = paste0(V1,":", V2,"-", V3), strand = "-")

tss_window.bed <- rbind.data.frame(tss_window_plus.bed, tss_window_minus.bed) %>%
    dplyr::select(group.var, strand) %>%
    left_join(offset.df) %>%
    drop_na() %>%
    separate(group.var, into=c("chr","start","end"), sep="[:-]", remove =FALSE) %>%
    mutate(offset = ifelse(strand =="+", offset, -offset)) %>% 
    mutate(end = as.numeric(end) + offset - 1500, start = as.numeric(start) + 1500 + offset)

peaks_for_stats <- file.path(out_dir, "peaks_for_stats.bed")
write.table(tss_window.bed %>%
            dplyr::select(chr, start, end), peaks_for_stats, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


    
