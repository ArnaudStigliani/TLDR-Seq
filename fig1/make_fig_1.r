library(ggplot2)
library(stringr)
library(magrittr)
library(tidyverse)
library(reshape2)
library(ggh4x)

arg  <-  commandArgs(trailingOnly=TRUE) 
lib1.dir <- "../mapping_preprocessing/results/format_adapter_information/lib1/" # first sequencing run = lib1, amplified libraries
lib2.dir <- "../mapping_preprocessing/results/format_adapter_information/lib2/" # second sequencing run = lib2, non amplified libraries

out_dir <-  "./results/" %>% file.path("make_fig_1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### read all files
lib1.trimmed.df <- read.table(file.path(lib1.dir, "trimmed.tsv"))  %>%
    dplyr::rename(barcode = V1, strand = V2, reads_trimmed = V3, percent_reads_trimmed = V4)
lib2.trimmed.df <- read.table(file.path(lib2.dir, "trimmed.tsv"))  %>%
    dplyr::rename(barcode = V1, strand = V2, reads_trimmed = V3, percent_reads_trimmed = V4)

trimmed <- bind_rows(lib1 = lib1.trimmed.df, lib2 = lib2.trimmed.df, .id="lib") %>%
    group_by(lib) %>%
    mutate(tot_reads_trimmed = sum(reads_trimmed)) %>%
    group_by(lib, barcode, tot_reads_trimmed ) %>%
    dplyr::summarize(reads_trimmed = sum(reads_trimmed)) %>%
    ungroup %>%
    dplyr::filter(reads_trimmed/tot_reads_trimmed > 0.01) 


### get total number of reads in that pass basecalling
log.file.path.lib1 <- lib1.dir %>% dirname %>% dirname %>% file.path("pipeline_map", "lib1"  "trimming", "bc1_s1.log")
command.line.lib1 <- paste("head -8 ", log.file.path.lib1, " | tail -1 | awk '{print $4}' | sed 's/,//g'",  sep="")
nread.lib1 <- system(command.line.lib1, intern = TRUE) %>% as.numeric

log.file.path.lib2 <- lib2.dir %>% dirname %>% dirname %>% file.path("pipeline_map", "lib2", "trimming", "bc1_s1.log")
command.line.lib2 <- paste("head -8 ", log.file.path.lib2, " | tail -1 | awk '{print $4}' | sed 's/,//g'",  sep="")
nread.lib2 <- system(command.line.lib2, intern = TRUE) %>% as.numeric

nread.df <- data.frame(lib=c("lib1", "lib2"), nread = c(nread.lib1/3, nread.lib2/12))

trimmed.nreads.df <- left_join(trimmed, nread.df) %>%
    dplyr::select(-tot_reads_trimmed) %>%
    mutate(percent_reads_trimmed = reads_trimmed/nread *100 ) %>%
    mutate(barcode = factor(barcode, levels = paste0("bc", 1:12))) 

#### plot read number (not average)
g <- ggplot(trimmed.nreads.df , aes(x = barcode, y = reads_trimmed, fill= lib)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("percentage of reads with 2 adapters")) +
    facet_grid2(.~ lib, scales="free", space="free_x", independent ="y") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust=0.5)) +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) 
ggpubr::ggexport(g, filename = file.path(out_dir, "nb_reads_2adapters.pdf"), width = 5, height =3 )


#### how many reads map
mapped.lib1.name <- file.path(out_dir, "mapped_lib1.tsv")
mapped.lib2.name <- file.path(out_dir, "mapped_lib2.tsv")
mapped.name <- c(mapped.lib1.name, mapped.lib2.name)


bam.lib1.name  <-  lib1.dir %>% dirname %>% dirname %>% file.path("pipeline_map", "lib1", "mapping", "trimmed_primary.bam")
bam.lib2.name  <-  lib2.dir %>% dirname %>% dirname %>% file.path("pipeline_map", "lib2", "mapping", "trimmed_primary.bam")
bam.name <- c(bam.lib1.name, bam.lib2.name)

command <- paste("module load samtools;samtools view", bam.name ,
                 " | awk -v OFS='\t' -v FS='\t' '{gsub( \"_\",\"\t\",$1)} ($2==0 || $2==16 ) {print $0}' | awk -v OFS='\t' '{print $2,$3}' | sort -k2,2 -u | awk '{print $1}' | sort | uniq -c | awk -v OFS='\t' '{print $1,$2}'", ">", mapped.name)
system(paste(paste(command, collapse=" & \n"), " & wait"))

mapped.df <- lapply(mapped.name, read.table, sep="\t") %>%
    mapply(mutate,. ,lib=c("lib1","lib2"), SIMPLIFY=FALSE) %>% 
    do.call(rbind.data.frame, .) %>%
    setNames(c("n.mapped", "barcode", "lib"))

######
mapped.percent.df <- trimmed.nreads.df %>%
    left_join(mapped.df, by=c("lib","barcode")) %>%
    mutate(percent_reads_mapped = n.mapped / reads_trimmed * 100)


g <- ggplot(mapped.percent.df, aes(x=lib, y=percent_reads_mapped,  colour=lib)) +
    geom_point(position = position_jitter(w = 0.2, h = 0), size=0.9) +
    theme_bw() +
    scale_y_continuous(breaks = c(90, 92, 94, 96, 98, 100), limits=c(90,100)) +
    theme(plot.title = element_text(size = 10, hjust=0.5)) +
    ggtitle(paste0("Mapping rate")) 
out_boxplot <- file.path(out_dir,"mapping_rate.pdf")
ggsave( filename = out_boxplot, g, width = 3, height = 5) 

#### How many reads with a polyA tail

pA.lib1.name <- lib1.dir %>% dirname %>% dirname %>% file.path("characterize_polyA", "lib1", "barcode_nreads_pA.tsv")
pA.lib2.name <- lib2.dir %>% dirname %>% dirname %>% file.path("characterize_polyA", "lib2", "barcode_nreads_pA.tsv")

pA.name <- c(pA.lib1.name, pA.lib2.name)

pA.df <- lapply(pA.name, read.table, sep="\t", header=TRUE) %>%
    mapply(mutate,. ,lib=c("lib1","lib2"), SIMPLIFY=FALSE) %>% 
    do.call(rbind.data.frame, .) 

pA_15.lib1.name <- lib1.dir %>% dirname %>% dirname %>% file.path("characterize_polyA", "lib1", "barcode_nreads_pA_15.tsv")
pA_15.lib2.name <- lib2.dir %>% dirname %>% dirname %>% file.path("characterize_polyA", "lib2", "barcode_nreads_pA_15.tsv")

pA_15.name <- c(pA_15.lib1.name, pA_15.lib2.name)

pA_15.df <- lapply(pA_15.name, read.table, sep="\t", header=TRUE) %>%
    mapply(mutate,. ,lib=c("lib1","lib2"), SIMPLIFY=FALSE) %>% 
    do.call(rbind.data.frame, .) 

pA_all.df <- bind_rows("reads with pA > 15 nt "=pA_15.df, "reads with pA "= pA.df, .id="wrap_pA")

pA.percent.df <- trimmed.nreads.df %>%
    left_join(pA_all.df, by=c("lib","barcode"), multiple= "all") %>%
    mutate(percent_reads_pA = nreads.polyA / reads_trimmed * 100)


g <- ggplot(pA.percent.df,  aes(x=lib, y=percent_reads_pA,  colour=lib)) +
    geom_point(position = position_jitter(w = 0.2, h = 0), size=0.9) +
    theme_bw() +
    scale_y_continuous(breaks = seq(0,20,2), limits=c(0,20)) +
    theme(plot.title = element_text(size = 10, hjust=0.5)) +
    ggtitle(paste0("%pA reads")) +
    facet_wrap(~wrap_pA)
out_boxplot <- file.path(out_dir,"pA_reads.pdf")
ggsave( filename = out_boxplot, g, width = 4, height = 5) 

