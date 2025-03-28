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

in_dir <- "./results/shared/eiffel_tower_plot/"


genos <- c("mES-WT1", "mES-RBM7", "mES-ZCCHC8")

## geno <- c("mES-WT1")

for (geno in genos)
{
    out_dir <- file.path("./results/shared/summary_stats_pAminus", geno)
    eiffel_dir  <-  file.path(in_dir, geno)
    ##
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )
    ##
   
### bam files
    reads.curated.name <- file.path("./results/", geno,
                                    "/get_primary_processed_bam/trimmed_primary_polyA_minus.bam")

    reads.curated_noAnnot.name <- file.path(out_dir, "all_reads_noAnnot.bam")
    reads.curated_Annot.name <- file.path(out_dir, "all_reads_Annot.bam")

    gencode_annot.name <- "../../../shared_data/Mus_musculus.GRCm39.103.bed6" ####  Update tha name

    command1 <- paste("module load bedtools;",
                  "bedtools intersect -v -abam",reads.curated.name,
                  "-b",gencode_annot.name, ">",reads.curated_noAnnot.name )
    command2 <- paste("module load bedtools;",
                  "bedtools intersect  -abam",reads.curated.name,
                  "-b",gencode_annot.name, ">",reads.curated_Annot.name )

    system(paste(command1,"&", command2))
    ## TSS
    tss_window_plus.bed.name <- file.path(in_dir, "peaks_for_stats_minus.bed")
    tss_window_minus.bed.name  <- file.path(in_dir, "peaks_for_stats_minus.bed")


###### get reads
#### get depth from bam file
    tss_window_plus_forward.bed12.name <-   file.path(out_dir, "tss_plus_window_forward.bed12")
    tss_window_plus_reverse.bed12.name <-   file.path(out_dir, "tss_plus_window_reverse.bed12")
    tss_window_minus_forward.bed12.name <-   file.path(out_dir, "tss_minus_window_forward.bed12")
    tss_window_minus_reverse.bed12.name <-   file.path(out_dir, "tss_minus_window_reverse.bed12")


    command.line.1 <- paste("module load samtools; module load bedtools;",
                            " samtools view -Sb -F 20", reads.curated_Annot.name,
                            "| bedtools intersect -abam - -b", tss_window_plus.bed.name, 
                            "| bedtools bamtobed -bed12 -i - >", 
                            tss_window_plus_forward.bed12.name)

    command.line.1bis <- paste("module load samtools; module load bedtools;",
                               " samtools view -Sb -f 16", reads.curated_noAnnot.name, 
                               "| bedtools intersect -abam - -b",tss_window_plus.bed.name,  
                               "| bedtools bamtobed -bed12 -i - >", 
                               tss_window_plus_reverse.bed12.name)

    command.line.2 <- paste("module load samtools; module load bedtools;",
                            " samtools view -Sb -f 16", reads.curated_Annot.name, 
                            "| bedtools intersect -abam - -b",tss_window_minus.bed.name,  
                            "| bedtools bamtobed -bed12 -i - >", 
                            tss_window_minus_forward.bed12.name)

    command.line.2bis <- paste("module load samtools; module load bedtools;",
                               " samtools view -Sb -F 20", reads.curated_noAnnot.name, 
                               "| bedtools intersect -abam - -b",tss_window_minus.bed.name, 
                               "| bedtools bamtobed -bed12 -i - >", 
                               tss_window_minus_reverse.bed12.name)


    system(paste(command.line.1, "&",command.line.1bis, "&", command.line.2,"&", command.line.2bis," & wait"))

    tss_window_plus_forward_wawb.bed12.name <-   file.path(out_dir, "tss_plus_window_forward_wawb.bed12")
    tss_window_plus_reverse_wawb.bed12.name <-   file.path(out_dir, "tss_plus_window_reverse_wawb.bed12")
    tss_window_minus_forward_wawb.bed12.name <-   file.path(out_dir, "tss_minus_window_forward_wawb.bed12")
    tss_window_minus_reverse_wawb.bed12.name <-   file.path(out_dir, "tss_minus_window_reverse_wawb.bed12")

    command.line.1 <- paste(" module load bedtools;",
                            "bedtools intersect -wa -wb -a",
                            tss_window_plus_forward.bed12.name," -b", tss_window_plus.bed.name, ">", 
                            tss_window_plus_forward_wawb.bed12.name)

    command.line.1bis <- paste(" module load bedtools;",
                            "bedtools intersect -wa -wb -a",
                            tss_window_plus_reverse.bed12.name," -b", tss_window_plus.bed.name, ">", 
                            tss_window_plus_reverse_wawb.bed12.name)

    command.line.2 <- paste(" module load bedtools;",
                            "bedtools intersect -wa -wb -a",
                            tss_window_minus_forward.bed12.name," -b", tss_window_minus.bed.name, ">", 
                            tss_window_minus_forward_wawb.bed12.name)

    command.line.2bis <- paste(" module load bedtools;",
                            "bedtools intersect -wa -wb -a",
                            tss_window_minus_reverse.bed12.name," -b", tss_window_minus.bed.name, ">", 
                            tss_window_minus_reverse_wawb.bed12.name)

    system(paste(command.line.1, "&",command.line.1bis, "&", command.line.2,"&", command.line.2bis," & wait"))

##### analyze data
}

for (geno in genos)
{
    out_dir <- file.path("./results/shared/summary_stats_pAminus", geno)
    tss_window_plus_forward_wawb.bed12.name <-   file.path(out_dir, "tss_plus_window_forward_wawb.bed12")
    tss_window_plus_reverse_wawb.bed12.name <-   file.path(out_dir, "tss_plus_window_reverse_wawb.bed12")
    tss_window_minus_forward_wawb.bed12.name <-   file.path(out_dir, "tss_minus_window_forward_wawb.bed12")
    tss_window_minus_reverse_wawb.bed12.name <-   file.path(out_dir, "tss_minus_window_reverse_wawb.bed12")

    prompts.df <- rbind.data.frame(read.table(tss_window_minus_reverse_wawb.bed12.name, sep="\t"),
                                   read.table(tss_window_plus_reverse_wawb.bed12.name, sep="\t")) %>%
        bind_rows(prompts = .,
                  mRNAs = rbind.data.frame(read.table(tss_window_minus_forward_wawb.bed12.name, sep="\t"),
                                           read.table(tss_window_plus_forward_wawb.bed12.name, sep="\t")),
                  .id = "transcript_type")

    prompts.summarized.df  <-  prompts.df %>%
        dplyr::select(transcript_type, nblock = V10) %>%
        mutate(read_type =  ifelse(nblock > 1, "spliced", "unspliced")) %>%
        group_by(transcript_type, read_type) %>%
        dplyr::summarize(n = n()) %>%
        group_by(transcript_type) %>%
        mutate(nreads = sum(n)) %>%
        ungroup %>%
        mutate(fraction = n/nreads, geno = geno) 
    write.table(prompts.summarized.df, file.path(out_dir, "summary_table.tsv"), sep="\t", quote =FALSE, row.names=FALSE )

    prompts.summarized_detail.df  <-  prompts.df %>%
        dplyr::select(transcript_type, nintrons = V10) %>%
        mutate(nintrons = nintrons - 1 ) %>% 
        group_by(transcript_type, nintrons) %>%
        dplyr::summarize(n = n()) %>%
        group_by(transcript_type) %>%
        mutate(nreads = sum(n)) %>%
        ungroup %>%
        mutate(fraction = n/nreads, geno = geno) 
    write.table(prompts.summarized_detail.df, file.path(out_dir, "summary_table_detail.tsv"), sep="\t", quote =FALSE, row.names=FALSE )
    
    prompts.summarized_length.df  <-  prompts.df %>%
        mutate(TSS=paste0(V13,":",V14,"-",V15), genomic.length = V3-V2) %>% 
        dplyr::select(transcript_type, blocks = V11, TSS, genomic.length) %>%
        mutate(length = blocks %>% str_split(",") %>% lapply(as.numeric) %>% sapply(sum, na.rm =TRUE)) %>% 
        group_by(transcript_type, TSS) %>%
        dplyr::summarize(average_length = mean(length), average_genomic_length = mean(genomic.length)) %>%
        ungroup %>%
        mutate(geno = geno)
    
    write.table(prompts.summarized_length.df, file.path(out_dir, "summary_table_length.tsv"), sep="\t", quote =FALSE, row.names=FALSE )
    
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
    theme_classic() +
    facet_grid(vars(transcript_type))
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
    theme_classic() +
    facet_grid(vars(transcript_type), scales = "free_x") ## otherwise, use faacet_grid2 with independand ="x"
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

g <- ggplot(all.summary_length.df , aes(colour = transcript_type, x = value)) +
    geom_density() +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2")[2:3])+
    theme_classic() +
    facet_grid(vars(geno), vars(variable), scales = "free_x") 
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types_length.pdf")
ggsave(out_eiffel, g, width=6, height=4)

g <- ggplot(all.summary_length.df , aes(colour = transcript_type, x = value, y=..ndensity..)) +
    geom_density() +
    ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values=RColorBrewer::brewer.pal(n=3, "Set2")[2:3])+
    theme_classic() +
    facet_grid(vars(geno), vars(variable), scales = "free_x") 
out_eiffel <- file.path(dirname(out_dir),  "summary_read_types_length_ndensity.pdf")
ggsave(out_eiffel, g, width=6, height=4)
