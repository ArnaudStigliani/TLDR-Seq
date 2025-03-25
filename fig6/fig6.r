library(ShortRead)
library(bambu)
library(tidyverse)
library(ggplot2)


out_dir <- "./results/make_fig6"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)


in.dir.RNASeq <- "./data/RNASeq/salmon_mapping"
gtf.name <- "../data/gencode.v39.primary_assembly.annotation_Salmon.gtf"
fas.name <- "../data/GRCh38.primary_assembly.genome.fa"

out_RNASeq <- "../fig1/annotate_reads/RNASeq"



### read RNASeq
RNASeq.name <- file.path(out_RNASeq, "quant.sf") 

RNASeq.df <-  read.table(RNASeq.name, header=TRUE)%>%
    dplyr::select(TXNAME=Name,  value = TPM) %>%
    mutate(variable = "RNASeq")


## Quantify all bambu based bam files

## lib1.bam and lib2.bam are symbolic links to the bam files in ../fig1/get_polyA_bam_files/trimmed_primary_renamed.bam

all_bam_files <- c(file.path(out_dir, c("lib1.bam","lib2.bam")),
                   file.path("../../data/nanopore/results_curated/nanopore_all_curated.bam"))
                   



annotations <- prepareAnnotations(gtf.name)
saveRDS(annotations, file.path(out_dir, "annotations.rds")) 
annotations <- readRDS(file.path(out_dir, "annotations.rds"))
LongReads.quant <- bambu(reads = all_bam_files,  annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 3)


TXinfos <- LongReads.quant %>% rowData %>% as.data.frame %>%
    dplyr::select(TXNAME, GENEID) %>%
    remove_rownames
write.table(TXinfos, file=file.path(out_dir, "txinfos.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
TXinfos <- read.table(file.path(out_dir, "txinfos.tsv"), sep="\t", header=TRUE) 
duplicates.df.name <-  "Path_tosalmon_index/duplicate_clusters.tsv"
duplicates.df <- duplicates.df.name %>% read.table(sep="\t", header=TRUE) %>%
    dplyr::rename(TXNAME=DuplicateRef)



LongReads.cpm.manual <-  assays(LongReads.quant)$counts %>%
                                               as.data.frame %>%
                                               rownames_to_column("TXNAME") %>%
                                               reshape2::melt() %>%
                                               left_join(duplicates.df) %>%
                                               mutate(TXNAME = ifelse(is.na(RetainedRef), TXNAME, RetainedRef)) %>%
                                               dplyr::select(-RetainedRef) %>%
                                               dplyr::filter(TXNAME %in% RNASeq.df$TXNAME) %>%
                                               group_by(TXNAME, variable) %>%
                                               dplyr::summarize(value = sum(value)) %>% 
                                               as.data.frame  %>%
                                               group_by(variable) %>%
                                               mutate(tot = sum(value)) %>%                                 
                                               mutate(value = 10**6 * value/tot) %>%
                                               dplyr::select(-tot)
write.table(LongReads.cpm.manual, file=file.path(out_dir, "tx.counts.gene.tpm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
LongReads.cpm.manual <- read.table(file.path(out_dir, "tx.counts.gene.tpm.tsv"), sep="\t", header=TRUE) 



all.tpm <- rbind.data.frame(LongReads.cpm.manual, RNASeq.df) %>%
    left_join(TXinfos) %>%
    group_by(GENEID, variable) %>%
    dplyr::summarize(value = sum(value))  %>%
    mutate(log.tpm = log10(value +1)) %>%
    reshape2::dcast(GENEID ~ variable, value.var="log.tpm") %>%
    dplyr::rename(TLDR.lib.Amp = lib1, TLDR.lib.NoAmp = lib2, ONT.Nanopore = nanopore_all_curated) 


### make figs

cor.stat <- cor.test(all.tpm$TLDR.lib.Amp, all.tpm$RNASeq, method="pearson") %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.tpm, aes(y = RNASeq, x = TLDR.lib.Amp)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_bw() +
    ggtitle(cor.stat) +
    xlab("TLDR-Seq  Amplified lib expression (log10 CPM)") +
    ylab("RNASeq expression (log10 TPM)")
out_dotplot <- file.path(out_dir,"comparison_RNASeq_libAmp_geneExp.png")
ggsave( filename = out_dotplot, g, width = 6, height = 6)


cor.stat <- cor.test(all.tpm$TLDR.lib.NoAmp, all.tpm$RNASeq, method="pearson") %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.tpm, aes(y = RNASeq, x = TLDR.lib.NoAmp)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_bw() +
    ggtitle(cor.stat) +
    xlab("TLDR-Seq  Non-amplified lib expression (log10 CPM)") +
    ylab("RNASeq expression (log10 TPM)")
out_dotplot <- file.path(out_dir,"comparison_RNASeq_lib_No_Amp_geneExp.png")
ggsave( filename = out_dotplot, g, width = 6, height = 6)


## ONT nanopore vs RNASeq  for reference
cor.stat <- cor.test(all.tpm$ONT.Nanopore, all.tpm$RNASeq, method="pearson") %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.tpm, aes(y = RNASeq, x = TLDR.lib.NoAmp)) +
    geom_point(alpha = 0.3, size = 0.5) +
    theme_bw() +
    ggtitle(cor.stat) +
    xlab("ONT-Nanopore (log10 CPM)") +
    ylab("RNASeq expression (log10 TPM)")
out_dotplot <- file.path(out_dir,"comparison_RNASeq_ONT_geneExp_FigS3.png")
ggsave( filename = out_dotplot, g, width = 6, height = 6)



