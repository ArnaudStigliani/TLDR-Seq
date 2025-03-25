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


in_dir <- "./results/shared/eiffel_tower_plot_enhancer/"

arg = commandArgs(trailingOnly=TRUE)

geno <- arg[1]
## geno <- "mES-RBM7"

out_dir <- file.path(in_dir, geno)
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )


reads.curated.name <- file.path("./results/", geno, 
                                "/get_primary_processed_bam/trimmed_primary_processed.bam")

reads.curated_pAp.name <- file.path("./results/", geno,
                                "/get_primary_processed_bam/trimmed_primary_polyA_plus.bam")

reads.curated_pAm.name <- file.path("./results/", geno,
                                "/get_primary_processed_bam/trimmed_primary_polyA_minus.bam")



gencode_annot.name <- "../../../shared_data/Mus_musculus.GRCm39.103.bed6"

#### make a bed file, not sure it is useful but maybe it will be so that we can calculate sequencing depth

read_curated_bed.name <- file.path("./results/shared/eiffel_tower_plot_enhancer",geno,"trimmed_primary.bed")
command1 <- paste("module load bedtools; bedtools bamtobed -i", reads.curated.name , ">", read_curated_bed.name)

system(command1)



############ generate non annotated  antisense reads ################
reads.curated_pAp_noAnnot.name <- file.path(out_dir, "all_reads_pAp_noAnnot.bam")
reads.curated_pAm_noAnnot.name <- file.path(out_dir, "all_reads_pAm_noAnnot.bam")


DHS.name <- "../../data/mES_histones_dnase/ENCFF048DWN_mm39.bed"


command1 <- paste("module load bedtools;", 
                  "bedtools intersect -v -abam ", reads.curated_pAp.name, 
                  "-b",gencode_annot.name, ">",reads.curated_pAp_noAnnot.name)

command2 <- paste("module load bedtools;", 
                  "bedtools intersect -v -abam ", reads.curated_pAm.name, 
                  "-b",gencode_annot.name, ">",reads.curated_pAm_noAnnot.name)

system(paste(command1, "&",command2," & wait"))


###### get primary alignment from mapped reads on transcriptome

tss_window_plus.bed.name <- file.path(in_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(in_dir, "tss_minus_window.bed")


tss_window_plus.bed <- read.table(tss_window_plus.bed.name, sep="\t")
tss_window_minus.bed <- read.table(tss_window_minus.bed.name, sep="\t")


################ pA plus ###################
############################################
tss_window_plus_forward_pAp.bdg.name <-   file.path(out_dir, "tss_plus_window_forward_pAp.bdg")
tss_window_plus_reverse_pAp.bdg.name <-   file.path(out_dir, "tss_plus_window_reverse_pAp.bdg")
tss_window_minus_forward_pAp.bdg.name <-   file.path(out_dir, "tss_minus_window_forward_pAp.bdg")
tss_window_minus_reverse_pAp.bdg.name <-   file.path(out_dir, "tss_minus_window_reverse_pAp.bdg")


command.line.1 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_pAp_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_forward_pAp.bdg.name)

command.line.1bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_pAp_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_reverse_pAp.bdg.name)

command.line.2 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_pAp_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_forward_pAp.bdg.name)

command.line.2bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_pAp_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_reverse_pAp.bdg.name)

################## pA minus ###########################
#######################################################

tss_window_plus_forward_pAm.bdg.name <-   file.path(out_dir, "tss_plus_window_forward_pAm.bdg")
tss_window_plus_reverse_pAm.bdg.name <-   file.path(out_dir, "tss_plus_window_reverse_pAm.bdg")
tss_window_minus_forward_pAm.bdg.name <-   file.path(out_dir, "tss_minus_window_forward_pAm.bdg")
tss_window_minus_reverse_pAm.bdg.name <-   file.path(out_dir, "tss_minus_window_reverse_pAm.bdg")


command.line.3 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_pAm_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_forward_pAm.bdg.name)

command.line.3bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_pAm_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_reverse_pAm.bdg.name)

command.line.4 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_pAm_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_forward_pAm.bdg.name)

command.line.4bis <- paste("module load samtools; module load bedtools;",
                           " samtools view -Sb -F 20", reads.curated_pAm_noAnnot.name, 
                           "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-2500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_reverse_pAm.bdg.name)

system(paste(command.line.1, "&",command.line.1bis, "&", command.line.2,"&", command.line.2bis, "&",
             command.line.3, "&",command.line.3bis, "&", command.line.4,"&", command.line.4bis,
             " & wait"))

## ########################## dnase 1
## ##################################

if (! is.na(arg[2]))
{
    dnase1.bdg.name <-   file.path("../../data/mES_histones_dnase/Dnase1_all-mm_norm_mm39s_split.bdg")
    dnase.tss_plus.bdg.name <- file.path(in_dir, "dnase1_tss_plus.bdg")
    dnase.tss_minus.bdg.name <- file.path(in_dir, "dnase1_tss_minus.bdg")

    command.line.1 <- paste("module load bedtools;",
                            "bedtools intersect -wa -wb  -a",tss_window_plus.bed.name,
                            "-b", dnase1.bdg.name, ">",    
                            dnase.tss_plus.bdg.name)

    command.line.2 <- paste("module load bedtools;",
                            "bedtools intersect -wa -wb -a",tss_window_minus.bed.name,
                            "-b", dnase1.bdg.name, ">",    
                            dnase.tss_minus.bdg.name)

    system(paste(command.line.1, "&",command.line.2," & wait"))
}
