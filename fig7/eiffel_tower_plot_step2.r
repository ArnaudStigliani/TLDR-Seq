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


in_dir <- "../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot/"

arg = commandArgs(trailingOnly=TRUE)

geno <- arg[1]
## geno <- "mES-ZCCHC8"

out_dir <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot/", geno)
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )


reads.curated.name <- file.path("../../results/sequencing_run_before_promethion2/", geno,
                                "/get_primary_processed_bam/trimmed_primary_polyA_minus.bam")
reads.curated_polyA.name <- file.path("../../results/sequencing_run_before_promethion2/",geno,
                                      "/get_primary_processed_bam/trimmed_primary_polyA_plus.bam")


gencode_annot.name <- "../../../shared_data/Mus_musculus.GRCm39.103.bed6" ####  Update tha name


#### make a bed file, not sure it is useful but maybe it will be so that we can calculate sequencing depth

read_curated_bed.name <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot",geno,"trimmed_primary.bed")
command1 <- paste("module load bedtools; bedtools bamtobed -i", reads.curated.name , ">", read_curated_bed.name)
system(command1)

read_curated_polyA.bed.name <- file.path("../../results/sequencing_run_before_promethion2/shared/eiffel_tower_plot",geno,"trimmed_primary_polyA_plus.bed")
command1 <- paste("module load bedtools; bedtools bamtobed -i", reads.curated_polyA.name , ">", read_curated_polyA.bed.name)
system(command1)


############ generate non annotated  antisense reads ################
reads.curated_noAnnot.name <- file.path(out_dir, "all_reads_noAnnot.bam")
reads.curated_polyA_noAnnot.name <- file.path(out_dir, "pA_reads_noAnnot.bam")

command1 <- paste("module load bedtools;",
                  "bedtools intersect -v -abam",reads.curated.name,
                  "-b",gencode_annot.name, ">",reads.curated_noAnnot.name )
command2 <- paste("module load bedtools;",
                  "bedtools intersect -v  -abam",reads.curated_polyA.name,
                  "-b",gencode_annot.name, ">",reads.curated_polyA_noAnnot.name )

system(paste(command1, "&", command2, "& wait"))


##### get  library n reads

nreads <- system(paste("module load samtools; samtools view", reads.curated.name, " | wc -l"), intern=TRUE)
norm.factor <- 10^6/as.numeric(nreads)

###### get primary alignment from mapped reads on transcriptome

tss_window_plus.bed.name <- file.path(in_dir, "tss_plus_window.bed")
tss_window_minus.bed.name  <- file.path(in_dir, "tss_minus_window.bed")


tss_window_plus.bed <- read.table(tss_window_plus.bed.name, sep="\t")
tss_window_minus.bed <- read.table(tss_window_minus.bed.name, sep="\t")

###### All reads
#### get depth from bam file
tss_window_plus_forward.bdg.name <-   file.path(out_dir, "tss_plus_window_forward.bdg")
tss_window_plus_reverse.bdg.name <-   file.path(out_dir, "tss_plus_window_reverse.bdg")
tss_window_minus_forward.bdg.name <-   file.path(out_dir, "tss_minus_window_forward.bdg")
tss_window_minus_reverse.bdg.name <-   file.path(out_dir, "tss_minus_window_reverse.bdg")


command.line.1 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_forward.bdg.name)

command.line.1bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_reverse.bdg.name)

command.line.2 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_forward.bdg.name)

command.line.2bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_reverse.bdg.name)


system(paste(command.line.1, "&",command.line.1bis, "&", command.line.2,"&", command.line.2bis," & wait"))


######## Same thing but with pA reads
###### pA reads
#### get depth from bam file
tss_window_plus_forward_polyA.bdg.name <-   file.path(out_dir, "tss_plus_window_forward_polyA.bdg")
tss_window_plus_reverse_polyA.bdg.name <-   file.path(out_dir, "tss_plus_window_reverse_polyA.bdg")
tss_window_minus_forward_polyA.bdg.name <-   file.path(out_dir, "tss_minus_window_forward_polyA.bdg")
tss_window_minus_reverse_polyA.bdg.name <-   file.path(out_dir, "tss_minus_window_reverse_polyA.bdg")


command.line.1 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_polyA.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_forward_polyA.bdg.name)

command.line.1bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_polyA_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_plus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_plus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_plus_reverse_polyA.bdg.name)


command.line.2 <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -f 16", reads.curated_polyA.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_forward_polyA.bdg.name)

command.line.2bis <- paste("module load samtools; module load bedtools;",
                       " samtools view -Sb -F 20", reads.curated_polyA_noAnnot.name, 
                       "| samtools depth -aa -b ",tss_window_minus.bed.name,"-",
                       "| awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}'",
                       "| bedtools intersect -wa -wb -a",tss_window_minus.bed.name , "-b -",
                       "| awk -v OFS='\t' '{print $1, $2, $3, $5-$2-1500, $7}'",
                       "| sed 's/\t/:/' | sed 's/\t/-/' > ",
                        tss_window_minus_reverse_polyA.bdg.name)
system(paste(command.line.1, "&",command.line.1bis, "&", command.line.2,"&", command.line.2bis," & wait"))



## Dnase1

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

