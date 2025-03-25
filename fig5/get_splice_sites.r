library(tidyverse)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(stringr)
library(data.table)
library(parallel)
library(plyranges)

out_dir <- "./results/get_splice_sites/"
dir.create(out_dir, showWarnings = FALSE, recursive=TRUE)

bam_tldr.name <- "./results/fig1/results/get_polyA_bam_files/lib1/trimmed_primary_renamed.bam"
bam_tldr.bed12.name <- file.path(out_dir, "tldr.bed12")

bam_tldr_lib2.name <- "./results/fig1/results/get_polyA_bam_files/lib2/trimmed_primary_renamed.bam"
bam_tldr_lib2.bed12.name <- file.path(out_dir, "tldr_lib2.bed12")

bam_nanopore1.name <- "../map_ONTnanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed"
bam_nanopore1.bed12.name <- file.path(out_dir, "nanopore1.bed12")
bam_nanopore2.name <- "../map_ONTnanopore/get_bedfiles/nanopore_wt2_curated_stranded.bed"
bam_nanopore2.bed12.name <- file.path(out_dir, "nanopore2.bed12")

command.line <- paste("module load bedtools;  bedtools bamtobed  -i", bam_tldr.name ," -bed12 > ", bam_tldr.bed12.name)
system(command.line) 
tldr.bed12.df <- read.table(bam_tldr.bed12.name, sep="\t") %>%
    dplyr::select(chr = V1, start.region = V2, strand=V6,  block.start = V12, block.size = V11, read.name =V4 )

command.line <- paste("module load bedtools;  bedtools bamtobed  -i", bam_tldr_lib2.name ," -bed12 > ", bam_tldr_lib2.bed12.name)
system(command.line) 
tldr_lib2.bed12.df <- read.table(bam_tldr_lib2.bed12.name, sep="\t") %>%
    dplyr::select(chr = V1, start.region = V2, strand=V6,  block.start = V12, block.size = V11, read.name =V4 )

## nanopore
command.line.1 <- paste("module load bedtools;  bedtools bamtobed  -i", bam_nanopore1.name ," -bed12 > ", bam_nanopore1.bed12.name)
command.line.2 <- paste("module load bedtools;  bedtools bamtobed  -i", bam_nanopore2.name ," -bed12 > ", bam_nanopore2.bed12.name)
all.command.lines <- paste(paste(command.line.1, command.line.2, sep=" & \n"), " & wait")
system(all.command.lines)

nanopore.bed12.df <- list(bam_nanopore1.bed12.name, bam_nanopore2.bed12.name) %>%
    lapply(read.table, sep="\t") %>% 
    do.call(rbind.data.frame, .) %>% 
    dplyr::select(chr = V1, start.region = V2, strand=V6,  block.start = V12, block.size = V11, read.name =V4 )

# gencode
gencode.bed12.name <- "../data/gencode.v39.primary_assembly.subset_selected_transcripts.bed12"
gencode.bed12.df <- read.table(gencode.bed12.name, sep="\t") %>%
    dplyr::select(chr = V1, start.region = V2, strand=V6, block.start = V12, block.size = V11) %>%
    mutate(block.size = block.size %>% str_replace(",$", ""), block.start = block.start %>% str_replace(",$", "")) 
    



gencode.splice.left <-  gencode.bed12.df$block.start %>%
    str_split(",")  %>%
    mclapply(as.numeric, mc.cores = 50) %>%
    mcmapply(function(x, y) x + y, . , gencode.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z) paste(x, y, y + 1, ".", ".", z, sep=":"),
             gencode.bed12.df$chr , .,
             gencode.bed12.df$strand,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-1], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") %>%
    mutate(start = as.numeric(start), end = as.numeric(end) -1 ) %>%
    mutate(name = ifelse(strand == "+", "exon.start", "exon.end"))
gencode.splice.right <-  mcmapply(function(x, y) x + y,
                                  gencode.bed12.df$block.start %>% str_split(",") %>% mclapply(as.numeric,mc.cores = 50),
                                  gencode.bed12.df$block.size %>% str_split(",") %>% mclapply(as.numeric, mc.cores = 50),
                                  SIMPLIFY=FALSE, mc.cores = 50)  %>%
    mcmapply(function(x, y) x + y, . , gencode.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z) paste(x, y, y + 1, ".", ".", z, sep=":"),
             gencode.bed12.df$chr , .,
             gencode.bed12.df$strand,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-length(x)], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") %>%
    mutate(start = as.numeric(start), end = as.numeric(end) ) %>%
    mutate(name = ifelse(strand == "+", "exon.end", "exon.start"))


gencode.exon.start  <-  rbind.data.frame(gencode.splice.left, gencode.splice.right) %>%
    dplyr::filter(name =="exon.start", grepl("chr", chr)) %>%
    distinct %>%
    arrange(chr, start, end)

gencode.exon.end  <-  rbind.data.frame(gencode.splice.left, gencode.splice.right) %>%
    dplyr::filter(name =="exon.end", grepl("chr", chr)) %>%
    distinct %>%
    arrange(chr, start, end) 


gencode.intron.bed <- cbind.data.frame(gencode.splice.left, gencode.splice.right) %>%
    setNames(paste0(names(.), c(".l", ".l", ".l", ".l", ".l", ".l", ".r", ".r", ".r", ".r", ".r", ".r"))) %>%
    dplyr::select(chr = chr.l, start = start.r, end = end.l, name = name.l, score = score.l, strand = strand.l ) %>%
    mutate(score=100, name = ".") %>%
    mutate(start = start %>% as.numeric(), end  = end %>% as.numeric) %>%
    distinct()
    


options(scipen=15)
write.table(gencode.exon.start, file.path(out_dir, "gencode_exon_start.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gencode.exon.end, file.path(out_dir, "gencode_exon_end.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gencode.intron.bed, file.path(out_dir, "gencode_intron.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

##### TLDR

tldr.splice.left <-  tldr.bed12.df$block.start %>%
    str_split(",")  %>%
    mclapply(as.numeric, mc.cores = 50) %>%
    mcmapply(function(x, y) x + y, . , tldr.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr.bed12.df$chr , .,
             tldr.bed12.df$strand,
             tldr.bed12.df$read.name,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-1], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") %>%
    mutate(start = as.numeric(start), end = as.numeric(end) -1 ) 
tldr.splice.right <-  mcmapply(function(x, y) x + y,
                                  tldr.bed12.df$block.start %>% str_split(",") %>% mclapply(as.numeric,mc.cores = 50),
                                  tldr.bed12.df$block.size %>% str_split(",") %>% mclapply(as.numeric, mc.cores = 50),
                                  SIMPLIFY=FALSE, mc.cores = 50)  %>%
    mcmapply(function(x, y) x + y, . , tldr.bed12.df$start.region, SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr.bed12.df$chr , .,
             tldr.bed12.df$strand,
             tldr.bed12.df$read.name,             
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-length(x)], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") 

tldr.intron.bed <- cbind.data.frame(tldr.splice.left, tldr.splice.right) %>%
    setNames(paste0(names(.), c(".l", ".l", ".l", ".l", ".l", ".l", ".r", ".r", ".r", ".r", ".r", ".r"))) %>%
    dplyr::select(chr = chr.l, start = start.r, end = end.l, name = name.l, score = score.l, strand = strand.l ) %>%
    mutate(score=1) %>%
    mutate(start = start %>% as.numeric(), end  = end %>% as.numeric)

options(scipen=15)
write.table(tldr.intron.bed, file.path(out_dir, "tldr_intron.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#### TLDR lib2

tldr_lib2.splice.left <-  tldr_lib2.bed12.df$block.start %>%
    str_split(",")  %>%
    mclapply(as.numeric, mc.cores = 50) %>%
    mcmapply(function(x, y) x + y, . , tldr_lib2.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr_lib2.bed12.df$chr , .,
             tldr_lib2.bed12.df$strand,
             tldr_lib2.bed12.df$read.name,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-1], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") %>% 
    mutate(start = as.numeric(start), end = as.numeric(end) -1 ) 
tldr_lib2.splice.right <-  mcmapply(function(x, y) x + y,
                                  tldr_lib2.bed12.df$block.start %>% str_split(",") %>% mclapply(as.numeric,mc.cores = 50),
                                  tldr_lib2.bed12.df$block.size %>% str_split(",") %>% mclapply(as.numeric, mc.cores = 50),
                                  SIMPLIFY=FALSE, mc.cores = 50)  %>%
    mcmapply(function(x, y) x + y, . , tldr_lib2.bed12.df$start.region, SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             tldr_lib2.bed12.df$chr , .,
             tldr_lib2.bed12.df$strand,
             tldr_lib2.bed12.df$read.name,             
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-length(x)], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") 

tldr_lib2.intron.bed <- cbind.data.frame(tldr_lib2.splice.left, tldr_lib2.splice.right) %>%
    setNames(paste0(names(.), c(".l", ".l", ".l", ".l", ".l", ".l", ".r", ".r", ".r", ".r", ".r", ".r"))) %>%
    dplyr::select(chr = chr.l, start = start.r, end = end.l, name = name.l, score = score.l, strand = strand.l ) %>%
    mutate(score=1) %>%
    mutate(start = start %>% as.numeric(), end  = end %>% as.numeric)

options(scipen=15)
write.table(tldr_lib2.intron.bed, file.path(out_dir, "tldr_lib2_intron.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


## ##### NANOPORE
nanopore.splice.left <-  nanopore.bed12.df$block.start %>%
    str_split(",")  %>%
    mclapply(as.numeric, mc.cores = 50) %>%
    mcmapply(function(x, y) x + y, . , nanopore.bed12.df$start.region ,  SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             nanopore.bed12.df$chr , .,
             nanopore.bed12.df$strand,
             nanopore.bed12.df$read.name,
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-1], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") 
nanopore.splice.right <-  mcmapply(function(x, y) x + y,
                                  nanopore.bed12.df$block.start %>% str_split(",") %>% mclapply(as.numeric,mc.cores = 50),
                                  nanopore.bed12.df$block.size %>% str_split(",") %>% mclapply(as.numeric, mc.cores = 50),
                                  SIMPLIFY=FALSE, mc.cores = 50)  %>%
    mcmapply(function(x, y) x + y, . , nanopore.bed12.df$start.region, SIMPLIFY=FALSE, mc.cores = 50) %>%
    mcmapply(function(x, y, z, z2) paste(x, y, y + 1, z2, ".", z, sep=":"),
             nanopore.bed12.df$chr , .,
             nanopore.bed12.df$strand,
             nanopore.bed12.df$read.name,             
             SIMPLIFY=FALSE, mc.cores = 50) %>%
    mclapply(function(x) x[-length(x)], mc.cores=50) %>%
    mclapply(function(x) seq=data.frame(x), mc.cores=50) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(x, into=c("chr","start","end","name","score","strand"), sep=":") 

nanopore.intron.bed <- cbind.data.frame(nanopore.splice.left, nanopore.splice.right) %>%
    setNames(paste0(names(.), c(".l", ".l", ".l", ".l", ".l", ".l", ".r", ".r", ".r", ".r", ".r", ".r"))) %>%
    dplyr::select(chr = chr.l, start = start.r, end = end.l, name = name.l, score = score.l, strand = strand.l ) %>%
    mutate(score=1) %>%
    mutate(start = start %>% as.numeric(), end  = end %>% as.numeric)


write.table(nanopore.intron.bed, file.path(out_dir, "nanopore_intron.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)






