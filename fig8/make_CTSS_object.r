## CRAN packages for data manipulation and plotting
rm(list=ls())
library(tidyverse)

# CAGEfightR and related packages
library(CAGEfightR)

# Bioconductor packages for differential expression
library(DESeq2)
library(limma)
library(edgeR)
library(sva)
library(BRGenomics)
library(ggpubr)
library(plyranges)
library(stringr)
library(data.table)
require(genomation)
require(GGally)
require(parallel)
## source("CAGEfightR_extensions/decompose.R")
## source("CAGEfightR_extensions/enhancers.R")
## source("CAGEfightR_extensions/normalize.R")
## source("CAGEfightR_extensions/heatmap.R")

# Bioconductor data packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
require(AnnotationDbi)
require(GenomicFeatures) 
bsg <- BSgenome.Hsapiens.UCSC.hg38
odb <- org.Hs.eg.db

chr.lengths <- seqlengths(bsg) %>% .[grepl("chr[0-9][0-9]*$", names(.)) | grepl("chr[XY]$", names(.))]


### code perso
out_dir <- "./results/make_CTSS_object/"
dir.create(out_dir, showWarnings=FALSE,recursive=TRUE)
reproducibility_dir <- "../fig2/results/reproducibility/"
polyA_dir <- "../fig1/results/get_polyA_bam_files/lib1"
#### read TLDR-Seq lib1   bed files
lib1.bam.bc.name <- file.path(reproducibility_dir, paste0("bc",1:3, ".bam" ))
lib1.bam.bc.bed <- file.path(out_dir, paste0("lib1_bc",1:3,".bed" ))


command.lines <- paste0("module load bedtools;  bedtools bamtobed -i ", lib1.bam.bc.name ," > ", 
                       lib1.bam.bc.bed)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)

### polyA files
polyA.list <- file.path(polyA_dir, "mapped_with_polyA.list")
no_polyA.list <- file.path(polyA_dir, "mapped_wo_polyA.list")

lib1.bam_pA.bc.bed <- file.path(out_dir, paste0("lib1_pA_bc",1:3,".bed" ))
lib1.bam_nopA.bc.bed <- file.path(out_dir, paste0("lib1_nopA_bc",1:3,".bed" ))


command.lines1 <- paste("grep -f", polyA.list, lib1.bam.bc.bed ," > ", 
                       lib1.bam_pA.bc.bed)

command.lines2 <- paste("grep -f", no_polyA.list, lib1.bam.bc.bed ," > ", 
                       lib1.bam_nopA.bc.bed)

command.lines <- paste(paste(c(command.lines1, command.lines2), collapse=" & \n"), " & wait") 
system(command.lines)


all.lib1.pA.bed <- file.path(out_dir, "all_lib1_pA.bed")
all.lib1.nopA.bed <- file.path(out_dir, "all_lib1_nopA.bed")

command.lines1 <- paste("cat",  paste(lib1.bam_pA.bc.bed, collapse=" "), ">", all.lib1.pA.bed)
command.lines2 <- paste("cat",  paste(lib1.bam_nopA.bc.bed, collapse=" "), ">", all.lib1.nopA.bed)
command.lines <- paste(paste(c(command.lines1, command.lines2), collapse=" & \n"), " & wait") 
system(command.lines)

all.bams <- c(all.lib1.pA.bed, all.lib1.nopA.bed)


### non polyA files


a <- mclapply(all.bams, readBed, mc.cores=2) %>% 
    lapply(function(x) plyranges::filter(x, grepl("chr[0-9][0-9]*$", seqnames) | grepl("chr[XY]$", seqnames))) %>%
    lapply(function(x) plyranges::filter(x, score > 0)) %>%
    lapply(anchor_5p) %>%
    lapply(mutate, width=1) %>%
    setNames(basename(all.bams) %>% str_replace(".bed", ""))  %>%
    lapply(as.data.frame) %>%
    mclapply(mutate, seqnames = seqnames %>% as.character ) %>% 
    lapply(mutate, strand = strand %>% as.character ) %>%
    lapply(function(x) split(x, f=x$strand))  %>%
    unlist(recursive=FALSE) %>%
    setNames(names(.) %>% str_replace("\\.","_")) %>%
    lapply(as_granges) %>%
    mclapply(compute_coverage, mc.cores = 6) %>%
    mclapply(function(x) plyranges::filter(x, score != 0)) %>%
    mclapply(convertGRanges2GPos, mc.cores=6)
for (i in 1:length(a))
{
    seqlengths(a[[i]]) <- chr.lengths[a[[i]] %>% seqnames %>% levels] 
}
mcmapply(export,a,  file.path(out_dir, paste0(names(a), ".bw")), SIMPLIFY=FALSE, mc.cores=6)




#### create design matrix
design  <-  list.files(out_dir) %>% .[grepl("bw",.)] %>% file.path(out_dir, .) %>%
    as.data.frame %>%
    setNames("file.name") %>%
    mutate(strand = ifelse(grepl("\\+",file.name), "BigWigPlus", "BigWigMinus")) %>%
    mutate(sampleIDs = file.name %>% basename() %>% str_replace("_[\\+\\-].*","")) %>%
    reshape2::dcast( sampleIDs ~ strand , value.var = "file.name") %>%
    mutate(condi = ifelse(grepl("nopA", sampleIDs), "nopA", "pA" )) %>%
    mutate(rownames = sampleIDs) %>%
    column_to_rownames("rownames")

bw_plus <- BigWigFileList(design$BigWigPlus)
bw_minus <- BigWigFileList(design$BigWigMinus) 
names(bw_plus) <- names(bw_minus) <- design$sampleIDs



#---Build CTSS object

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome = seqinfo(bsg)[seqnames(bsg) %>% .[grepl("chr[0-9][0-9]*$", .)| grepl("chr[XY]$",.)]],
                       design = design)

CTSSs <- calcPooled(CTSSs, inputAssay="counts")
CTSSs <- calcTPM(CTSSs, inputAssay="counts",
                 outputAssay = "TPM",
                 totalTags = NULL, 
                 outputColumn = 'totalTags')


supported.CTSSs <- subsetBySupport(CTSSs, 
                                   inputAssay="counts", 
                                   outputColumn = "support",
                                   unexpressed=0,
                                   minSamples=0)  

supported.CTSSs <- calcPooled(supported.CTSSs, inputAssay="counts")
saveRDS(supported.CTSSs, file=file.path(out_dir, "CTSS.rds"))
