library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(CAGEfightR)

# Bioconductor packages for differential expression
library(BRGenomics)
library(plyranges)
library(stringr)
library(data.table)
require(genomation)
require(GGally)
require(parallel)

# Bioconductor data packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
require(AnnotationDbi)
require(GenomicFeatures) 
bsg <- BSgenome.Hsapiens.UCSC.hg38
odb <- org.Hs.eg.db
chr.lengths <- seqlengths(bsg) %>% .[grepl("chr[0-9][0-9]*$", names(.)) | grepl("chr[XY]$", names(.))]


out_dir <- "./results/compare_TLDR_SLIC-CAGE/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )

sliccage_out_dir <- "./data/sliccage_bedfiles"
slic_cage_1_name <- file.path(sliccage_out_dir,"hg38.SLICCAGE_CPH4_SiGFB_1.bed")
slic_cage_2_name <- file.path(sliccage_out_dir,"hg38.SLICCAGE_CPH4_SiGFB_2.bed")
slic_cage_3_name <- file.path(sliccage_out_dir,"hg38.SLICCAGE_CPH4_SiGFB_3.bed")

lib1_tss_name="../fig3/results/bedfile_from_bam/lib1.bed"
lib2_tss_name="../fig3/results/bedfile_from_bam/lib2.bed"

nanopore1_name <- "../map_ONTnanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed"
nanopore2_name <- "../map_ONTnanopore/get_bedfiles/nanopore_wt2_curated_stranded.bed"

all.beds <- c(slic_cage_1_name, slic_cage_2_name, slic_cage_3_name,
              lib1_tss_name, lib2_tss_name,
              nanopore1_name, nanopore2_name)
names(all.beds) <- c("CAGE","CAGE","CAGE","tldr_lib1", "tldr_lib2", "ONT", "ONT")

### read beds

a <- mclapply(all.beds, readBed, mc.cores=7) %>% 
    lapply(function(x) plyranges::filter(x, grepl("chr[0-9][0-9]*$", seqnames) | grepl("chr[XY]$", seqnames))) %>%
    lapply(function(x) plyranges::filter(x, score > 0)) %>%
    lapply(anchor_5p) %>%
    lapply(mutate, width=1) %>%
    setNames(names(all.beds))

b <- list(c(a[[1]],  a[[2]],  a[[3]]),
          c(a[[4]]),
          c(a[[5]]),
          c(a[[6]], a[[7]]))

names(b) <- c("CAGE","tldr_lib1", "tldr_lib2", "ONT")


c <- b %>%
    setNames(names(b)) %>% 
    lapply(as.data.frame) %>%
    mclapply(mutate, seqnames = seqnames %>% as.character ) %>% 
    lapply(mutate, strand = strand %>% as.character ) %>%
    lapply(function(x) split(x, f=x$strand))  %>%
    unlist(recursive=FALSE) %>%
    setNames(names(.) %>% str_replace("\\.","_")) %>%
    lapply(as_granges) %>%
    mclapply(compute_coverage, mc.cores = 7) %>%
    mclapply(function(x) plyranges::filter(x, score != 0)) %>%
    mclapply(convertGRanges2GPos, mc.cores=6)
for (i in 1:length(c))
{
    seqlengths(c[[i]]) <- chr.lengths[c[[i]] %>% seqnames %>% levels] 
}
mcmapply(export,c,  file.path(out_dir, paste0(names(c), ".bw")), SIMPLIFY=FALSE, mc.cores=6)




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
                                   minSamples=1)  

supported.CTSSs <- calcPooled(supported.CTSSs, inputAssay="counts")
saveRDS(supported.CTSSs, file=file.path(out_dir, "CTSS.rds"))

rowRanges(CTSSs)

orig.TCs <- clusterUnidirectionally(supported.CTSSs)


#---TC expression

TC.expr.TCs <- quantifyClusters(CTSSs, orig.TCs)
TC.expr.TCs <- calcTPM(TC.expr.TCs, 
                      inputAssay = "counts",
                      outputAssay = "TPM",
                      totalTags = NULL)   




txdb <- loadDb("../data/gencode.no_scaffolds.v39.primary_assembly.annotation.sqlite")

genomeInfo <- bwCommonGenome(plusStrand=bw_plus, minusStrand=bw_minus, method='intersect')
seqlevels(txdb) <- seqlevels(genomeInfo)
supported.TCs <- assignTxType(TC.expr.TCs, txdb, tssUpstream = 100, tssDownstream = 100)
supported.TCs <- assignGeneID(supported.TCs, geneModels = txdb, outputColumn = "geneID")


### max per supported TC ###

all.df.dcast <- as.data.frame(as.matrix(assays(supported.TCs)$TPM)) %>%
    mutate(TC = rownames(supported.TCs)) %>%
    mutate(gene.id = mcols(supported.TCs)$geneID) %>%
    group_by(gene.id) %>%
    dplyr::summarize(CAGE = log(1 + sum(CAGE)), ONT.Nanopore = log(1 + sum(ONT)),
                     TLDR.lib_amp = log(1 + sum(tldr_lib1)),   TLDR.lib_noamp = log(1 + sum(tldr_lib2)))  


cor.stat <- cor.test(all.df.dcast$TLDR.lib_amp, all.df.dcast$CAGE) %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.df.dcast , aes(x=TLDR.lib_amp, y = CAGE)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    xlab("TLDR-Seq amplified library, log10 Tags per million") +
    ylab("SLIC-CAGE, log10 Tags per million") +
    ggtitle(cor.stat)
out_dotplot <- file.path(out_dir,"fig6B_left_panel_CAGE_vs_TLDR_amplifified_library.png")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 


all.df.dcast <- as.data.frame(as.matrix(assays(supported.TCs)$TPM)) %>%
    mutate(TC = rownames(supported.TCs)) %>%
    mutate(gene.id = mcols(supported.TCs)$geneID) %>%
    group_by(gene.id) %>%
    dplyr::summarize(CAGE = log(1 + sum(CAGE)), ONT.Nanopore = log(1 + sum(ONT)),
                     TLDR.lib_amp = log(1 + sum(tldr_lib1)),   TLDR.lib_noamp = log(1 + sum(tldr_lib2)))  



cor.stat <- cor.test(all.df.dcast$TLDR.lib_noamp, all.df.dcast$CAGE) %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.df.dcast , aes(x=TLDR.lib_noamp, y = CAGE)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    xlab("TLDR-Seq non-amplified library, log10 Tags per million") +
    ylab("SLIC-CAGE, log10 Tags per million") +
    ggtitle(cor.stat)
out_dotplot <- file.path(out_dir,"fig6B_right_panel_CAGE_vs_TLDR_non_amplifified_library.png")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 



######

all.df.dcast <- as.data.frame(as.matrix(assays(supported.TCs)$TPM)) %>%
    mutate(TC = rownames(supported.TCs)) %>%
    mutate(gene.id = mcols(supported.TCs)$geneID) %>%
    group_by(gene.id) %>%
    dplyr::summarize(CAGE = log(1 + sum(CAGE)), ONT.Nanopore = log(1 + sum(ONT)),
                     TLDR.lib_amp = log(1 + sum(tldr_lib1)),   TLDR.lib_noamp = log(1 + sum(tldr_lib2)))  



cor.stat <- cor.test(all.df.dcast$ONT.Nanopore, all.df.dcast$CAGE) %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(all.df.dcast , aes(x=ONT.Nanopore, y = CAGE)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    xlab("ONT-Nanopore, log10 Tags per million") +
    ylab("SLIC-CAGE, log10 Tags per million") +
    ggtitle(cor.stat)
out_dotplot <- file.path(out_dir,"figS4_CAGE_vs_ONT_library.png")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 


