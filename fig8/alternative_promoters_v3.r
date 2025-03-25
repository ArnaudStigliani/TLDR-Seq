## CRAN packages for data manipulation and plotting
rm(list=ls())

library(plyranges)
library(tidyverse)

# CAGEfightR and related packages
library(CAGEfightR)
library(ggseqlogo)
# Bioconductor packages for differential expression
library(BRGenomics)
library(data.table)
library(ggpubr)
library(stringr)
library(data.table)
library(GenomicRanges)
require(genomation)
require(GGally)
require(parallel)
library(BiocParallel)
library(memes)
library(entropy)
# Bioconductor data packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
require(AnnotationDbi)
require(GenomicFeatures) 
bsg <- BSgenome.Hsapiens.UCSC.hg38
odb <- org.Hs.eg.db
chr.lengths <- seqlengths(bsg) %>% .[grepl("chr[0-9][0-9]*$", names(.)) | grepl("chr[XY]$", names(.))]


in_dir <- "./results/make_CTSS_object/"
out_dir <- "./results/alternative_promoters_v3/"
dir.create(out_dir, showWarnings=FALSE,recursive=TRUE)

#### read the gene expression files
expression.dir <- "../fig6/results/make_fig6"
LongReads.cpm.manual <- read.table(fgile.path(expression.dir, "tx.counts.gene.tpm.tsv"), sep="\t", header=TRUE) 


Nanocount.placeholder.name <- "../data/Nanocount_placeholder.tsv" ## can be done from a gtf file instead
tr.type <- read.table(Nanocount.placeholder.name, header=TRUE) %>%
    dplyr::select(transcript_name) %>%
    mutate(tr.type = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,(ncol(.)-1 )]) %>%
    mutate(gene.name = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,(ncol(.)-3 )]) %>%
    mutate(gene.id = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,2 ]) %>%
    mutate(transcript.name = transcript_name %>% str_replace("\\|.*",""))  %>%
    dplyr::select(-transcript_name) 

### read annotation
gtf.df <- readGFF("../data/gencode.no_scaffolds.v39.primary_assembly.annotation.sorted.gtf")
gene.annot.df <- gtf.df %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(gene_name, gene_type, geneID=gene_id)


## extract first exons
first.exons.plus.df <- gtf.df %>%
    lapply(unlist) %>%
    data.frame %>%
    dplyr::filter(type=="exon") %>%
    dplyr::filter(strand =="+") %>%
    dplyr::arrange(seqid, start) %>% 
    group_by(transcript_id) %>%
    add_count() %>%
    dplyr::filter(n!=1) %>%
    dplyr::slice(1)


first.exons.minus.df <- gtf.df %>%
    lapply(unlist) %>%
    data.frame %>%
    dplyr::filter(type=="exon") %>%
    dplyr::filter(strand =="-") %>%
    dplyr::arrange(seqid, desc(end))  %>% 
    group_by(transcript_id) %>%
    add_count() %>%
    dplyr::filter(n!=1) %>%
    dplyr::slice(1)

first.exons.df <- bind_rows(first.exons.plus.df, first.exons.minus.df) %>%
    as.data.frame %>%
    arrange(strand, seqid, start) 

#### work on expression file  (maybe useless, see later.)

expression.gene.name.df <-  LongReads.cpm.manual %>%
    dplyr::filter(variable =="lib1") %>%
    dplyr::select(transcript.name = TXNAME, cpm = value) %>%
    left_join(tr.type)  %>%
    arrange(gene.id, desc(cpm)) %>%
    group_by(gene.id) %>%
    dplyr::slice(1) 


#---Read CTSS object
CTSS <- readRDS(file.path(in_dir, "CTSS.rds"))
orig.CTSSs <- calcTotalTags(CTSS)
orig.CTSSs <- calcPooled(orig.CTSSs, inputAssay="counts")
orig.CTSSs <- calcTPM(orig.CTSSs, inputAssay="counts")

TC <-  first.exons.df %>%
    dplyr::filter(transcript_id %in% expression.gene.name.df$transcript.name, seqid != "chrM") %>% 
    dplyr::rename(seqnames = seqid) %>%
    mutate( seqnames = factor(as.character(seqnames))) %>% 
    as_granges %>%
    mutate(start = start -100, end = end + 100)




overlaps_list <- findOverlaps(TC, drop.self=TRUE, select="first")
overlap.all <- !is.na(overlaps_list)
overlap.shift <- c(FALSE, overlap.all[-length(overlap.all)])
overlap.rm <- which(overlap.all & overlap.shift)



TC.keep <- TC[-c(overlap.rm)] 
seqlevels(TC.keep) <-  seqlevels(orig.CTSSs)
seqinfo(TC.keep) <-  seqinfo(orig.CTSSs)

TC.expr.TCs <- quantifyClusters(orig.CTSSs, TC.keep)
TC.expr.TCs <- calcTPM(TC.expr.TCs,  inputAssay = "counts",  outputAssay = "TPM",  totalTags = NULL)   
rownames(TC.expr.TCs) <-  paste0(seqnames(TC.expr.TCs), ":", start(TC.expr.TCs), "-" ,end(TC.expr.TCs), ";", strand(TC.expr.TCs))


a <- rowRanges(orig.CTSSs) %>%
    mutate(pA = data.frame(as.matrix(assays(orig.CTSSs)$count))$all_lib1_pA,
           nopA = data.frame(as.matrix(assays(orig.CTSSs)$count))$all_lib1_nopA)



b <- TC.expr.TCs




create_distributions <- function(stitched)
{
    stitched.df <- stitched %>%
        as.data.frame %>%
        dplyr::select(pos, pA, nopA)
    res <- list(pA=NULL, nopA=NULL)
    pA <- unlist(mapply(function(x, y) rep(x, y), stitched.df$pos, stitched.df$pA, SIMPLIFY=FALSE))
    nopA <- unlist(mapply(function(x, y) rep(x, y), stitched.df$pos, stitched.df$nopA, SIMPLIFY=FALSE))
    res[["pA"]] <- pA
    res[["nopA"]] <- nopA
    return(res)
}




ov <- findOverlaps(b, a)
q.list <- split(subjectHits(ov), queryHits(ov))




overlaps_list <- mclapply(q.list, function(x,y) y[x], a, mc.cores=50)
names(overlaps_list) <- as.numeric(names(overlaps_list)) %>% 
    {paste0(rownames(b)[.], "#",  mcols(b[.])$gene_name,"#",  mcols(b[.])$gene_type)}

index <- grep("S100A16", names(overlaps_list))
## v <- overlaps_list[[7136]]
name.test <- names(overlaps_list[7136])

overlaps_list.distrib <- mclapply(overlaps_list, create_distributions, mc.cores=50)

width.range <- overlaps_list.distrib %>% mclapply(function(x) max(x$pA, x$nopA) - min(x$pA, x$nopA) , mc.cores=50) %>% unlist
exp.range <- overlaps_list.distrib %>% mclapply(function(x) (length(x$pA) > 10) & (length(x$nopA) > 10) , mc.cores=50) %>% unlist







overlaps_list.distrib.filtered <- overlaps_list.distrib[(width.range > 20) & exp.range]##  %>%
    ## .[(mclapply(., function(x) abs(median(x$pA)  - median(x$nopA)) , mc.cores=50) %>% unlist) >10 ]

wilcox.list <- mclapply(overlaps_list.distrib.filtered,
                        function(x) return(c(wilcox.test(x$pA, x$nopA)$p.value,
                                             wilcox.test(x$pA, x$nopA)$statistic,
                                             wilcox.test(x$pA, x$nopA)$statistic/(length(x$pA) * length(x$nopA)))),
                        mc.cores=50)

wilcox.df <- cbind.data.frame( names(wilcox.list),
                              do.call(rbind.data.frame, wilcox.list)) %>%
    setNames(c("TC", "p.value", "D", "D.norm")) %>%
    arrange(desc(abs(D.norm - 0.5))) %>%
    remove_rownames() %>% 
    mutate(padj=p.adjust(p.value, method="fdr")) %>%
    separate(TC, into=c("TC.2","gene.name", "biotype"), remove=TRUE,  sep="#")  %>%
    mutate(group = "All") %>%
    mutate(biotype =str_replace_all(biotype, "_", "\n"))
wilcox.df.significant <-  wilcox.df %>% dplyr::filter(padj < 0.05) %>%  mutate(group = "FDR < 0.05")

l1 <- wilcox.df.significant$gene.name


write.table(wilcox.df.significant %>% mutate(biotype = str_replace_all(biotype, "\n", "_")), file=file.path(out_dir, "wilcox_significant.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

write.table(wilcox.df %>% mutate(biotype = str_replace_all(biotype, "\n", "_")), file=file.path(out_dir, "wilcox_bg.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


wilcox.df.all <- bind_rows(wilcox.df, wilcox.df.significant)



g <- ggplot(wilcox.df.all, 
            aes(x = biotype, y=..count..))+
    scale_fill_brewer(palette = "Set2") +
    geom_bar(position ="dodge") +
    theme_bw() +
    facet_wrap(~group, ncol=1, scales="free") 
ggsave( filename = file.path(out_dir, "TSS_distribution.pdf"), g, width = 5, height = 6)


g <- ggplot(wilcox.df.all %>% dplyr::filter(group !="All"),  
            aes(x = biotype, y=..count..))+
    scale_fill_brewer(palette = "Set2") +
    geom_bar(position ="dodge") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw() 
ggsave( filename = file.path(out_dir, "TSS_distribution_significant.pdf"), g, width = 4, height = 5)





###### entropy !!

list.to.boostrap <- overlaps_list.distrib.filtered %>%
    setNames(names(.) %>% str_replace("#.*","")) %>% 
    .[wilcox.df.significant$TC.2]

popA_data <- list.to.boostrap %>% lapply(function(x) x$pA)
popB_data <- list.to.boostrap %>% lapply(function(x) x$nopA)

n.c <- list.to.boostrap %>% lapply(function(x) min(length(x$pA), length(x$nopA))) %>% unlist
n.c2 <- n.c/2


bootstrap_entropy <- function(data, R = 1000, n) {
  ret <- numeric(R)
  sample.size <- length(data)
  if(sample.size == n)
  {
      freqs <- table(data)
      entropy_values <- entropy.empirical(freqs, unit = "log2")  # Compute entropy
      H <- entropy_values 
      H_max <- 1  ## log2(length(freqs)) 
      ret <- H / H_max
  }else{
      for (i in 1:R) {
          resampled_data <- sample(data, size = n, replace = FALSE)  # Bootstrap sample
          freqs <- table(resampled_data)
          entropy_values <- entropy.empirical(freqs, unit = "log2")  # Compute entropy
          H <- entropy_values 
          H_max <- 1 # log2(length(freqs))
          ret[i] <- H / H_max
      }
      ret <- mean(ret)
  }
  return(ret)  # Return mean bootstrap entropy
}

# Compute entropy for each individual in both populations with bootstrapping
entropy_A_boot <- mcmapply(bootstrap_entropy, popA_data, R = 1000, n.c2, SIMPLIFY = FALSE, mc.cores=50) %>%
    unlist
entropy_B_boot <- mcmapply(bootstrap_entropy, popB_data, R = 1000, n.c2, SIMPLIFY = FALSE, mc.cores=50) %>%
    unlist

# Print results
entropy_boot_results <- data.frame( pA = entropy_A_boot, nopA = entropy_B_boot) %>%
    rownames_to_column("TSS")

entropy_boot_results.melt <- entropy_boot_results %>% reshape2::melt() 
## print(entropy_boot_results)

pval.wilcox <- wilcox.test(entropy_boot_results$pA, entropy_boot_results$nopA, paired = TRUE)$p.value %>%
                                                                                            formatC(1,format="e") %>%
                                                                                            paste("p =",.)

g <- ggplot(entropy_boot_results.melt, aes( x = value, colour = variable)) +
    geom_density() +
    theme_bw() +
    scale_colour_brewer(palette = "Set2") +
    geom_text(data=data.frame(TSS ="", variable = ".", value = ".", label = pval.wilcox),
              aes(label = label), x=2.5, y=0.38, color = "black")
ggsave( filename = file.path(out_dir, "entropy_2_distributions.pdf"), g, width = 6, height = 4)

g <- ggplot(entropy_boot_results %>%
            mutate(diff = pA - nopA),
            aes( x = diff)) +
    xlab("Paired entropy difference: entropy(pA+) - entropy(pA-)") +
    geom_density() +
    theme_bw() 
ggsave( filename = file.path(out_dir, "entropy_difference.pdf"), g, width = 6, height = 4)



g <- ggplot(entropy_boot_results.to_plot, aes(y = ))

### map reads to TSS
old_dir <- "../../results/reviews/alternative_promoters/"

all.lib1.pA.bed <- file.path(old_dir, "all_lib1_pA.bed")
all.lib1.nopA.bed <- file.path(old_dir, "all_lib1_nopA.bed")
all.bed6 <- c(all.lib1.pA.bed, all.lib1.nopA.bed)

bed6.df <- mclapply(all.bed6, read.table, header=FALSE, sep="\t", mc.cores=2) %>%
    mapply(mutate, . , tail =c("pA", "nopA"), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .)
read.pA.df <- bed6.df %>% dplyr::select(name = V4, tail) %>%
    distinct


bed12.all.name <- file.path(old_dir, "all_reads.bed12")


bed12 <- readBed(bed12.all.name) %>% 
    plyranges::filter( grepl("chr[0-9][0-9]*$", seqnames) | grepl("chr[XY]$", seqnames)) %>% 
    plyranges::filter( score > 0) %>%
    mutate(genomic.length = end - start) %>%
    anchor_5p %>% 
    mutate(width=1) 



TC.to.overlap <-  wilcox.df %>%
    dplyr::select(TC.2) %>%
    separate(TC.2, into=c("seqnames", "range", "strand"), remove=TRUE,  sep="[:;]")  %>%
    separate(range, into=c("start", "end"), remove=TRUE,  sep="-")   %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    as_granges()
seqlevels(TC.to.overlap) <- seqlevels(bed12)
seqinfo(TC.to.overlap) <- seqinfo(bed12)

ov <- findOverlaps(TC.to.overlap, bed12)
q.list <- split(subjectHits(ov), queryHits(ov))
overlaps_list <- mclapply(q.list, function(x,y) y[x], bed12, mc.cores=50)

overlaps_list.df <- overlaps_list %>%
    lapply(as.data.frame) %>%
    mapply(mutate, ., TC=wilcox.df$TC.2, SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    mutate(exon.length.list = str_split(blockSizes,",") %>% lapply(as.numeric) ) %>%
    mutate(readLength = sapply(exon.length.list, sum)) %>%
    dplyr::select(-exon.length.list) %>%
    left_join(read.pA.df)    %>%
    drop_na() %>% 
    dplyr::select(genomic.length, readLength, tail, n.exons = blockCount, TC, name) 



overlaps_list.df.summarized.by.tail <- overlaps_list.df %>%
    group_by(TC, tail) %>%
    summarize(readLength = median(readLength), genomic.length = median(genomic.length), n.exon = median(n.exons))  %>%
    melt() %>%
    reshape2::dcast(TC + variable ~ tail, value.var ="value")  
    ## mutate(shift = ifelse(shift >0 , "positive\nshift", "negative\nshift")) %>%


diff.distrib.median.df <- overlaps_list.df.summarized.by.tail %>%
    mutate(diff = pA - nopA) %>%
    dplyr::filter(TC %in% wilcox.df.significant$TC.2) %>% 
    group_by(variable) %>%
    dplyr::summarize(median = as.integer(median(diff)))


x1 <- overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "readLength",
                                                            TC %in% wilcox.df.significant$TC.2) %>%
    mutate(diff=pA- nopA) %>%
    .$diff
g <- ggplot(overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "readLength",
                                                                  TC %in% wilcox.df.significant$TC.2),
            aes(x= pA - nopA)) +
    geom_density() +
    theme_bw() +
    ggtitle("Read Length difference pA+ - pA-") +
    geom_vline(xintercept = diff.distrib.median.df %>% dplyr::filter(variable =="readLength") %>% .$median,
               colour ="blue", linetype ="dashed") +
    scale_x_continuous(breaks = c(pretty(x1),
                                  diff.distrib.median.df %>% dplyr::filter(variable =="readLength") %>% .$median ))
ggsave( filename = file.path(out_dir, "read_length_density_diff.pdf"), g, width = 4, height = 4)


x1 <- overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "genomic.length",
                                                            TC %in% wilcox.df.significant$TC.2) %>%
    mutate(diff=pA- nopA) %>%
    .$diff
labs <-  c(pretty(x1), diff.distrib.median.df %>% dplyr::filter(variable =="genomic.length") %>%
                       .$median)  %>% as.character
labs[1] <- ""
g <- ggplot(overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "genomic.length",
                                                                  TC %in% wilcox.df.significant$TC.2),
            aes(x= pA - nopA)) +
    geom_density() +
    theme_bw() +
    ggtitle("Read Length difference pA+ - pA-") +
    geom_vline(xintercept = diff.distrib.median.df %>% dplyr::filter(variable =="genomic.length") %>% .$median,
               colour ="blue", linetype ="dashed") +
    scale_x_continuous(breaks = c(pretty(x1),
                                  diff.distrib.median.df %>% dplyr::filter(variable =="genomic.length") %>% .$median ),
                       labels = labs )
ggsave( filename = file.path(out_dir, "genomic_length_density_diff.pdf"), g, width = 4, height = 4)



x1 <- overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "n.exon",
                                                            TC %in% wilcox.df.significant$TC.2) %>%
    mutate(diff=pA- nopA) %>%
    .$diff
g <- ggplot(overlaps_list.df.summarized.by.tail %>% dplyr::filter(variable == "n.exon",
                                                                  TC %in% wilcox.df.significant$TC.2),
            aes(x= pA - nopA)) +
    geom_density() +
    theme_bw() +
    ggtitle("n exon difference pA+ - pA-") +
    geom_vline(xintercept = diff.distrib.median.df %>% dplyr::filter(variable =="n.exon") %>% .$median,
               colour ="blue", linetype ="dashed") +
    scale_x_continuous(breaks = c(pretty(x1),
                                  diff.distrib.median.df %>% dplyr::filter(variable =="n.exon") %>% .$median ))
ggsave( filename = file.path(out_dir, "n_exon_density_diff.pdf"), g, width = 4, height = 4)


### get readdy to create motif
genome.name <- "../../../shared_data/human_genome/GRCh38.primary_assembly.genome.fa"
genome <- readDNAStringSet(genome.name) %>%
    setNames( names(.) %>% str_replace(" .*", ""))


ribo.TC.list <- wilcox.df.significant %>%
    dplyr::filter(grepl( "(^RPS)|(^RPL)", gene.name))  %>%
    .$TC.2

overlaps_list.ribo <- overlaps_list.distrib.filtered %>%
    setNames(names(.) %>% str_replace("#.*", "")) %>%
    .[ribo.TC.list] %>%
    lapply(unlist, recursive = TRUE) %>%
    lapply(function(x) data.frame(center = x, tail = names(x))) %>%
    mapply(mutate, ., TC=names(.), SIMPLIFY=FALSE)  %>%
    do.call(rbind.data.frame, .) %>%
    mutate(tail = tail %>% str_replace("pA.*", "pA")) %>%
    mutate(start = center - 15,
           end = center + 15,
           seqnames = TC %>% str_replace(":.*",""),
           strand = TC %>% str_replace(".*;","")) %>%
    dplyr::select(seqnames, start, end, TC, strand, tail)  %>%
    remove_rownames()  %>%
    split(f = paste(.$TC, .$tail, sep="_")) %>% 
    lapply(as_granges)




seq.logo <- overlaps_list.ribo %>%
    mclapply(get_sequence, genome = genome, mc.cores=50) 

nt.freq <- function(seq.logo)
{
    freq <- seq.logo %>%
        as.character()  %>%
        str_split(., "") %>%
        do.call(rbind, .) %>%
        as.data.frame %>%
        as.list()  %>%
        lapply(paste, collapse = "") %>%
        lapply(function(x) sapply(list(A="A",C="C",G="G",T="T"), function(y) str_count( x ,y ))) %>%
        do.call(cbind.data.frame, .) %>%
        t %>%
        as.data.frame %>%
        remove_rownames %>%
        mutate(position = 1:31) %>%
        reshape2::melt(id.vars=c("position"))
    return(freq)
}

freq.logo <- mclapply(seq.logo, nt.freq, mc.cores=50)


all.freq <- mapply(mutate, freq.logo, seq_group = names(freq.logo), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    dplyr::rename(nt=variable) %>%
    group_by(seq_group, position) %>%
    mutate(n.seq = sum(value)) %>% 
    group_by(nt, seq_group) %>%
    mutate(nt.abundance = sum(value)) %>%
    ungroup %>%
    group_by(seq_group) %>%
    mutate(nt.tot=sum(value)) %>%
    ungroup %>% 
    mutate(freq=value/n.seq, freq.norm = value/(n.seq * (nt.abundance/nt.tot))) %>%
    mutate(position = position - 16) 

## write.table(all.freq, file.path(out_dir, "all_freq_unique_3support.tsv"), row.names=FALSE, sep="\t", quote=FALSE)        

## all.freq <- read.table(file.path(out_dir, "all_freq_unique_3support.tsv"), sep="\t", header=TRUE)


all.logo.towrite <- all.freq %>%
    mutate(freq = (freq.norm))  %>%
    dplyr::select(seq_group, position, nt, freq) %>%
    group_by(position, seq_group) %>%
    mutate(freq = freq + 0.01) %>%
    mutate(tot = sum(freq), freq = (freq/tot) ) %>%
    ungroup %>% 
    reshape2::dcast(seq_group  + position ~ nt, value.var = "freq") %>% 
    remove_rownames %>% 
    split(f=.$seq_group) %>%
    lapply(dplyr::select, A,C,G,T) 


p_list <- all.logo.towrite %>%
    lapply(function(x) t(as.matrix(x) *  matrix(rep(2 + rowSums(( as.matrix(x)) *log (as.matrix(x), 2)),4),ncol=4)))


## colnames(all.logo.towrite$SLIC_CAGE.TSS) <- c(2:10) %>% as.character

g <- ggseqlogo(p_list,
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=2, scales='free_x') +
    scale_x_continuous(breaks=seq(1,31, 2), labels=seq(-15,15, 2))  
ggpubr::ggexport(g, filename = file.path(out_dir, "logo.pdf"), width = 10, height =50 )

##### do an averaged version
###### start from the granges object an reduce it

overlaps_list.ribo.summarized <- overlaps_list.ribo %>%
    lapply(as.data.frame) %>%
    do.call(rbind.data.frame, .) %>%
    group_by(TC, tail, start, strand) %>%
    add_count()  %>%
    ungroup %>%
    group_by(TC, tail) %>%
    mutate(max.count = max(n)) %>%
    ungroup() %>%
    dplyr::filter(n > 0.4 * max.count)  %>%
    distinct() %>% 
    split(f =  .$tail)  %>%
    lapply(function(x) mutate(x, min.dist = c(diff(start)[1],diff(start)))) %>%
    lapply(function(x) mutate(x, cs=cumsum(as.numeric(abs(min.dist)>5)))) %>%
    lapply(group_by, cs) %>%
    lapply( arrange, cs, desc(n)) %>%
    lapply(dplyr::slice, 1) %>% 
    lapply(as_granges)

    



seq.logo.summarized <- overlaps_list.ribo.summarized %>%
    mclapply(get_sequence, genome = genome, mc.cores=50) 


freq.logo.summarized <- mclapply(seq.logo.summarized, nt.freq, mc.cores=50)


all.freq.summarized <- mapply(mutate, freq.logo.summarized, seq_group = names(freq.logo.summarized), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    dplyr::rename(nt=variable) %>%
    group_by(seq_group, position) %>%
    mutate(n.seq = sum(value)) %>% 
    group_by(nt, seq_group) %>%
    mutate(nt.abundance = sum(value)) %>%
    ungroup %>%
    group_by(seq_group) %>%
    mutate(nt.tot=sum(value)) %>%
    ungroup %>% 
    mutate(freq=value/n.seq, freq.norm = value/(n.seq * (nt.abundance/nt.tot))) %>%
    mutate(position = position - 16) 

## write.table(all.freq, file.path(out_dir, "all_freq_unique_3support.tsv"), row.names=FALSE, sep="\t", quote=FALSE)        

## all.freq <- read.table(file.path(out_dir, "all_freq_unique_3support.tsv"), sep="\t", header=TRUE)


all.logo.towrite.summarized <- all.freq.summarized %>%
    mutate(freq = (freq.norm))  %>%
    dplyr::select(seq_group, position, nt, freq) %>%
    group_by(position, seq_group) %>%
    mutate(freq = freq + 0.01) %>%
    mutate(tot = sum(freq), freq = (freq/tot) ) %>%
    ungroup %>% 
    reshape2::dcast(seq_group  + position ~ nt, value.var = "freq") %>% 
    remove_rownames %>% 
    split(f=.$seq_group) %>%
    lapply(dplyr::select, A,C,G,T) 


p_list.summarized <- all.logo.towrite.summarized %>%
    lapply(function(x) t(as.matrix(x) *  matrix(rep(2 + rowSums(( as.matrix(x)) *log (as.matrix(x), 2)),4),ncol=4)))


## colnames(all.logo.towrite$SLIC_CAGE.TSS) <- c(2:10) %>% as.character

g <- ggseqlogo(p_list.summarized,
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=2, scales='free_x') +
    scale_x_continuous(breaks=seq(1,31, 2), labels=seq(-15,15, 2))  
ggpubr::ggexport(g, filename = file.path(out_dir, "logo_summarized.pdf"), width = 10, height =3 )
