rm(list=ls())
library(optparse)

option_list <- list(make_option(c("-f","--foreground"), 
                                  help="path to the foreground file(s) separated by a comma"), 
                    make_option(c("-b","--background"), 
                                 help="path to the background file(s) separated by a comma"),
                    make_option(c("-s", "--score"), type="integer", default=NA,
                                 help="minimum score to define a TSS in foreground"),
                    make_option(c("-t", "--score_bg"), type="integer", default=NA,
                                 help="minimum score to define a TSS in the background"),
                    make_option(c("-o", "--output_table1"), type="character", 
                                 help="path output table fg in bg"),
                    make_option(c("-u", "--output_table2"), type="character",  
                                 help="path output table bg in fg"),
                    make_option(c("-r", "--randomize"), type="character", default="1000",
                                help="sample the TSS in an inerval [p-r,p+r] where p is the TSS position"),
                    make_option(c("-a", "--auto_threshold"), action="store_true",default=FALSE,
                                help="decides the thresholds according to sequencing depth")
                    )



parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args


fg.tables <- unlist(strsplit(opt$foreground, ","))
bg.tables <- unlist(strsplit(opt$background, ","))
score.min <- as.numeric(opt$score)
score.min.bg <- as.numeric(opt$score_bg)
output_fg_in_bg <- opt$output_table1
output_bg_in_fg <- opt$output_table2
window.sample <- as.numeric(unlist(strsplit(opt$randomize, ",")))
auto_threshold <- opt$auto_threshold


## debug
## out_dir <- "../../results/sequencing_july_2023/shared/reviews/write_distance_tables"
## fg.tables <- "../../results/sequencing_july_2023/sequencing_july_2023_lib1/bedfile_from_bam/09-08-2023/reads_genome_stranded.bed"
## fg.tables <- c("/maps/projects/scarball/people/fzh976/binf-isilon/shared_data/gencode_v39_transcripts.bed")
## fg.tables <- "/maps/projects/sandelin_main/scratch/TERA_Seq_data/data/process_bam_files/reads.1.sanitize.toGenome.sorted_5prime.bed"
## bg.tables <- c("/maps/projects/scarball/people/fzh976/binf-isilon/shared_data/gencode_v39_transcripts.bed", "/maps/projects/scarball/people/fzh976/binf-isilon/shared_data/gencode_v39_transcripts.bed")
## bg.tables <- file.path("../../../polyA_CAGE/results/20200315_no_trimming/compare_TSSs/reanalysis_20230110/2023-04-19", c("hg38.SLICCAGE_CPH4_SiGFB_1.bed", "hg38.SLICCAGE_CPH4_SiGFB_2.bed", "hg38.SLICCAGE_CPH4_SiGFB_3.bed"))
## bg.tables <- "../../data/nanopore/get_bedfiles/nanopore_wt1_curated_stranded.bed"
## bg.tables <- "../../../shared_data/human_genome/gencode_v39_albin_selected_transcripts_for_TLDR.bed"
## out_dir <- "../../results/sequencing_july_2023/shared/write_distance_tables_splice"
## score.min <- NA
## score.min.bg <- 1
## window.sample <- 1000
## output_fg_in_bg <- file.path(out_dir, "test1.tsv")
## output_bg_in_fg <- file.path(out_dir, "test2.tsv")
## fg.tables <- "../../results/sequencing_july_2023/sequencing_july_2023_lib1/get_splice_sites/09-08-2023/tldr_intron.bed"
## bg.tables <- "../../results/sequencing_july_2023/sequencing_july_2023_lib1/get_splice_sites/09-08-2023/nanopore_intron.bed"
## auto_threshold <- TRUE



library(tidyverse)
library(stringr)
require(genomation)
library(parallel)

fg.tables <- lapply( fg.tables, readBed) %>% do.call(what = c, args = .) %>% as.data.frame
bg.tables <- lapply(bg.tables, readBed) %>% do.call(what = c, args = .) %>% as.data.frame

if(auto_threshold)
{
    if(is.na(score.min))
    {
        score.min <- round(nrow(fg.tables)/10**6) %>% ifelse(. < 3, 3, .)
    }
    if(is.na(score.min.bg))
    {
        score.min.bg <- round(nrow(bg.tables)/10**6) %>% ifelse(. < 3, 3, .)
    }
}

set.seed(1)

format_bed <- function(bed, cov) 
{
    bed.formatted <- bed %>%
        mutate(start.tss = ifelse(strand=="+", start, end)) %>%
        mutate(seq=paste0(seqnames,":",  start.tss, ":", strand)) %>%
        dplyr::rename(chr=seqnames) %>% 
        dplyr::select( -end, -width, -score, -start) %>% 
        group_by(seq) %>%
        add_count() %>%
        ungroup %>%
        dplyr::filter(n >= cov )  %>%
        distinct(seq, .keep_all=TRUE) %>%
        dplyr::select(-seq) %>% 
        arrange(chr, start.tss)
    return(bed.formatted)
}

distance.2beds <- function(bed1, bed2)
{
    a <- bed1 %>%
        mutate(name.fg = paste0("fg_", 1:nrow(.) %>% as.character)) %>%
        dplyr::rename(start.tss.fg = start.tss)
    b <- bed2 %>% 
        dplyr::rename(start.tss.bg = start.tss) %>% 
        mutate(name.bg=paste0("bg", 1:nrow(.) %>% as.character))
    ### files need to be sorted !!
    bed1.plus <- a %>% dplyr::filter(strand=="+")
    bed2.plus <- b %>% dplyr::filter(strand=="+")
    bed1.minus <- a %>% dplyr::filter(strand=="-")
    bed2.minus <- b %>% dplyr::filter(strand=="-")
    ## + strand ##
    if(nrow(bed1.plus)==0 | nrow(bed2.plus)==0)
    {
        bed.plus.closest <- data.frame()
    } else {
    fg.start.plus <- bed1.plus$start.tss.fg
    bg.start.plus <- bed2.plus$start.tss.bg
    cuts.plus <- c(-Inf, bg.start.plus[-1]-diff(bg.start.plus)/2, Inf)
    bed.plus.closest <- bed1.plus %>%
        mutate(cut1 = cut(fg.start.plus, breaks=cuts.plus,labels=bed2.plus$name.bg, right=TRUE) %>% as.character,
               cut2 = cut(fg.start.plus, breaks=cuts.plus, labels=bed2.plus$name.bg, right=FALSE) %>% as.character) %>%
        mutate(ran = sample(2, nrow(.), replace = TRUE)) %>%
        mutate(name.bg = ifelse(ran==1, cut1, cut2)) %>%
        left_join(bed2.plus %>% dplyr::select(name.bg, start.tss.bg), by="name.bg") %>%
        mutate(distance = start.tss.fg - start.tss.bg)
    }
    ## - strand ##
    if(nrow(bed1.minus)==0 | nrow(bed2.minus)==0)
    {
        bed.minus.closest <- data.frame()
    } else {
        fg.start.minus <- bed1.minus$start.tss.fg
    bg.start.minus <- bed2.minus$start.tss.bg
    cuts.minus <- c(-Inf, bg.start.minus[-1]-diff(bg.start.minus)/2, Inf)
    bed.minus.closest <- bed1.minus %>%
        mutate(cut1 = cut(fg.start.minus, breaks=cuts.minus,labels=bed2.minus$name.bg, right=TRUE) %>% as.character,
               cut2 = cut(fg.start.minus, breaks=cuts.minus, labels=bed2.minus$name.bg, right=FALSE) %>% as.character) %>%
        mutate(ran = sample(2, nrow(.), replace = TRUE)) %>%
        mutate(name.bg = ifelse(ran==1, cut1, cut2)) %>%
        left_join(bed2.minus %>% dplyr::select(name.bg, start.tss.bg), by="name.bg") %>%
        mutate(distance = -(start.tss.fg - start.tss.bg))
    }
    return(rbind.data.frame(bed.plus.closest, bed.minus.closest))
}


parallel_distance <- function(fg, bg)
{
    fg.list <- fg %>%
        as.data.frame %>%
        mutate(chr=as.character(chr)) %>% 
        split(f = .$chr) 
    bg.list <- bg %>%
        as.data.frame %>%
        mutate(chr=as.character(chr)) %>% 
        split(f = .$chr)
    name.lists <- intersect(names(fg.list), names(bg.list)) %>% .[grepl("chr[0-9]",.)] %>% .[! grepl("random",.)]
    fg.list <- fg.list[name.lists]
    bg.list <- bg.list[name.lists]
    dist.list <- mcmapply(distance.2beds, fg.list, bg.list, SIMPLIFY=FALSE, mc.cores=25)
    return(dist.list)
}



bed1 <- format_bed(fg.tables, score.min)
bed2 <- format_bed(bg.tables, score.min.bg)

test <-  parallel_distance(bed1, bed2)  


fg_in_bg.all_chrom <- parallel_distance(bed1, bed2)  %>% do.call(rbind.data.frame, .) %>% mutate(n.TSS=nrow(.))  %>% dplyr::filter(abs(distance) < 100)
bg_in_fg.all_chrom <- parallel_distance(bed2, bed1)  %>% do.call(rbind.data.frame, .)  %>% mutate(n.TSS=nrow(.)) %>% dplyr::filter(abs(distance) < 100)

write.table(fg_in_bg.all_chrom, output_fg_in_bg, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(bg_in_fg.all_chrom, output_bg_in_fg, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


print("generate null distributions")
for (i in window.sample)
{
    bed1.neg <- bed1 %>% mutate(start.tss = start.tss + sample(-i:i, nrow(.), replace=TRUE))
    bed2.neg <- bed2 %>% mutate(start.tss = start.tss + sample(-i:i, nrow(.), replace=TRUE))   
    fg.neg_in_bg.all_chrom <- parallel_distance(bed1.neg, bed2)  %>% do.call(rbind.data.frame, .) %>% mutate(n.TSS=nrow(.))  %>% dplyr::filter(abs(distance) < 100)
    bg.neg_in_fg.all_chrom <- parallel_distance(bed2.neg, bed1)  %>% do.call(rbind.data.frame, .)  %>% mutate(n.TSS=nrow(.)) %>% dplyr::filter(abs(distance) < 100)
    write.table(fg.neg_in_bg.all_chrom, output_fg_in_bg %>% str_replace("\\.tsv",paste0("_",i %>% as.character,"_neg.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(bg.neg_in_fg.all_chrom, output_bg_in_fg %>% str_replace("\\.tsv",paste0("_",i %>% as.character,"_neg.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}


quit()
