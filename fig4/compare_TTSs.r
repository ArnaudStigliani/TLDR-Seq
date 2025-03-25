options(warn=-1)
rm(list=ls())
suppressMessages(library(optparse))

option_list <- list(make_option(c("-f","--foreground"), 
                                  help="path to the foreground file(s) separated by a comma"), 
                    make_option(c("-b","--background"), 
                                 help="path to the background file(s) separated by a comma"),
                    make_option(c("-s", "--score"), type="integer", default=NA,
                                 help="minimum score to define a TTS in foreground"),
                    make_option(c("-t", "--score_bg"), type="integer", default=NA,
                                 help="minimum score to define a TTS in the background"),
                    make_option(c("-o", "--output_table1"), type="character", 
                                 help="path output table fg in bg"),
                    make_option(c("-u", "--output_table2"), type="character",  
                                 help="path output table bg in fg"),
                    make_option(c("-r", "--randomize"), type="character", default="1000",
                                help="sample the TTS in an inerval [p-r,p+r] where p is the TTS position"),
                    make_option(c("-a", "--auto_threshold"), action="store_true",
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

## ## debug
## out_dir <- "../../results/sequencing_july_2023/shared/write_distance_tables_splice"
## ## out_dir <- "../../results/sequencing_july_2023/shared/write_distance_tables_TTS"
## ## fg.tables <- "../../../shared_data/gencode_v39_transcripts.bed"
## ## fg.tables <- "../../data/quantSeq_files/download/mapping_STAR/bed/siGFP_noPAP_in_1.bed"
## fg.tables <- "../../results/sequencing_july_2023/shared/write_distance_tables_splice/tldr_intron_polyA.bed"
## ## bg.tables <- "../../results/sequencing_july_2023/shared/write_distance_tables_TTS/lib1_polyA_reads_stranded.bed"
## bg.tables <- "../../results/sequencing_july_2023/sequencing_july_2023_lib1/get_splice_sites/09-08-2023/nanopore_intron.bed"
## score.min <- NA
## score.min.bg <- NA
## window.sample <- 1000
## output_fg_in_bg <- file.path(out_dir, "test1.tsv")
## output_bg_in_fg <- file.path(out_dir, "test2.tsv")
## auto_threshold <- TRUE


suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(require(genomation))
suppressMessages(library(parallel))

fg.tables <- mclapply(fg.tables, readBed, mc.cores=5) %>% do.call(what = c, args = .) %>% as.data.frame
bg.tables <- mclapply(bg.tables, readBed, mc.cores=5) %>% do.call(what = c, args = .) %>% as.data.frame



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
        mutate(start.tts = ifelse(strand=="+", end, start)) %>%
        mutate(seq=paste0(seqnames,":",  start.tts, ":", strand)) %>%
        dplyr::rename(chr=seqnames) %>% 
        dplyr::select( -end, -width, -score, -start) %>% 
        group_by(seq) %>%
        add_count() %>%
        ungroup %>%
        dplyr::filter(n >= cov )  %>%
        distinct(seq, .keep_all=TRUE) %>%
        dplyr::select(-seq) %>% 
        arrange(chr, start.tts)
    return(bed.formatted)
}

distance.2beds <- function(bed1, bed2)
{
    a <- bed1 %>%
        mutate(name.fg = paste0("fg_", 1:nrow(.) %>% as.character)) %>%
        dplyr::rename(start.tts.fg = start.tts)
    b <- bed2 %>% 
        dplyr::rename(start.tts.bg = start.tts) %>% 
        mutate(name.bg=paste0("bg", 1:nrow(.) %>% as.character))
    ### files need to be sorted !!
    bed1.plus <- a %>% dplyr::filter(strand=="+")
    bed2.plus <- b %>% dplyr::filter(strand=="+")
    bed1.minus <- a %>% dplyr::filter(strand=="-")
    bed2.minus <- b %>% dplyr::filter(strand=="-")
    ## + strand ##
    fg.start.plus <- bed1.plus$start.tts.fg
    bg.start.plus <- bed2.plus$start.tts.bg
    cuts.plus <- c(-Inf, bg.start.plus[-1]-diff(bg.start.plus)/2, Inf)
    bed.plus.closest <- bed1.plus %>%
        mutate(cut1 = cut(fg.start.plus, breaks=cuts.plus,labels=bed2.plus$name.bg, right=TRUE) %>% as.character,
               cut2 = cut(fg.start.plus, breaks=cuts.plus, labels=bed2.plus$name.bg, right=FALSE) %>% as.character) %>%
        mutate(ran = sample(2, nrow(.), replace = TRUE)) %>%
        mutate(name.bg = ifelse(ran==1, cut1, cut2)) %>%
        left_join(bed2.plus %>% dplyr::select(name.bg, start.tts.bg), by="name.bg") %>%
        mutate(distance = start.tts.fg - start.tts.bg)
    ## - strand ##
    fg.start.minus <- bed1.minus$start.tts.fg
    bg.start.minus <- bed2.minus$start.tts.bg
    cuts.minus <- c(-Inf, bg.start.minus[-1]-diff(bg.start.minus)/2, Inf)
    bed.minus.closest <- bed1.minus %>%
        mutate(cut1 = cut(fg.start.minus, breaks=cuts.minus,labels=bed2.minus$name.bg, right=TRUE) %>% as.character,
               cut2 = cut(fg.start.minus, breaks=cuts.minus, labels=bed2.minus$name.bg, right=FALSE) %>% as.character) %>%
        mutate(ran = sample(2, nrow(.), replace = TRUE)) %>%
        mutate(name.bg = ifelse(ran==1, cut1, cut2)) %>%
        left_join(bed2.minus %>% dplyr::select(name.bg, start.tts.bg), by="name.bg") %>%
        mutate(distance = -(start.tts.fg - start.tts.bg))
    return(rbind.data.frame(bed.plus.closest, bed.minus.closest))
}


parallel_distance <- function(fg, bg)
{
    fg.list <- fg %>%
        as.data.frame %>%
        mutate(chr=as.character(chr)) %>%  ### look here if there is a bug
        split(f = .$chr)
    bg.list <- bg %>%
        as.data.frame %>%
        mutate(chr=as.character(chr)) %>% 
        split(f = .$chr)
    name.lists <- intersect(names(fg.list), names(bg.list)) %>% .[grepl("chr[0-9]",.)] %>% .[! grepl("random",.)]
    fg.list <- fg.list[name.lists]
    bg.list <- bg.list[name.lists]
    dist.list <- mcmapply(distance.2beds, fg.list, bg.list, SIMPLIFY=FALSE, mc.cores=50)
    return(dist.list)
}

bed1 <- format_bed(fg.tables, score.min)
bed2 <- format_bed(bg.tables, score.min.bg)


## we keep "n.TSS" to conserve compatibilitty with plot_distance_table_v2.r
## test <-  parallel_distance(bed1, bed2)#  %>% do.call(rbind.data.frame, .) %>% mutate(n.TSS=nrow(.)) 


fg_in_bg.all_chrom <- parallel_distance(bed1, bed2)  %>% do.call(rbind.data.frame, .) %>% mutate(n.TSS=nrow(.))  %>% dplyr::filter(abs(distance) < 100)
bg_in_fg.all_chrom <- parallel_distance(bed2, bed1)  %>% do.call(rbind.data.frame, .)  %>% mutate(n.TSS=nrow(.)) %>% dplyr::filter(abs(distance) < 100)


write.table(fg_in_bg.all_chrom, output_fg_in_bg, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(bg_in_fg.all_chrom, output_bg_in_fg, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


print("generate null distributions")
for (i in window.sample)
{
    bed1.neg <- bed1 %>% mutate(start.tts = start.tts + sample(-i:i, nrow(.), replace=TRUE))
    bed2.neg <- bed2 %>% mutate(start.tts = start.tts + sample(-i:i, nrow(.), replace=TRUE))   
    fg.neg_in_bg.all_chrom <- parallel_distance(bed1.neg, bed2)  %>% do.call(rbind.data.frame, .) %>% mutate(n.TSS=nrow(.))  %>% dplyr::filter(abs(distance) < 100)
    bg.neg_in_fg.all_chrom <- parallel_distance(bed2.neg, bed1)  %>% do.call(rbind.data.frame, .)  %>% mutate(n.TSS=nrow(.)) %>% dplyr::filter(abs(distance) < 100)
    write.table(fg.neg_in_bg.all_chrom, output_fg_in_bg %>% str_replace("\\.tsv",paste0("_",i %>% as.character,"_neg.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(bg.neg_in_fg.all_chrom, output_bg_in_fg %>% str_replace("\\.tsv",paste0("_",i %>% as.character,"_neg.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}


quit()
