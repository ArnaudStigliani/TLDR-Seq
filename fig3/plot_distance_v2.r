rm(list=ls())
library(optparse)

option_list <- list(make_option(c("-p","--pos"), type="character", 
                                help="path to the foreground table"), 
                    make_option(c("-n","--neg"), type="character" , default="",
                                help="path to the background table(s) separated by a comma"),
                    make_option(c("-o", "--output_dir"), type="character", default="",
                                help="path output directory"),
                    make_option(c("-y", "--ymax"), type="integer", default=NA,
                                 help="Set manually the y axis max coordinate")
                    )


parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args




pos.tables <- opt$pos
neg.tables <- unlist(strsplit(opt$neg, ","))
out_dir <- opt$output_dir
ymax <- opt$ymax

## ## debug
## out_dir <- "../../results/sequencing_july_2023/shared/plot_distance_v2"
## pos.tables <- "../../results/sequencing_july_2023/shared/write_distance_tables/lib1_in_gencode_a.tsv"
## neg.tables <- c("../../results/sequencing_july_2023/shared/write_distance_tables/lib1_in_gencode_a_50_neg.tsv",
##                 "../../results/sequencing_july_2023/shared/write_distance_tables/lib1_in_gencode_a_100_neg.tsv",
##                 "../../results/sequencing_july_2023/shared/write_distance_tables/lib1_in_gencode_a_1000_neg.tsv")



dir.create(out_dir, showWarnings=FALSE)

library(tidyverse)
library(stringr)
require(genomation)
library(parallel)

basename_out_file <- basename(out_dir) %>% str_replace("\\.tsv", "")


calc.frac <- function(object.dist, nTSS)
{
    object.dist.th <- object.dist %>%
        abs %>% 
        table %>% 
        as.data.frame %>%
        dplyr::rename(dist=".", n=Freq) 
    object.dist.th <- transform(object.dist.th, cumFreq = cumsum(n)) %>%
        mutate(frac=cumFreq/nTSS)
    return(object.dist.th)
}


### plot distributions

pos.df <- read.table(pos.tables, sep="\t", header=TRUE) %>%
    mutate(id="pos")

## neg.id <- neg.tables %>% str_replace(".*_(\\d*_neg).*","\\1")
## neg.df <- lapply(neg.tables, read.table,  sep="\t", header=TRUE) %>%
##     setNames(neg.id) %>% 
##     {do.call(bind_rows, args=list(.,  .id="id"))}

## tab.all <- rbind.data.frame( pos.df, neg.df)
tab.all <- rbind.data.frame( pos.df)


## g <- ggplot(tab.all %>% dplyr::filter(id == "pos"), aes(x=distance, colour=id, y=after_stat(count))) +
##     geom_density()  +
##     theme_bw()
## ggpubr::ggexport(g, filename = file.path(out_dir, "distance_distribution.pdf"), width = 5, height =5 )


## g <- ggplot(tab.all %>% dplyr::filter(id == "pos"), aes(x=distance, colour=id)) +
##     geom_density(bw=0.5)  +
##     theme_bw()
## ggpubr::ggexport(g, filename = file.path(out_dir, "distance_distribution_density.pdf"), width = 5, height =5 )


tab.histo.pos  <- tab.all  %>% dplyr::filter(id =="pos", abs(distance) <=50 )  
## tab.histo.neg  <- tab.all  %>% dplyr::filter(id =="50_neg", abs(distance) <=50)  


if(is.na(ymax))
{
    g <- ggplot(tab.histo.pos , aes(x=distance, y=after_stat(count))) +
        geom_histogram( binwidth=1)  +
        theme_bw() +
        xlim(-50,50)
    ggpubr::ggexport(g, filename = file.path(out_dir, paste0(basename_out_file, "_distHist.pdf")), width = 5, height =5 )
}else
{
    g <- ggplot(tab.histo.pos , aes(x=distance, y=after_stat(count))) +
        geom_histogram( binwidth=1)  +
        theme_bw() +
        xlim(-50,50) +
        ylim(0, ymax)
    ggpubr::ggexport(g, filename = file.path(out_dir, paste0(basename_out_file, "_distHist.pdf")), width = 5, height =5 )
}


## g <- ggplot(tab.histo.neg , aes(x=distance, y=after_stat(count))) +
##     geom_histogram( binwidth=1)  +
##     theme_bw() +
##     xlim(-50,50)
## ggpubr::ggexport(g, filename = file.path(out_dir, "distance_histogram_neg.pdf"), width = 5, height =5 )


## tab.all.counts <- tab.all %>%
##     group_by(id, distance, n.TSS) %>%
##     summarize(n=n()) %>%
##     reshape2::dcast( distance + n.TSS  ~ id , value.var="n") %>%
##     reshape2::melt(id.vars=c("distance", "pos", "n.TSS" ), value.name="neg", variable.name="background") %>%
##     mutate(ratio=pos/neg)

## g <- ggplot(tab.all.counts %>% dplyr::filter(abs(distance) <100 ), aes(x=distance, y=ratio)) +
##     geom_histogram(binwidth=1, stat="identity") +
##     facet_wrap(~ background , ncol=4,  scales="free_y")
## ggpubr::ggexport(g, filename = file.path(out_dir, "distance_distribution_ratio.pdf"), width = 15, height =8 )



## tab.all.list <- tab.all %>% split(f=.$id)
## n.TSS.c <- sapply(tab.all.list, function(x) x$n.TSS[1])

## tab.cumdist.onesided <-  tab.all.list %>%
##     mapply(function(x, y, z) calc.frac(x$distance,y) %>% mutate(id=z), .,  n.TSS.c, names(.),  SIMPLIFY=FALSE) %>%
##     do.call(rbind.data.frame, .) %>%
##     remove_rownames %>%
##     mutate(dist = dist %>% as.character %>% as.numeric)



## g <- ggplot(tab.cumdist.onesided %>% dplyr::filter(abs(dist) <30), aes(x=dist, y = frac)) +
##     geom_line(aes(colour=id))  +
##     theme_bw()
## out_dotplot <- file.path(out_dir,"all_dist_cumDensity_one_sided.png")
## ggsave( filename = out_dotplot, g, width = 5, height = 5)
