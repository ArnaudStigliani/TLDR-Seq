library(ggplot2)
library(stringr)
library(magrittr)
library(tidyverse)
library(reshape2)

arg  <-  commandArgs(trailingOnly=TRUE)

## arg[1] <- "./results/{name sequening run}/format_adapter_information/"


if (is.null(arg[1]))
{
    print("provide arguments in such a manner\nRscript ./results/{name sequening run}/format_adapter_information/")
    quit()
}

in_dir <- arg[1]
date <- in_dir %>% basename
out_dir <- in_dir %>% dirname %>% dirname %>% file.path("analyze_adapter_information", date)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### read all files
adapter_details.df <- read.table(file.path(in_dir, "adapter_details.tsv")) %>%
    dplyr::rename(barcode = V1, strand = V2, adapter_5prime = V3, adapter_3prime = V4)

trimmed.df <- read.table(file.path(in_dir, "trimmed.tsv"))  %>%
    dplyr::rename(barcode = V1, strand = V2, reads_trimmed = V3, percent_reads_trimmed = V4)

trimmed_one_sided.df <- read.table(file.path(in_dir, "trimmed_one_sided.tsv")) %>%
        dplyr::rename(barcode = V1, strand = V2, reads_trimmed = V3, percent_reads_trimmed = V4)

tsv_files.dir <- file.path(in_dir, "tsv_files")
tsv_files.list <- list.files(tsv_files.dir) %>% file.path(tsv_files.dir, .)

tsv_files.df <- mapply(read.table, tsv_files.list, MoreArgs=list(sep = "\t", header=TRUE), SIMPLIFY = FALSE) %>%
    setNames(names(.) %>% basename %>% str_replace(".tsv", "")) %>%
    mapply(mutate, ., adapter = names(.), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    separate(adapter, into=c("side", "barcode", "strand")) %>%
    remove_rownames()


### get total number of reads in that pass basecalling
log.file.path <- in_dir %>% dirname %>% dirname %>% file.path("pipeline_map", date, "trimming", "bc1_s1.log")
command.line <- paste("head -8 ", log.file.path, " | tail -1 | awk '{print $4}' | sed 's/,//g'",  sep="")
nread <- system(command.line, intern = TRUE) %>% as.numeric


##### look at adapter proportion
adapter_details.to_plot <- adapter_details.df %>%
    mutate(A2 = ifelse(strand=="s1", adapter_5prime, adapter_3prime),
           A1 = ifelse(strand != "s1", adapter_5prime, adapter_3prime)) %>%
    dplyr::select(-adapter_5prime, -adapter_3prime) %>%
    melt(variable.name="adapter", value.name = "reads_with_adapters") %>%
    arrange(adapter,barcode,  strand)   %>%
    head(grep("A1", .$adapter)[2]) %>%
    mutate(percent_reads_with_adapters = reads_with_adapters/nread * 100) %>%
    group_by(adapter) %>%
    mutate(total_percent = sum(percent_reads_with_adapters) %>% round(0) %>% as.character) %>%
    ungroup() %>%
    mutate(adapter=paste0(adapter, " (", total_percent, "%)"))

g <- ggplot(adapter_details.to_plot, aes(x=barcode, y=reads_with_adapters)) +
    geom_bar(stat='identity') +
    facet_grid(vars(strand),vars(adapter), scales="free_x", space="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "adapter_details.pdf"), width = 8, height = 6)

g <- ggplot(adapter_details.to_plot, aes(x=barcode, y=percent_reads_with_adapters)) +
    geom_bar(stat='identity') +
    facet_grid(vars(strand),vars(adapter), scales="free_x", space="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "adapter_details_percent.pdf"), width = 8, height = 6)


## look at detected adapters
trimmed_one_sided.to_plot <- trimmed_one_sided.df %>%
    dplyr::rename(percent_reads_with_1_ormore_adapter = percent_reads_trimmed) %>%
    dplyr::rename(reads_with_1_ormore_adapter = reads_trimmed)


g <- ggplot(trimmed_one_sided.to_plot, aes(x=barcode, y=percent_reads_with_1_ormore_adapter)) +
    geom_bar(stat='identity') +
    facet_wrap(~ strand, scales="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_one_ormore_adapters_percent.pdf"), width = 12, height = 6)

g <- ggplot(trimmed_one_sided.to_plot, aes(x=barcode, y=reads_with_1_ormore_adapter)) +
    geom_bar(stat='identity') +
    facet_wrap(~ strand, scales="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_one_ormore_adapters.pdf"), width = 12, height = 6)

### 2 adapters
trimmed.to_plot <- trimmed.df %>%
    dplyr::rename(percent_reads_with_2_adapter = percent_reads_trimmed) %>%
    dplyr::rename(reads_with_2_adapter = reads_trimmed)  %>%
    group_by(strand) %>%
    mutate(total_percent = sum(percent_reads_with_2_adapter) %>% round(0) %>% as.character) %>%
    ungroup() %>%
    mutate(strand=paste0(strand, " (", total_percent, "%)"))



g <- ggplot(trimmed.to_plot, aes(x=barcode, y=percent_reads_with_2_adapter)) +
    geom_bar(stat='identity') +
    facet_wrap(~ strand, scales="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_two_adapters_percent.pdf"), width = 12, height = 6)

g <- ggplot(trimmed.to_plot, aes(x=barcode, y=reads_with_2_adapter)) +
    geom_bar(stat='identity') +
    facet_wrap(~ strand, scales="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_two_adapters.pdf"), width = 12, height = 6)

### Exclusively one adapter
one_adapter_exclusive.df <- left_join(adapter_details.to_plot %>%
                                      mutate(adapter=adapter %>% str_replace(" .*", "")) %>%
                                      dplyr::select(-total_percent),
                                      trimmed.to_plot %>%
                                      mutate(strand=strand %>% str_replace(" .*", "")) %>%
                                      dplyr::select(-total_percent),
                                      by=c("barcode","strand")) %>%
    dplyr::select(-starts_with("percent")) %>%
    group_by(adapter, strand) %>%
    mutate(sum_reads_2_adapters = sum(reads_with_2_adapter)) %>%
    ungroup() %>%
    arrange(adapter, barcode, strand) %>%
    mutate(sum_reads_2_adapters = rep(.$sum_reads_2_adapters[.$adapter=="A2"][1:2], nrow(.)/2)) %>% 
    mutate(reads_with_only_1_adapter = ifelse(adapter=="A2",
                                              reads_with_adapters - reads_with_2_adapter,
                                              reads_with_adapters - sum_reads_2_adapters)) %>%
    as.data.frame  %>%
    dplyr::select(- reads_with_adapters,  - reads_with_2_adapter, - sum_reads_2_adapters) %>%
    mutate(percent_reads_with_only_1_adapter = reads_with_only_1_adapter/nread * 100) %>%
    group_by(adapter) %>%
    mutate(total_percent = sum(percent_reads_with_only_1_adapter) %>% round(0) %>% as.character) %>%
    ungroup() %>%
    mutate(adapter=paste0(adapter, " (", total_percent, "%)"))
    

g <- ggplot(one_adapter_exclusive.df, aes(x=barcode, y=reads_with_only_1_adapter)) +
    geom_bar(stat='identity') +
    facet_grid(vars(strand),vars(adapter), scales="free_x", space="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_only_one_adapter.pdf"), width = 12, height = 6)

g <- ggplot(one_adapter_exclusive.df, aes(x=barcode, y=percent_reads_with_only_1_adapter)) +
    geom_bar(stat='identity') +
    facet_grid(vars(strand),vars(adapter), scales="free_x", space="free_x") +
    ggtitle(paste0("Total read number: ", formatC(nread, format='e', digits=3)))
ggpubr::ggexport(g, filename = file.path(out_dir, "reads_with_only_one_adapter_percent.pdf"), width = 12, height = 6)

### plot adapter position

adapter_position.df <- tsv_files.df  %>%
    group_by(barcode) %>%
    mutate(n=sum(count)) %>%
    dplyr::filter(n > 10000) %>%
    ungroup %>%
    dplyr::filter(count > 200) %>%
    mutate(barcode = ifelse((strand == "s1" & side == "3") | (strand == "s2" & side == "5"),
                            paste0(barcode, " (A1)"), paste0(barcode, " (A2)")))


g <- ggplot(adapter_position.df, aes(x=length, y=count)) +
    geom_density(stat='identity') +
    scale_x_log10() +
    facet_wrap(~ barcode + strand + side, ncol=4, scale='free_y')
ggpubr::ggexport(g, filename = file.path(out_dir, "density_length_trimmed.pdf"), width = 12,
                 height = 1.5 * length(adapter_position.df$barcode %>% unique))


g <- ggplot(adapter_position.df, aes(x=length, y=count)) +
    geom_histogram(stat='identity') +
    xlim(c(40,110)) +
    ## scale_x_log10() +
    facet_wrap(~ barcode + strand + side, ncol=4, scale='free_y')
ggpubr::ggexport(g, filename = file.path(out_dir, "histogram_length_trimmed.pdf"), width = 12,
                 height = 1.5 * length(adapter_position.df$barcode %>% unique))

