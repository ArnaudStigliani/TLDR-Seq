library(tidyverse)
library(ggplot2)
library(Biostrings)
library(stringr)
library(data.table)
library(parallel)
library(ggseqlogo)

out_dir <- "./results/reviews/plot_motifs_TTS"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##### get sequences for the gencode TTS motifs

genome_fas <- "../data/GRCh38.primary_assembly.genome.fa"
gencode_bed.name    <- "../data/gencode_v39_prot_coding_transcripts_for_TLDR.bed")

gencode_tts_bed.name    <- file.path(out_dir, "gencode_tts.bed")
gencode_fas.name    <- file.path(out_dir ,"gencode_tts.fas")

awk_command <- paste("awk -v OFS='\\t' '$1~/chr/ && $6==\"-\" && $2-50 >0 {print $1,$2-50,$2+50,\".\",\".\",$6} $1~/chr/ && $6==\"+\" && $3-50 >0 {print $1,$3-50,$3+50,\".\",\".\",$6}'",  gencode_bed.name  , ">",  gencode_tts_bed.name )

system(awk_command)

bedtools.command <- paste("module load bedtools; bedtools getfasta -fi",  genome_fas, "-fo", gencode_fas.name,
      "-bed", gencode_tts_bed.name, " -s 2> /dev/null")
system(bedtools.command)


#########
gencode_fas <- readDNAStringSet(gencode_fas.name)



determine.windows <- function(df)
{
    w.size <- 20
    a <- lapply(0:(w.size-1), function(x, y) (y[["start"]] + x) %/% w.size , df ) %>%
        setNames(paste0("window.", 0:(w.size-1))) %>%
        do.call(cbind.data.frame, .) %>%
        cbind.data.frame(df, .)
    for(i in 0:(w.size-1))
    {
        oldname <- paste0("window.",i)
        a <- a %>%
            mutate(w := get(!!oldname)) %>%
            arrange(strand, chr, w, desc(score)) %>% 
            group_by(chr, w, strand) %>% 
            dplyr::slice(1)
    }
    a  <- a %>% ungroup
    return(a)
}


gencode_fas.keep  <- gencode_fas %>%
    names %>%
    table %>%
    sort(decreasing = TRUE) %>% 
    head(10000) %>% 
    as.data.frame %>%
    dplyr::select(score=Freq, seqnames=".") %>%
    mutate(seqnames = seqnames %>% str_replace("\\(-\\)","(.)")) %>% 
    separate(seqnames, into=c("chr","start","end","strand"), sep="[:\\-\\(]") %>%
    mutate(strand = ifelse(strand ==".)", "-", "+")) %>%
    dplyr::filter(strand == "+") %>%
    mutate( start = start %>% as.numeric, end  = end %>% as.numeric ) %>% 
    determine.windows  %>%
    dplyr::select(!starts_with("window")) %>%
    dplyr::select(-w) %>%
    mutate(name.seq=paste0(chr,":", start, "-", end, "(",strand , ")")) %>%
    arrange(desc(score)) %>%
    head(750) %>%
    .$name.seq


gencode_fas.nt.freq <- gencode_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% gencode_fas.keep)] %>% 
    as.character()  %>%
    str_split(., "") %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    as.list()  %>%
    mclapply(paste, collapse = "", mc.cores=100) %>%
    mclapply(function(x) sapply(list(A="A",C="C",G="G",T="T"), function(y) str_count( x ,y )),mc.cores=100) %>%
    do.call(cbind.data.frame, .) %>%
    t %>%
    as.data.frame %>%
    remove_rownames %>%
    mutate(position = 1:100) %>%
    reshape2::melt(id.vars=c("position"))
 




all.freq <- bind_rows(Gencode = gencode_fas.nt.freq, 
                      .id="Method") %>%
    dplyr::rename(nt=variable) %>%
    group_by(Method, position) %>%
    mutate(n.seq = sum(value)) %>% 
    group_by(nt, Method) %>%
    mutate(nt.abundance = sum(value)) %>%
    ungroup %>%
    group_by(Method) %>%
    mutate(nt.tot=sum(value)) %>%
    ungroup %>% 
    mutate(freq=value/n.seq, freq.norm = value/(n.seq * (nt.abundance/nt.tot))) %>%
    mutate(position = position - 50) %>%
    mutate(Method=factor(Method, levels = c( "Gencode")))

write.table(all.freq, file.path(out_dir, "all_freq_unique_3support.tsv"), row.names=FALSE, sep="\t", quote=FALSE)


g <- ggplot(all.freq, aes(x=position, y = freq, colour = nt)) +
    geom_density(stat="identity")  +
    facet_wrap( ~ Method, nrow = 1) +
    theme_bw()
ggpubr::ggexport(g, filename = file.path(out_dir, "nucleotide_frequency_position_unique_3support.pdf"), width = 15, height = 5 )

g <- ggplot(all.freq, aes(x=position, y = freq.norm, colour = nt)) +
    geom_density(stat="identity")  +
    facet_wrap( ~ Method, nrow = 1) +
    theme_bw()
ggpubr::ggexport(g, filename = file.path(out_dir, "nucleotide_frequency_norm_position_unique_3support.pdf"), width = 15, height = 5 )
