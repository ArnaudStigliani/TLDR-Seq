library(tidyverse)
library(ggplot2)
library(Biostrings)
library(stringr)
library(data.table)
library(parallel)
library(ggseqlogo)

in_dir <- "./results/get_motifs_from_read_end/"
out_dir <- "./results/get_motifs_TTSs/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##### Same thing but with unique sequences


nanopore_fas.name    <- file.path(in_dir ,"nanopore_tts.fas")
quantSeq_nopap_fas.name <- file.path(in_dir , "quantSeq_nopap_tts.fas")
quantSeq_xpap_fas.name <- file.path(in_dir ,  "quantSeq_xpap_tts.fas")
tldr_lib1_fas.name <- file.path(in_dir ,   "tldr_lib1_tts.fas")
tldr_lib2_fas.name <- file.path(in_dir ,   "tldr_lib2_tts.fas")
tldrpA_lib1_fas.name <- file.path(in_dir , "tldrpA_lib1_tts.fas")
tldrpA_lib2_fas.name <- file.path(in_dir , "tldrpA_lib2_tts.fas" )                   

#########
nanopore_fas <- readDNAStringSet(nanopore_fas.name)
quantSeq_nopap_fas <- readDNAStringSet(quantSeq_nopap_fas.name)
quantSeq_xpap_fas <- readDNAStringSet(quantSeq_xpap_fas.name)
tldr_lib1_fas <- readDNAStringSet(tldr_lib1_fas.name)
tldr_lib2_fas <- readDNAStringSet(tldr_lib2_fas.name)
tldrpA_lib1_fas <- readDNAStringSet(tldrpA_lib1_fas.name)
tldrpA_lib2_fas <- readDNAStringSet(tldrpA_lib2_fas.name)



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

tldr_lib1_fas.keep  <- tldr_lib1_fas %>%
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

tldr_lib2_fas.keep  <- tldr_lib2_fas %>%
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

tldrpA_lib1_fas.keep  <- tldrpA_lib1_fas %>%
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

tldrpA_lib2_fas.keep  <- tldrpA_lib2_fas %>%
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


nanopore_fas.keep  <- nanopore_fas %>%
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

quantSeq_nopap_fas.keep  <- quantSeq_nopap_fas %>%
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

quantSeq_xpap_fas.keep  <- quantSeq_xpap_fas %>%
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



tldr_lib1_fas.nt.freq <- tldr_lib1_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% tldr_lib1_fas.keep)] %>% 
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
tldr_lib2_fas.nt.freq <- tldr_lib2_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% tldr_lib2_fas.keep)] %>% 
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
tldrpA_lib1_fas.nt.freq <- tldrpA_lib1_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% tldrpA_lib1_fas.keep)] %>% 
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
tldrpA_lib2_fas.nt.freq <- tldrpA_lib2_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% tldrpA_lib2_fas.keep)] %>% 
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
nanopore_fas.nt.freq <- nanopore_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% nanopore_fas.keep)] %>% 
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
quantSeq_nopap_fas.nt.freq <- quantSeq_nopap_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% quantSeq_nopap_fas.keep)] %>% 
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
quantSeq_xpap_fas.nt.freq <- quantSeq_xpap_fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% quantSeq_xpap_fas.keep)] %>% 
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
 




all.freq <- bind_rows(TLDR.lib1=tldr_lib1_fas.nt.freq, TLDR.lib2=tldr_lib2_fas.nt.freq,
                      TLDR_polyA.lib1=tldrpA_lib1_fas.nt.freq, TLDR_polyA.lib2=tldrpA_lib2_fas.nt.freq,
                      Nanopore = nanopore_fas.nt.freq, quantSeq.nopap = quantSeq_nopap_fas.nt.freq,
                      quantSeq.xpap = quantSeq_xpap_fas.nt.freq,
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
    mutate(Method=factor(Method, levels = c( "TLDR.lib1", "TLDR.lib2", "TLDR_polyA.lib1", "TLDR_polyA.lib2",
                                            "quantSeq.nopap",  "quantSeq.xpap",  "Nanopore")))

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



#### write logo
all.logo.towrite <- all.freq %>%
    mutate(freq = (freq.norm))  %>%
    dplyr::select(Method, position, nt, freq) %>%
    group_by(position, Method) %>%
    mutate(tot = sum(freq), freq = freq/tot) %>%
    ungroup %>% 
    reshape2::dcast(Method  + position ~ nt, value.var = "freq") %>% 
    mutate(group = ifelse(position > -25 &  position < -12, paste0(Method, ".preTTS"), paste0(Method,".other"))) %>%
    mutate(group = ifelse(grepl("other", group) & position > -5 &  position < 5, paste0(Method, ".TTS"), group)) %>%
    dplyr::filter(!grepl("other", group)) %>%
    split(f=.$group) %>%
    lapply(dplyr::select, A,C,G,T) %>%
    lapply(remove_rownames) %>% 
    lapply(function(x) t(as.matrix(x) *  matrix(rep(2 + rowSums(( as.matrix(x)) *log (as.matrix(x), 2)),4),ncol=4)))





g <- ggseqlogo(all.logo.towrite[grepl("pre", names(all.logo.towrite))],
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=7, scales='free_x', dir="v") +
    scale_x_continuous(breaks=seq(1,11,3), labels=seq(-24,-13, 3)) 
ggpubr::ggexport(g, filename = file.path(out_dir, "ggseqlogo_preTTS.pdf"), width = 12, height =3 )


g <- ggseqlogo(all.logo.towrite[!grepl("pre", names(all.logo.towrite))],
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=7, scales='free_x', dir="v") +
    scale_x_continuous(breaks=seq(1,9,2), labels=seq(-4,4, 2)) 
ggpubr::ggexport(g, filename = file.path(out_dir, "ggseqlogo_TTS.pdf"), width = 12, height =3 )
