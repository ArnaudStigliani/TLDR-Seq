library(tidyverse)
library(ggplot2)
library(Biostrings)
library(stringr)
library(plyranges)
library(data.table)
library(parallel)
library(ggseqlogo)

in_dir <- "./results/get_motifs_from_read_starts/"
in_dir <- "./results/plot_motifs"

n.seq.tldr <- system(paste("wc -l" ,file.path(in_dir , "tldr_tss.fas")), intern=TRUE ) %>%
    str_replace(" .*","") %>%
    as.numeric %>%
    "/" (2)

n.seq.tldrlib2 <- system(paste("wc -l" ,file.path(in_dir , "tldrlib2_tss.fas")), intern=TRUE ) %>%
    str_replace(" .*","") %>%
    as.numeric %>%
    "/" (2)

n.seq.nanopore <- system(paste("wc -l" ,file.path(in_dir , "nanopore_tss.fas")), intern=TRUE ) %>%
    str_replace(" .*","") %>%
    as.numeric %>%
    "/" (2)

n.seq.slic <- system(paste("wc -l" ,file.path(in_dir , "slic_cage_tss.fas")), intern=TRUE ) %>%
    str_replace(" .*","") %>%
    as.numeric %>%
    "/" (2)

n.seq.df <- list(Nanopore = n.seq.nanopore,
                 TLDR.lib1 = n.seq.tldr,
                 TLDR.lib2 = n.seq.tldrlib2,
                 SLIC_CAGE = n.seq.slic) %>%
    lapply(as.data.frame) %>%
    lapply(setNames, "n.seq") %>% 
    mapply(mutate, ., Method = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames






##### look at the average nucleotide frequency

tldr_fas.name <- file.path(in_dir , "tldr_tss.fas")
tldrlib2_fas.name <- file.path(in_dir , "tldrlib2_tss.fas")
nanopore_fas.name <- file.path(in_dir , "nanopore_tss.fas")
slic_cage_fas.name <- file.path(in_dir , "slic_cage_tss.fas")

tldr.fas <- readDNAStringSet(tldr_fas.name)
tldrlib2.fas <- readDNAStringSet(tldrlib2_fas.name)
nanopore.fas <- readDNAStringSet(nanopore_fas.name)
slic_cage.fas <- readDNAStringSet(slic_cage_fas.name)



##### Same thing but with unique sequences
nTSS.tldr.fas <- round((tldr.fas %>% length())/(5*10**6))
nTSS.tldrlib2.fas <- round((tldrlib2.fas %>% length())/(5*10**6))
nTSS.nanopore.fas <- round((nanopore.fas %>% length())/(5*10**6))
nTSS.slic_cage.fas <- round((slic_cage.fas %>% length())/(5*10**6))

determine.windows <- function(df)
{
    w.size <- 40
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


    
tldr.fas.keep  <- tldr.fas %>%
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



tldrlib2.fas.keep  <- tldrlib2.fas %>%
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

nanopore.fas.keep  <- nanopore.fas %>%
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


slic_cage.fas.keep  <- slic_cage.fas %>%
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





 

tldr.nt.freq <- tldr.fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% tldr.fas.keep)] %>% 
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

nanopore.nt.freq <- nanopore.fas %>%
    .[!duplicated(names(.))] %>% 
    .[which( names(.) %in% nanopore.fas.keep)] %>% 
    as.character() %>%
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


slic_cage.nt.freq <- slic_cage.fas %>%
    .[!duplicated(names(.))] %>% 
    .[which(names(.) %in% slic_cage.fas.keep )] %>% 
    as.character() %>%
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

tldrlib2.nt.freq <- tldrlib2.fas %>%
    .[!duplicated(names(.))] %>% 
    .[which(names(.) %in% tldrlib2.fas.keep )] %>% 
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



all.freq <- bind_rows(TLDR.lib1=tldr.nt.freq, TLDR.lib2=tldrlib2.nt.freq,
                      Nanopore = nanopore.nt.freq, SLIC_CAGE = slic_cage.nt.freq,  .id="Method") %>%
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
    mutate(Method = factor(Method, levels=c("TLDR.lib1",  "TLDR.lib2" , "SLIC_CAGE", "Nanopore" )))

write.table(all.freq, file.path(out_dir, "all_freq_unique_3support.tsv"), row.names=FALSE, sep="\t", quote=FALSE)        


g <- ggplot(all.freq, aes(x=position, y = freq, colour = nt)) +
    geom_density(stat="identity")  +
    facet_wrap( ~ Method, nrow=1) +
    theme_bw()
ggpubr::ggexport(g, filename = file.path(out_dir, "nucleotide_frequency_position_unique_3support.pdf"), width = 15, height =5 )

g <- ggplot(all.freq, aes(x=position, y = freq.norm, colour = nt)) +
    geom_density(stat="identity")  +
    facet_wrap( ~ Method, nrow=1) +
    theme_bw()
ggpubr::ggexport(g, filename = file.path(out_dir, "nucleotide_frequency_norm_position_unique_3support.pdf"), width = 15, height =5 )




### prepare to write sites


all.freq.towrite <- all.freq %>%
    mutate(freq = (freq.norm))  %>%
    dplyr::select(Method, position, nt, freq) %>%
    group_by(position, Method) %>%
    mutate(tot = sum(freq), freq = freq/tot) %>%
    ungroup %>% 
    reshape2::dcast(Method  + position ~ nt, value.var = "freq") %>% 
    mutate(group = ifelse(position > -35 &  position < -20, paste0(Method, ".TATA"), paste0(Method,".other"))) %>%
    mutate(group = ifelse(grepl("other", group) & position > -5 &  position < 5,
                          paste0(Method, ".TSS"), group)) %>%
    dplyr::filter(!grepl("other", group)) %>%
    split(f=.$group) %>%
    lapply(dplyr::select, A,C,G,T) %>%
    lapply(remove_rownames) %>%
    mapply(write.table,.,  file.path(out_pfm, paste0(names(.),".pfm")),
           quote=FALSE,
           col.names=FALSE,
           row.names=FALSE,
           sep="  ",
           SIMPLIFY=FALSE)


all.logo.towrite <- all.freq %>%
    mutate(freq = (freq.norm))  %>%
    dplyr::select(Method, position, nt, freq) %>%
    group_by(position, Method) %>%
    mutate(tot = sum(freq), freq = freq/tot) %>%
    ungroup %>% 
    reshape2::dcast(Method  + position ~ nt, value.var = "freq") %>% 
    mutate(group = ifelse(position > -35 &  position < -20, paste0(Method, ".TATA"), paste0(Method,".other"))) %>%
    mutate(group = ifelse(grepl("other", group) & position > -5 &  position < 5,
                          paste0(Method, ".TSS"), group)) %>%
    dplyr::filter(!grepl("other", group)) %>%
    split(f=.$group) %>%
    lapply(dplyr::select, A,C,G,T) %>%
    lapply(remove_rownames) %>% 
    lapply(function(x) t(as.matrix(x) *  matrix(rep(2 + rowSums(( as.matrix(x)) *log (as.matrix(x), 2)),4),ncol=4)))


p_list <- lapply(all.logo.towrite$SLIC_CAGE.TSS)


g <- ggseqlogo(all.logo.towrite[grepl("TATA", names(all.logo.towrite))],
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=4, scales='free_x', dir="v") +
    scale_x_continuous(breaks=seq(1,14,2), labels=seq(-34,-21, 2)) 
    ggpubr::ggexport(g, filename = file.path(out_dir, "ggseqlogo_TATA.pdf"), width = 10, height =3 )


g <- ggseqlogo(all.logo.towrite[grepl("TSS", names(all.logo.towrite))],
               method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=4, scales='free_x', dir="v") +
    scale_x_continuous(breaks=seq(1,9,2), labels=seq(-4,4, 2)) 
ggpubr::ggexport(g, filename = file.path(out_dir, "ggseqlogo_TSS.pdf"), width = 10, height =3 )


