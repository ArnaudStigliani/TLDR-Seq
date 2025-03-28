rm(list=ls())
library(tidyverse)
library(parallel)
library(ggplot2)



lib1.dir <- "../mapping_preprocessing/results/get_polyA_tailfindr/lib1"
lib2.dir <- "../mapping_preprocessing/results/get_polyA_tailfindr/lib2"
names.lib <- c("lib1", "lib2")

out_dir.shared <- out_dir <- "./results/characterize_polyA"

i <- 0
for (lib.dir in c(lib1.dir, lib2.dir))
{
    i <- i+1
    in_dir <- lib.dir
    out_dir <- file.path(out_dir.shared, names.lib[i])
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    list_barcodes <- list.files(in_dir) %>% .[grepl("^bc[0-9]+", .)]
    print(list_barcodes)
    print("reading files")
    bc_polyA.list <- list()
    for (bc in list_barcodes)
    {
        bc_polyA.name <- file.path(in_dir, bc, "cDNA_tails.csv")
        if(! file.exists(bc_polyA.name))
            next
        bc_polyA.list[[bc]] <- read.table(bc_polyA.name, sep=",", header=TRUE)  %>% dplyr::select(-file_path)
    }
    bc_polyA.list <- bc_polyA.list %>% .[sapply(., nrow)> 1000]

    print("files read")
    pa.df <- bc_polyA.list %>%
        do.call(rbind.data.frame, args=.) %>%
        remove_rownames() %>%
        dplyr::filter(tail_length < 15 & ! is.na(tail_length)) %>%
        drop_na
    
    
    bc_poly.df <- sapply(bc_polyA.list, function(a)sum( !(a$tail_length == 0 | is.na(a$tail_length))))%>%
        data.frame(nreads.polyA = .) %>%
        rownames_to_column("barcode") %>%
        write.table(file.path(out_dir, "barcode_nreads_pA.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    
    bc_poly.df <- sapply(bc_polyA.list, function(a)sum( !(a$tail_length < 15 | is.na(a$tail_length))))%>%
        data.frame(nreads.polyA = .) %>%
        rownames_to_column("barcode") %>%
        write.table(file.path(out_dir, "barcode_nreads_pA_15.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    bc_poly.df <- sapply(bc_polyA.list,
                         function(a)sum(a$tail_length < 15 & a$tail_length > 0 & (!is.na(a$tail_length)))) %>%
        data.frame(nreads.polyA = .) %>%
        rownames_to_column("barcode") %>%
        write.table(file.path(out_dir, "barcode_nreads_pA_0_15.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    

    bc_poly.prop <- sapply(bc_polyA.list, function(a)sum( !(a$tail_length == 0 | is.na(a$tail_length)))/ nrow(a)) %>% round(2) %>% as.character
    bc_polyA.prop <- sapply(bc_polyA.list, function(a)sum( a$read_type == 'contains_a_polyA_tail')/ sum( !(a$tail_length == 0 | is.na(a$tail_length)))) %>% round(2) %>% as.character
    bc_polyT.prop <- sapply(bc_polyA.list, function(a)sum( a$read_type == 'contains_a_polyT_tail')/ sum( !(a$tail_length == 0 | is.na(a$tail_length)))) %>% round(2) %>% as.character


#### to generate bam files (with get_bam_polyA.sh)
    ## 15 nt threshold
    reads_with_polyA <- lapply(bc_polyA.list, function(a)(a  %>%
                                                          dplyr::filter(!(tail_length == 0 |
                                                                          is.na(a$tail_length))) %>%
                                                          dplyr::filter(tail_length > 15))) %>%
        do.call(rbind.data.frame, args=.) 
    write.table(reads_with_polyA %>% dplyr::select(read_id), file.path(out_dir, "reads_with_polyA.list"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    write.table(reads_with_polyA  %>% dplyr::select(read_id, read_type, tail_length), file.path(out_dir, "reads_with_polyA.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
    ## 0 nt threshold    
    reads_with_polyA_all <- lapply(bc_polyA.list, function(a)(a  %>%
                                                              dplyr::filter(!(tail_length == 0 |
                                                                              is.na(a$tail_length))) %>%
                                                              dplyr::filter(tail_length > 0))) %>%
        do.call(rbind.data.frame, args=.) 
    write.table(reads_with_polyA_all%>% dplyr::select(read_id), file.path(out_dir, "reads_with_polyA_all.list"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    write.table(reads_with_polyA_all  %>% dplyr::select(read_id, read_type, tail_length), file.path(out_dir, "reads_with_polyA_all.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
    ## 0-15 nt 
    reads_with_polyA_0_15 <- lapply(bc_polyA.list, function(a)(a  %>%
                                                               dplyr::filter(!(tail_length == 0 |
                                                                               is.na(a$tail_length))) %>%
                                                               dplyr::filter(tail_length <15 ))) %>%
        do.call(rbind.data.frame, args=.) 
    write.table(reads_with_polyA_0_15 %>% dplyr::select(read_id), file.path(out_dir, "reads_with_polyA_0-15.list"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    write.table(reads_with_polyA_0_15  %>% dplyr::select(read_id, read_type, tail_length), file.path(out_dir, "reads_with_polyA_0-15.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

########

plot_polyA_density <- function(polya.df, prop_tot, prop_A, prop_T, bc)
{
    g <- ggplot(polya.df, aes(x=tail_length)) +
        geom_density() +
        ggtitle(paste(bc, ", fraction reads with tail: ", prop_tot, "% polyA: ", prop_A, "% polyT: ", prop_T )) +
        facet_wrap(~read_type)
    return(g)
}

g.list <- mapply(plot_polyA_density, bc_polyA.list %>% lapply(drop_na),
                 bc_poly.prop, bc_polyA.prop, bc_polyT.prop, names(bc_polyA.list), ## arguments
                 SIMPLIFY=FALSE)


ggpubr::ggexport(g.list, filename = file.path(out_dir, "polyA_density.pdf"),
                 width = 12,  height =7 )
