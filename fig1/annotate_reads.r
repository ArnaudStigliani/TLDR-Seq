library(tidyverse)
library(parallel)
library(stringr)
library(data.table)
library(Biostrings)
library(ShortRead)
library(bambu)

shared.data.dir <- "../../../shared_data/human_genome"

out_dir <- "./annotate_reads"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

in_dir.lib1 <- "./get_polyA_bam_files/lib1"
in_dir.R10 <- "./get_polyA_bam_files/R10"
in_dir.lib2 <- "./get_polyA_bam_files/lib2"
in.dir.ont <- "../map_ONTnanopore/results_curated"
in.dir.RNASeq <- "../data/RNASeq/salmon_mapping" ## path to RNASeq fastq files
out_RNASeq <- file.path(out_dir,"./RNASeq/" )


gtf.name <- "../data/gencode.v39.primary_assembly.annotation.sorted.gtf" # path to gtf file
fas.name <- "../data/GRCh38.primary_assembly.genome.fa" # path to fasta genome  file



## create bambu annotation file
annotations <- prepareAnnotations(gtf.name)
saveRDS(annotations, file.path(out_dir, "annotations.rds")) 
annotations <- readRDS(file.path(out_dir, "annotations.rds"))


#### for ont files, create a merged bam file
bam_ont_merged <- file.path(out_dir, "ont_merged.bam")

command <- paste(" module load samtools; samtools merge /dev/stdout",
                 file.path(in.dir.ont, c("nanopore_wt1_curated.bam")),
                 file.path(in.dir.ont, c("nanopore_wt2_curated.bam")),                 
                 "| samtools sort -@ 20 > ", bam_ont_merged ,
                 "; samtools index ", bam_ont_merged) 
system(command)


### for R10 data, separate pA reads from non pA  reads (diffeent from lib1 and lib2 because of some symbolic link silly issue)

bam_R10_merged <- file.path(in_dir.R10, "trimmed_primmary_renamed.bam")
bam_R10_pA <- file.path(out_dir, "R10_polyA_plus.bam")
bam_R10_npA <- file.path(out_dir, "R10_polyA_minus.bam")

R10_pA_list.name <- file.path( "../mapping_preprocessing/get_pA_reads_R10", "reads_pA.txt")
R10_npA_list.name <- file.path( "../mapping_preprocessing/get_pA_reads_R10", "no_polyA.list")

command1 <- paste("module load samtools; samtools view -b -N",  R10_pA_list.name, bam_R10_merged, ">", bam_R10_pA , "&\n" )
command2 <- paste("module load samtools; samtools view -b -N",  R10_npA_list.name, bam_R10_merged, ">", bam_R10_npA)

system(paste(command1 ,command2, "& wait" ))


## Quantify RNASeq files using salmon

salmon.index <-  "Path_to_salmon_index/salmon_index"

RNASeq.pair1.name <- file.path(in.dir.RNASeq, "HeLA_1.fastq")
RNASeq.pair2.name <- file.path(in.dir.RNASeq, "HeLA_2.fastq")

command <- paste("module unload R-packages; module unload R; module unload gcc; module load gcc/13.2.0; salmon quant -i",
                 salmon.index," -l IU -1",
                 RNASeq.pair1.name, "-2",  RNASeq.pair2.name ,
                 "--allowDovetail -o", out_RNASeq,
                 "--validateMappings -p 20")
system(command)


### create symbolic links because bambu only looks at the basename, making it crash if 2 files gave the same basename
lib1.links <- list.files(in_dir.lib1, pattern="bam$", full.names=TRUE) %>%
    {file.path(out_dir, paste0("lib1_",
                               basename(.) %>%
                               str_replace("trimmed_primmary_", "")))}
command.lib1 <- paste("ln -s", list.files(in_dir.lib1, pattern="bam$", full.names=TRUE), lib1.links)

lib2.links <- list.files(in_dir.lib2, pattern="bam$", full.names=TRUE) %>%
        {file.path(out_dir, paste0("lib2_",
                                   basename(.) %>%
                                   str_replace("trimmed_primmary_", "")))}
command.lib2 <- paste("ln -s", list.files(in_dir.lib2, pattern="bam$", full.names=TRUE), lib2.links)

commands <- c(command.lib1, command.lib2) %>%
    str_replace_all("../../", "/maps/projects/scarball/people/fzh976/binf-isilon/polyA_CAGE/") %>% paste(collapse="; ")
system(commands)

####

all_bam_files <- c(list.files(out_dir, pattern="^lib.*bam$", full.names=TRUE),
                   c(bam_R10_pA, bam_R10_npA), 
                   bam_ont_merged) %>% .[!grepl("renamed",.)] %>%
    setNames(basename(.) %>% str_replace( "trimmed_primary_polyA_", "") %>% str_replace(".bam", "") )



LongReads.quant <- bambu(reads = all_bam_files,  annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 20)

# LongReads.gene.tpm <- transcriptToGeneExpression(LongReads.quant)


df.quant <- LongReads.quant %>%
    assays() %>%
    .$CPM %>%
    as.data.frame %>%
    rownames_to_column("tr.name")
write.table(df.quant, file=file.path(out_dir, "quant.tsv"), row.names=FALSE, sep="\t", quote=FALSE) 
df.quant <- read.table(file=file.path(out_dir, "quant.tsv"), header=TRUE, sep="\t")

######## read nanocounts files

RNASeq.name <- file.path(out_RNASeq, "quant.sf")

tr.types.c <- c("lncRNA" ,"nonsense_mediated_decay" ,"non_stop_decay" ,"retained_intron" ,"protein_coding" ,"protein_coding_LoF" ,"protein_coding_CDS_not_defined" ,"processed_transcript" ,"non_coding" ,"ambiguous_orf" ,"sense_intronic" ,"sense_overlapping" ,"antisense/antisense_RNA" ,"transcribed_processed_pseudogene" ,"transcribed_unprocessed_pseudogene" ,"transcribed_unitary_pseudogene" ,"translated_processed_pseudogene" ,"translated_unprocessed_pseudogene" ,"lincRNA" ,"macro_lncRNA" ,"3prime_overlapping_ncRNA" ,"disrupted_domain" ,"bidirectional_promoter_lncRNA")


RNASeq.name <- file.path(out_RNASeq, "quant.sf") 
RNASeq.df <-  read.table(RNASeq.name, header=TRUE) %>%
    dplyr::select(Name, RNASeq = TPM)

Nanocount.placeholder.name <- "../../results/sequencing_july_2023/shared/annotate_reads_v2/tx_lib1_polyA_minus.tsv" ## Use old Nanocount output to get transcript types.  One can also use the gtf file to obtain the same information

tr.type <- read.table(Nanocount.placeholder.name, header=TRUE) %>%
    dplyr::select(transcript_name) %>%
    mutate(tr.type = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,(ncol(.)-1 )]) %>%
    mutate(transcript_name = transcript_name %>% str_replace("\\|.*","")) 
    



names.LongReads <- all_bam_files %>%
    basename() %>%
    str_replace( "_trimmed_primary_polyA", "_pA") %>% 
    str_replace(".bam", "") %>%
    str_replace("annotate_reads_v2_ont_merged", "ONT") %>%
    str_replace("lib1", "Lib_Amp") %>%
    str_replace("lib2", "Lib_NoAmp") %>%
    str_replace("all", "0") %>%
    str_replace("0_15", "1_14") 
    
    
    

y.levels <- c("Lib_Amp_pA_plus_1_14",
              "Lib_NoAmp_pA_plus_1_14",
              "Lib_Amp_pA_minus_0",
              "Lib_NoAmp_pA_minus_0",
              "Lib_Amp_pA_plus",
              "Lib_NoAmp_pA_plus",
              "ONT",
              "RNA-Seq")               

LongReads.gene.tpm.df <- LongReads.quant %>%
    assays() %>%
    .$CPM %>%
    as.data.frame %>% 
    setNames(names.LongReads) %>%
    rownames_to_column("transcript_name") %>%
    left_join(tr.type, multiple="all") %>%
    dplyr::filter(tr.type %in% tr.types.c)  %>%
    dplyr::rename(Name = transcript_name) %>%
    left_join(RNASeq.df) %>%
    reshape2::melt(variable.name="Library") %>%
    drop_na() %>% 
    group_by(tr.type, Library) %>%
    summarize(tpm.tr.type = sum(value)) %>%
    dplyr::filter(tpm.tr.type > 20) %>%
    mutate(Lib = Library %>% as.character  %>% ifelse(. =="RNASeq", "RNA-Seq", .)) %>%
    dplyr::filter( Lib %in% y.levels) %>%
    mutate(Lib = factor(Lib, y.levels %>% rev))


#### to get LNC RNA (snapshots)
tr.type2 <- read.table(Nanocount.placeholder.name, header=TRUE) %>%
    dplyr::select(transcript_name) %>%
    mutate(tr.type = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,(ncol(.)-1 )]) %>%
    mutate(gene.name = str_split(transcript_name, "\\|", simplify=TRUE) %>% .[,(ncol(.)-3 )])  %>% 
    mutate(transcript_name = transcript_name %>% str_replace("\\|.*","")) 

df.annot <- df.quant %>%
    dplyr::select(transcript_name = tr.name, lib1_trimmed_primary_polyA_minus, lib1_trimmed_primary_polyA_plus) %>% 
    left_join(tr.type2, multiple="all") %>%
    dplyr::rename(Name = transcript_name) %>%
    dplyr::filter(tr.type == "lncRNA") %>%
    dplyr::filter(lib1_trimmed_primary_polyA_minus > 10 & lib1_trimmed_primary_polyA_minus > 10 )
#######

 
   
g <- ggplot(LongReads.gene.tpm.df, aes(x=tr.type, y=tpm.tr.type, fill=Lib))+
    geom_bar(stat="identity", position=position_dodge(preserve = "single"))  +
    ## scale_fill_manual(values=RColorBrewer::brewer.pal(n=12, "Set2"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_bw() +
    coord_flip()
out_dotplot <- file.path(out_dir, "annotation_allSeq_TPM.pdf")
ggsave( filename = out_dotplot, g, width = 7, height = 10) 



################ make paper figure


names.LongReads.figure <- all_bam_files %>%
    basename() %>%
    str_replace( "_trimmed_primary_polyA", "_pA") %>% 
    str_replace(".bam", "") %>%
    str_replace("ont_merged", "ONT") %>%
    str_replace("lib1", "Lib_Amp") %>%
    str_replace("lib2", "Lib_NoAmp")  %>%
    str_replace("Lib_Amp_pA_plus", "TLDR lib amplified pA+") %>%
    str_replace("Lib_NoAmp_pA_plus", "TLDR lib non-amplified pA+") %>%
    str_replace("Lib_Amp_pA_minus", "TLDR lib amplified pA-") %>%
    str_replace("Lib_NoAmp_pA_minus", "TLDR lib non-amplified pA-") %>%
    str_replace("ONT", "ONT-Nanopore") 
    

y.levels.paper <- c("TLDR lib amplified pA+",
              "TLDR lib amplified pA-",
              "TLDR lib non-amplified pA+",
              "TLDR lib non-amplified pA-",
              "ONT-Nanopore",
              "RNA-Seq")

    


LongReads.gene.tpm.paper.df <- LongReads.quant %>%
    assays() %>%
    .$CPM %>%
    as.data.frame %>%
    setNames(names.LongReads.figure) %>%
    rownames_to_column("transcript_name") %>%
    left_join(tr.type, multiple="all") %>%
    dplyr::filter(tr.type %in% tr.types.c)  %>%
    dplyr::rename(Name = transcript_name) %>%
    left_join(RNASeq.df) %>%
    reshape2::melt(variable.name="Library") %>%
    drop_na() %>% 
    group_by(tr.type, Library) %>%
    summarize(tpm.tr.type = sum(value)) %>%
    dplyr::filter(tpm.tr.type > 20) %>%
    mutate(Lib = Library %>% as.character  %>% ifelse(. =="RNASeq", "RNA-Seq", .)) %>%
    dplyr::filter( Lib %in% y.levels.paper) %>%
    mutate(Lib = factor(Lib, y.levels.paper))


g <- ggplot(LongReads.gene.tpm.paper.df, aes(x=tr.type, y=tpm.tr.type, fill=Lib))+
    geom_bar(stat="identity", position=position_dodge(preserve = "single"))  +
    ## scale_fill_manual(values=RColorBrewer::brewer.pal(n=12, "Set2"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_bw() +
    coord_flip()
out_dotplot <- file.path(out_dir, "annotation_allSeq_TPM_paper.pdf")
ggsave( filename = out_dotplot, g, width = 7, height = 5) 


########################################
########### make R10 figure



names.LongReads.figure <- all_bam_files %>%
    basename() %>%
    str_replace( "_trimmed_primary_polyA", "_pA") %>% 
    str_replace(".bam", "") %>%
    str_replace("ont_merged", "ONT") %>%
    str_replace("lib1", "Lib_Amp") %>%
    str_replace("lib2", "Lib_NoAmp")  %>%
    str_replace("Lib_Amp_pA_plus", "TLDR lib amplified pA+") %>%
    str_replace("Lib_NoAmp_pA_plus", "TLDR lib non-amplified pA+") %>%
    str_replace("Lib_Amp_pA_minus", "TLDR lib amplified pA-") %>%
    str_replace("Lib_NoAmp_pA_minus", "TLDR lib non-amplified pA-") %>%
    str_replace("R10_polyA_plus", "R10 TLDR lib amplified pA+")  %>%
    str_replace("R10_polyA_minus", "R10 TLDR lib amplified pA-")  %>% 
    str_replace("ONT", "ONT-Nanopore")  

    

y.levels.paper <- c(
    "TLDR lib amplified pA+",
    "TLDR lib amplified pA-",
    "R10 TLDR lib amplified pA+",
    "R10 TLDR lib amplified pA-",
    "ONT-Nanopore",
    "RNA-Seq")

    


LongReads.gene.tpm.paper.df <- LongReads.quant %>%
    assays() %>%
    .$CPM %>%
    as.data.frame %>%
    setNames(names.LongReads.figure) %>%
    rownames_to_column("transcript_name") %>%
    left_join(tr.type, multiple="all") %>%
    dplyr::filter(tr.type %in% tr.types.c)  %>%
    dplyr::rename(Name = transcript_name) %>%
    left_join(RNASeq.df) %>%
    reshape2::melt(variable.name="Library") %>%
    drop_na() %>% 
    group_by(tr.type, Library) %>%
    summarize(tpm.tr.type = sum(value)) %>%
    dplyr::filter(tpm.tr.type > 20) %>%
    mutate(Lib = Library %>% as.character  %>% ifelse(. =="RNASeq", "RNA-Seq", .)) %>%
    dplyr::filter( Lib %in% y.levels.paper) %>%
    mutate(Lib = factor(Lib, y.levels.paper))


g <- ggplot(LongReads.gene.tpm.paper.df, aes(x=tr.type, y=tpm.tr.type, fill=Lib))+
    geom_bar(stat="identity", position=position_dodge(preserve = "single"))  +
    ## scale_fill_manual(values=RColorBrewer::brewer.pal(n=12, "Set2"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_bw() +
    coord_flip()
out_dotplot <- file.path(out_dir, "annotation_allSeq_TPM_R10.pdf")
ggsave( filename = out_dotplot, g, width = 7, height = 5) 
