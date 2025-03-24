library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2)
library(GGally)
library(data.table)
library(parallel)
library(ggpubr)
library(ggcorrplot)
library(bambu)

out_dir <- "./results/reproducibility_R10"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )




reads.curated.dir <- "../../results/get_processed_bam_files_R10/"

gtf.name <- "../data/gencode.v39.primary_assembly.annotation.sorted.gtf" # path to gtf file
fas.name <- "../data/GRCh38.primary_assembly.genome.fa" # path to fasta genome  file

### bambu #####
##############

annotations <- prepareAnnotations(gtf.name)
saveRDS(annotations, file.path(out_dir, "annotations.rds")) 
annotations <- readRDS(file.path(out_dir, "annotations.rds"))

### demultiplex reads


bc.bam.name <- file.path(reads.curated.dir, paste0("bc",1:3, ".bam"))



se.quant <- bambu(reads = bc.bam.name,  annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 3)


tx.counts.all.tpm <- assays(se.quant)$CPM %>%
                                    as.data.frame %>%
                                    rownames_to_column("transcript_name") %>%
                                    reshape2::melt() %>% 
                                    mutate(log.tpm=log10(1 + value))   %>%
                                    reshape2::dcast(transcript_name ~ variable, value.var="log.tpm") 
write.table(tx.counts.all.tpm, file=file.path(out_dir, "tx.counts.all.tpm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
tx.counts.all.tpm <- read.table(file.path(out_dir, "tx.counts.all.tpm.tsv"), sep="\t", header=TRUE) 


g <- ggpairs(tx.counts.all.tpm, columns= 2:4, title ="log10 TPM", 
             lower = list(continuous = wrap("points", size=0.2, alpha = 0.2)),
             ## upper = list(continuous = wrap(ggally_cor, stars = FALSE))) +
             upper = list(continuous = wrap("cor", method = "pearson", stars = FALSE)))+
    theme_bw() +
    scale_x_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib3.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 

M <- tx.counts.all.tpm %>% column_to_rownames("transcript_name") %>% cor(method="pearson") 
M.lib1.counts <- M[lower.tri(M)]


#######################"
#### binned genome style
#######################

genome.size.name <- file.path(out_dir, "hg38.size")
genome.binned.name <- file.path(out_dir, "genome_binned_500.bed")
system(paste("fetchChromSizes hg38 | grep -v '_*_' > ", genome.size.name ))

command <- paste("module load bedops; awk -v FS='\t' -v OFS='\t' '{ print $1, \"0\", $2 }'",  genome.size.name, "| sort-bed - | bedops --chop  500 - >", genome.binned.name )
system(command)


######
######
bc.bed.binned.name <- file.path(out_dir, paste0("bc",1:3, "_binned.bed"))

command.lines <- paste("module load samtools; module load bedtools; module load bedops;
  samtools view -b", bc.bam.name,
  " | samtools depth - | awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}' | sort-bed - | ",
  " bedtools map -a", genome.binned.name, " -b - -o sum -c 4 > ",
  bc.bed.binned.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)

######
######

bc.bin.list <- mclapply(bc.bed.binned.name, fread, mc.cores=3)  %>%
    setNames(paste0("bc",1:3))
bc.all.df <- rbindlist(bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]
bc.all.df <- bc.all.df[,nna:=sum(V4=="." | V4 == 0 ), by=.(region)][nna != 3,][,nna:=NULL][,V4:=ifelse(V4==".", "0", V4 )]

fwrite(bc.all.df, file.path(out_dir, "all_lib3_binned.tsv"), sep="\t", quote=FALSE)
bc.all.df <- fread(file.path(out_dir, "all_lib3_binned.tsv"), sep="\t") %>%
    as.data.frame


bc.all.df.dcast <- bc.all.df %>% mutate(V4=(V4+1)/500) %>% 
    reshape2::dcast(region ~ id , value.var ="V4", fill=1/500) %>%
    dplyr::filter((bc1+bc2+bc3)> 3) %>% 
    mutate(bc1=log10(1+ bc1)) %>% 
    mutate(bc2=log10(1+ bc2)) %>%
    mutate(bc3=log10(1+ bc3))


g <- ggpairs(bc.all.df.dcast, columns= 2:4, title ="log10 coverage in binned genome",
             lower = list(continuous = wrap("points", size=0.2, alpha = 0.2)),
             ## upper = list(continuous = wrap(ggally_cor, stars = FALSE))) +
             upper = list(continuous = wrap("cor", method = "pearson", stars = FALSE))) +
    theme_bw() +
    scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib3.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 



M <- bc.all.df.dcast %>% column_to_rownames("region") %>% cor(method="pearson")
M.lib3.binned <- M[lower.tri(M)]
