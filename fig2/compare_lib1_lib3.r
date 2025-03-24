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

out_dir <- "./results/compare_lib1_lib3"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )




dir.lib1 <- "./results/reproducibility/"
dir.lib3 <- "./results/reproducibility_R10/"

### bambu #####
##############



tx.counts.all.tpm.lib1 <- read.table(file.path(dir.lib1, "tx.counts.all.tpm.tsv"), sep="\t", header=TRUE) 
tx.counts.all.tpm.lib3 <- read.table(file.path(dir.lib3, "tx.counts.all.tpm.tsv"), sep="\t", header=TRUE) 


tx.counts.all.tpm.lib1.melt <- tx.counts.all.tpm.lib1 %>%
    reshape2::melt(variable.name="barcode", value.name ="log.CPM") %>%
    mutate(CPM=(10**log.CPM) -1)
tx.counts.all.tpm.lib3.melt <- tx.counts.all.tpm.lib3 %>%
    reshape2::melt(variable.name="barcode", value.name ="log.CPM") %>%
    mutate(CPM=(10**log.CPM) -1)

tx.counts.all.tpm <- full_join(tx.counts.all.tpm.lib1.melt, tx.counts.all.tpm.lib3.melt,
                               by=c("transcript_name","barcode"), suffix=c(".R9",".R10"))



g <- ggplot(tx.counts.all.tpm , aes(x=log.CPM.R9, y=log.CPM.R10)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    facet_grid(vars(barcode)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white")) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) +
    theme(strip.text = element_text(size = 7))
out_dotplot <- file.path(out_dir,"reproducibility_pairwise_TPM_lib1_lib3.png")
ggsave( filename = out_dotplot, g, width = 2.5, height = 6) 




#######################"
#### binned genome style
#######################



dir.lib1.binned <- "./results/reproducibility/"
dir.lib3.binned <- "./results/reproducibility_R10/"

bc.lib1.df <- fread(file.path(dir.lib1.binned, "all_lib1_binned.tsv"), sep="\t") %>%
    as.data.frame %>%
    dplyr::rename(barcode=id, coverage=V4)
bc.lib3.df <- fread(file.path(dir.lib3.binned, "all_lib3_binned.tsv"), sep="\t") %>%
    as.data.frame %>%
    dplyr::rename(barcode=id, coverage=V4)


bc.binned.all <- full_join(bc.lib1.df, bc.lib3.df, by=c("region", "barcode"), suffix=c(".R9",".R10"))  %>%
    mutate(coverage.R9  = ifelse(is.na(coverage.R9), 0, coverage.R9),
           coverage.R10  = ifelse(is.na(coverage.R10), 0, coverage.R10)) %>%
    group_by(barcode) %>% 
    mutate(coverage.R9.tot = sum(coverage.R9, na.rm=TRUE),
           coverage.R10.tot = sum(coverage.R10, na.rm=TRUE)) %>%
    mutate(coverage.R10.s = coverage.R10*(10**6)/coverage.R10.tot,
           coverage.R9.s = coverage.R9*(10**6)/coverage.R9.tot) %>%
    mutate(coverage.R9=log10(1+coverage.R9.s),
           coverage.R10=log10(1+coverage.R10.s)) %>%
    ungroup

g <- ggplot(bc.binned.all, aes(x=coverage.R9, y=coverage.R10)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    facet_grid(vars(barcode)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white")) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) +
    theme(strip.text = element_text(size = 7))
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib1_lib3_bcall.png")
ggsave( filename = out_dotplot, g, width = 2.5, height = 6) 






#####################
#### pooled libraries

R9.bam <- "../mapping_preprocessing/results/get_polyA_bam_files/lib2/trimmed_primary_renamed.bam"

## format R10
R10.bams <- list.files("./results/get_processed_bam_files_R10/", full.names=TRUE, pattern="bam$")  %>%
    paste(collapse=" ")

R10.bam <- file.path(out_dir,"R10.bam" )
command.line <- paste("module load samtools; samtools merge /dev/stdout ", R10.bams,
                      "| samtools sort -@45 >", R10.bam)
system(command.line)

### start
bam.names <- c(R9.bam, R10.bam)

## TPM

se.quant <- bambu(reads = bam.names,  annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 2)


tx.counts.pooled.all.tpm <- assays(se.quant)$CPM %>%
                                    as.data.frame %>%
                                    rownames_to_column("transcript_name") %>%
                                    reshape2::melt() %>% 
                                    mutate(log.tpm=log10(1 + value))   %>%
                                    reshape2::dcast(transcript_name ~ variable, value.var="log.tpm") 
write.table(tx.counts.pooled.all.tpm, file=file.path(out_dir, "tx.counts.pooled.all.tpm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
tx.counts.pooled.all.tpm <- read.table(file.path(out_dir, "tx.counts.pooled.all.tpm.tsv"), sep="\t", header=TRUE) %>%
    dplyr::rename(R9.TPM=trimmed_primary_renamed, R10.TPM=R10)


    

g <- ggplot(tx.counts.pooled.all.tpm, aes(x=R9.TPM, y=R10.TPM)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) +
    theme(strip.text = element_text(size = 7))
out_dotplot <- file.path(out_dir,"reproducibility_pooled_tpm_lib1_lib3_bcall.png")
ggsave( filename = out_dotplot, g, width = 3, height = 3) 

#################
### binned style


genome.size.name <- file.path(out_dir, "hg38.size")
genome.binned.name <- file.path(out_dir, "genome_binned_500.bed")
system(paste("fetchChromSizes hg38 | grep -v '_*_' > ", genome.size.name ))

command <- paste("module load bedops; awk -v FS='\t' -v OFS='\t' '{ print $1, \"0\", $2 }'",  genome.size.name, "| sort-bed - | bedops --chop  500 - >", genome.binned.name )
system(command)


######
######
bed.binned.name <- file.path(out_dir, paste0("R",9:10, "_binned.bed"))

command.lines <- paste("module load samtools; module load bedtools; module load bedops;
  samtools view -b", bam.names,
  " | samtools depth - | awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}' | sort-bed - | ",
  " bedtools map -a", genome.binned.name, " -b - -o sum -c 4 > ",
  bed.binned.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)



bc.bin.list <- mclapply(bed.binned.name, fread, mc.cores=3)  %>%
    setNames(paste0("R",9:10))
bc.all.df <- rbindlist(bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]
bc.all.df <- bc.all.df[,nna:=sum(V4=="." | V4 == 0 ), by=.(region)][nna != 3,][,nna:=NULL][,V4:=ifelse(V4==".", "0", V4 )]

fwrite(bc.all.df, file.path(out_dir, "lib1_lib3_pooled_binned.tsv"), sep="\t", quote=FALSE)
bc.all.df <- fread(file.path(out_dir, "lib1_lib3_pooled_binned.tsv"), sep="\t") %>%
    as.data.frame %>%
    dplyr::rename(lib=id, coverage=V4)


bc.binned.all <- bc.all.df %>%
    reshape2::dcast(region ~ lib, value.var="coverage", fill=0) %>% 
    mutate(R9  = ifelse(is.na(R9), 0, R9),
           R10  = ifelse(is.na(R10), 0, R10)) %>%
    dplyr::filter(!(R9==0 & R10==0)) %>% 
    mutate(R9.tot = sum(R9, na.rm=TRUE),
           R10.tot = sum(R10, na.rm=TRUE)) %>%
    mutate(R10.s = R10*(10**6)/R10.tot,
           R9.s = R9*(10**6)/R9.tot) %>%
    mutate(coverage.R9=log10(1+R9.s),
           coverage.R10=log10(1+R10.s))  

g <- ggplot(bc.binned.all, aes(x=coverage.R9, y=coverage.R10)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) 
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_pooled_lib1_lib3_bcall.png")
ggsave( filename = out_dotplot, g, width = 3, height = 3) 




