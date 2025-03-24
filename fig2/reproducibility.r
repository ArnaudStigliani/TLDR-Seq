rm(list=ls())
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

out_dir <- "./results/reproducibility/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE )

gtf.name <- "../data/gencode.v39.primary_assembly.annotation.sorted.gtf" # path to gtf file
fas.name <- "../data/GRCh38.primary_assembly.genome.fa" # path to fasta genome  file


bc.tsv.name <- "../mapping_preprocessing/results/demultiplex_read_names/lib1/reads_specifications.tsv"
reads.curated.name <- "../fig1/results/get_polyA_bam_files/lib1/trimmed_primary_renamed.bam"



### bambu #####
##############

annotations <- prepareAnnotations(gtf.name)
saveRDS(annotations, file.path(out_dir, "annotations.rds")) 
annotations <- readRDS(file.path(out_dir, "annotations.rds"))

### demultiplex reads


bc.tsv.list <- read.table(bc.tsv.name, sep="\t", header=FALSE) %>%
    dplyr::select(bc=V2, read.name=V3) %>% 
    dplyr::filter(bc %in% c("bc1", "bc2", "bc3")) %>%
    split(f=.$bc) %>%
    lapply(function(x) dplyr::select(x, read.name)) %>%
    mapply(write.table, .,file.path(out_dir, paste0(names(.),".list")), MoreArgs=list(col.names=FALSE, row.names=FALSE, quote=FALSE),  SIMPLIFY=FALSE)



bc.list.name <-  file.path(out_dir, paste0("bc",1:3,".list"))
bc.bam.name <- file.path(out_dir, paste0("bc",1:3, ".bam"))

command.lines <- paste("module load samtools;
  samtools view -b", reads.curated.name,
  "-N", bc.list.name, " > ",
  bc.bam.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)




se.quant <- bambu(reads = bc.bam.name,  annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 20)


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
             upper = list(continuous = wrap("cor", method = "pearson", stars = FALSE)))+
    theme_bw() +
    scale_x_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib1.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 

M <- tx.counts.all.tpm %>% column_to_rownames("transcript_name") %>% cor(method="pearson")
M.lib1.counts <- M[lower.tri(M)]





########## binning method
#########################


####### ####### ####### ######
#### generate genome binned

genome.size.name <- file.path(out_dir, "hg38.size")
genome.binned.name <- file.path(out_dir, "genome_binned_500.bed")
system(paste("fetchChromSizes hg38 | grep -v '_*_' > ", genome.size.name ))

command <- paste("module load bedops; awk -v FS='\t' -v OFS='\t' '{ print $1, \"0\", $2 }'",  genome.size.name, "| sort-bed - | bedops --chop 500 - >", genome.binned.name )
system(command)




######
######

bc.tsv.name <-  file.path(out_dir, paste0("bc",1:3,".list"))

bc.bed.binned.name <- file.path(out_dir, paste0("bc",1:3, "_binned.bed"))

command.lines <- paste("module load samtools; module load bedtools; module load bedops;
  samtools view -b", reads.curated.name,
  "-N", bc.tsv.name,
  " | samtools depth - | awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}' | sort-bed - | ",
  " bedtools map -a", genome.binned.name, " -b - -o sum -c 4 > ",
  bc.bed.binned.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)

######
######

bc.bin.list <- mclapply(bc.bed.binned.name, fread, mc.cores=12)  %>%
    setNames(paste0("bc",1:3))
bc.all.df <- rbindlist(bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]
bc.all.df <- bc.all.df[,nna:=sum(V4=="." | V4 == 0 ), by=.(region)][nna != 3,][,nna:=NULL][,V4:=ifelse(V4==".", "0", V4 )]

fwrite(bc.all.df, file.path(out_dir, "all_lib1_binned.tsv"), sep="\t", quote=FALSE)
bc.all.df <- fread(file.path(out_dir, "all_lib1_binned.tsv"), sep="\t") %>%
    as.data.frame


bc.all.df.dcast <- bc.all.df %>% mutate(V4=(V4+1)/500) %>% 
    reshape2::dcast(region ~ id , value.var ="V4", fill=1/500) %>%
    dplyr::filter((bc1+bc2+bc3)> 3) %>% 
    mutate(bc1=log10(1+ bc1)) %>% 
    mutate(bc2=log10(1+ bc2)) %>%
    mutate(bc3=log10(1+ bc3))

g <- ggpairs(bc.all.df.dcast, columns= 2:4, title ="log10 coverage in binned genome",
             lower = list(continuous = wrap("points", size=0.2, alpha = 0.2)),
             upper = list(continuous = wrap(ggally_cor, stars = FALSE))) +
    theme_bw() +
    scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib1.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 

######################################
#####################################
### restart for  lib2 ########"


bc_lib2.tsv.name <- "../mapping_preprocessing/results/demultiplex_read_names/lib1/reads_specifications.tsv"
reads_lib2.curated.name <-  "../fig1/resullts/get_polyA_bam_files/lib1/trimmed_primary_renamed.bam"


bc_lib2.tsv.list <- read.table(bc_lib2.tsv.name, sep="\t", header=FALSE) %>%
    dplyr::select(bc=V2, read.name=V3) %>% 
    dplyr::filter(bc %in% paste0("bc",1:12)) %>%
    split(f=.$bc) %>%
    mclapply(function(x) dplyr::select(x, read.name), mc.cores=12) %>%
    mcmapply(write.table, .,file.path(out_dir, paste0(names(.),"_lib2.list")), MoreArgs=list(col.names=FALSE, row.names=FALSE, quote=FALSE),  SIMPLIFY=FALSE, mc.cores=12)



bc_lib2.list.name <-  file.path(out_dir, paste0("bc",1:12,"_lib2.list"))
bc_lib2.bam.name <- file.path(out_dir, paste0("bc",1:12, "_lib2.bam"))

command.lines <- paste("module load samtools;
  samtools view -b", reads_lib2.curated.name,
  "-N", bc_lib2.list.name, " > ",
  bc_lib2.bam.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)


se_lib2.quant <- bambu(reads = bc_lib2.bam.name,   annotations = annotations,
                          genome = fas.name, discovery = FALSE, ncore = 20)


tx_lib2.counts.all.tpm <- assays(se_lib2.quant)$CPM %>%
                                    as.data.frame %>%
                                    rownames_to_column("transcript_name") %>%
                                    reshape2::melt() %>% 
                                    mutate(log.tpm=log10(1 + value))   %>%
                                    reshape2::dcast(transcript_name ~ variable, value.var="log.tpm") 
write.table(tx_lib2.counts.all.tpm, file=file.path(out_dir, "tx_lib2.counts.all.tpm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
tx_lib2.counts.all.tpm <- read.table(file.path(out_dir, "tx_lib2.counts.all.tpm.tsv"), sep="\t", header=TRUE) 



g <- ggpairs(tx_lib2.counts.all.tpm, columns= 2:13, title ="log10 TPM", 
             lower = list(continuous = wrap("points", size=0.2, alpha = 0.2)),
             upper = list(continuous = wrap("cor", method = "pearson", stars = FALSE)))+
    theme_bw() +
    scale_x_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib2.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 

M <- tx_lib2.counts.all.tpm %>% column_to_rownames("transcript_name") %>% cor(method="pearson")
M.lib2.counts <- M[lower.tri(M)]



###### binning methods
######

lib2.bc.tsv.name <- file.path(out_dir, paste0("bc",1:12, "_lib2.list"))


####
####
lib2.bc.bed.binned.name <- file.path(out_dir, paste0("bc",1:12, "_lib2_binned.bed"))

command.lines <- paste("module load samtools; module load bedtools; module load bedops;
  samtools view -b", reads_lib2.curated.name,
  "-N", lib2.bc.tsv.name,
  " | samtools depth - | awk -v OFS='\t' '$1~/^chr/{print $1, $2, $2+1, $3}' | sort-bed - | ",
  " bedtools map -a", genome.binned.name, " -b - -o sum -c 4 > ",
  lib2.bc.bed.binned.name)

all.command.lines <- paste(paste(command.lines, collapse=" & \n"), " & wait")
system(all.command.lines)

## read files

lib2.bc.bin.list <- mclapply(lib2.bc.bed.binned.name, fread, mc.cores=12)  %>%
    setNames(paste0("bc",1:12))
lib2.bc.all.df <- rbindlist(lib2.bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]
lib2.bc.all.df <- lib2.bc.all.df[,nna:=sum(V4=="." | V4 == 0 ), by=.(region)][nna != 12,][,nna:=NULL][,V4:=ifelse(V4==".", "0", V4 )]

fwrite(lib2.bc.all.df, file.path(out_dir, "all_lib2_binned.tsv"), sep="\t", quote=FALSE)
lib2.bc.all.df <- fread(file.path(out_dir, "all_lib2_binned.tsv"), sep="\t") %>%
    as.data.frame


lib2.bc.all.df.dcast <- lib2.bc.all.df %>% mutate(V4=(V4+1)/500) %>% 
    reshape2::dcast(region ~ id , value.var ="V4", fill=1/500) %>%
    dplyr::filter((bc1+bc2+bc3+bc4+bc5+bc6+bc7+bc8+bc9+bc10+bc11+bc12)> 12) %>% 
    mutate(bc1=log10(1+ bc1)) %>% 
    mutate(bc2=log10(1+ bc2)) %>%
    mutate(bc3=log10(1+ bc3)) %>%
    mutate(bc4=log10(1+ bc4)) %>% 
    mutate(bc5=log10(1+ bc5)) %>%
    mutate(bc6=log10(1+ bc6)) %>%
    mutate(bc7=log10(1+ bc7)) %>% 
    mutate(bc8=log10(1+ bc8)) %>%
    mutate(bc9=log10(1+ bc9)) %>%
    mutate(bc10=log10(1+ bc10)) %>% 
    mutate(bc11=log10(1+ bc11)) %>%
    mutate(bc12=log10(1+ bc12))


g <- ggpairs(lib2.bc.all.df.dcast, columns= 2:13, title ="log10 coverage in binned genome",
             lower = list(continuous = wrap("points", size=0.2, alpha = 0.2)),
             upper = list(continuous = wrap(ggally_cor, stars = FALSE))) +
    theme_bw() +
    scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white"))
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib2.png")
ggsave( filename = out_dotplot, g, width = 10, height = 10) 

M <- lib2.bc.all.df.dcast %>% column_to_rownames("region") %>% cor %>% round(.,2)
cols <- c("#BB4444",  "#FFFFFF", "#4477AA") %>% rev
g <- ggcorrplot(M, type = "upper", 
                lab = TRUE, colors=cols, lab_size=2)
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_corrplot_lib2.png")
ggsave(filename = out_dotplot, g, width = 5, height = 5)
M <- lib2.bc.all.df.dcast %>% column_to_rownames("region") %>% cor 
M.lib2.binned <- M[lower.tri(M)]




##############################
##### All correlations
### read old binned files

### Compile all correlations

M.df <- list(lib1.TPM = M.lib1.counts, lib1.binned = M.lib1.binned,
             lib2.TPM = M.lib2.counts, lib2.binned = M.lib2.binned)  %>%
    lapply(data.frame) %>% 
    lapply(setNames, "pearson.corr") %>%
    mapply(mutate, ., id = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    separate(id, into = c("lib", "metric"), sep="\\.") 


g <- ggplot(M.df, aes(x=lib, y=pearson.corr,  colour=lib)) +
    theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), size=2) +
    facet_wrap(~ metric, ncol=2) +
    theme(strip.background = element_rect(fill="white")) +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
    ylab("Pearson pairwise correlations")
out_boxplot <- file.path(out_dir,"correlation_dotplot.pdf")
ggsave( filename = out_boxplot, g, width = 5, height = 5) 

g <- ggplot(M.df, aes(x=lib, y=pearson.corr,  colour=lib)) +
    theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), size=2) +
    facet_wrap(~ metric, ncol=2) +
    theme(strip.background = element_rect(fill="white")) +
    ylab("Pearson pairwise correlations")
out_boxplot <- file.path(out_dir,"correlation_dotplot_default.pdf")
ggsave( filename = out_boxplot, g, width = 5, height = 5) 


#### Compare lib1 and lib2
### counts
tx.counts_lib1.all.df.melt <- tx.counts.all.tpm %>%
    reshape2::melt()
tx.counts_lib2.all.df.melt <- tx_lib2.counts.all.tpm %>%
    reshape2::melt()



tx.counts_lib1_2.all <- bind_rows(lib1 = tx.counts_lib1.all.df.melt, lib2 = tx.counts_lib2.all.df.melt, .id="lib")  %>%
    dplyr::rename(tpm=value, bc=variable) %>%
    mutate(tpm=(10**tpm)-1)



tx.counts_joined <-  full_join(tx.counts_lib1_2.all %>% dplyr::filter(lib=="lib1"),
                               tx.counts_lib1_2.all %>% dplyr::filter(lib=="lib2"),
                               by="transcript_name", multiple = "all") 


tx.counts_joined.plot  <- tx.counts_joined %>% 
    dplyr::select(transcript_name,  lib1 = tpm.x, lib2 =  tpm.y, bc.lib1 = bc.x, bc.lib2 = bc.y)  %>%
    group_by(bc.lib1, bc.lib2 ) %>%
    group_by(transcript_name) %>%
    mutate(lib1.tr = sum(lib1)/12, lib2.tr = sum(lib2)/3) %>%
    ungroup  %>%
    dplyr::filter(lib1.tr + lib2.tr > 10) %>% 
    mutate(lib1=log10(1 + lib1), lib2=log10(1 + lib2)) 
    

g <- ggplot(tx.counts_joined.plot , aes(x=lib1, y=lib2)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    facet_grid(vars(bc.lib1), vars(bc.lib2)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white")) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) +
    theme(strip.text = element_text(size = 7))
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib1_lib2.png")
ggsave( filename = out_dotplot, g, width = 15, height = 5) 


#### pool libraries

tx.counts_lib1_2.pooled <- tx.counts_lib1_2.all %>%
    group_by(lib, transcript_name) %>%
    summarize(mean.tpm= mean(tpm)) %>%
    mutate(log.mean.tpm=log10(mean.tpm + 1)) %>%
    reshape2::dcast(transcript_name ~ lib, value.var="log.mean.tpm" )


cor.stat <- cor.test(tx.counts_lib1_2.pooled$lib1, tx.counts_lib1_2.pooled$lib2, method="pearson") %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(tx.counts_lib1_2.pooled , aes(x=lib1, y=lib2)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    ggtitle(cor.stat) +
    xlab("Mean expression of amplified libraries (log10 TPM)") +
    ylab("Mean expression of non-amplified libraries (log10 TPM)")
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib1_lib2_pooled.png")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 
out_dotplot <- file.path(out_dir,"reproducibility_tpm_lib1_lib2_pooled.pdf")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 



#### binned genome


bc.bin.list <- mclapply(bc.bed.binned.name, fread, mc.cores=12)  %>%
    setNames(paste0("bc",1:3))
bc.all.df <- rbindlist(bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]


lib2.bc.bin.list <- mclapply(lib2.bc.bed.binned.name, fread, mc.cores=12)  %>%
    setNames(paste0("bc",1:12))
lib2.bc.all.df <- rbindlist(lib2.bc.bin.list , idcol="id")[,region:=paste0(V1,":",V2,"-",V3)][, c("V1", "V2", "V3"):=NULL]

lib1_2.df <- rbindlist(list(lib1=bc.all.df, lib2=lib2.bc.all.df), idcol="lib")[,nna:=sum(V4=="." | V4 == 0 ), by=.(region)][nna != 15,][,nna:=NULL][,V4:=ifelse(V4==".", "0", V4 )]


fwrite(lib1_2.df, file.path(out_dir, "lib1_lib2_binned.tsv"), sep="\t", quote=FALSE)
lib1_2.df <- fread(file.path(out_dir, "lib1_lib2_binned.tsv"), sep="\t") %>%
    as.data.frame




lib1_2.df.joined.2 <-  full_join(lib1_2.df %>% dplyr::filter(lib=="lib1"),
                               lib1_2.df %>% dplyr::filter(lib=="lib2"),
                               by="region", multiple = "all")


lib1_2.df.plot  <- lib1_2.df.joined.2 %>% 
    dplyr::select(region,lib1 = V4.x, lib2 = V4.y, bc.lib1=id.x, bc.lib2=id.y)  %>%
    dplyr::filter(!(((lib2 == 0 ) & is.na(lib1)) | ((lib1 == 0 ) & is.na(lib2)))) %>%    
    group_by(bc.lib1, bc.lib2 ) %>%
    mutate(lib1.tot = sum(lib1, na.rm=TRUE), lib2.tot = sum(lib2, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(lib1 = ifelse(is.na(lib1), 0 , lib1), lib2 = ifelse(is.na(lib2), 0, lib2)) %>% 
    mutate(lib2.s = lib2*(10**6)/lib2.tot, lib1.s = lib1*(10**6)/lib1.tot) %>%
    mutate(lib1=lib1.s, lib2=lib2.s) %>%
    group_by(region) %>%
    mutate(lib1.region = sum(lib1)/12, lib2.region = sum(lib2)/3) %>%
    ungroup  %>%
    dplyr::filter(lib1.region + lib2.region > 10) %>% 
    mutate(lib1=log10(1 + lib1.s), lib2=log10(1 + lib2.s)) 
    


g <- ggplot(lib1_2.df.plot , aes(x=lib1, y=lib2)) +
    theme_bw() +
    geom_point(size=0.2, alpha=0.2) +
    facet_grid(vars(bc.lib1), vars(bc.lib2)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(strip.background = element_rect(fill="white")) +
    stat_cor(method="pearson", size=2, aes(label = ..r.label..)) +
    theme(strip.text = element_text(size = 7))
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib1_lib2.png")
ggsave( filename = out_dotplot, g, width = 12, height = 4) 



######## do the same for the pooled library

lib1_2.df.pooled <- lib1_2.df %>%
    group_by(lib, region) %>%
    summarize(s.V4= mean(V4))


lib1_2.df.pooled.dcast <- lib1_2.df.pooled %>%
    mutate(s.V4=s.V4/500) %>%
    reshape2::dcast(region ~ lib, value.var="s.V4") %>%
    mutate(lib2.tot = sum(lib2, na.rm=TRUE), lib1.tot = sum(lib1, na.rm=TRUE)) %>%
    mutate(lib2.s = lib2*(10**6)/lib2.tot, lib1.s = lib1*(10**6)/lib1.tot) %>%
    mutate(lib1.log = log(1+lib1.s), lib2.log = log(1+lib2.s)) 

cor.stat <- cor.test(lib1_2.df.pooled.dcast$lib1.log, lib1_2.df.pooled.dcast$lib2.log, method="pearson") %>% 
    {paste0("R = ", round(.$estimate, 2),
            ifelse(.$p.value<2.2e-16, ", p < 2.2e-16",
                   paste0(", p = ", formatC(.$p.value, format = "e", digits = 1))))}
g <- ggplot(lib1_2.df.pooled.dcast, aes(x = lib1.log, y = lib2.log)) +
    theme_bw() +
    geom_point(alpha=0.1, size=0.1) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10),labels= c(0,2,4,6,8,10)) +
    scale_y_continuous(breaks = c(0,2,4,6,8,10),labels= c(0,2,4,6,8,10)) +
    ggtitle(cor.stat) +
    xlab("Mean coverage of amplified libraries (log10)") +
    ylab("Mean coverage of non-amplified libraries (log10)")
out_dotplot <- file.path(out_dir,"reproducibility_binned_genome_lib1_lib2_pooled.png")
ggsave( filename = out_dotplot, g, width = 5, height = 5) 
