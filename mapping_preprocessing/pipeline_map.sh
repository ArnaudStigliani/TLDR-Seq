#!/bin/bash
# set -euo pipefail

module load samtools
module load cutadapt
module load minimap2

export TPMDIR=/maps/projects/sandelin_main/scratch

sequencing_run=$1

results=../../results/TLDRSeq/pipeline_map/$sequencing_run
mkdir -p  $results

fastq_in=$data/all_pass.fastq.gz
echo -e "results in $results\n"

results_trimmed=$results/trimming
results_trimmed_onesided=$results/trimming_onesided
results_mapped=$results/mapping
mkdir -p $results $results_trimmed_onesided $results_trimmed $results_mapped 




#################  trimming

## adapters/barcodes
A2_bc=("TGGCTTGATCCCTCATCTTGTGAAGTTGTTTCGGTCGATTCCGTTTGTAGTCGTCTGT"
       "ACTGGTGCAGCTTTGAACATCTAGAGAGGGTACTATGTGCCTCAGCACCGTACAGCAA"
       "ATGGACTTTGGTAACTTCCTGCGTCACCCACACTTACTTCAGGACGTACCAGTAGAAG"
       "GTTGAATGAGCCTACTGGGTCCTCTTCTGAAGTTCCTGGGTCTTGAACCAGACTTGGT"
       "TGAGAGACAAGATTGTTCGTGGACGACAGACACCGTTCATCGACTTTCGGACGAAGAA"
       "AGATTCAGACCGTCTCATGCAAAGTTCTCAGTCTTCCTCCAGACAAGGCTTACGAAGC"
       "CAAGAGCTTTGACTAAGGAGCATGCCGATCCTTGTGGCTTCTAACTTCATGTCCCAGT"
       "TGGAAGATGAGACCCTGATCTACGGTTTGTCATACTCGTGTGCTCACCCTAGACACCT"
       "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"
       "GCTAGGTCAATCTCCTTCGGAAGTTACAGTCCGAGCCTCATGTGATCTACAACGTCAT"
       "CAGGTTACTCCTCCGTGAGTCTGAACCGAGATCCTACGAATGGAGTGTGTTGGGTAAC"
       "TCAATCAAGAAGGGAAAGCAAGGTCCTGGGAGCATCAGGTAGTAACAGGGAGGAAACA"
       "CATGTTCAACCAAGGCTTCTATGGTAGCTGACTGTCTTCCATACCGACAAGAAGGCAC"
      )

A1rev=AAGAAAGTTGTCGGTGTCTTTGTGCGGATGAACATAGGATAGCGATTCACGGTATGTC

A2_bc_rev=($(printf "%s\n"  "${A2_bc[@]}" | xargs -I {} bash -c " echo {} | tr '[ACGT]' '[TGCA]' | rev"))
A1=$(echo -n $A1rev | tr '[ACGT]' '[TGCA]' | rev)
ix=($(seq 1 13))


### run cutadapt
if [[ $2 =~ "1" ]]
then
    echo "cuting adapters 2 sided"
    parallel -j 13 --link  "cutadapt  -g {1}...{2}  -o $results_trimmed/bc{3}_s1.fastq.gz   -O 8 -m 1 --discard-untrimmed -x s1_bc{3}_  --buffer-size=16000000 $fastq_in > $results_trimmed/bc{3}_s1.log" ::: "${A2_bc[@]}" :::  "$A1rev"  :::  "${ix[@]}" &
    parallel -j 13 --link  "cutadapt  -g {1}...{2}  -o $results_trimmed/bc{3}_s2.fastq.gz   -O 8 -m 1 --discard-untrimmed -x s2_bc{3}_  --buffer-size=16000000 $fastq_in > $results_trimmed/bc{3}_s2.log" ::: "$A1"  ::: "${A2_bc_rev[@]}" :::  "${ix[@]}"
    wait
    parallel 'gzip -dc {}'  ::: $results_trimmed/*.gz |  pigz  > $results/trimmed.fastq.gz
    nerr=$(wc -l  $results_trimmed/*\.log | awk '$1==3{print $2}' | wc -l)
    echo $nerr
    while [ ! $nerr -eq 0 ]
    do
    	cmd=$(wc -l  $results_trimmed/*\.log | awk '$1==3{print $2}' | parallel grep "parameter" {}| sed 's/Command line parameters:/cutadapt/' | awk '{print $0, "> "$5}' | sed 's/fastq.gz/log/3' | sed 's/-o/-j 50 -o/'| tr "[\n]" "[&]" | sed 's/&/ & /g')
    	echo "running cutadapt again to correct few  errors"
    	bash -c "$cmd wait"
    	wait
    	nerr=$(wc -l  $results_trimmed/*\.log | awk '$1==3{print $2}' | wc -l)
    	echo "$nerr starting egain"
    done
    cat  $results_trimmed/*.gz > $results/trimmed.fastq.gz
fi



if [[ $2 =~ "2" ]]
then
    echo "mapping"
    genome_index=../../../shared_data/human_genome/minimap2_index/GRCh38.mmi # path of human genome index
    genome_bedgff=../../../shared_data/human_genome/minimap2_index/anno.bed # path of splice junctions

    echo "minimap2 -ax splice --junc-bed -un $genome_bedgff -t 5  $genome_index $results_trimmed/trimmed.fastq.gz"
    minimap2 -ax splice -un  -t 50  $genome_index $results/trimmed.fastq.gz   | samtools view  -S -b  > $results_mapped/trimmed.bam  2> $results_mapped/trimmed.log 

    # get primary alignment 

    bam=$results_mapped/trimmed.bam
    bam_H=$results_mapped/trimmed.headers
    bam_primmary=$results_mapped/trimmed_primary.bam

    samtools view -H $bam > $bam_H
    samtools view $bam  -F 256 -F 4 | awk -v OFS="\t" '{gsub("AS:i:","",$14)}{print $0}' |\
	awk  -v OFS="\t" '(!($1 in max) || ($14 > max[$1])) {max[$1] = $14; line[$1] = $0} END{for(i in line) print line[i]}' |  awk -v OFS="\t" '{gsub("^","AS:i:",$14)} {print $0}' | cat $bam_H  - | samtools view -S -b   | samtools sort -@ 45 | tee $bam_primmary  | samtools index  - $bam_primmary.bai  
fi

### run cutadapt for one sided adpters
if [[ $2 =~ "3" ]]
then
    echo "cuting adapters 1 sided for log files"
    parallel -j 13 --link  "cutadapt  -a {1}...{2}  -O 8 -m 1 --discard-untrimmed -x s1_bc{3}_ --buffer-size=16000000 $fastq_in  2> $results_trimmed_onesided/bc{3}_s1.log | cat >  /dev/null" ::: "${A2_bc[@]}" ::: "$A1rev" ::: "${ix[@]}" &
    parallel -j 13 --link  "cutadapt  -a {1}...{2}   -O 8 -m 1 --discard-untrimmed -x s2_bc{3}_ --buffer-size=16000000 $fastq_in 2> $results_trimmed_onesided/bc{3}_s2.log | cat > /dev/null" ::: "$A1" ::: "${A2_bc_rev[@]}" ::: "${ix[@]}"
    wait
fi

exit 0
