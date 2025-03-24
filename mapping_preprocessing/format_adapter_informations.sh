#!/bin/bash


sequencing_run=$1
data=./results/pipeline_map/$sequencing_run

results=./results/format_adapter_information/$sequencing_run
echo $results
tsv_files=$results/tsv_files
mkdir -p $results $tsv_files

grep "Reads with adapters" $data/trimming/*.log | sed 's/.*bc\(.*\)\.log.*:/bc\1/' | sed 's/(\|)\|%\|,//g' | sed 's/_/\t/' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | sort -k2,2 -k1,1g > $results/trimmed.tsv
grep "Reads with adapters" $data/trimming_onesided/*.log | sed 's/.*bc\(.*\)\.log.*:/bc\1/' | sed 's/(\|)\|%\|,//g' | sed 's/_/\t/' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | sort -k2,2 -k1,1g > $results/trimmed_one_sided.tsv

grep " Type: linked" $data/trimming_onesided/*.log | sed 's/.*bc\(.*\)\.log:/bc\1/' | sed 's/Sequence.*; 5..trimmed:/\t/' | sed 's/times; 3..trimmed:/\t/' | sed 's/times//' | sed 's/_/\t/' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | sort -k2,2 -k1,1g > $results/adapter_details.tsv


### get 5' and 3' details
logs=$(find $data/trimming  -name "*.log" )

for file in $logs
do
    file_basename=$(basename $file .log).tsv
    # echo "$file_basename"
    start_line=$(grep -n "Overview of removed sequences at 5' end" $file | awk -v FS=":" '{print $1+1}')
    end_line=$(grep -n "Overview of removed sequences at 3' end" $file | awk -v FS=":" '{print $1}')
    diff=$(($end_line - $start_line - 3))
    tail -n +$start_line $file | head -n $diff | awk -v OFS="\t" '{print $1,$2}' > $tsv_files/5_$file_basename
    tail -n +$(($end_line + 1))  $file | awk -v OFS="\t" '{print $1,$2}' > $tsv_files/3_$file_basename
done


