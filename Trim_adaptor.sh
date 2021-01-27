#!/bin/bash


for file in $(ls ~/Documents/NGS/RNA_analysis/*.fastq.gz | grep "R1"); do 
echo $file; 
FILENAME=$(basename $file | cut -d "_" -f 1); 
echo $FILENAME; 

java -jar ~/Software/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -threads 4 -trimlog trim.log -summary trim.summary \
$file ${file//R1/R2} \
$(dirname $file)/${FILENAME}_trimmed_paired_R1.fastq \
$(dirname $file)/${FILENAME}_trimmed_unpaired_R1.fastq \
$(dirname $file)/${FILENAME}_trimmed_paired_R2.fastq \
$(dirname $file)/${FILENAME}_trimmed_unpaired_R2.fastq \
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:true;
#LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;
            
#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>

done

gzip ~/Documents/NGS/RNA_analysis/*.fastq
