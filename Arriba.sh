#!/bin/bash

REFERENCE=~/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
ANNOTATION=~/Reference/Gencode/gencode.v36.primary_assembly.annotation.gtf
BLACKLIST=~/Software/arriba_v2.1.0/database/blacklist_hg38_GRCh38_v2.1.0.tsv.gz
KNOWN_FUSION=~/Software/arriba_v2.1.0/database/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz
PROTEIN_DOMAINS=~/Software/arriba_v2.1.0/database/protein_domains_hg38_GRCh38_v2.1.0.gff3
STARINDEX=~/Reference/Gencode/STAR_index/


for file in $(ls ~/Documents/NGS/RNA_analysis/*.fastq.gz | grep "R1"); do 
echo $file; 
FILENAME=$(basename $file | cut -d "_" -f 1); 
echo $FILENAME; 

mkdir -p $(dirname $file)/Output/$FILENAME/;

STAR \
    --runThreadN 8 \
    --genomeDir $STARINDEX \
    --readFilesIn ${file} ${file//R1/R2} \
    --readFilesCommand "gunzip -c" \
    --outFileNamePrefix $(dirname $file)/Output/$FILENAME/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outBAMcompression 0 \
    --outFilterMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50;
   
#    --twopassMode Basic \
    
arriba \
    -x $(dirname $file)/Output/$FILENAME/Aligned.sortedByCoord.out.bam \
    -o $(dirname $file)/Output/$FILENAME/fusions.tsv \
    -O $(dirname $file)/Output/$FILENAME/fusions.discarded.tsv \
    -a $REFERENCE -g $ANNOTATION \
    -b $BLACKLIST -k $KNOWN_FUSION -t $KNOWN_FUSION -p $PROTEIN_DOMAINS;
            
done
