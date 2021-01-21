#!/bin/bash

STARINDEX=~/Reference/Gencode/STAR_index/
CTAT_RESOURCE_LIB=~/Reference/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/

for file in $(ls ~/Documents/NGS/RNA_analysis/*.fastq.gz | grep "R1"); do 
echo $file; 
FILENAME=$(basename $file | cut -d "_" -f 1); 
echo $FILENAME; 

# Use STAR-Aligner to align paired reads to reference genome
STAR --genomeDir $STARINDEX \
--readFilesIn $file ${file//R1/R2} \
--outReadsUnmapped None \
--runThreadN 4 \
--twopassMode Basic \
--readFilesCommand "gunzip -c" \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 8 \
--chimOutJunctionFormat 1 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 \
--alignInsertionFlush Right \
--alignSplicedMateMapLminOverLmate 0 \
--alignSplicedMateMapLmin 30 \
--outFileNamePrefix $(dirname $file)/Output/$FILENAME/ \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts;

### --outSAMstrandField intronMotif \ # include for potential use with StringTie for assembly
### --chimSegmentMin 12 \ # ** essential to invoke chimeric read detection & reporting **
### --chimOutJunctionFormat 1  \ # **essential** includes required metadata in Chimeric.junction.out file.
### --alignMatesGapMax 100000 \ # avoid readthru fusions within 100k
### --alignSJstitchMismatchNmax 5 -1 5 5 \  # settings improved certain chimera 

# Use STAR-Fusion to call any fusion
STAR-Fusion --genome_lib_dir $CTAT_RESOURCE_LIB \
            -J $(dirname $file)/Output/$FILENAME/Chimeric.out.junction \
            --output_dir $(dirname $file)/Output/$FILENAME/;
            
## Set minimum FFPM to filter; Default 0.1; set to 0 to turn off

mkdir -p $(dirname $file)/Output/$FILENAME/Splicing;

samtools index $(dirname $file)/Output/$FILENAME/Aligned.sortedByCoord.out.bam;

# Use CTAT-splicing to call any alternative splicing
~/Software/CTAT-SPLICING.v0.0.1/STAR_to_cancer_introns.py \
--SJ_tab_file $(dirname $file)/Output/$FILENAME/SJ.out.tab \
--chimJ_file $(dirname $file)/Output/$FILENAME/Chimeric.out.junction \
--vis \
--bam_file $(dirname $file)/Output/$FILENAME/Aligned.sortedByCoord.out.bam \
--output_prefix $(dirname $file)/Output/$FILENAME/Splicing/${filename} \
--sample_name ${filename} \
--ctat_genome_lib $CTAT_RESOURCE_LIB;

# --min_total_reads, default 5 minimum read alignments required 
# Cases which MET exon14 skipping was not expected, were called with total reads of 5-16, though these may be real cases that have not been confirmed with other methods.
# On the other hand, EGFRVIII variants were called at read of 20
# So a cut-off of 20 may be appropriate

cat $(dirname $file)/Output/$FILENAME/Splicing/${filename}.cancer.introns > ${FILEBASE}/${filename}/Splicing/${filename}.cancer.introns.tsv;

done
