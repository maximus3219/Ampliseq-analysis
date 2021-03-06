#!/bin/bash

REFERENCE=~/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=~/Reference/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=~/Reference/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/

for file in $(ls ~/Documents/NGS/DNA_analysis/*.fastq.gz | grep "R1"); do 
echo $file; 
FILENAME=$(basename $file | cut -d "_" -f 1); 

mkdir -p ~/Documents/NGS/DNA_analysis/Output/${FILENAME};

#Alignment by bwa
echo "###Aligning '$FILENAME' to Reference genome................";
bwa mem -t 8 \
-R "@RG\tID:$FILENAME\tLB:lib1\tPL:ILLUMINA\tPM:MINISEQ\tSM:$FILENAME\tPU:unit1" \
$REFERENCE $file ${file//R1/R2} | samtools view -bS > ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.bam;

#Sort by coordinate
echo "###Sort by coordinate...........................";
samtools sort ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.bam -o ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_sorted.bam;

#Mark duplicates with Picard tool
echo "####Mark Duplicates..............................";
gatk MarkDuplicates \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_sorted.bam \
-M ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_metrics.txt \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_sorted_marked.bam;

#CREATE_INDEX=true 

#Base recalibration
gatk BaseRecalibrator \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_sorted_marked.bam \
-R $REFERENCE \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recal_table.table \
--known-sites ~/Reference/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ctat_mutation_lib/dbsnp.vcf.gz \
--disable-read-filter NotDuplicateReadFilter;

gatk ApplyBQSR \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_sorted_marked.bam \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recalibrated.bam \
-bqsr ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recal_table.table;

# Create index for BAM file
samtools index ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recalibrated.bam;

# Get Pileup Summaries
echo "#############Getting Pileup Summaries..............";
gatk GetPileupSummaries \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recalibrated.bam \
-V ~/Reference/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
-L ~/Reference/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.getpileupsummaries.table;

# Calculate contamination
echo "##############Calculating contamination..............";
gatk CalculateContamination \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.getpileupsummaries.table \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.contamination.table \
-segments ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.segments.table;

# Call variant
echo "#############Calling somatic variants.................";
gatk Mutect2 \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_recalibrated.bam \
-R $REFERENCE \
-L ~/Reference/Ampliseq_focus_interval_list/Ampliseq_focus_interval_list_hg38.bed \
--germline-resource $GERMLINE  \
--panel-of-normals ~/Reference/PON/pon.vcf.gz \
--f1r2-tar-gz ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.f1r2.tar.gz \
--minimum-allele-fraction 0.01 \
--native-pair-hmm-threads 8 \
--base-quality-score-threshold 18 \
--disable-read-filter NotDuplicateReadFilter \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_unfiltered.vcf.gz \
-bamout ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_realigned.bam;

# --callable-depth 100 \
# --disable-read-filter NotDuplicateReadFilter \
# --af-of-alleles-not-in-resource 0.0000025

echo "###############Learning Read Orientation Model...........";
gatk LearnReadOrientationModel \
-I ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.f1r2.tar.gz \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.read-orientation-model.tar.gz;

#Filter Mutect calls
echo "################Filtering Mutect calls....................";
gatk FilterMutectCalls \
-V ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_unfiltered.vcf.gz \
-R $REFERENCE \
--ob-priors ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.read-orientation-model.tar.gz \
--contamination-table ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.contamination.table \
--tumor-segmentation ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}.segments.table \
--min-allele-fraction 0.01 \
--max-events-in-region 5 \
--max-alt-allele-count 2 \
--min-slippage-length 13 \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_filtered.vcf.gz;

# --unique-alt-read-count 100 \

#Annotate variants
echo "################Annotating variants.......................";
gatk Funcotator \
-R $REFERENCE \
-V ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_filtered.vcf.gz \
-O ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_annotated.vcf.gz \
--output-file-format VCF \
--data-sources-path $ANNOTATION_FILE \
--ref-version hg38;

#Output important informations to tsv files
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t [%AD\t%AF\t%DP]\n' ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_annotated.vcf.gz > ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_core_temp.tsv; 

bcftools query -f '%INFO\n' ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_annotated.vcf.gz | cut -d ";" -f 5 | cut -d "[" -f 2 | awk -F '|' -v OFS='\t' '{print $1,$6,$8,$17,$19,$49}' > ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_annotation_temp.tsv;

paste ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_core_temp.tsv \
~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_annotation_temp.tsv |\
awk -v OFS='\t' '{print$1,$2,$3,$4,$5,$9,$12,$13,$14,$7,$6,$8,$10,$11}'> \
~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_combined_temp.tsv;

cat ~/Reference/Heading.tsv ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_combined_temp.tsv > ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/${FILENAME}_Final.tsv;

rm -f ~/Documents/NGS/DNA_analysis/Output/${FILENAME}/*_temp.tsv;

done
