#!/bin/bash

REFERENCE=~/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
BASEPATH=~/Documents/NGS/RNA_analysis/Mutation2
INTERVAL_LIST=~/Reference/Ampliseq_focus_interval_list/Ampliseq_focus_interval_list_hg38.bed

gatk AddOrReplaceReadGroups \
-I ~/Documents/NGS/RNA_analysis/Mutation/work_dir/Aligned.sortedByCoord.out.bam \
-O ${BASEPATH}/sorted.bam \
-SO coordinate \
-RGID id \
-RGLB library \
-RGPL ILLUMINA \
-RGPU machine \
-RGSM E6146R

gatk MarkDuplicates \
-I ${BASEPATH}/sorted.bam \
-M ${BASEPATH}/marked_dup_metrics.txt \
-O ${BASEPATH}/marked.bam
#CREATE_INDEX=true

gatk SplitNCigarReads \
-R $REFERENCE \
-I ${BASEPATH}/marked.bam \
-O ${BASEPATH}/split.bam \
--read-validation-stringency LENIENT

gatk BaseRecalibrator \
-I ${BASEPATH}/split.bam \
-R $REFERENCE \
--known-sites ~/Reference/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ctat_mutation_lib/dbsnp.vcf.gz \
-O ${BASEPATH}/recal_data.table

gatk ApplyBQSR \
-R $REFERENCE \
-I ${BASEPATH}/split.bam \
--bqsr-recal-file ${BASEPATH}/recal_data.table \
-O ${BASEPATH}/recalibrated.bam

# Variant calling
gatk HaplotypeCaller \
-R $REFERENCE \
-I ${BASEPATH}/recalibrated.bam \
-O ${BASEPATH}/Output_unfiltered.vcf.gz \
--recover-dangling-heads true \
--dont-use-soft-clipped-bases \
-stand-call-conf 20.0

#-L $INTERVAL_LIST

# Variant filtration
gatk VariantFiltration \
-R $REFERENCE \
-V ${BASEPATH}/Output_unfiltered.vcf.gz \
-O ${BASEPATH}/Output_filtered.vcf.gz \
-window 35 \
-cluster 3 \
--filter-name "FS" \
-filter "FS > 30.0" \
--filter-name "QD" \
-filter "QD < 2.0" \
--filter-name "SPLICEADJ" \
-filter "SPLICEADJ < 3"

# Annotate variant using open cravat
oc run ./Output_filtered.vcf.gz -l hg38 -t text excel -a cosmic clinvar
