# Ampliseq-analysis

This repository contains a collection of bash scripts written by Maximus Yeung as bioinformatics pipeline for DNA and RNA analysis on Ampliseq data.
This is intended for internal use at Department of Pathology, Queen Mary Hospital, and the scripts are provided without being actively supported or maintained.

# DNA somatic short variant calling pipeline

![DNA analysis workflow](https://github.com/maximus3219/Ampliseq-analysis/blob/main/images/DNA.jpg)


The DNA analysis pipeline requires the following software:
- BWA MEM version 0.7.17
- samtools version 1.11
- GATK version 4.1.9

This pipeline follows the GATK best practice for somatic short variants (SNP + indel) calling.
https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-

One point to emphasize is that, since an amplicon-based method was used for library preparation, the --disable-read-filter NotDuplicateReadFilter should be enabled to avoid removing those amplified duplicated reads. However, rarely in one case, disabling the filter may lead to failure to call for some important variants e.g. ERBB2 p.770_771insYVMA, for unknown reason.


# RNA analysis pipeline

![RNA analysis workflow](https://github.com/maximus3219/Ampliseq-analysis/blob/main/images/RNA_analysis.jpg)

The RNA analysis pipeline requires the following software:
- STAR version 2.7.7a
- STAR-Fusion version 1.9.1
- CTAT-Splicing version 0.0.1
- samtools version 1.11
- GATK version 4.1.9
- Open Cravat version 2.2
- Arriba version 2.1.0


STAR is the most widely used splice-aware aligner for RNA-seq data, and is shown to be of high accuracy and ultra-fast (outperforms other aligners by more than a factor of 50 in mapping speed), but it is memory intensive (a minimum of 32GB RAM is recommended).
Different output files generated from the STAR aligner can be utilized for different downstream analysis.
- The chimeric.out.junction contains information about the reads and location at the breakpoint, and can be utilized to find out any fusion
- The bam file generated can be integrated with GATK Best Practices as downstream analysis to call, filter and annotate variants, as well as to prioritize variants that are more likely somatic mutations. However, this may not be practically useful, as the RNA panel covered in Ampliseq does not include all the clinically relevant genes e.g. KRAS. On the other hand, this information about the covered gene panel can supplement those from the DNA-analysis.
- The SJ.out.tab contains information about the junctional reads at the exon-intron boundaries, which together with information from chimeric.out.junction and bam file, could determine any alternative splicing event e.g. MET exon14 skipping, EGFR VIII variant.

STAR aligner:
https://github.com/alexdobin/STAR

STAR-Fusion is one of the most accurate and fastest method for fusion transcript detection on cancer transcriptomes. It can bridge nicely with the output from STAR aligner.
Original paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9
Github for STAR-Fusion:
https://github.com/STAR-Fusion/STAR-Fusion/wiki

The reference datasource for STAR-Fusion can be downloaded via:
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/

Download the one of GRCh38 with latest gencode version - GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play.tar.gz

Reference for CTAT-splicing:
https://github.com/NCIP/CTAT-SPLICING/wiki

Arriba is an award-winning command-line tool for the detection of gene fusions from RNA-Seq data in a clinical research setting. Apart from gene fusions, it can detect other structural rearrangements with potential clinical relevance, such as internal tandem duplications, whole exon duplications, truncations of genes (i.e., breakpoints in introns and intergenic regions). However, it cannot detect intragenic deletions, as they are difficult to distinguish from ordinary splicing.
https://arriba.readthedocs.io/en/latest/


Some of the above softwares depends on Python and/or R scripts and some of their related packages or modules. They should be downloaded according to instructions.
