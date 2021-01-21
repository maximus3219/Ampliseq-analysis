# Ampliseq-analysis

This repository contains a collection of bash scripts written by Maximus Yeung as bioinformatics pipeline for DNA and RNA analysis on Ampliseq data.
This is intended for internal use at Department of Pathology, Queen Mary Hospital, and the scripts are provided without being actively supported or maintained.

[My image](https://github.com/maximus3219/Ampliseq-analysis/blob/main/DNA.jpg)


<img
src=“https://github.com/maximus3219/Ampliseq-analysis/blob/main/DNA.jpg”
raw=true
alt=“Subject Pronouns”
style=“margin-right: 10px;”
/>

<img width=“964” alt=“DNA workflow” src=“https://github.com/maximus3219/Ampliseq-analysis/blob/main/DNA.jpg?raw=true”>


The DNA analysis pipeline requires the following software:
- BWA MEM version 0.7.17
- samtools version 1.11
- GATK version 4.1.9

This pipeline follows the GATK best practice for somatic short variants (SNP + indel) calling.
https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-

One point to emphasize is that, since an amplicon-based method was used for library preparation, the --disable-read-filter NotDuplicateReadFilter should be enabled to avoid removing those amplififed duplicated reads. However, rarely in one case, disabling the filter may lead to failure to call for some important variants e.g. ERBB2 p.770_771insYVMA.

The RNA analysis pipeline requires the following software:
- STAR version 2.7.7a
- STAR-Fusion version 1.9.1
- CTAT-Splicing version 0.0.1
- samtools version 1.11
- GATK version 4.1.9
- Open Cravat version 2.2

STAR aligner:
https://github.com/alexdobin/STAR

Reference for STAR-Fusion:
https://github.com/STAR-Fusion/STAR-Fusion/wiki

The reference datasource for STAR-Fusion can be downloaded via:
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/

Download the one of GRCh38 with latest gencode version - GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play.tar.gz

Reference for CTAT-splicing:
https://github.com/NCIP/CTAT-SPLICING/wiki


Some of the above softwares depends on Python and/or R scripts and some of their related packages or modules. They should be downloaded according to instructions.
