# Parasite_screening_code

# Overview
-This repository contains a collection of Bash scripts designed to map unmapped reads to identify parasite reads from WGS and RNA-seq data, including map unmmaped reads to parasite genome/transcriptome, filter rRNA, and generate various statistics on mapped reads. 

-The scripts utilize tools such as bwa, samtools, and awk for the processing of genomic data. The project involves several workflows for analyzing TopMed DNA and RNA-seq data, including single-cell RNA samples.

## The scripts rely on several bioinformatics tools that must be installed in your environment:
bwa: For aligning sequences to the Loa loa genome.
samtools: For manipulating SAM/BAM files and performing various operations on mapped reads.
awk and sed: For text processing and filtering.
