#!/bin/bash

# Author: Yuanqing Feng <fengyq>
# Email: fengyq2018@gmail.com

# Variables for paths
DATA_DIR=~/data/parasite                       # Base directory for data
GENOME_DIR=~/genome/parasite/Plasmodium      # Genome directory for P. falciparum
RRNA_INDEX_DIR=~/genome/rRNA/rna_bwaidx        # Directory for rRNA index
TOOL_DIR=~/tools/prinseq-lite-0.20.4           # Directory for Prinseq tool
PROJECT_DIR=/YourPath  # Project directory for Chao's unmapped data

# Activate conda environment
source activate snakemake

# Split 496 samples into 10 files for parallel processing
split -d -l 50 Chao_sample.ids Chao_S

# ---------------------- filter all unmapped fastq ------------------------
for x in $(cat ${DATA_DIR}/sample.id); do
  # Change to the input directory
  cd ${DATA_DIR}/Chao_unmap_fq || exit
  # Create a directory for the sample
  echo "Sample ID: $x"
  mkdir "${x}" || exit
  cd "${x}" || exit
  # Generate fastq from bam, filter out read-pairs with alignment score [AS]<100 and low complexity prinseq filter
  samtools view -e '[AS]<100' ${PROJECT_DIR}/${x}/${x}.fully.unaligned.cram | samtools fastq | \
    perl ${TOOL_DIR}/prinseq-lite.pl \
    -fastq stdin \
    -lc_method dust -lc_threshold 7 -min_len 50 \
    -graph_data "${x}"_unaligned.stat.gd -log "${x}"_unaligned.stat.log \
    -out_good "${x}"_unaligned.good -out_bad "${x}"_unaligned.bad \
    -trim_ns_left 10 -trim_ns_left 10
  
  # Compress the output fastq files
  bgzip --threads 10 "${x}"_unaligned.good.fastq
  bgzip --threads 10 "${x}"_unaligned.bad.fastq
  
  # Optionally generate the HTML report of reads
  # perl ${TOOL_DIR}/prinseq-graphs.pl -i "${x}"_unaligned.stat.gd -html_all -o "${x}"_unaligned.stat
done

# ---------------------- remap unmapped RNA reads to P. falciparum genome ------------------------
for x in $(cat ${DATA_DIR}/sample.id); do
  # Change to the P. falciparum input directory
  cd ${DATA_DIR}/Plasmodium || exit
  # Create a directory for the sample
  echo "Sample ID: $x"
  mkdir "${x}" || exit
  cd "${x}" || exit
  # Map reads to P. falciparum genome using STAR
  STAR --genomeDir ${GENOME_DIR}/STAR_index/ \
    --runThreadN 10 \
    --readFilesCommand zcat \
    --readFilesIn ${DATA_DIR}/Chao_unmap_fq/${x}/${x}_unaligned.good.fastq.gz \
    --outFileNamePrefix ${x}_Plasmodium_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --outSAMattributes Standard \
    --twopassMode Basic \
    --outFilterMatchNmin 30 \
    --outFilterScoreMinOverLread 0.6 \
    --outFilterMatchNminOverLread 0.6 \
    --limitBAMsortRAM=30000000000
done

# ---------------------- remove ribosomal rRNA and count reads ------------------------
for x in $(cat ${DATA_DIR}/sample.id); do
  cd ${DATA_DIR}/Plasmodium/$x || exit
  # Create a directory for the sample statistics
  echo "Sample ID: $x" >> ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.tsv
  
  # Map reads to the rRNA reference
  samtools fastq *out.bam | bwa mem -t4 -p ${RRNA_INDEX_DIR}/rRNA - | samtools view -hb > rRNA.bam
  # Get statistics of mapped reads
  samtools flagstat rRNA.bam | grep "primary"
  # Extract the read IDs of non-rRNA reads
  samtools view -f 4 rRNA.bam | cut -f 1 > non_rRNA.reads
  # Extract non-rRNA reads mapped to P. falciparum, keeping only primary mapped reads
  samtools view -F 256 -N non_rRNA.reads *out.bam > non_rRNA.sam
  # Count Plasmodium reads and contigs
  wc -l non_rRNA.sam >> ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.tsv
  cut -f 3 non_rRNA.sam | sort | uniq -c >> ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.tsv
done

# Clean up and format the read statistics
grep -E 'ID|rRNA.sam' ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.tsv | \
  sed 'N;s/\n/\t/;s/non_rRNA.sam//; s/Sample ID://' | \
  awk '{print $1,$2}' OFS="\t" > ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.reads

# ---------------------- count contigs ------------------------
for x in $(cat ${DATA_DIR}/sample.id); do
  cd ${DATA_DIR}/Plasmodium/$x || exit
  # Count Plasmodium contigs for each sample
  echo "Sample ID: $x" >> ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.contigs.tsv
  cut -f 3 non_rRNA.sam | sort | uniq -c | wc -l >> ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.contigs.tsv
done

# Clean up and format the contig statistics
sed 'N;s/\n/\t/;s/non_rRNA.sam//; s/Sample ID://' ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.contigs.tsv | \
  awk '{print $1,$2}' OFS="\t" | sort > ${DATA_DIR}/Plasmodium/Plasmodium.map.stat.contigs.tsv2
