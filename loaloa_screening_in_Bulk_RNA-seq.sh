#!/bin/bash
# Author: Yuanqing Feng <fengyq>
# Email: fengyq2018@gmail.com

# Variables for paths
DATA_DIR=/YOURPATH           # Base directory for data
TOOL_DIR=/YOURPATH  # Directory for Prinseq tool
GENOME_DIR=~/genome/parasite/loaloa   # Directory for Loa loa genome
RRNA_INDEX_DIR=~/genome/rRNA/rna_bwaidx  # Directory for rRNA index
PROJECT_DIR=/YOURPATH  # Project directory

# Activate conda environment
source activate snakemake

# ---------------------- filter all unmapped fastq ------------------------
for x in $(cat ${DATA_DIR}/sample.ids); do
  # Change to the input directory
  cd ${DATA_DIR}/unmap_fq || exit
  # Create a directory for the sample
  echo "Sample ID: $x"
  mkdir "${x}" || exit
  cd "${x}" || exit
  # Generate fastq from bam, filter out read-pairs with alignment score [AS]<100 and low complexity prinseq filter
  samtools view -e '[AS]<100' ${PROJECT_DIR}/${x}/${x}.fully.unaligned.cram | samtools fastq | perl ${TOOL_DIR}/prinseq-lite.pl \
    -fastq stdin \
    -lc_method dust -lc_threshold 7 -min_len 50 \
    -graph_data "${x}"_unaligned.stat.gd -log "${x}"_unaligned.stat.log \
    -out_good "${x}"_unaligned.good -out_bad "${x}"_unaligned.bad \
    -trim_ns_left 10 -trim_ns_left 10
  
  # Compress the output fastq files
  bgzip --threads 10 "${x}"_unaligned.good.fastq
  bgzip --threads 10 "${x}"_unaligned.bad.fastq
  
  # Optionally generate the HTML report for the reads
  # perl ${TOOL_DIR}/prinseq-graphs.pl -i "${x}"_unaligned.stat.gd -html_all -o "${x}"_unaligned.stat
done

# ----------------------remap unmapped RNA reads to loaloa genome ------------------------
for x in $(cat ${DATA_DIR}/sample.ids); do
  # Change to the loaloa input directory
  cd ${DATA_DIR}/loaloa || exit
  # Create a directory for the sample
  echo "Sample ID: $x"
  mkdir "${x}" || exit
  cd "${x}" || exit
  # Map reads to Loa loa genome using STAR
  STAR --genomeDir ${GENOME_DIR}/STAR_index/ \
    --runThreadN 10 \
    --readFilesCommand zcat \
    --readFilesIn ${DATA_DIR}/unmap_fq/${x}/${x}_unaligned.good.fastq.gz \
    --outFileNamePrefix ${x}_loaloa_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --outSAMattributes Standard \
    --twopassMode Basic \
    --outFilterMatchNmin 30 \
    --outFilterScoreMinOverLread 0.6 \
    --outFilterMatchNminOverLread 0.6 \
    --limitBAMsortRAM=30000000000
done

# ----------------------remove ribosomal rRNA and count reads ------------------------
for x in $(cat ${DATA_DIR}/sample.ids); do
  cd ${DATA_DIR}/loaloa/$x || exit
  # Create a directory for the sample
  echo "Sample ID: $x" >> ${DATA_DIR}/loaloa/loaloa.map.stat.tsv
  
  # Map reads to the rRNA reference
  samtools fastq *out.bam | bwa mem -t4 -p ${RRNA_INDEX_DIR}/rRNA - | samtools view -hb > rRNA.bam
  # Get statistics of mapped reads
  samtools flagstat rRNA.bam | grep "primary"
  # Extract the read IDs of non-rRNA reads
  samtools view -f 4 rRNA.bam | cut -f 1 > non_rRNA.reads
  # Extract non-rRNA reads mapped to Loa loa, keeping only primary mapped reads
  samtools view -F 256 -N non_rRNA.reads *out.bam > non_rRNA.sam
  # Count Plasmodium reads and contigs
  wc -l non_rRNA.sam >> ${DATA_DIR}/loaloa/loaloa.map.stat.tsv
  cut -f 3 non_rRNA.sam | sort | uniq -c >> ${DATA_DIR}/loaloa/loaloa.map.stat.tsv
done

grep -E 'ID|rRNA.sam' ${DATA_DIR}/loaloa/loaloa.map.stat.tsv | sed 'N;s/\n/\t/;s/non_rRNA.sam//; s/Sample ID://' | awk '{print $1,$2}' OFS="\t" > ${DATA_DIR}/loaloa/loaloa.map.stat.reads

# ----------------------count contigs ------------------------
for x in $(cat ${DATA_DIR}/sample.ids); do
  cd ${DATA_DIR}/loaloa/$x || exit
  # Count Plasmodium contigs for each sample
  echo "Sample ID: $x" >> ${DATA_DIR}/loaloa/loaloa.map.stat.contigs.tsv
  cut -f 3 non_rRNA.sam | sort | uniq -c | wc -l >> ${DATA_DIR}/loaloa/loaloa.map.stat.contigs.tsv
done
