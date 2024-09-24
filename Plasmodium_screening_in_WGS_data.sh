#!/bin/bash

# Author: Yuanqing Feng <fengyq>
# Email: fengyuanqing2010@gmail.com

# Variables for directories and paths
BASE_DIR=~/genome/parasite/Plasmodium/topmed_DNA  # Base directory for topmed DNA samples
GENOME_DIR=~/genome/parasite/Plasmodium/bwa_index  # Directory for P. falciparum genome index
CHM13_INDEX_DIR=~/genome/CHM13                      # CHM13 genome index directory
STAT_DIR=${BASE_DIR}                          # Directory to store statistics

# ------------ 1. Map unmapped reads to P. falciparum genome for all topmed_DNA samples -----------------
# Loop through two columns of `ghana.bwa_00` file, where $x = sample ID and $y = sample directory
while read -r x y; do
  echo "Sample ID: $x"
  echo "Sample Dir: $y"
  
  # Change to the output directory and create a new directory for each sample
  cd ${BASE_DIR}/ghana || exit
  mkdir -p $x
  cd ${BASE_DIR}/ghana/$x || exit
  
  echo -e "Starting mapping $x to P. falciparum"

  # Map unmapped reads to P. falciparum genome using bwa and samtools
  samtools fastq -n $y | bwa mem -t10 -p ${GENOME_DIR}/Plasmodium.idx - | \
    samtools view -@4 -hb -q 20 | samtools sort -@4 - > ${x}.Plasmodium.q20.bam
  
  # Index the BAM file
  samtools index ${x}.Plasmodium.q20.bam
  
  # Filter BAM and produce statistics for M100 reads (not useful in this version)
  samtools view ${x}.Plasmodium.q20.bam | \
    awk -v OFS="\t" '{$16=$6;gsub(/[0-9]*[S,D,I,H]/,"",$6); split($6,ma,"M"); if(ma[1]+ma[2]+ma[3]+ma[4] > 100) print $0}' | \
    egrep -v "A{20,}|G{20,}|C{20,}|T{20,}" > ${x}.Plasmodium.q20.M100.sam
  
  # Count the number of unique contigs
  cut -f 3 ${x}.Plasmodium.q20.M100.sam | sort | uniq -c | sed s'/refseq/\t/' | sort -k1,1n > ${x}.Plasmodium.q20.M100.stat.tsv
done < ghana.bwa_00

# ------------ 2. Statistic of M151 reads after removing CHM13 mappings ---------------------
for x in $(ls -d N* | grep "NWD"); do

  cd ${BASE_DIR}/ghana/$x || exit

  # 1.1 Filter M151 reads, including ambiguous bases but excluding clipping
  samtools view -h ${x}.Plasmodium.q20.bam | awk '$1 ~"^@" || $6=="151M"' > 151M.sam

  # 1.2 Map M151 reads to CHM13 genome
  samtools fastq 151M.sam | bwa mem -t18 -S -P ${CHM13_INDEX_DIR}/CHM13.idx - | samtools view -h -b > 151M.CHM13.bam

  # 1.3 Filter out M151 reads that map to CHM13
  samtools view -f 4 151M.CHM13.bam | cut -f 1 > 151M.nonCHM13.names
  samtools view -N 151M.nonCHM13.names ${x}.Plasmodium.q20.bam | awk '$6=="151M"' > 151M.rm_CHM13.sam

  # 2. Count M151 reads with 0-9 mismatches
  grep -w -E "NM:i:[0]" 151M.rm_CHM13.sam > 151M.rm_CHM13.nm0.sam
  grep -w -E "NM:i:[0123]" 151M.rm_CHM13.sam > 151M.rm_CHM13.nm3.sam
  grep -w -E "NM:i:[0123456789]" 151M.rm_CHM13.sam > 151M.rm_CHM13.nm9.sam

  # 3. Find paired-end M151 reads
  grep -w 'MC:Z:151M' 151M.rm_CHM13.sam > 151M.rm_CHM13.PE.sam

  # Append statistics to summary files
  wc -l 151M.rm_CHM13.sam >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.stat
  echo "$x" >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.stat

  wc -l 151M.rm_CHM13.nm0.sam >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.nm0.stat
  echo "$x" >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.nm0.stat

  wc -l 151M.rm_CHM13.PE.sam >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.PE.stat
  echo "$x" >> ${STAT_DIR}/Plasmodium_151M.rm_CHM13.PE.stat

  # Clean up temporary files
  rm -f 151M.sam 151M.nonCHM13.names

done
