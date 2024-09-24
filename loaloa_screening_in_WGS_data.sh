#!/bin/bash
# Author: Yuanqing Feng <fengyq>
# Email: fengyuanqing2010@gmail.com

# Variables for directories and paths
GENOME_DIR=~/genome/parasite/loaloa          # Base directory for Loa loa genome
BWA_INDEX_DIR=${GENOME_DIR}/bwa_index        # Directory for BWA index
TOPMED_DNA_DIR=${GENOME_DIR}/topmed_DNA      # Base directory for TopMed DNA samples
RRNA_INDEX_DIR=~/genome/rRNA/rna_bwaidx      # rRNA index directory

# ---------------------- 2. Generate bwa index for parasites ------------------------------------
cd ${BWA_INDEX_DIR} || exit
bwa index -p loaloa.idx -a rb2 ./loaloa.fa

# ---------------------- 3. Map unmapped reads to Loa loa genome --------------------------------
# Loop through the sample IDs and directories provided in `bwa_00` file
while read -r x y; do
  echo "Sample ID: $x"
  echo "Sample Dir: $y"

  # Create directories for each sample
  cd ${TOPMED_DNA_DIR}/tzetbt || exit
  mkdir -p $x
  cd ${TOPMED_DNA_DIR}/tzetbt/$x || exit

  echo "Start mapping $x to Loa loa"

  # Map unmapped reads to Loa loa genome
  samtools fastq -n $y | bwa mem -t10 -p ${BWA_INDEX_DIR}/loaloa.idx - | \
    samtools view -@4 -hb -q 20 | samtools sort -@4 - > ${x}.loaloa.q20.bam
  samtools index ${x}.loaloa.q20.bam

  # Filter BAM and generate stats for reads with >100M
  samtools view ${x}.loaloa.q20.bam | \
    awk -v OFS="\t" '{$16=$6; gsub(/[0-9]*[S,D,I,H]/, "", $6); split($6, ma, "M"); if(ma[1]+ma[2]+ma[3]+ma[4] > 100) print $0}' | \
    egrep -v "A{20,}|G{20,}|C{20,}|T{20,}" > ${x}.loaloa.q20.M100.sam
  cut -f 3 ${x}.loaloa.q20.M100.sam | sort | uniq -c | sed 's/refseq/\t/' | sort -k1,1n > ${x}.loaloa.q20.M100.stat.tsv
  wc -l ${x}.loaloa.q20.M100.sam

  # Filter reads with >150M and exclude specific contigs
  awk '$6=="151M" && $3!="NW_018108197.1"' ${x}.loaloa.q20.M100.sam > ${x}.loaloa.q20.M150.sam
  wc -l ${x}.loaloa.q20.M150.sam
  cut -f 3 ${x}.loaloa.q20.M150.sam | sort | uniq -c | sed 's/refseq/\t/' | sort -k1,1n > ${x}.loaloa.q20.M150.stat.tsv
done < bwa_00

# ---------------------- 4. Remove ribosomal rRNA and count reads ------------------------
for x in $(cat ${TOPMED_DNA_DIR}/sample.ids); do
  cd ${TOPMED_DNA_DIR}/$x || exit
  echo "Sample ID: $x" >> ${TOPMED_DNA_DIR}/loa.map.stat.tsv
  
  # Map reads to the rRNA reference
  samtools fastq *q20.bam | bwa mem -t4 -p ${RRNA_INDEX_DIR}/rRNA - | samtools view -hb > rRNA.bam
  
  # Count primary mapped reads
  samtools flagstat rRNA.bam | grep "primary"
  
  # Extract non-rRNA reads
  samtools view -f 4 rRNA.bam | cut -f 1 > non_rRNA.reads
  samtools view -F 256 -N non_rRNA.reads *q20.bam > non_rRNA.sam

  # Count Loa loa reads and contigs
  wc -l non_rRNA.sam >> ${TOPMED_DNA_DIR}/loa.map.stat.tsv
  cut -f 3 non_rRNA.sam | sort | uniq -c >> ${TOPMED_DNA_DIR}/loa.map.stat.tsv
done

grep -E 'ID|rRNA.sam' loa.map.stat.tsv | sed 'N;s/\n/\t/;s/non_rRNA.sam//;s/Sample ID://' | awk '{print $1,$2}' OFS="\t" > loa.map.stat.reads

# ------------------------- 5. Statistic of 150M and 151M reads ----------------------------------------
cd ${TOPMED_DNA_DIR}
for x in $(ls | grep "NWD"); do
  cd ${TOPMED_DNA_DIR}/$x || exit
  echo "Processing sample: $x"
  
  # Keep 151M reads and remove unwanted contig
  awk '$6=="151M" && $3!="NW_018108197.1"' non_rRNA.sam > non_rRNA.M150.sam
  cut -f 3 non_rRNA.M150.sam | sort | uniq -c | sed 's/refseq/\t/' | sort -k1,1n > non_rRNA.M150.stat.tsv
done

for x in $(ls | grep "NWD"); do
  cd ${TOPMED_DNA_DIR}/$x || exit
  echo "${x}" >> ../non_rRNA.M150.reads.tsv
  wc -l non_rRNA.M150.sam >> ../non_rRNA.M150.reads.tsv
  echo "${x}" >> ../non_rRNA.M150.contigs.tsv
  wc -l non_rRNA.M150.stat.tsv >> ../non_rRNA.M150.contigs.tsv
done

# Generate summary statistics
sed 'N;s/\n/\t/;s/non_rRNA.*sam//;' non_rRNA.M150.reads.tsv | awk '{print $1,$2}' OFS="\t" | sort > non_rRNA.M150.reads.stat
sed 'N;s/\n/\t/;s/non_rRNA.*tsv//;' non_rRNA.M150.contigs.tsv | awk '{print $1,$2}' OFS="\t" | sort > non_rRNA.M150.contigs.stat

# Combine counts from logs
awk -v FS="\t" -v OFS="\t" 'NR==FNR {a[$1]=$0; next} {if($1 in a) print $0, a[$1]}' M150_Chrs_counts.log.sorted M150_reads_counts.log.sorted | cut -f 1,2,4 > M150_Chrs_counts.log.tsv


