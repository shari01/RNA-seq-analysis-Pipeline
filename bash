#!/bin/bash
# use linux to run the file Execute the Bash script to process RNA-seq data and generate counts:
#bash rna_seq_pipeline.sh

# Set up directories
mkdir -p fastq raw_data trimmed_data alignment counts results figures

# 1. Quality Control with FastQC
fastqc raw_data/*.fastq.gz -o fastq/

# 2. Trimming low-quality bases using Trimmomatic
for file in raw_data/*.fastq.gz
do
  base=$(basename ${file} .fastq.gz)
  java -jar trimmomatic.jar SE -phred33 \
      raw_data/${base}.fastq.gz \
      trimmed_data/${base}_trimmed.fastq.gz \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# 3. Align reads to reference genome using HISAT2
hisat2-build reference_genome.fa reference_index
for file in trimmed_data/*_trimmed.fastq.gz
do
  base=$(basename ${file} _trimmed.fastq.gz)
  hisat2 -x reference_index -U ${file} -S alignment/${base}.sam
done

# 4. Convert SAM to BAM, sort, and index
for file in alignment/*.sam
do
  base=$(basename ${file} .sam)
  samtools view -bS ${file} | samtools sort -o alignment/${base}.sorted.bam
  samtools index alignment/${base}.sorted.bam
done

# 5. Count features using featureCounts
featureCounts -T 4 -a genes.gtf -o counts/counts.txt alignment/*.sorted.bam
