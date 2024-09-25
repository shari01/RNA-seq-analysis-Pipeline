<h1> RNA-Seq Analysis Pipeline </h1>

This repository contains a complete pipeline for RNA-seq analysis, starting from raw FASTQ files and ending with differential expression analysis and data visualization.

<h1> Workflow Overview </h1>
The pipeline consists of the following steps:

Quality Control: Assess raw sequencing data quality using FastQC.
Read Trimming: Remove low-quality bases and adapter sequences using Trimmomatic.
Alignment: Map reads to a reference genome using HISAT2.
Quantification: Count the number of reads aligned to each gene using featureCounts.
Differential Gene Expression: Perform differential gene expression analysis with DESeq2 in R.
Visualization: Generate PCA plots, volcano plots, and heatmaps to visualize the results.


<h2> How to Use This Pipeline</h2>
1. Prerequisites
You will need the following tools installed to run this pipeline:

FastQC: Quality control for raw reads.
Trimmomatic: Trimming of low-quality reads and adapters.
HISAT2: Read alignment to a reference genome.
SAMtools: For converting, sorting, and indexing SAM/BAM files.
featureCounts: Counting reads mapped to genes.
R and DESeq2: For differential gene expression analysis.
Additionally, ensure that you have:

A reference genome and annotation file (FASTA and GTF).
Your raw sequencing files (FASTQ format).

<h2> 2. Running the Pipeline</h2>
You can run each step of the pipeline manually or automate it with a shell script. Below are the individual steps.


# Step 1: Quality Control (FastQC)
Run FastQC on your raw sequencing data:

#bash
fastqc raw_data/*.fastq.gz -o fastqc_output/
Check the fastqc_output/ folder for the HTML reports.

Step 2: Trimming Reads (Trimmomatic)
Use Trimmomatic to clean your raw reads:

# bash
for file in raw_data/*.fastq.gz
do
  base=$(basename ${file} .fastq.gz)
  java -jar trimmomatic.jar SE -phred33 \
      raw_data/${base}.fastq.gz \
      trimmed_data/${base}_trimmed.fastq.gz \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Aligning Reads (HISAT2)
Align the trimmed reads to the reference genome:

# bash
for file in trimmed_data/*_trimmed.fastq.gz
do
  base=$(basename ${file} _trimmed.fastq.gz)
  hisat2 -x reference_index -U ${file} -S alignment/${base}.sam

# Convert SAM to BAM and sort:

# bash
for file in alignment/*.sam
do
  base=$(basename ${file} .sam)
  samtools view -bS ${file} | samtools sort -o alignment/${base}.sorted.bam
  samtools index alignment/${base}.sorted.bam


# Step 4: Gene Counting (featureCounts)
Count reads mapped to genes:

# bash
featureCounts -T 4 -a genes.gtf -o counts/counts.txt alignment/*.sorted.bam
Step 5: Differential Gene Expression Analysis (DESeq2)
Run DESeq2 in R to perform differential expression analysis. Use the provided R script located in scripts/deseq2_analysis.R.



# Run the script inside R
source("scripts/deseq2_analysis.R")
The differential expression results will be saved as results/differential_expression_results.csv.

# Step 6: Visualization
You can visualize the results using PCA, volcano plots, and heatmaps. These are included in the provided R scripts under the scripts/ folder.

# Run PCA plot script
source("scripts/plot_pca.R")

# Run Volcano plot script
source("scripts/plot_volcano.R")

# Run Heatmap script
source("scripts/plot_heatmap.R")
The plots will be saved in the figures/ directory.

# Contributing
Feel free to open an issue or submit a pull request if you want to contribute to improving this pipeline.

# License
This project is licensed under the MIT License - see the LICENSE file for details.

# Acknowledgments
This pipeline is based on standard RNA-seq analysis practices, and tools used include FastQC, Trimmomatic, HISAT2, featureCounts, and DESeq2. Many thanks to the open-source bioinformatics community for developing these tools.

