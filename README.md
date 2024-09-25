<h1> RNA-Seq Analysis Pipeline </h1>

This repository contains a complete pipeline for RNA-seq analysis, starting from raw FASTQ files and ending with differential expression analysis and data visualization.

<h1> Workflow Overview </h1>
The pipeline consists of the following steps:

<p>Quality Control: Assess raw sequencing data quality using FastQC.</p>
<p>Read Trimming: Remove low-quality bases and adapter sequences using Trimmomatic.</p>
<p>Alignment: Map reads to a reference genome using HISAT2.</p>
<p>Quantification: Count the number of reads aligned to each gene using featureCounts.</p>
<p>Differential Gene Expression: Perform differential gene expression analysis with DESeq2 in R.</p>
<p>Visualization: Generate PCA plots, volcano plots, and heatmaps to visualize the results.</p>

<pre><code>
# 1. Quality Control (FastQC)
# Create directory for FastQC output
mkdir -p fastqc_output

# Run FastQC on all FASTQ files
fastqc raw_data/*.fastq.gz -o fastqc_output/
</code>
<p>Check the fastqc_output/ folder for reports and ensure your data has acceptable quality.
</p>
</pre>

# 2. Trimming Reads (Trimmomatic)
Use Trimmomatic to clean and trim the raw reads.
<pre><code>
# Create directory for trimmed FASTQ files
mkdir -p trimmed_data

# Loop through all FASTQ files and trim them using Trimmomatic
for file in raw_data/*.fastq.gz
do
  base=$(basename ${file} .fastq.gz)
  java -jar trimmomatic.jar SE -phred33 \
      raw_data/${base}.fastq.gz \
      trimmed_data/${base}_trimmed.fastq.gz \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
  <p>This will trim poor-quality sequences and remove adapters.

</p>
</code></pre>

