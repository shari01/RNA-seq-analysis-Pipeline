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

