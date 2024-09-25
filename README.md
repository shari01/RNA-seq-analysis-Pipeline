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

<ul>
├── fastqc_output/            # FastQC output files
├── trimmed_data/             # Trimmed FASTQ files
├── alignment/                # BAM files and their indices
├── counts/                   # Gene count matrix
├── results/                  # Differential expression analysis results
├── figures/                  # PCA plots, volcano plots, and heatmaps
├── scripts/                  # Bash and R scripts for each step
└── README.md                 # This README file
</ul>



