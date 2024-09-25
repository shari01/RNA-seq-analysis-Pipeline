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

# 3. Aligning Reads (HISAT2)
Next, align the trimmed reads to a reference genome using HISAT2.

<pre><code>
# Create directory for alignment results
mkdir -p alignment

# Run HISAT2 on all trimmed FASTQ files
for file in trimmed_data/*_trimmed.fastq.gz
do
  base=$(basename ${file} _trimmed.fastq.gz)
  hisat2 -x reference_index -U ${file} -S alignment/${base}.sam
done

# Convert SAM to BAM, sort, and index the files
for file in alignment/*.sam
do
  base=$(basename ${file} .sam)
  samtools view -bS ${file} | samtools sort -o alignment/${base}.sorted.bam
  samtools index alignment/${base}.sorted.bam
done
</code></pre>

# 4. Gene Counting (featureCounts)
Once alignment is done, count the number of reads mapped to each gene using featureCounts.
<pre><code>
# Create directory for counts
mkdir -p counts

# Use featureCounts to count reads aligned to genes
featureCounts -T 4 -a reference_annotation.gtf -o counts/gene_counts.txt alignment/*.sorted.bam
  <p>The output gene_counts.txt contains the count matrix for all your samples.

</p>
</code></pre>

# 5. Differential Gene Expression (DESeq2)
Now, we move to R for differential expression analysis using DESeq2.

<pre><code>
# Load necessary libraries
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)

# Load the counts data (from featureCounts) and remove the first few metadata rows
countData <- read.table("counts/gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")
countData <- countData[ , -c(1:5)] # Remove non-sample columns

# Prepare metadata (sample information)
# You need to modify this according to your experimental design (e.g., condition information)
sampleNames <- colnames(countData)
condition <- c("control", "treated", "control", "treated") # Modify according to your samples
colData <- data.frame(row.names = sampleNames, condition = factor(condition))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results (differential expression)
res <- results(dds)

# Sort results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Save the results
write.csv(as.data.frame(resOrdered), file="results/differential_expression_results.csv")

# Visualizations

# 1. PCA Plot
rld <- rlog(dds) # Regularized-log transformation for PCA
plotPCA(rld, intgroup = "condition")

# 2. Volcano Plot
EnhancedVolcano(resOrdered,
                lab = rownames(resOrdered),
                x = "log2FoldChange",
                y = "pvalue",
                title = 'Volcano Plot',
                pCutoff = 0.05, # Adjust p-value cutoff
                FCcutoff = 1.5) # Adjust log2 fold change cutoff

# 3. Heatmap for top 30 genes
top30genes <- head(rownames(resOrdered[order(resOrdered$padj), ]), 30)
heatmapData <- assay(rld)[top30genes, ]
pheatmap(heatmapData, cluster_rows = TRUE, cluster_cols = TRUE)
</code></pre>

# 6. Running DESeq2 and Visualizing Results
PCA Plot: This will help you assess how samples cluster based on gene expression.
Volcano Plot: Shows log2 fold change vs. p-value, highlighting significant genes.
Heatmap: Displays the top 30 most differentially expressed genes, with clustering of samples and genes.

# 7. Final Structure of Files
The pipeline is complete, and your folder should look something like this:

<pre><code>
├── raw_data/                 # Your raw FASTQ files
├── fastqc_output/            # FastQC output files
├── trimmed_data/             # Trimmed FASTQ files
├── alignment/                # SAM/BAM files after alignment
├── counts/                   # FeatureCounts output (gene count matrix)
├── results/                  # DESeq2 results and visualizations
├── scripts/                  # Scripts for trimming, alignment, counting, and DESeq2 analysis
├── figures/                  # PCA, volcano, and heatmap images
└── README.md                 # Instructions (provided earlier)
</code></pre>
