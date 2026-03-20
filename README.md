# transcriptomics-pipeline

RNA-Seq data analysis pipeline covering read alignment and differential expression analysis, developed as part of a bioinformatics course project.

---
 
## Overview
 
This project consists of two parts:
 
1. **Part 1 – Read Alignment**: Align a RNA-Seq FASTQ file to the hg19 reference genome using HISAT2.
2. **Part 2 – Differential Expression Analysis**: Perform Differential expression analysis on a prostate cancer RNA-Seq dataset (pre- vs post-docetaxel treatment) using metaseqR2.
 
---

## Part 1

### Requirements
- HISAT2 — spliced-aware RNA-Seq aligner
- samtools — for SAM/BAM conversion and indexing
- hg19 HISAT2 genome index

### Step 1: Set up and data acquisition
Create a working directory and retrieve the required files to perform further analysis.
```
# Create a working directory
mkdir align && cd align

# Get the FASTQ file
wget http://epigenomics.fleming.gr/~panos/appbio/human.fastq.gz

# Get and uncompress the index
wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
tar -xvf hg19_genome.tar.gz

# Get and unzip the aligner itself
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1-Linux_x86_64.zip
unzip hisat2-2.2.1-Linux_x86_64.zip

# Test the installation
./hisat2-2.2.1/hisat2 --help

# Get samtools (if not already in the system you are working at)
wget https://github.com/samtools/samtools/releases/download/1.23/samtools-1.23.tar.bz2
tar -xvf samtools-1.23.tar.bz2
cd samtools-1.23
./configure --without-curses --disable-bz2 --disable-lzma && make
cd ..

# Test the installaiton
./samtools-1.23/samtools --help

```

### Step 2: Alignment with HISAT2
HISAT2 is a fast and sensitive alignment program for mapping RNA-Sequencing reads to a reference genome.

```
# Align the FASTQ file to hg19 with HISAT2
./hisat2-2.2.1/hisat2 -x genome -U human.fastq.gz -S HISAT_OUT.sam
```
Flags:
- -x: the path to the genome index
- -U: the input FASTQ file 
- -S: outfile in SAM format

### Step 3: Convert SAM to BAM
```
# Convert the SAM file to BAM file
samtools view -bS HISAT_OUT.sam > HISAT_OUT.bam
```
- -b: output in BAM format
- -S: input in SAM format

### Step 4: Sorting and Indexing
```
# Sort BAM
samtools sort human.bam -o human.sorted.bam

# Index BAM
samtools index human.sorted.bam
```

## Part 2

### Requirements
 
- R 
- metaseqR2

Metaseq2 is an R package for analyzing RNA-seq data and producing interactive reports with results and quality control metrics. It integrates multiple normalization methods and statistical algorithms into a single pipeline.

The specific dataset used in this part consists transcriptomics data of pre- and post-docetaxel treatment biopsies from patients with advanced hormone-naive prostate cancer treated with docetaxel chemotherapy.

### Step 1: Create the targets file
The targets file is a simple text tab-delimited file describing the samples and experimental conditions. It is required by the metaseqR2 pipeline to point out which samples correspond to which BAM files.
A targets file has the following columns: 1) Sample name, 2) BAM file path, 3) Condition, 4) Paired (if TRUE), 5) Strand direction.

```
samplename	filename	condition	paired	stranded
PT_1	[path]/Post_1.bam	Post	paired	forward
PT_2	[path]/Post_2.bam	Post	paired	forward
PT_3	[path]/Post_3.bam	Post	paired	forward
PR_1	[path]/Pre_1.bam	Pre	paired	forward
PR_2	[path]/Pre_2.bam	Pre	paired	forward
PR_3	[path]/Pre_3.bam	Pre	paired	forward
```

### Step 2: DEA with metaseqR2
Differential gene expression analysis aims to identify genes that exhibit significant changes in expression levels between two or more conditions.

In R, run the following command.
```
# example command
library(metaseqR2)

targetsFile <- read.delim("[path]/targets.txt", header=TRUE, stringsAsFactors=FALSE)

# Define the contrast: Post vs Pre
theContrasts <- c("Post_vs_Pre")

metaseqr2(
    sampleList=targetsFile,    # sample metadata
    contrast=theContrasts,     # comparison to perform
    org="hg19",                # reference genome
    countType="exon",          # count reads per exon, summarize to gene level
    normalization="deseq2",    # DESeq2 normalisation
    statistics="deseq2",       # DESeq2 statistical testing
    figFormat="png",
    qcPlots=c(
        "mds","biodetection","countsbio","saturation",
        "correl","boxplot","meandiff","meanvar","volcano","mastat"
    ),
    exportWhere="./metaseqr2_output",       # output folder
    pcut=0.05,                              # adjusted p-value cutoff for DE calls
    restrictCores=0.25,
    exportWhat=c("annotation","p_value","adj_p_value","fold_change",
        "counts","flags"),
    exportScale=c("natural","log2","rpgm"),     # export counts in multiple scales
    exportValues="normalized",                  # export normalised counts 
    saveGeneModel=TRUE,
    reportTop=0.1                               # show top 10% DE genes in the HTML report
)

```

## Acknowledgements
 
This project was developed as part of the Applied Bioinformatics course in the MSc programme. The project design, dataset, tutorial materials, and HPC environment were provided by the course instructor.
