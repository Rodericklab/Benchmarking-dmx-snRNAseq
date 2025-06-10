# Benchmarking Genetic Demultiplexing Tools for Pooled snRNA-Seq Data
This repository provides the code used for benchmarking genetic demultiplexing methods in single-nucleus RNA sequencing (snRNA-Seq), with a focus on comparing different tools and SNP sources (SNP array vs. RNA-Seq-derived) in real-world and simulated datasets. The benchmarking evaluates both accuracy and computational efficiency, offering guidance for researchers seeking optimal demultiplexing strategies.

## Overview
Demultiplexing genetically pooled snRNA-Seq data allows increased sample throughput and reduces batch effects. This project evaluates four widely used tools:

Vireo
Souporcell
Freemuxlet
scSplit

These tools are tested using SNPs derived from:
SNP arrays (genomic DNA)
Matched bulk RNA-Seq (variant calling via BCFtools, cellSNP, and FreeBayes)

### Repository Structure

.
├── ApplicationDataset.R         # Applies demultiplexing pipeline to real-world datasets
├── Concordance.R                # Concordance analysis between SNP array and RNA-Seq derived genotypes
├── SpeciesDataset.R            # Pipeline for species-mixed sample demultiplexing (e.g., human/mouse)
├── multiplex_6samples/         # Directory containing simulated 6-sample data and configurations
├── LICENSE                     # Licensing information (MIT)
└── README.md                   # This file

### Requirements
R (≥ 4.0.0)
Bioconductor packages: GenomicRanges, VariantAnnotation, AnnotationHub, etc.
VCF and annotation files (from SNP arrays or RNA-Seq)
Tools installed: bcftools, cellSNP-lite, FreeBayes, Vireo, Souporcell, Freemuxlet, scSplit

### Usage
1. Run demultiplexing on simulated or real pooled datasets
source("ApplicationDataset.R")
2. Run concordance analysis between RNA-Seq and SNP-array genotypes
source("Concordance.R")
3. Run species-mixing demultiplexing benchmark
source("SpeciesDataset.R")

Concordance Analysis
The script Concordance.R includes functions to:
Parse VCF files from both SNP array and RNA-Seq
Annotate SNPs with genic features using AnnotationHub
Compute SNP-level and per-sample concordance statistics
Filter for genic SNPs (exonic, CDS, etc.) for targeted benchmarking

Citation
If you use this code, please cite our manuscript:

License
This project is licensed under the MIT License.

